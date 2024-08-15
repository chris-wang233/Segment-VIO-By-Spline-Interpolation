#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>

#include "estimator.h"
#include "parameters.h"
#include "utility/visualization.h"
//addition
#include "./spline_algo.cpp"
#include "../../pose_graph/src/so3.h"

vector<Estimator> Spline_estimator;
bool pose_op_lable=false;
string spline_op_lable;
Estimator estimator;

std::condition_variable con;
double current_time = -1;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::PointCloudConstPtr> relo_buf;
int sum_of_wait = 0;

double sum_of_acc_Norm = 0;
int sum_of_acc = 0;
double sum_of_gy_Norm = 0;
int sum_of_gy = 0;
int sum_of_all = 0;
bool spline_acc_adjust = true;
bool spline_adjust_ = false;
double chazhi_sum = 0;

bool start_lable_=false;
Vector3d start_pose_test;
Vector3d end_pose_test;
Vector3d start_velocity_test;
Vector3d end_velocity_test;
double start_time;
double end_time;
Spline spline_test;

vector<Vector3d> velocity_buffer_tmp;
vector<Matrix3d> Rotation_buffer_tmp;

Vector3d Head_pose;
Vector3d Mid_pose;
Vector3d Tail_pose;

double Head_time_;
double Mid_time_;
double Tail_time_;

bool acc_mid_change = true;

std::mutex m_buf;
std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Gy;
Eigen::Vector3d tmp_acc;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = 0;

void predict(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};

    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};

    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - estimator.g;

    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);

    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - estimator.g;

    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);

    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void update()
{
    TicToc t_predict;
    latest_time = current_time;
    tmp_P = estimator.Ps[WINDOW_SIZE];
    tmp_Q = estimator.Rs[WINDOW_SIZE];
    tmp_V = estimator.Vs[WINDOW_SIZE];
    tmp_Ba = estimator.Bas[WINDOW_SIZE];
    tmp_Bg = estimator.Bgs[WINDOW_SIZE];
    acc_0 = estimator.acc_0;
    gyr_0 = estimator.gyr_0;

    queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;
    for (sensor_msgs::ImuConstPtr tmp_imu_msg; !tmp_imu_buf.empty(); tmp_imu_buf.pop())
        predict(tmp_imu_buf.front());

}

std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>>
getMeasurements()
{
    std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;

    while (true)
    {
        if (imu_buf.empty() || feature_buf.empty())
            return measurements;

        if (!(imu_buf.back()->header.stamp.toSec() > feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            //ROS_WARN("wait for imu, only should happen at the beginning");
            sum_of_wait++;
            return measurements;
        }

        if (!(imu_buf.front()->header.stamp.toSec() < feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            ROS_WARN("throw img, only should happen at the beginning");
            feature_buf.pop();
            continue;
        }
        sensor_msgs::PointCloudConstPtr img_msg = feature_buf.front();
        feature_buf.pop();

        std::vector<sensor_msgs::ImuConstPtr> IMUs;
        while (imu_buf.front()->header.stamp.toSec() < img_msg->header.stamp.toSec() + estimator.td)
        {
            IMUs.emplace_back(imu_buf.front());
            imu_buf.pop();
        }
        IMUs.emplace_back(imu_buf.front());
        if (IMUs.empty())
            ROS_WARN("no imu between two image");
        measurements.emplace_back(IMUs, img_msg);
    }
    return measurements;
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();

    m_buf.lock();
    imu_buf.push(imu_msg);
    m_buf.unlock();
    con.notify_one();

    last_imu_t = imu_msg->header.stamp.toSec();

    {
        std::lock_guard<std::mutex> lg(m_state);
        predict(imu_msg);
        std_msgs::Header header = imu_msg->header;
        header.frame_id = "world";
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            pubLatestOdometry(tmp_P, tmp_Q, tmp_V, header);
    }
}


void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    if (!init_feature)
    {
        //skip the first detected feature, which doesn't contain optical flow speed
        init_feature = 1;
        return;
    }
    m_buf.lock();
    feature_buf.push(feature_msg);
    m_buf.unlock();
    con.notify_one();
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while(!feature_buf.empty())
            feature_buf.pop();
        while(!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        m_estimator.lock();
        estimator.clearState();
        estimator.setParameter();
        m_estimator.unlock();
        current_time = -1;
        last_imu_t = 0;
    }
    return;
}

void relocalization_callback(const sensor_msgs::PointCloudConstPtr &points_msg)
{
    //printf("relocalization callback! \n");
    m_buf.lock();
    relo_buf.push(points_msg);
    m_buf.unlock();
}

// thread: visual-inertial odometry JB_融合代码
void process()
{
    Spline Spline;//JB:Spline对象
    double sum_reproject=0;
    double mid_sum_reproject=0;
    int num_reproject = 0;
    double sum_V=0;

    int num_reproject_test = 0;
    double sum_V_test=0;

    double sum_gyc = 0;
    double sum_acc = 0;


    int buffer_sount = 0;
    while (true)
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;
        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
                 {
            return (measurements = getMeasurements()).size() != 0;
                 });
        lk.unlock();
        m_estimator.lock();
        
        for (auto &measurement : measurements)
        {
            auto img_msg = measurement.second;
            
            double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;

            int i = 1;

            double tmp_acc_bias=0,tmp_gyc_bias=0;
            for (auto &imu_msg : measurement.first)
            {
                double t = imu_msg->header.stamp.toSec();
                            
                double img_t = img_msg->header.stamp.toSec() + estimator.td;
                if (t <= img_t)
                { 
                    
                    if (current_time < 0)
                        current_time = t;
                    double dt = t - current_time;
                    ROS_ASSERT(dt >= 0);
                    current_time = t;
                    dx = imu_msg->linear_acceleration.x;
                    dy = imu_msg->linear_acceleration.y;
                    dz = imu_msg->linear_acceleration.z;
                    rx = imu_msg->angular_velocity.x;
                    ry = imu_msg->angular_velocity.y;
                    rz = imu_msg->angular_velocity.z;
                    Vector3d tmp_acc_;
                    estimator.processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz),tmp_acc_bias,tmp_gyc_bias);
                    //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
                    //printf("dimu:  w: %f\n", Vector3d(rx, ry, rz).squaredNorm());
                    //tmp_Gy=Vector3d(rx, ry, rz);
                    //tmp_acc=Vector3d(dx, dy, dz);
                    if(i==1){
                        tmp_acc = tmp_acc_;
                    }
                   // printf("*******%ffront acc is:%f\n",current_time,tmp_gyc_bias);
                }
                else
                {
                    //printf("inside,%f\n",img_msg->header.stamp.toSec());
                    double dt_1 = img_t - current_time;
                    double dt_2 = t - img_t;
                    current_time = img_t;
                    ROS_ASSERT(dt_1 >= 0);
                    ROS_ASSERT(dt_2 >= 0);
                    ROS_ASSERT(dt_1 + dt_2 > 0);
                    double w1 = dt_2 / (dt_1 + dt_2);
                    double w2 = dt_1 / (dt_1 + dt_2);
                    dx = w1 * dx + w2 * imu_msg->linear_acceleration.x;
                    dy = w1 * dy + w2 * imu_msg->linear_acceleration.y;
                    dz = w1 * dz + w2 * imu_msg->linear_acceleration.z;
                    rx = w1 * rx + w2 * imu_msg->angular_velocity.x;
                    ry = w1 * ry + w2 * imu_msg->angular_velocity.y;
                    rz = w1 * rz + w2 * imu_msg->angular_velocity.z;
                    //IMU数据更新！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
                    //IMU预积分计算出大概的pose,并且计算出加速度
                    
                    estimator.processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz),tmp_acc_bias,tmp_gyc_bias);
                    //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
                    //printf("dimu:  w: %f\n", Vector3d(rx, ry, rz).squaredNorm());
                    //tmp_Gy=Vector3d(rx, ry, rz);
                    //tmp_acc=Vector3d(dx, dy, dz);
                    //printf("*******%ffront acc is:%f\n",current_time,tmp_gyc_bias);
                }
               i++;
            }
            // printf("***************************acc is%f\n",tmp_acc_bias);
           //  printf("***************************gyc is%f\n\n",tmp_gyc_bias);
            // set relocalization frame
            sensor_msgs::PointCloudConstPtr relo_msg = NULL;
            while (!relo_buf.empty())
            {
                relo_msg = relo_buf.front();
                relo_buf.pop();
            }
            if (relo_msg != NULL)
            {
                vector<Vector3d> match_points;
                double frame_stamp = relo_msg->header.stamp.toSec();
                for (unsigned int i = 0; i < relo_msg->points.size(); i++)
                {
                    Vector3d u_v_id;
                    u_v_id.x() = relo_msg->points[i].x;
                    u_v_id.y() = relo_msg->points[i].y;
                    u_v_id.z() = relo_msg->points[i].z;
                    match_points.push_back(u_v_id);
                }
                Vector3d relo_t(relo_msg->channels[0].values[0], relo_msg->channels[0].values[1], relo_msg->channels[0].values[2]);
                Quaterniond relo_q(relo_msg->channels[0].values[3], relo_msg->channels[0].values[4], relo_msg->channels[0].values[5], relo_msg->channels[0].values[6]);
                Matrix3d relo_r = relo_q.toRotationMatrix();
                int frame_index;
                frame_index = relo_msg->channels[0].values[7];
                estimator.setReloFrame(frame_stamp, frame_index, match_points, relo_t, relo_r);
            }

            ROS_DEBUG("processing vision data with stamp %f \n", img_msg->header.stamp.toSec());

            TicToc t_s;
            map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> image;
            for (unsigned int i = 0; i < img_msg->points.size(); i++)
            {       
                int v = img_msg->channels[0].values[i] + 0.5;
                int feature_id = v / NUM_OF_CAM;
                int camera_id = v % NUM_OF_CAM;
                double x = img_msg->points[i].x;
                double y = img_msg->points[i].y;
                double z = img_msg->points[i].z;
                double p_u = img_msg->channels[1].values[i];
                double p_v = img_msg->channels[2].values[i];
                double velocity_x = img_msg->channels[3].values[i];
                double velocity_y = img_msg->channels[4].values[i];
                ROS_ASSERT(z == 1);
                Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
                xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
                image[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
            }

/*             if(estimator.Spline_segment_lable_==false && estimator.frame_count==10&&Spline.compute_spline_lable(tmp_V,tmp_Gy))
            {   //如果目前没有待优化的段落，先记录段首数据,直到出现待优化段落
                    estimator.Spline_Head_Position_ = estimator.Ps[WINDOW_SIZE];
                    estimator.Spline_Head_time_ = current_time;
                    estimator.Spline_Head_Velocity_ = estimator.Vs[WINDOW_SIZE];
                    estimator.Spline_Head_header_ = img_msg->header;
                    
            }
            if(!Spline.compute_spline_lable(tmp_V,tmp_Gy)&& estimator.frame_count==10)
            {
                
                    estimator.Spline_segment_lable_ = true;//true 代表存在一个队列等待插值
                    estimator.Spline_mid_Position_.push_back(estimator.Ps[WINDOW_SIZE]);//将中间帧传入
                    estimator.Spline_mid_time_.push_back(current_time);
                    estimator.Spline_mid_header_.push_back(img_msg->header);
                    estimator.Spline_mid_Velocity_.push_back(estimator.Vs[WINDOW_SIZE]);
            }; */
            
            //传入Spline数据
            //estimator.Sping_lable_estimatior_ = Spline.compute_spline_lable(tmp_V,tmp_Gy);
            //printf("%f Spline Befor v is:%f,%f,%f\n",current_time,estimator.Vs[WINDOW_SIZE].x(),estimator.Vs[WINDOW_SIZE].y(),estimator.Vs[WINDOW_SIZE].z());
            //printf("outside,%f\n",img_msg->header.stamp.toSec());
            double current_rej = 0;
            estimator.current_bias = 0;
            double delta_time = 0;
            estimator.IMU_bias_sum = 0;
            estimator.processImage(image, img_msg->header,sum_reproject,current_rej,delta_time);
           // printf("******mamba out is:%f\n",tmp_gyc_bias);
            estimator.num_bias++;
            //printf("dt is:%f\n",tmp_gyc_bias);

/*              ofstream loop_path_file("/home/wjb/output/reprejectrion.txt", ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(9);
            loop_path_file << current_time << " ";
            loop_path_file.precision(5);
            loop_path_file <<estimator.IMU_bias_sum<< endl;
            loop_path_file.close();
            double correction_time  = current_time+0.1*delta_time;

            ofstream loop_path_file_("/home/wjb/output/VW.txt", ios::app);
            loop_path_file_.setf(ios::fixed, ios::floatfield);
            loop_path_file_.precision(9);
            loop_path_file_ << correction_time << " ";
            loop_path_file_.precision(5);
            loop_path_file_ << estimator.Vs[WINDOW_SIZE].squaredNorm() << " ";
            loop_path_file_.precision(5);
            loop_path_file_ << tmp_gyc_bias << endl;
            loop_path_file_.close();   */
            
            //利用重投影误差和速度的加权平均作为分段依据
            //if(!spline_adjust_){
                
                sum_V+=(estimator.Vs[WINDOW_SIZE]).squaredNorm();
                mid_sum_reproject = sum_reproject;
            //}
                if(!spline_adjust_){
                   // num_reproject_test++;
                    sum_V_test+=(estimator.Vs[WINDOW_SIZE]).squaredNorm();
                }
            
            if(tmp_acc_bias>50){
                sum_acc+=tmp_acc_bias;
                num_reproject_test++;
            } 

            if(tmp_gyc_bias>0.02){
                sum_gyc+=tmp_gyc_bias;
                num_reproject++;
            }
            

            estimator.theta_gyc = 0.8*estimator.theta_gyc+0.2*tmp_gyc_bias;

            double test_sum_reprj = num_reproject*abs(current_rej-mid_sum_reproject/num_reproject)/mid_sum_reproject;
            double test_sum_V = num_reproject*abs((estimator.Vs[WINDOW_SIZE]).squaredNorm()-sum_V/num_reproject)/sum_V;
            double test_sum_gyc = num_reproject*abs(tmp_gyc_bias-sum_gyc/num_reproject)/sum_gyc;
            double test_sum_IMU = num_reproject*abs(estimator.current_bias-estimator.IMU_bias_sum/num_reproject)/estimator.IMU_bias_sum;


            
            //printf("outside V is:%f\n",tmp_acc_bias/sum_acc*num_reproject);
            //if((0.2*test_sum_V+0.8*test_sum_reprj)>0.5)
            //printf("man is:%f\n",tmp_gyc_bias);
            if(tmp_gyc_bias>sum_gyc/num_reproject||tmp_acc_bias >sum_acc/num_reproject_test)
            //if(tmp_gyc_bias>2*sum_gyc/num_reproject ||tmp_acc_bias >200)
            //if(tmp_gyc_bias>sum_gyc/num_reproject)
            {
                //printf("\ninside V is:%f\n",tmp_gyc_bias);
                //printf("test V is:%f\n",estimator.current_bias);
                //num_reproject = 0;
                //sum_reproject = 0;
                sum_V =0;
                spline_adjust_ = true;
                //buffer_sount = 2;
            }
            else
            {
              //printf("%facc bias is:%f\n",current_time,tmp_acc_bias);
                spline_adjust_ = false;
            } 
/* 
             if((0.5*test_sum_V+0.5*test_sum_gyc)>0.4){
                if(!spline_adjust_){
                    num_reproject = 0;
                    sum_gyc = 0;
                    sum_V = 0;
                }
                buffer_sount = 3;

                spline_adjust_ = true;
                
            }
            else{
                spline_adjust_ = false;
            }          
          //  printf("pose is:%f,%f,%f\n",pose_camera.x(),pose_camera.y(),pose_camera.z());
           // if(current_time>=1403715564.512143 && current_time <=1403715564.812143)
            //{
                //printf("%f a is:%f,%f,%f\n",current_time,tmp_acc.x(),tmp_acc.y(),tmp_acc.z());
                //printf("%f Spline -2 op is:%f,%f,%f\n",current_time,estimator.Vs[WINDOW_SIZE-2].x(),estimator.Vs[WINDOW_SIZE-2].y(),estimator.Vs[WINDOW_SIZE-2].z());
                //printf("%f Spline -1 op is:%f,%f,%f\n",current_time,estimator.Vs[WINDOW_SIZE-1].x(),estimator.Vs[WINDOW_SIZE-1].y(),estimator.Vs[WINDOW_SIZE-1].z());
               // printf("%f v is:%f,%f,%f\n",current_time,estimator.Vs[WINDOW_SIZE].x(),estimator.Vs[WINDOW_SIZE].y(),estimator.Vs[WINDOW_SIZE].z());
                //printf("%f Spline After r is:%f,%f,%f\n\n",current_time,estimator.Ps[WINDOW_SIZE].x(),estimator.Ps[WINDOW_SIZE].y(),estimator.Ps[WINDOW_SIZE].z());
               //printf("%f,front end current acc is:%f,%f,%f\n\n",current_time,tmp_acc.x(),tmp_acc.y(),tmp_acc.z());
                //printf("current acc is:%f\n",tmp_acc.squaredNorm());
            //}
                sum_of_acc++;
                sum_of_gy++;
                sum_of_acc_Norm += tmp_acc.squaredNorm();
                sum_of_gy_Norm += tmp_Gy.squaredNorm();

                double k = ((estimator.Vs[WINDOW_SIZE].cross(tmp_acc).norm())/(estimator.Vs[WINDOW_SIZE].norm()*estimator.Vs[WINDOW_SIZE].norm()*estimator.Vs[WINDOW_SIZE].norm()));
                //printf("k is:%f\n\n",k);
                   // if(k>5||tmp_acc.squaredNorm()>(sum_of_acc_Norm/sum_of_acc)){
                    if(k>1.5){
                    //if(tmp_Gy.squaredNorm()>(sum_of_gy_Norm/sum_of_gy)){
                    //sum_of_all++;
                    spline_acc_adjust = false;
                    sum_of_acc_Norm = 0;
                    sum_of_acc = 0;
                    sum_of_gy_Norm = 0;
                    sum_of_gy = 0;
                    if(!start_lable_){
                        start_pose_test = estimator.Ps[WINDOW_SIZE];
                        start_velocity_test = estimator.Vs[WINDOW_SIZE];
                        start_time = current_time;
                        start_lable_ = true;
                    }
                    else{
                        end_pose_test = estimator.Ps[WINDOW_SIZE];
                        end_velocity_test = estimator.Vs[WINDOW_SIZE];
                        end_time = current_time;
                        start_pose_test = end_pose_test;
                        start_velocity_test = end_velocity_test;
                        start_time = end_time;
                    }
                    //printf("acc count is:%f\n\n",sum_of_acc_Norm/sum_of_acc);
                    sum_of_acc_Norm = 0;
                    sum_of_acc =0; 
                }
                else{
                    spline_acc_adjust = true;
                }                 
            //}
            
            //acc_adjustment_2
/*             if(Head_pose.isZero()){
                Head_pose = estimator.Ps[WINDOW_SIZE];
                Head_time_ = current_time;
            }
            else if(Mid_pose.isZero()){
                Mid_pose = estimator.Ps[WINDOW_SIZE];
                Mid_time_ = current_time;
            }
            else if(Tail_pose.isZero()){
                Tail_pose = estimator.Ps[WINDOW_SIZE];
                Tail_time_ = current_time;
            }
            else{
                if(acc_mid_change){
                    Head_pose = Mid_pose;
                    Mid_pose = Tail_pose;
                    Head_time_ = Mid_time_;
                    Mid_time_ = Tail_time_;
                }
                Tail_pose = estimator.Ps[WINDOW_SIZE];
                Tail_time_ = current_time;
                //printf("Head pose is:%f,%f,%f\n",Head_pose.x(),Head_pose.y(),Head_pose.z());
                //printf("Mid pose is:%f,%f,%f\n",Head_pose.x(),Head_pose.y(),Head_pose.z());
                //printf("Tail pose is:%f,%f,%f\n",Head_pose.x(),Head_pose.y(),Head_pose.z());
                //printf("Head_Time is:%f\n",Head_time_);
                //printf("Mid_Time is:%f\n",Mid_time_);
                //printf("Tail_Time is:%f\n\n",Tail_time_);

                Vector3d vector_front = Head_pose - Mid_pose;
                Vector3d vector_back =  Mid_pose - Tail_pose;
                double alpha = acos( vector_front.dot(vector_back) / ( (vector_front.norm())*(vector_back.norm()) ) );
                //printf("alpha is %f",alpha);
                if(alpha<0.5){
                    acc_mid_change = false;
                }
                else{
                    acc_mid_change = true;
                }
            } */

            //JB add
            
            velocity_buffer_tmp.push_back(estimator.Vs[WINDOW_SIZE]);
            Rotation_buffer_tmp.push_back(estimator.Rs[WINDOW_SIZE]);
            if(velocity_buffer_tmp.size()>=10){
              //  double current_Velocity = velocity_buffer_tmp[velocity_buffer_tmp.size()-1].squaredNorm();
               // double history_Velocity = velocity_buffer_tmp[estimator.buff_count].squaredNorm();

                Matrix3d current_Rotation = Rotation_buffer_tmp[Rotation_buffer_tmp.size()-1];
                Matrix3d history_Rotation = Rotation_buffer_tmp[estimator.buff_count];
                Matrix3d diff_Rotation = current_Rotation.transpose()*history_Rotation;

                double trace = diff_Rotation.trace();
                //double theta = std::acos((trace-1)/2);

                Eigen::Quaterniond q(diff_Rotation);

                Eigen::Vector3d rotation_vector = q.normalized().vec();

               // Vector3d omega_c =SO3_::anti_hat( SO3_::log(current_Rotation));
               // Vector3d omega_h = SO3_::anti_hat( SO3_::log(history_Rotation));

                //double diff = (omega_c - omega_h).norm();
                double diff = rotation_vector.norm();
                //printf("diff is: %f\n theta is :%f\n\n",diff,theta);

               // Vector3d current_Velocity_ = velocity_buffer_tmp[velocity_buffer_tmp.size()-1];
               // Vector3d history_Velocity_ = velocity_buffer_tmp[velocity_buffer_tmp.size()-2];
                
                //double velocity_angle = acos( history_Velocity_.dot(current_Velocity_) / ( (history_Velocity_.norm())*(current_Velocity_.norm() )));

                //double diff = current_Velocity - history_Velocity;
                //printf("diff_R is :%f\n\n",diff);
                
                /* bool v_adj=true;
                if(abs(current_Velocity_.x()>0.5*( abs(history_Velocity_.x()-current_Velocity_.x()) ))||
                abs(current_Velocity_.y()>0.5*( abs(history_Velocity_.y()-current_Velocity_.y()) ))||
                abs(current_Velocity_.z()>0.5*( abs(history_Velocity_.z()-current_Velocity_.z()) ))){
                    v_adj = false;
                }; */
                //printf("chahzi is:%f\n",chazhi_sum/sum_of_all);            
                //if(abs(velocity_angle)>0.3)
                
                /*
                if(abs(diff)>0.02||abs(velocity_angle)>0.01)
                {   
                    //printf("current time is: %f\n",current_time);
                    estimator.Velocity_buff.push_back(estimator.Vs[WINDOW_SIZE]);
                    estimator.Sping_lable_estimatior_ = false;
                    estimator.buff_count=velocity_buffer_tmp.size()-1;
                    estimator.spline_count = 0;
                    pose_op_lable = true;
                }
                */
                
               // if(abs(diff)>0.05)
              // if(tmp_acc_bias>50&&tmp_gyc_bias>0.15)
                //if(abs(diff)>0.1||std::isnan(diff)||estimator.spline_count>5)
                //if(!v_adj)
               //if((!spline_acc_adjust))
            if(spline_adjust_||buffer_sount > 0)    
               //if(estimator.spline_count>4)
                //if(estimator.IMU_bias_sum>10)
                {   
                    //printf("bias is:%f\n",tmp_acc_bias);
                    estimator.buff_count=velocity_buffer_tmp.size()-1;
                    estimator.Velocity_buff.push_back(estimator.Vs[WINDOW_SIZE]);
                    estimator.Sping_lable_estimatior_ = false;
                    estimator.spline_count = 0; 
                    pose_op_lable = true;
                    buffer_sount --;
                    
                   // printf("\n%fsuper pose is:%f,%f,%f\n",current_time,estimator.Ps[WINDOW_SIZE].x(),estimator.Ps[WINDOW_SIZE].y(),estimator.Ps[WINDOW_SIZE].z());
                   // printf("%fsuper velocity is:%f,%f,%f\n\n",current_time,estimator.Vs[WINDOW_SIZE].x(),estimator.Vs[WINDOW_SIZE].y(),estimator.Vs[WINDOW_SIZE].z());
                } 
                else
                {   
                    
                    //此处改为false就是不进行分段的意思，分段就改成true
                    estimator.Sping_lable_estimatior_ = true;
                    estimator.spline_count++;
/*                     if(pose_op_lable){
                       estimator.Sping_lable_estimatior_ = false;
                    } */
                    pose_op_lable = true;
                 //   printf("%fmid pose is:%f,%f,%f\n",current_time,estimator.Ps[WINDOW_SIZE].x(),estimator.Ps[WINDOW_SIZE].y(),estimator.Ps[WINDOW_SIZE].z());
                 //   printf("%fmid velocity is:%f,%f,%f\n",current_time,estimator.Vs[WINDOW_SIZE].x(),estimator.Vs[WINDOW_SIZE].y(),estimator.Vs[WINDOW_SIZE].z()); 
                }
            }
            //end


            //视觉数据包含BA运算,如果变化不大,便可以跳过视觉优化（BA）
/*            if(Spline.compute_spline_lable(tmp_V,tmp_Gy)&&estimator.frame_count==10)
            {
                if(estimator.Spline_segment_lable_ ==true)
                {   
                    estimator.Spline_Tail_Position_ = estimator.Ps[WINDOW_SIZE];
                    estimator.Spline_Tail_time_ = current_time;
                    for(int i = 0;i<estimator.Spline_mid_Position_.size();i++)
                    {   
                        estimator.processSpline(estimator.Spline_Head_Position_,estimator.Spline_Tail_Position_,
                                                estimator.Spline_Head_Velocity_,estimator.Spline_Tail_Velocity_,
                                                estimator.Spline_Head_time_,estimator.Spline_Tail_time_,estimator.Spline_mid_time_[i]);
                   // printf("not spline is:%F,%f,%f\n\n",estimator.Spline_mid_Position_[i].x(),estimator.Spline_mid_Position_[i].y(),estimator.Spline_mid_Position_[i].z());
                    }
                    estimator.Spline_mid_Position_.clear();
                    estimator.Spline_mid_time_.clear();
                    estimator.Spline_mid_header_.clear();
                    estimator.Spline_mid_Velocity_.clear();
                    estimator.Spline_segment_lable_ = false;//结束优化
                } 
            } */
            
            //输出信息
            //if(Spline.compute_spline_lable(tmp_V)&&estimator.frame_count==10&&estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
              // printf("R is:\n%f,%f,%f,%f\n",
               // estimator.pre_integrations[9]->delta_q.x(),estimator.pre_integrations[9]->delta_q.y(),
               // estimator.pre_integrations[9]->delta_q.z(),estimator.pre_integrations[9]->delta_q.w());

           // double whole_t = t_s.toc();
            //printStatistics(estimator, whole_t);
            std_msgs::Header header = img_msg->header;
            header.frame_id = "world";

           // printf("pose camara is:%f,%f,%f\n",pose_camera.x(),pose_camera.y(),pose_camera.z());
           // printf(" Ps is:%f,%f,%f\n\n",estimator.Ps[WINDOW_SIZE].x(),estimator.Ps[WINDOW_SIZE].y(),estimator.Ps[WINDOW_SIZE].z());
                    pubOdometry(estimator, header);
                    pubKeyPoses(estimator, header);
                    pubCameraPose(estimator, header);
                    pubPointCloud(estimator, header);
                    pubTF(estimator, header);
                    pubKeyframe(estimator,tmp_acc,tmp_Gy);

            if (relo_msg != NULL)
                pubRelocalization(estimator);
            //ROS_ERROR("end: %f, at %f", img_msg->header.stamp.toSec(), ros::Time::now().toSec());
        }
        m_estimator.unlock();
        m_buf.lock();
        m_state.lock();
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            update();
        m_state.unlock();
        m_buf.unlock();
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "vins_estimator");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n);
    estimator.setParameter();
#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif
    ROS_WARN("waiting for image and imu...");

    registerPub(n);

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_image = n.subscribe("/feature_tracker/feature", 2000, feature_callback);
    ros::Subscriber sub_restart = n.subscribe("/feature_tracker/restart", 2000, restart_callback);
    ros::Subscriber sub_relo_points = n.subscribe("/pose_graph/match_points", 2000, relocalization_callback);

    std::thread measurement_process{process};
    ros::spin();

    return 0;
}
