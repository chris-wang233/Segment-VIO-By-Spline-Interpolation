#include "./pose_graph.h"
// #include "../../vins_estimator/src/utility/visualization.h"

double op_time_sum = 0;
int loop_times = 0;
double ceres_op_time_sum = 0;
double LERP_op_time_sum = 0;
double SE_op_time_sum = 0;
int path_change = 2;
double time_stap_path = -1;

double total_ceres_t_sum = 0;
double total_lerp_t_sum = 0;

int V_sum =0;
int noV_sum = 0;
std::vector<geometry_msgs::Point32> p_v_SEspline;
sensor_msgs::PointCloud point_cloud_SE;
PoseGraph::PoseGraph()
{
    posegraph_visualization = new CameraPoseVisualization(1.0, 0.0, 1.0, 1.0);
    posegraph_visualization->setScale(0.1);
    posegraph_visualization->setLineWidth(0.01);
    t_optimization = std::thread(&PoseGraph::optimize4DoF, this);
    // t_part_optimization = std::thread(&PoseGraph::optimize_part, this);
    earliest_loop_index = -1;
    t_drift = Eigen::Vector3d(0, 0, 0);
    yaw_drift = 0;
    r_drift = Eigen::Matrix3d::Identity();
    w_t_vio = Eigen::Vector3d(0, 0, 0);
    w_r_vio = Eigen::Matrix3d::Identity();
    global_index = 0;
    sequence_cnt = 0;
    sequence_loop.push_back(0);
    base_sequence = 1;
}

PoseGraph::~PoseGraph()
{
    t_optimization.join();
    // t_part_optimization.join();
}

void PoseGraph::registerPub(ros::NodeHandle &n)
{
    pub_point_cloud_SE = n.advertise<sensor_msgs::PointCloud>("point_cloud_SE", 1000);
    pub_pg_path = n.advertise<nav_msgs::Path>("pose_graph_path", 1000);
    pub_base_path = n.advertise<nav_msgs::Path>("base_path", 1000);
    pub_pose_graph = n.advertise<visualization_msgs::MarkerArray>("pose_graph", 1000);
    for (int i = 1; i < 10; i++)
        pub_path[i] = n.advertise<nav_msgs::Path>("path_" + to_string(i), 1000);
}

void PoseGraph::pubPointCloud_backend(const Vector3d pose)
{

    // point_cloud_SE.header = header;
    geometry_msgs::Point32 p_SEspline;
    point_cloud_SE.header.frame_id = "world";
    p_SEspline.x = pose.x();
    p_SEspline.y = pose.y();
    p_SEspline.z = pose.z();

    p_v_SEspline.push_back(p_SEspline);
    for (int i = 0; i < p_v_SEspline.size(); i++)
    {
        point_cloud_SE.points.push_back(p_v_SEspline[i]);
    }
    pub_point_cloud_SE.publish(point_cloud_SE);
}

void PoseGraph::loadVocabulary(std::string voc_path)
{
    voc = new BriefVocabulary(voc_path);
    db.setVocabulary(*voc, false, 0);
}

void PoseGraph::addKeyFrame(KeyFrame *cur_kf, bool flag_detect_loop)
{
    // shift to base frame
    Vector3d vio_P_cur;
    Matrix3d vio_R_cur;
    if (sequence_cnt != cur_kf->sequence)
    {
        sequence_cnt++;
        sequence_loop.push_back(0);
        w_t_vio = Eigen::Vector3d(0, 0, 0);
        w_r_vio = Eigen::Matrix3d::Identity();
        m_drift.lock();
        t_drift = Eigen::Vector3d(0, 0, 0);
        r_drift = Eigen::Matrix3d::Identity();
        m_drift.unlock();
    }

    cur_kf->getVioPose(vio_P_cur, vio_R_cur);
    vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
    vio_R_cur = w_r_vio * vio_R_cur;
    cur_kf->updateVioPose(vio_P_cur, vio_R_cur);
    cur_kf->index = global_index;
    global_index++;
    int loop_index = -1;
    if (flag_detect_loop)
    {
        TicToc tmp_t;
        loop_index = detectLoop(cur_kf, cur_kf->index);
    }
    else
    {
        addKeyFrameIntoVoc(cur_kf);
    }
    if (loop_index != -1) // 如果发现产生回环
    {
        // printf(" %d detect loop with %d \n", cur_kf->index, loop_index);
        KeyFrame *old_kf = getKeyFrame(loop_index);

        if (cur_kf->findConnection(old_kf))
        {
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
                earliest_loop_index = loop_index;
            earliest_loop_index = -1; // JBJB 这个地方有两处需要修改，这里仅对一处进行批注

            Vector3d w_P_old, w_P_cur, vio_P_cur;
            Matrix3d w_R_old, w_R_cur, vio_R_cur;
            old_kf->getVioPose(w_P_old, w_R_old);
            cur_kf->getVioPose(vio_P_cur, vio_R_cur);

            Vector3d relative_t;
            Quaterniond relative_q;
            relative_t = cur_kf->getLoopRelativeT();
            relative_q = (cur_kf->getLoopRelativeQ()).toRotationMatrix();
            w_P_cur = w_R_old * relative_t + w_P_old;
            w_R_cur = w_R_old * relative_q;
            double shift_yaw;
            Matrix3d shift_r;
            Vector3d shift_t;
            shift_yaw = Utility::R2ypr(w_R_cur).x() - Utility::R2ypr(vio_R_cur).x();
            shift_r = Utility::ypr2R(Vector3d(shift_yaw, 0, 0));
            shift_t = w_P_cur - w_R_cur * vio_R_cur.transpose() * vio_P_cur;
            // shift vio pose of whole sequence to the world frame
            if (old_kf->sequence != cur_kf->sequence && sequence_loop[cur_kf->sequence] == 0)
            {
                w_r_vio = shift_r;
                w_t_vio = shift_t;
                vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
                vio_R_cur = w_r_vio * vio_R_cur;
                cur_kf->updateVioPose(vio_P_cur, vio_R_cur);
                list<KeyFrame *>::iterator it = keyframelist.begin();
                for (; it != keyframelist.end(); it++)
                {
                    if ((*it)->sequence == cur_kf->sequence)
                    {
                        Vector3d vio_P_cur;
                        Matrix3d vio_R_cur;
                        (*it)->getVioPose(vio_P_cur, vio_R_cur);
                        vio_P_cur = w_r_vio * vio_P_cur + w_t_vio;
                        vio_R_cur = w_r_vio * vio_R_cur;
                        (*it)->updateVioPose(vio_P_cur, vio_R_cur);
                    }
                }
                sequence_loop[cur_kf->sequence] = 1;
            }
            m_optimize_buf.lock();
            optimize_buf.push(cur_kf->index);
            m_optimize_buf.unlock();
        }
    }
    m_keyframelist.lock();
    Vector3d P;
    Matrix3d R;
    cur_kf->getVioPose(P, R);
    P = r_drift * P + t_drift;
    R = r_drift * R;
    cur_kf->updatePose(P, R);
    Quaterniond Q{R};
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.stamp = ros::Time(cur_kf->time_stamp);
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
    pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
    pose_stamped.pose.position.z = P.z();
    pose_stamped.pose.orientation.x = Q.x();
    pose_stamped.pose.orientation.y = Q.y();
    pose_stamped.pose.orientation.z = Q.z();
    pose_stamped.pose.orientation.w = Q.w();
    path[sequence_cnt].poses.push_back(pose_stamped);
    path[sequence_cnt].header = pose_stamped.header;
    /*######################################################################################
        //orientation
        if (SAVE_LOOP_PATH)
        {
            ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(0);
            loop_path_file << cur_kf->time_stamp * 1e9 << ",";
            loop_path_file.precision(5);
            loop_path_file  << P.x() << ","
                  << P.y() << ","
                  << P.z() << ","
                  << Q.w() << ","
                  << Q.x() << ","
                  << Q.y() << ","
                  << Q.z() << ","
                  << endl;
            loop_path_file.close();
        }
     ########################################################################################   */
    // euroc
    if (SAVE_LOOP_PATH)
    {
        ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
        loop_path_file.setf(ios::fixed, ios::floatfield);
        loop_path_file.precision(9);
        loop_path_file << cur_kf->time_stamp << " ";
        loop_path_file.precision(5);
        loop_path_file << P.x() << " "
                       << P.y() << " "
                       << P.z() << " "
                       << Q.x() << " "
                       << Q.y() << " "
                       << Q.z() << " "
                       << Q.w() << endl;
        loop_path_file.close();
    }
    // TUM
    /*######################################################################################
    if (SAVE_LOOP_PATH)
   {
           ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
           loop_path_file.setf(ios::fixed, ios::floatfield);
           loop_path_file.precision(9);
           loop_path_file << cur_kf->time_stamp << " ";
           loop_path_file.precision(5);
           loop_path_file  << P.x() << " "
                           << P.y() << " "
                           << P.z() << " "
                           << Q.x() << " "
                           << Q.y() << " "
                           << Q.z() << " "
                           << Q.w() << endl;
           loop_path_file.close();
    }
   ########################################################################################   */
    // draw local connection
    if (SHOW_S_EDGE)
    {
        list<KeyFrame *>::reverse_iterator rit = keyframelist.rbegin();
        for (int i = 0; i < 4; i++)
        {
            if (rit == keyframelist.rend())
                break;
            Vector3d conncected_P;
            Matrix3d connected_R;
            if ((*rit)->sequence == cur_kf->sequence)
            {
                (*rit)->getPose(conncected_P, connected_R);
                posegraph_visualization->add_edge(P, conncected_P);
            }
            rit++;
        }
    }
    if (SHOW_L_EDGE)
    {
        if (cur_kf->has_loop)
        {
            // printf("has loop \n");
            KeyFrame *connected_KF = getKeyFrame(cur_kf->loop_index);
            Vector3d connected_P, P0;
            Matrix3d connected_R, R0;
            connected_KF->getPose(connected_P, connected_R);
            // cur_kf->getVioPose(P0, R0);
            cur_kf->getPose(P0, R0);
            if (cur_kf->sequence > 0)
            {
                // printf("add loop into visual \n");
                posegraph_visualization->add_loopedge(P0, connected_P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0));
            }
        }
    }
    // posegraph_visualization->add_pose(P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0), Q);

    keyframelist.push_back(cur_kf);
    KeyFrame_buf.push(cur_kf->index); // add
    publish();
    m_keyframelist.unlock();
}

void PoseGraph::loadKeyFrame(KeyFrame *cur_kf, bool flag_detect_loop)
{
    cur_kf->index = global_index;
    global_index++;
    int loop_index = -1;
    if (flag_detect_loop)
        loop_index = detectLoop(cur_kf, cur_kf->index);
    else
    {
        addKeyFrameIntoVoc(cur_kf);
    }
    if (loop_index != -1)
    {
        printf(" %d detect loop with %d \n", cur_kf->index, loop_index);
        KeyFrame *old_kf = getKeyFrame(loop_index);
        if (cur_kf->findConnection(old_kf))
        {
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
                earliest_loop_index = loop_index;
            m_optimize_buf.lock();
            optimize_buf.push(cur_kf->index);
            m_optimize_buf.unlock();
        }
    }
    m_keyframelist.lock();
    Vector3d P;
    Matrix3d R;
    cur_kf->getPose(P, R);
    Quaterniond Q{R};
    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.stamp = ros::Time(cur_kf->time_stamp);
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
    pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
    pose_stamped.pose.position.z = P.z();
    pose_stamped.pose.orientation.x = Q.x();
    pose_stamped.pose.orientation.y = Q.y();
    pose_stamped.pose.orientation.z = Q.z();
    pose_stamped.pose.orientation.w = Q.w();
    base_path.poses.push_back(pose_stamped);
    base_path.header = pose_stamped.header;

    // draw local connection
    if (SHOW_S_EDGE)
    {
        list<KeyFrame *>::reverse_iterator rit = keyframelist.rbegin();
        for (int i = 0; i < 1; i++)
        {
            if (rit == keyframelist.rend())
                break;
            Vector3d conncected_P;
            Matrix3d connected_R;
            if ((*rit)->sequence == cur_kf->sequence)
            {
                (*rit)->getPose(conncected_P, connected_R);
                posegraph_visualization->add_edge(P, conncected_P);
            }
            rit++;
        }
    }
    /*
    if (cur_kf->has_loop)
    {
        KeyFrame* connected_KF = getKeyFrame(cur_kf->loop_index);
        Vector3d connected_P;
        Matrix3d connected_R;
        connected_KF->getPose(connected_P,  connected_R);
        posegraph_visualization->add_loopedge(P, connected_P, SHIFT);
    }
    */

    keyframelist.push_back(cur_kf);
    // publish();
    m_keyframelist.unlock();
}

KeyFrame *PoseGraph::getKeyFrame(int index)
{
    //    unique_lock<mutex> lock(m_keyframelist);
    list<KeyFrame *>::iterator it = keyframelist.begin();
    for (; it != keyframelist.end(); it++)
    {
        if ((*it)->index == index)
            break;
    }
    if (it != keyframelist.end())
        return *it;
    else
        return NULL;
}

int PoseGraph::detectLoop(KeyFrame *keyframe, int frame_index)
{
    // put image into image_pool; for visualization
    cv::Mat compressed_image;
    if (DEBUG_IMAGE)
    {
        int feature_num = keyframe->keypoints.size();
        cv::resize(keyframe->image, compressed_image, cv::Size(376, 240));
        putText(compressed_image, "feature_num:" + to_string(feature_num), cv::Point2f(10, 10), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
        image_pool[frame_index] = compressed_image;
    }
    TicToc tmp_t;
    // first query; then add this frame into database!
    QueryResults ret;
    TicToc t_query;
    db.query(keyframe->brief_descriptors, ret, 4, frame_index - 50);
    // printf("query time: %f", t_query.toc());
    // cout << "Searching for Image " << frame_index << ". " << ret << endl;

    TicToc t_add;
    db.add(keyframe->brief_descriptors);
    // printf("add feature time: %f", t_add.toc());
    //  ret[0] is the nearest neighbour's score. threshold change with neighour score
    bool find_loop = false;
    cv::Mat loop_result;
    if (DEBUG_IMAGE)
    {
        loop_result = compressed_image.clone();
        if (ret.size() > 0)
            putText(loop_result, "neighbour score:" + to_string(ret[0].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255));
    }
    // visual loop result
    if (DEBUG_IMAGE)
    {
        for (unsigned int i = 0; i < ret.size(); i++)
        {
            int tmp_index = ret[i].Id;
            auto it = image_pool.find(tmp_index);
            cv::Mat tmp_image = (it->second).clone();
            putText(tmp_image, "index:  " + to_string(tmp_index) + "loop score:" + to_string(ret[i].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255));
            cv::hconcat(loop_result, tmp_image, loop_result);
        }
    }
    // a good match with its nerghbour
    if (ret.size() >= 1 && ret[0].Score > 0.05)
        for (unsigned int i = 1; i < ret.size(); i++)
        {
            // if (ret[i].Score > ret[0].Score * 0.3)
            if (ret[i].Score > 0.015)
            {
                find_loop = true;
                int tmp_index = ret[i].Id;
                if (DEBUG_IMAGE && 0)
                {
                    auto it = image_pool.find(tmp_index);
                    cv::Mat tmp_image = (it->second).clone();
                    putText(tmp_image, "loop score:" + to_string(ret[i].Score), cv::Point2f(10, 50), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
                    cv::hconcat(loop_result, tmp_image, loop_result);
                }
            }
        }
    /*
        if (DEBUG_IMAGE)
        {
            cv::imshow("loop_result", loop_result);
            cv::waitKey(20);
        }
    */
    if (find_loop && frame_index > 50)
    {
        int min_index = -1;
        for (unsigned int i = 0; i < ret.size(); i++)
        {
            if (min_index == -1 || (ret[i].Id < min_index && ret[i].Score > 0.015))
                min_index = ret[i].Id;
        }
        return min_index;
    }
    else
        return -1;
}

void PoseGraph::addKeyFrameIntoVoc(KeyFrame *keyframe)
{
    // put image into image_pool; for visualization
    cv::Mat compressed_image;
    if (DEBUG_IMAGE)
    {
        int feature_num = keyframe->keypoints.size();
        cv::resize(keyframe->image, compressed_image, cv::Size(376, 240));
        putText(compressed_image, "feature_num:" + to_string(feature_num), cv::Point2f(10, 10), CV_FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(255));
        image_pool[keyframe->index] = compressed_image;
    }

    db.add(keyframe->brief_descriptors);
}


// JB_会一直在后台运行
int test_sum = 0;
void PoseGraph::optimize4DoF()
{
    // int current_keyframe=0;
    // list<KeyFrame*>::iterator current_now=keyframelist.begin();
    // list<KeyFrame*>::iterator current_old=keyframelist.begin();

    while (true)
    {
        int cur_index = -1;
        int first_looped_index = -1;
        // int part_op_index = -1;//add

        m_optimize_buf.lock();
        while (!optimize_buf.empty())
        {
            cur_index = optimize_buf.front();
            first_looped_index = earliest_loop_index;
            optimize_buf.pop();
        }
        m_optimize_buf.unlock();
        if (cur_index != -1)
        {
            printf("optimize pose graph \n");

           // TicToc tmp_t;
            m_keyframelist.lock();
            KeyFrame *cur_kf = getKeyFrame(cur_index);

            int max_length = cur_index + 1;

            // w^t_i   w^q_i
            double t_array[max_length][3];
            double t_array_lerp[max_length][3];
            Quaterniond q_array[max_length];
            double euler_array[max_length][3];
            double sequence_array[max_length];

            ceres::Problem problem;
            ceres::Solver::Options options;
            options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
            // options.minimizer_progress_to_stdout = true;
            // options.max_solver_time_in_seconds = SOLVER_TIME * 3;
            options.max_num_iterations = 5;
            ceres::Solver::Summary summary;
            ceres::LossFunction *loss_function;
            loss_function = new ceres::HuberLoss(0.1);
            // loss_function = new ceres::CauchyLoss(1.0);
            ceres::LocalParameterization *angle_local_parameterization =
                AngleLocalParameterization::Create();
            list<KeyFrame *>::iterator it;
            int skip_optimistic_count = 0;
            int i = 0;

            // count 记录了对应关键帧的编号，两个lable记录spline开始与结束的关键帧编号
            int test_count = 0, start_lable = -1, end_Lable = -1;
            vector<KeyFrame *> keyframe_buffer;

            int ceres_num = 0;
            int spline_num = 0;

            for (it = keyframelist.begin(); it != keyframelist.end(); it++) // 全局优化？？？
            {
                std::string Spline_lable((*it)->Spline_frame_);
                if ((*it)->index < first_looped_index)
                {
                    continue;
                }

               // if (Spline_lable != "PoseGraph" && !(*it)->has_loop)
               if (Spline_lable != "PoseGraph" && !(*it)->has_loop)
                {
                    spline_num++;
                    if (start_lable < 0)
                        start_lable = test_count - 1;
                    (*it)->local_index = i;
                    Quaterniond tmp_q;
                    Matrix3d tmp_r;
                    Vector3d tmp_t;
                    (*it)->getVioPose(tmp_t, tmp_r);
                    tmp_q = tmp_r;

                    t_array[i][0] += tmp_t(0);
                    t_array[i][1] += tmp_t(1);
                    t_array[i][2] += tmp_t(2);

                    q_array[i] *= tmp_q;

                    Vector3d euler_angle = Utility::R2ypr(tmp_q.toRotationMatrix());
                    euler_array[i][0] += euler_angle.x();
                    euler_array[i][1] += euler_angle.y();
                    euler_array[i][2] += euler_angle.z();

                    sequence_array[i] = (*it)->sequence;

                    continue;
                }
                ceres_num++;

                (*it)->local_index = i;
                Quaterniond tmp_q;
                Matrix3d tmp_r;
                Vector3d tmp_t;

                (*it)->getVioPose(tmp_t, tmp_r);
                tmp_q = tmp_r;

                t_array[i][0] = tmp_t(0);
                t_array[i][1] = tmp_t(1);
                t_array[i][2] = tmp_t(2);

                t_array_lerp[i][0] = tmp_t(0);
                t_array_lerp[i][1] = tmp_t(1);
                t_array_lerp[i][2] = tmp_t(2);

                q_array[i] = tmp_q;

                Vector3d euler_angle = Utility::R2ypr(tmp_q.toRotationMatrix());
                euler_array[i][0] = euler_angle.x();
                euler_array[i][1] = euler_angle.y();
                euler_array[i][2] = euler_angle.z();

                sequence_array[i] = (*it)->sequence;
                // if(Spline_lable == "Spline")
                //     printf("position is:%f,%f,%f\n",(*it)->point_3d.x(),(*it)->point_3d.y(),(*it)->point_3d.z());
                Vector3d vio_t;
                Matrix3d vio_r;
                (*it)->getVioPose(vio_t, vio_r);
                // if(Spline_lable == "Spline")
                // printf("position %d is:%f,%f,%f\n",p,vio_t.x(),vio_t.y(),vio_t.z());

                problem.AddParameterBlock(euler_array[i], 1, angle_local_parameterization);
                problem.AddParameterBlock(t_array[i], 3);

                if ((*it)->index == first_looped_index || (*it)->sequence == 0)
                {
                    problem.SetParameterBlockConstant(euler_array[i]);
                    problem.SetParameterBlockConstant(t_array[i]);
                }

                // add edge
                for (int j = 1; j < 5; j++)
                {
                    if (i - j >= 0 && sequence_array[i] == sequence_array[i-j])
                    {
                        Vector3d euler_conncected = Utility::R2ypr(q_array[i-j].toRotationMatrix());
                        Vector3d relative_t(t_array[i][0] - t_array[i-j][0], t_array[i][1] - t_array[i - j][1], t_array[i][2] - t_array[i-j][2]);
                        relative_t = q_array[i-j].inverse() * relative_t;
                        double relative_yaw = euler_array[i][0] - euler_array[i-j][0];
                        ceres::CostFunction *cost_function = FourDOFError::Create(relative_t.x(), relative_t.y(), relative_t.z(),
                                                                                  relative_yaw, euler_conncected.y(), euler_conncected.z());
                        problem.AddResidualBlock(cost_function, NULL, euler_array[i-j],
                                                 t_array[i-j],
                                                 euler_array[i],
                                                 t_array[i]);
                    }
                }
                // add loop edge

                if ((*it)->has_loop)
                {
                    assert((*it)->loop_index >= first_looped_index);
                    int connected_index = getKeyFrame((*it)->loop_index)->local_index;
                    Vector3d euler_conncected = Utility::R2ypr(q_array[connected_index].toRotationMatrix());
                    Vector3d relative_t;
                    relative_t = (*it)->getLoopRelativeT();//并不是接近0向量
                    double relative_yaw = (*it)->getLoopRelativeYaw();
                    ceres::CostFunction *cost_function = FourDOFWeightError::Create(relative_t.x(), relative_t.y(), relative_t.z(),
                                                                                    relative_yaw, euler_conncected.y(), euler_conncected.z());
                    problem.AddResidualBlock(cost_function, loss_function, euler_array[connected_index],
                                             t_array[connected_index],
                                             euler_array[i],
                                             t_array[i]);
                }

                if ((*it)->index == cur_index)
                    break;
                i++;
            }

            m_keyframelist.unlock();
            TicToc tmp_t;
            ceres::Solve(options, &problem, &summary);
            op_time_sum += tmp_t.toc();
            loop_times ++;
            //printf("Pose optimization time: %f \n", tmp_t.toc());
            printf("%d Pose optimization time: %f \n",loop_times,op_time_sum);
            printf("ceres num is:%d,spline num is:%d\n",ceres_num,spline_num);

            m_keyframelist.lock();
            i = 0;
            bool Spline_start_lable = false;
            bool Spline_start_head = false;
            // 更新之前的威姿 P3
            Vector3d no_temp_start_pose;
            Vector3d no_temp_end_pose;
            // 更新之前的威姿 R3
            Matrix3d no_temp_start_rotation;
            Matrix3d no_temp_end_rotation;
            // 更新之后的威姿 R3
            Vector3d temp_start_pose;
            Vector3d temp_end_pose;
            Matrix3d temp_start_rotation;
            Matrix3d temp_end_rotation;
            // 更新之后的威姿 SE3
            Sophus::SE3d se_start_pose;
            Sophus::SE3d se_end_pose;

            double temp_start_time;
            double temp_end_time;

            Vector3d temp_start_velocity;
            Vector3d temp_end_velocity;
            Vector3d temp_start_acc;
            Vector3d temp_end_acc;
            Vector3d temp_start_gyc;
            Vector3d temp_end_gyc;
            // 单数存开头，复数存结尾

            //vector<Matrix3d> temp_rotation_buffer;
            vector<Vector3d> temp_no_op_pose_buffer;
            vector<Matrix3d> temp_no_op_rotation_buffer;
            vector<Vector3d> temp_velocity_buffer;
            vector<Vector3d> temp_acc_buffer;
            temp_gyc_buffer.clear();
            temp_time_buffer.clear();
            temp_pose_buffer.clear();
            temp_rotation_buffer.clear();
            vector<Sophus::SE3d> temp_se_pose_buffer;
            vector<Sophus::SE3d> temp_se_after_pose_buffer;
            // Bspline
            vector<Vector3d> bspline_pose;
            
            TicToc firstLoop;
            //************************************第一次loop，记录pose和velocity等信息**************************************
            for (it = keyframelist.begin(); it != keyframelist.end(); it++)
            {
                if ((*it)->index < first_looped_index)
                    continue;
                std::string Spline_lable((*it)->Spline_frame_);

                if (Spline_lable != "PoseGraph" && !(*it)->has_loop)
                {
                    Matrix3d bspline_r;
                    Vector3d bspline_t;
                    (*it)->getVioPose(bspline_t, bspline_r);
                    bspline_pose.push_back(bspline_t);
                    Spline_start_lable = true;
                    continue;
                }

                Matrix3d no_op_r;
                Vector3d no_op_t;
                (*it)->getVioPose(no_op_t, no_op_r);
                Quaterniond tmp_q;
                tmp_q = Utility::ypr2R(Vector3d(euler_array[i][0], euler_array[i][1], euler_array[i][2]));
                Vector3d tmp_t = Vector3d(t_array[i][0], t_array[i][1], t_array[i][2]);
                Vector3d tmp_t_lerp = Vector3d(t_array_lerp[i][0], t_array_lerp[i][1], t_array_lerp[i][2]);
                Matrix3d tmp_r = tmp_q.toRotationMatrix();
                Vector3d delta_t = tmp_t - no_op_t;
                // printf("change is:%f,%f,%f\n\n",delta_t.x(),delta_t.y(),delta_t.z());
                (*it)->updatePose(tmp_t, tmp_r); //**************************************************************
                Quaterniond tmp_r_Q = tmp_q.normalized();

                // Matrix3d tmp_r_input = tmp_r_Q.toRotationMatrix();
                // add
                // if(Spline_lable == "Spline" || (*it)->has_loop){
                if (!Spline_start_lable)
                {
                    no_temp_start_pose = no_op_t;
                    no_temp_start_rotation = no_op_r;
                    temp_start_pose = tmp_t;
                    temp_start_rotation = tmp_r;
                    temp_start_time = (*it)->time_stamp;
                    Sophus::SE3d tmp_se_start_pose(tmp_r_Q, tmp_t);
                    se_start_pose = tmp_se_start_pose;
                    (*it)->getVelocity(temp_start_velocity);
                    (*it)->getIMU(temp_start_acc, temp_start_gyc);
                    if (Spline_start_head)
                    {
                        temp_se_after_pose_buffer.push_back(se_start_pose);
                    }
                }
                else
                {
                    no_temp_end_pose = no_op_t;
                    no_temp_end_rotation = no_op_r;
                    temp_end_pose = tmp_t;
                    temp_end_rotation = tmp_r;
                    temp_end_time = (*it)->time_stamp;
                    Sophus::SE3d tmp_se_end_pose(tmp_r_Q, tmp_t);
                    se_end_pose = tmp_se_end_pose;
                    (*it)->getVelocity(temp_end_velocity);
                    (*it)->getIMU(temp_end_acc, temp_end_gyc);
                    // 传入段首信息
                    temp_pose_buffer.push_back(temp_start_pose);
                    temp_rotation_buffer.push_back(temp_start_rotation);

                    temp_no_op_pose_buffer.push_back(no_temp_start_pose);
                    temp_no_op_rotation_buffer.push_back(no_temp_start_rotation);
                    
                    temp_velocity_buffer.push_back(temp_start_velocity);
                    temp_acc_buffer.push_back(temp_start_acc);
                    temp_gyc_buffer.push_back(temp_start_gyc);
                    temp_time_buffer.push_back(temp_start_time);
                    temp_se_pose_buffer.push_back(se_start_pose);
                    // 传入段尾信息
                    temp_pose_buffer.push_back(temp_end_pose);
                    temp_rotation_buffer.push_back(temp_end_rotation);
                    temp_no_op_pose_buffer.push_back(no_temp_end_pose);
                    temp_no_op_rotation_buffer.push_back(no_temp_end_rotation);
                    temp_velocity_buffer.push_back(temp_end_velocity);
                    temp_acc_buffer.push_back(temp_end_acc);
                    temp_gyc_buffer.push_back(temp_end_gyc);
                    temp_time_buffer.push_back(temp_end_time);
                    temp_se_pose_buffer.push_back(se_end_pose);
                    // 同步信息
                    temp_start_pose = temp_end_pose;
                    temp_start_rotation = temp_end_rotation;
                    no_temp_start_pose = no_temp_end_pose;
                    no_temp_start_rotation = no_temp_end_rotation;
                    temp_start_velocity = temp_end_velocity;
                    temp_start_acc = temp_end_acc;
                    temp_start_gyc = temp_end_gyc;
                    temp_start_time = temp_end_time;
                    se_start_pose = se_end_pose;
                    Spline_start_head = true;
                }
                Spline_start_lable = false;
                //}
                // end
                if ((*it)->index == cur_index)
                    break;
                i++;
            }
            printf("first loop time is:%f\n",firstLoop.toc());

            TicToc test_algs;
            // add loop again start
            list<KeyFrame *>::iterator it_;
            list<KeyFrame *>::iterator it_P;
            list<KeyFrame *>::iterator it_P_his;
            list<KeyFrame *>::iterator it_his;

            list<KeyFrame *>::iterator it_B;

            list<KeyFrame *>::iterator _it;
            list<KeyFrame *>::iterator P_it;
            segment_index = 0;
            int segment_index_inside = 0;
            int pose_count = 0;
            bool switch_lable = false; // 确定是否切换segement
            int SE_count_A = 0;
            int SE_count_B = 0;
            int R3_count = 0;
            int time_test_lable_ = 0;
            int time_test_lable = 0;
            vector<Eigen::Vector3d> spline_pose_;
            vector<Eigen::Vector3d> spline_velocity_;
            vector<Eigen::Vector3d> spline_omega_;
            vector<Eigen::Matrix3d> spline_Rotaion_;

            std::map<double,Eigen::Vector3d> pose_map;
            std::map<double,Eigen::Vector3d> velocity_map;
            std::map<double,Eigen::Vector3d> rotation_map;

            vector<Eigen::Vector3d> spline_acc_;
            vector<Eigen::Vector3d> spline_gyc_;
            vector<Eigen::Vector3d> _spline_pose;
            vector<Eigen::Vector3d> _spline_velocity;
            vector<double> spline_pose_time;
            vector<double> _spline_pose_time;

            Spline spline;

            
            
            SO3Spline SO3Spline;
            SE3spline SE3spline;

            Eigen::Vector3d start_pose_test;
            Vector3d end_pose_test;
            Vector3d start_velocity_test;
            Vector3d end_velocity_test;
            Vector3d start_acc_test;
            Vector3d end_acc_test;
            Vector3d velocity_coutinue;

            //*******************************************
            /* for (_it = keyframelist.begin(); _it != keyframelist.end(); _it++)
            {
                if ((*_it)->index < first_looped_index)
                    continue;
                std::string Spline_lable_((*_it)->Spline_frame_);
                Vector3d start_V;
                Vector3d end_V;
                if (Spline_lable_ != "PoseGraph" && !(*_it)->has_loop)
                {
                    if (!switch_lable)
                        P_it = _it;
                    // printf("start time is:%f,end time is:%f\n\n",temp_start_time,temp_end_time);

                    // 进入样条段内
                    switch_lable = true;
                    SESpline SEspline;

                    // Bspline testB;

                    // 未更新的差值帧位姿
                    Matrix3d tmp_r;
                    Vector3d tmp_t;

                    double current_time = (*_it)->time_stamp;
                    time_test_lable_++;

                    Vector3d current_V;
                    Vector3d current_Acc;
                    Vector3d current_Gyc;
                    (*_it)->getVelocity(current_V);
                    (*_it)->getIMU(current_Acc, current_Gyc);
                    double start_time = temp_time_buffer[2 * segment_index];
                    double end_time = temp_time_buffer[2 * segment_index + 1];
                    double time_persent = (current_time - start_time) / (end_time - start_time);

                    if (abs(time_persent) > 1)
                    {
                        continue;
                    }
                    // pose
                    Vector3d start_pose = temp_pose_buffer[2 * segment_index];
                    Vector3d end_pose = temp_pose_buffer[2 * segment_index + 1];
                    // rotation
                    Matrix3d start_rotation = temp_rotation_buffer[2 * segment_index];
                    Matrix3d end_rotation = temp_rotation_buffer[2 * segment_index + 1];

                    // no optimistic pose
                    Vector3d no_start_pose = temp_no_op_pose_buffer[2 * segment_index];
                    Vector3d no_end_pose = temp_no_op_pose_buffer[2 * segment_index + 1];
                    start_pose_test = no_start_pose;
                    end_pose_test = no_end_pose;
                    start_V = temp_no_op_pose_buffer[2 * segment_index];
                    end_V = temp_no_op_pose_buffer[2 * segment_index + 1];
                    // velocity
                    Vector3d start_velocity = temp_velocity_buffer[2 * segment_index];
                    Vector3d end_velocity = temp_velocity_buffer[2 * segment_index + 1];

                    start_velocity_test = start_velocity;
                    end_velocity_test = end_velocity;
                    // IMU
                    Vector3d start_acc_ = temp_acc_buffer[2 * segment_index];
                    Vector3d end_acc_ = temp_acc_buffer[2 * segment_index + 1];

                    start_acc_test = start_acc_;
                    end_acc_test = end_acc_;

                    Vector3d start_gyc_ = temp_gyc_buffer[2 * segment_index];
                    Vector3d end_gyc_ = temp_gyc_buffer[2 * segment_index + 1];

                    // 获取位置信息
                    (*_it)->getVioPose(tmp_t, tmp_r);
                    _spline_pose.push_back(tmp_t);
                    _spline_velocity.push_back(current_V);
                    _spline_pose_time.push_back(time_persent);

                    // ###################################################################################################
                    spline.compute_spline_argument(no_start_pose, no_end_pose, start_velocity, end_velocity);
                    segment_index_inside = segment_index;
                } */
              //  else if (switch_lable)
                //{
                    //****************************************************************
                    /*                     //ceres优化需要找18维数组指针，vector pose，vector timepersent
                                        Eigen::Vector3d P0_,P1_,P2_,P3_,P4_,P5_;
                                        spline.get_spline_para(P0_,P1_,P2_,P3_);
                                        // double splinetest[18]={P0_.x(),P0_.y(),P0_.z(),P1_.x(),P1_.y(),P1_.z(),P2_.x(),P2_.y(),P2_.z(),P3_.x(),P3_.y(),P3_.z(),
                                         //                       P4_.x(),P4_.y(),P4_.z(),P5_.x(),P5_.y(),P5_.z()};
                                        double splinetest[12]={P0_.x(),P0_.y(),P0_.z(),P1_.x(),P1_.y(),P1_.z(),P2_.x(),P2_.y(),P2_.z(),P3_.x(),P3_.y(),P3_.z()};
                                        SolveCeres solve_problem(end_pose_test,splinetest,_spline_pose,_spline_pose_time);
                                       // printf("test:%d,%d\n\n",spline_pose_.size(),spline_pose_time.size());
                                        //printf("time is:%f,%f,%f\n",spline_pose_time[0]);
                                        //printf("out spline P0 is:%f,%f.%f\n\n",splinetest[0],splinetest[1],splinetest[2]);
                                        //solve_problem.get_para()[2]);
                                       // printf("v0 befor is:%f,%f,%f\n",spline.C1.x(),spline.C1.y(),spline.C1.z());
                                       solve_problem.Solve();
                                       //if(spline_pose_.size()!=0)
                                            //solve_problem.Solve();
                                        //printf("origin p0 is:%f,%f,%f\n",spline.C1.x(),spline.C1.y(),spline.C1.z());
                                        //更新数据
                                        //spline.C0 = solve_problem.P0_op;
                                        spline.C1 = solve_problem.P1_op;
                                        spline.C2 = solve_problem.P2_op;
                                        spline.C3 = solve_problem.P3_op;
                                        //spline.C1.x()+=10;
                                        //spline.C2.y()+=10;
                                       //spline.C3.z()+=10;
                                        //spline.P0 = solve_problem.P0_op;
                                        //spline.P1 = solve_problem.P1_op;
                                        //spline.P2 = solve_problem.P2_op;
                                        //spline.P3 = solve_problem.P3_op;
                                        //spline.P4 = solve_problem.P4_op;
                                        //spline.P5 = solve_problem.P5_op;
                                        //printf("v0 after is:%f,%f,%f\n\n",spline.C1.x(),spline.C1.y(),spline.C1.z());

                                        //Update Velocity************************************************************************
                                        if(segment_index <= segment_index_inside){
                                            temp_velocity_buffer[2*segment_index] = spline.C1;
                                       // Vector3d est =  3*spline.C3+2*spline.C2+spline.C1;
                                        //printf("velocity start: is:%f,%f,%f\n",spline.C1.x(),spline.C1.y(),spline.C1.z());
                                        //printf("velocity end: is:%f,%f,%f\n\n",est.x(),est.y(),est.z());

                                            //printf("buffer size is:%d,and segment is :%d\n\n",temp_velocity_buffer.size(),segment_index);
                                            temp_velocity_buffer[2*segment_index+1] = 3*spline.C3+2*spline.C2+spline.C1;
                                           //velocity_coutinue = 3*spline.C3+2*spline.C2+spline.C1;
                                        }  */

                    /*                          int i =0;
                                        for(P_it;P_it!=it_&&i<_spline_pose_time.size();P_it++){
                                            Vector3d tmp_v_ ;
                                            (*P_it)->getVelocity(tmp_v_);
                                            //printf("%f current v is:%f,%f,%f\n",(*P_it)->time_stamp,tmp_v_.x(),tmp_v_.y(),tmp_v_.z());
                                            Vector3d tmp_v = 3*spline.C3*_spline_pose_time[i]*_spline_pose_time[i]+
                                            2*spline.C2*_spline_pose_time[i]+spline.C1;
                                            Vector3d tmp_a = 6*spline.C3*_spline_pose_time[i]+2*spline.C2;
                                            //printf("%f com a is:%f,%f,%f\n\n",(*P_it)->time_stamp,tmp_a.x(),tmp_a.y(),tmp_a.z());
                                            Vector3d tmp_a_ ;
                                            Vector3d tmp_w_;
                                            //(*P_it)->getIMU(tmp_a_,tmp_w_);no_start_pose

                                            (*P_it)->updateVelocity(tmp_v);
                                             //(*P_it)->updateAcc(tmp_a);
                                                i++;
                                        }   */
                    // else{
                    // temp_velocity_buffer[2*segment_index] = velocity_coutinue;
                    // temp_velocity_buffer[2*segment_index+1] = 3*spline.C3+2*spline.C2+spline.C1;
                    // velocity_coutinue = 3*spline.C3+2*spline.C2+spline.C1;
                    //  }
/*                     segment_index++;
                    switch_lable = false;
                    time_test_lable_ = 0;
                    _spline_pose.clear();
                    _spline_pose_time.clear();
                    _spline_velocity.clear(); */
               // }
            //}
            //***************************************************************************************
            //segment_index = 0;

            // Bspline
            //      for (it_B = keyframelist.begin(); it_B != keyframelist.end(); it_B++){
            // 实际上 cp_2才是current_pose
            //  printf("man!!!!!!\n");
            /*                list<KeyFrame*>::iterator contral_poins_1 = it_B;
                           // '++it' 是前置自增运算符，会先增加再使用迭代器
                           list<KeyFrame*>::iterator contral_poins_2 = contral_poins_1;

                           if(contral_poins_1 != keyframelist.end()) contral_poins_2++;
                        //   printf("manba1 out!\n");
                           list<KeyFrame*>::iterator contral_poins_3 = contral_poins_2;
                           if(contral_poins_2 != keyframelist.end()) contral_poins_3++;
                        //   printf("manba2 out!\n");
                           list<KeyFrame*>::iterator contral_poins_4 = contral_poins_3;
                           if(contral_poins_3 != keyframelist.end()) contral_poins_4++;
                           if(contral_poins_4 == keyframelist.end())
                               break; */
            TicToc Bspline_t;
            const int step = 5; // 定义步长为5
            bool ifpolo = true;
            double total_time_Bspline = 0;
            Vector3d alg_R_;
            Matrix3d deltaR;
            
/*             for (it_B = keyframelist.begin(); it_B != keyframelist.end();)
            {
                // 获取当前迭代器指向的控制点
                list<KeyFrame *>::iterator contral_poins_1 = it_B;

                // 递进迭代器到下一个控制点
                for (int i = 0; i < step && it_B != keyframelist.end(); ++i, ++it_B)
                    ;
                // 检查是否有足够的关键帧
                if (it_B == keyframelist.end())
                    break;
                list<KeyFrame *>::iterator contral_poins_2 = it_B;

                for (int i = 0; i < step && it_B != keyframelist.end(); ++i, ++it_B)
                    ;
                if (it_B == keyframelist.end())
                    break;
                list<KeyFrame *>::iterator contral_poins_3 = it_B;

                for (int i = 0; i < step && it_B != keyframelist.end(); ++i, ++it_B)
                    ;
                if (it_B == keyframelist.end())
                    break;
                list<KeyFrame *>::iterator contral_poins_4 = it_B;

                Bspline_ Bspline;
                Spline spline_test;
                //  printf("manba3 out!\n\n");
                Vector3d pose1;
                Matrix3d rotation1;
                Vector3d pose2;
                Matrix3d rotation2;
                Vector3d pose3;
                Matrix3d rotation3;
                Vector3d pose4;
                Matrix3d rotation4;

                Vector3d V1;
                Vector3d V2;
                Vector3d V3;
                Vector3d V4;

                list<KeyFrame *>::iterator contral_poins_1_ = contral_poins_1;
                list<KeyFrame *>::iterator contral_poins_2_ = contral_poins_2;
                list<KeyFrame *>::iterator contral_poins_3_ = contral_poins_3;
                list<KeyFrame *>::iterator contral_poins_4_ = contral_poins_4;

                vector<double *> para_buffer;
                // 获得各控制点的位姿信息
                (*contral_poins_1)->getPose(pose1, rotation1);
                // printf("manba25 out!\n");
                (*contral_poins_2)->getPose(pose2, rotation2);
                //   printf("manba26 out!\n");
                (*contral_poins_3)->getPose(pose3, rotation3);
                //  printf("manba27 out!\n");
                (*contral_poins_4)->getPose(pose4, rotation4);

                (*contral_poins_1)->getVelocity(V1);
                (*contral_poins_2)->getVelocity(V2);
                (*contral_poins_3)->getVelocity(V3);
                (*contral_poins_4)->getVelocity(V4);

                // 1 is start*******************************************************************************
                spline_test.compute_spline_argument(pose1, pose2, V1, V2);
                double head_time = (*contral_poins_1)->time_stamp;
                double tail_time = (*contral_poins_2)->time_stamp;
                vector<Vector3d> test_spline_velocity_1;
                vector<Vector3d> test_spline_pose_1;
                vector<double> test_spline_pose_time_1;
                // 3 is end*******************************************************************************
                //    }
                //    else{ *************************************************************
                                        it_B = contral_poins_2;
                                        double head_time_ = (*contral_poins_2)->time_stamp;
                                        //double current_time1_ = (*contral_poins_2)->time_stamp;
                                       //double current_time2_ = (*contral_poins_3)->time_stamp;
                                        double tail_time_ = (*contral_poins_3) ->time_stamp;

                                        double bspline_time=0;

                                        for (int i = 0; i < step ; ++i, contral_poins_2++){
                                        Matrix3d Bspline_Rotation_C;
                                        Vector3d Bspline_t_C;
                                        //double time_persent2_ = (current_time2_-head_time)/(tail_time-head_time);
                                        list<KeyFrame*>::iterator current_poins = contral_poins_2;

                                        (*current_poins)->getPose(Bspline_t_C,Bspline_Rotation_C);

                                        double current_time = (*current_poins)->time_stamp;
                                        double time_persent = (current_time-head_time_)/(tail_time_-head_time_);
                                        
                                        

                                        Vector3d Bspline_pose = Bspline.Pose_Bspline(time_persent,pose1,pose2,pose3,pose4);
                                        
                                        Matrix3d Bspline_Rotation = Bspline.Rotation_Bspline(time_persent,rotation1,rotation2,rotation3,rotation4);

                                        Vector3d Bspline_pose2_ = Bspline.Pose_Bspline(1,pose1,pose2,pose3,pose4);

                                        //Bspline_pose.z()+=5;
                                        //Bspline_pose2_.z()+=5;

                                       // (*contral_poins_2)->updatePose(Bspline_pose, Bspline_Rotation); 
                                       // (*contral_poins_3)->updatePose(Bspline_pose2_, rotation3);
                                       bspline_time += Bspline_t.toc();

                }
                
                total_time_Bspline +=bspline_time;
            } */
            printf("Bspline time is:%f\n",total_time_Bspline);

            //}
            //**********************************************************************************************

            double ceres_t_sum = 0;

            bool need_spline = false;
            Vector3d pose_continue;
            bool need_add_position = true;
            TicToc SecondLoop;

            std::ofstream clear_file("/home/wjb/output/omega.txt", std::ios::trunc);
            clear_file.close();
            //*********************************第二次loop，差值计算，补全路径***************************************************
            for (it_ = keyframelist.begin(), it_P = keyframelist.begin(),it_his = keyframelist.begin(); it_ != keyframelist.end();it_++)
            {
               /*  Vector3d k_v, k_a, k_w, k_t,kk_t;
                Matrix3d k_r,kk_r;
                (*it_)->getPose(k_t, k_r);
                (*it_his)->getPose(kk_t, kk_r);
                (*it_)->getVelocity(k_v);
                (*it_)->getIMU(k_a, k_w);
                double k = ((k_v.cross(k_a).norm())/(k_v.norm()*k_v.norm()*k_v.norm())); */
                //如果此时路径弯曲程度过大就应加入中间值
/*                 if(k>0.4){
                    need_add_position = true;
                    //pubPointCloud_backend(k_t);
                    //pubPointCloud_backend(kk_t);
                }
                else{
                    need_add_position = false; 
                } */
/*                 if((*it_)->has_loop){
                    pubPointCloud_backend(k_t);
                    pubPointCloud_backend(kk_t);
                } */
                
                //  printf("lable is: %s\n\n",(*it_)->Spline_frame_.c_str());
                //  printf("%fpose is:%f,%f,%f\n",(*it_)->time_stamp,k_t.x(),k_t.y(),k_t.z());
                if ((*it_)->index < first_looped_index)
                    continue;
                std::string Spline_lable_((*it_)->Spline_frame_);
                Vector3d start_V;
                Vector3d end_V;
                // 用于存储最小二乘所需的结尾帧位置数据
                Vector3d end_Pose;

                Vector3d start_pose ;
                Vector3d end_pose ;

                Vector3d start_velocity ;
                Vector3d end_velocity ;

                if (Spline_lable_ != "PoseGraph" && !(*it_)->has_loop)
                {
                    
                    if (!switch_lable)
                        it_P = it_;
                    // printf("start time is:%f,end time is:%f\n\n",temp_start_time,temp_end_time);

                    // 进入样条段内
                    switch_lable = true;
                    SESpline SEspline;
                    Bspline_ Bspline;

                    // 未更新的差值帧位姿
                    Matrix3d tmp_r;
                    Vector3d tmp_t;

                    double current_time = (*it_)->time_stamp;
                    // double current_time_t = spline_pose_time_buffer[segment_index][time_test_lable];
                    // double a_t = current_time_t+current_time;
                    // printf("current time is%f\n",current_time);
                    time_test_lable++;

                    Vector3d current_V;
                    Vector3d current_Acc;
                    Vector3d current_Gyc;
                    (*it_)->getVelocity(current_V);
                    (*it_)->getIMU(current_Acc, current_Gyc);
                    // printf("acc is:%f,%f,%f\n",current_Acc.x(),current_Acc.y(),current_Acc.z());
                    // printf("gyc is:%f,%f,%f\n\n",current_Gyc.x(),current_Gyc.y(),current_Gyc.z());
                    double start_time = temp_time_buffer[2 * segment_index];
                    double end_time = temp_time_buffer[2 * segment_index + 1];
                    double time_persent = (current_time - start_time);

                    if (abs(time_persent/(end_time-start_time)) > 1)
                    {
                       // printf("error time!\n\n");
                        continue;
                    }

                    // pose
                    start_pose = temp_pose_buffer[2 * segment_index];
                    end_pose = temp_pose_buffer[2 * segment_index + 1];

                    Vector3d pose1 = temp_pose_buffer[2 * segment_index];
                    Vector3d pose2 = bspline_pose[2 * segment_index];
                    Vector3d pose3 = bspline_pose[2 * segment_index + 1];
                    Vector3d pose4 = temp_pose_buffer[2 * segment_index + 1];

                    // rotation
                    Matrix3d start_rotation = temp_rotation_buffer[2 * segment_index];
                    Matrix3d end_rotation = temp_rotation_buffer[2 * segment_index + 1];

                    start_pose_test = start_pose;
                    end_pose_test = end_pose;
                    // no optimistic pose
                    Vector3d no_start_pose = temp_no_op_pose_buffer[2 * segment_index];
                    Vector3d no_end_pose = temp_no_op_pose_buffer[2 * segment_index + 1];
                    Matrix3d no_start_rotation = temp_no_op_rotation_buffer[2 * segment_index];
                    Matrix3d no_end_rotation = temp_no_op_rotation_buffer[2 * segment_index + 1];
                    end_Pose = no_end_pose;
                    start_V = temp_no_op_pose_buffer[2 * segment_index];
                    end_V = temp_no_op_pose_buffer[2 * segment_index + 1];
                    // velocity
                    start_velocity = temp_velocity_buffer[2 * segment_index];
                    end_velocity = temp_velocity_buffer[2 * segment_index + 1];

                    start_velocity_test = start_velocity;
                    end_velocity_test = end_velocity;
                    // IMU
                    Vector3d start_acc_ = temp_acc_buffer[2 * segment_index];
                    Vector3d end_acc_ = temp_acc_buffer[2 * segment_index + 1];

                    start_acc_test = start_acc_;
                    end_acc_test = end_acc_;

                    Vector3d start_gyc_ = temp_gyc_buffer[2 * segment_index];
                    Vector3d end_gyc_ = temp_gyc_buffer[2 * segment_index + 1];

                    // translation
                    Matrix3d delta_Q = start_rotation.transpose() * end_rotation;
                    Vector3d delta_angle = SO3_::anti_hat(SO3_::log(delta_Q));
                    Vector3d delta_pose_Q = start_rotation.transpose() * (end_pose - start_pose);

                    Matrix4d start_T = SE3_::setSE3(start_pose, start_rotation);
                    //  Matrix4d end_T = SE3_::setSE3(end_pose_Q,end_rotation);
                    // Matrix4d start_T = SE3_::setSE3(start_pose,Matrix3d::Identity());
                    // Matrix4d end_T = SE3_::setSE3(end_pose,Matrix3d::Identity());
                    Matrix4d delta_T = SE3_::setSE3(delta_pose_Q, delta_Q);
                    VectorXd Lie_alg_start(6);
                    VectorXd Lie_alg_end(6);

                    VectorXd Lie_alg_delta(6);
                    //  Vector3d rou_start = start_velocity-SO3_::hat(start_gyc_)*start_pose;
                    //  Vector3d rou_end = end_velocity-SO3_::hat(end_gyc_)*end_pose;
                    Vector3d rou_start = start_velocity;
                    Vector3d rou_end = end_velocity;
                    Lie_alg_start << rou_start.x(), rou_start.y(), rou_start.z(),
                        start_gyc_.x(), start_gyc_.y(), start_gyc_.z();
                    Lie_alg_end << rou_end.x(), rou_end.y(), rou_end.z(),
                        end_gyc_.x(), end_gyc_.y(), end_gyc_.z();

                    Lie_alg_delta << delta_pose_Q.x(), delta_pose_Q.y(), delta_pose_Q.z(),
                        delta_angle.x(), delta_angle.y(), delta_angle.z();

                    // SE3
                    Sophus::SE3d se_start_pose_ = temp_se_pose_buffer[2 * segment_index];
                    Sophus::SE3d se_after_start_pose_ = temp_se_after_pose_buffer[2 * segment_index];
                    Sophus::SE3d se_end_pose_ = temp_se_pose_buffer[2 * segment_index + 1];
                    // 计算“速度”
                    Eigen::Matrix<double, 6, 1> se_start_velocity;
                    Eigen::Matrix<double, 6, 1> se_end_velocity;

                    if (segment_index != 0 && (2 * segment_index + 1 != temp_se_pose_buffer.size() - 1))
                    {

                        se_start_velocity = 0.5 * (temp_se_pose_buffer[2 * segment_index + 1].log() -
                                                   temp_se_pose_buffer[2 * segment_index - 1].log());
                        se_end_velocity = 0.5 * (temp_se_pose_buffer[2 * segment_index + 2].log() -
                                                 temp_se_pose_buffer[2 * segment_index].log());
                        // se_start_velocity =   (temp_se_pose_buffer[2*segment_index+1].log()-
                        // temp_se_pose_buffer[2*segment_index].log());
                        // se_end_velocity   =   (temp_se_pose_buffer[2*segment_index+2].log()-
                        // temp_se_pose_buffer[2*segment_index+1].log());
                        // se_start_velocity =   temp_se_pose_buffer[2*segment_index].log();
                        // se_end_velocity   =   temp_se_pose_buffer[2*segment_index+1].log();
                    }
                    else if (segment_index = 0)
                    {
                        se_start_velocity = se_end_pose_.log() - se_start_pose_.log();
                        se_end_velocity = 0.5 * (temp_se_pose_buffer[2 * segment_index + 2].log() -
                                                 temp_se_pose_buffer[2 * segment_index].log());
                    }
                    else
                    {
                        se_start_velocity = 0.5 * (temp_se_pose_buffer[2 * segment_index + 1].log() -
                                                   temp_se_pose_buffer[2 * segment_index - 1].log());
                        se_end_velocity = se_end_pose_.log() - se_start_pose_.log();
                    }

                    // 获取位置信息
                    Vector3d tmp_a, tmp_w;
                   // (*it_)->getIMU(k_a, k_w);
                    (*it_)->getVioPose(tmp_t, tmp_r);
                    //(*it_)->getIMU()
                    // spline_pose_.push_back(tmp_t);//**********************
                    spline_velocity_.push_back(current_V);
                    spline_pose_time.push_back(time_persent);
                    spline_acc_.push_back(current_Acc);
                    spline_gyc_.push_back(current_Gyc);
                    // spline_omega_.push_back(tmp_w);
                   // spline_Rotaion_.push_back(start_rotation.transpose() * tmp_r);
                    // spline_Rotaion_.push_back(tmp_r);

                    TicToc spline_tmp_t;
                    double factor_up = (end_pose-start_pose).squaredNorm();
                    double factor_down = (no_end_pose-no_start_pose).squaredNorm();
                    double factor = factor_up/factor_down;

                    // 计算未更新前的delta
                    Vector3d delta_pose_head = no_start_pose - tmp_t;
                    Vector3d delta_pose_tail = no_end_pose - tmp_t;
                    Matrix3d delta_rotation_head = no_start_rotation.transpose()*tmp_r;
                    Matrix3d delta_rotation_tail = no_end_rotation.transpose()*tmp_r;
                    
                    Vector3d t_HC = start_pose - factor*delta_pose_head;
                    Vector3d t_TC = end_pose - factor*delta_pose_tail;

                    Matrix3d r_HC = start_rotation*delta_rotation_head;
                    Matrix3d r_TC = end_rotation*delta_rotation_tail;

                    Matrix3d rHC ;

                    

                    // 利用LERP加权平均计算更新后的威姿
                    Matrix3d SLERP_R = SO3_::SLERP(time_persent, r_HC, r_TC);
                    Eigen::Quaterniond mamba = Eigen::Quaterniond(start_rotation.transpose()*end_rotation);
                    Matrix3d Rotation_mamba = start_rotation*Eigen::Quaterniond::Identity().slerp(time_persent,mamba).toRotationMatrix();
                    Matrix3d final_R = SLERP_R;
                    spline_Rotaion_.push_back(final_R);

                    //printf("favtor is%f\n\n",factor);
                    
                    Vector3d LERP_test_t = (t_HC + ((t_HC - t_TC) * time_persent));

                   // printf("%d delata t is:%f,%f,%f\n",segment_index,(t_HC - t_TC).x(),(t_HC - t_TC).y(),(t_HC - t_TC).z());
                   // printf("%d delata t about H2T is:%f,%f,%f\n",segment_index,(end_pose - start_pose).x(),(end_pose - start_pose).y(),(end_pose - start_pose).z());
                    //Vector3d LERP_test_t = t_HC ;
                    total_lerp_t_sum+=spline_tmp_t.toc();

                    spline_pose_.push_back(LERP_test_t);
                    
                    //if(need_add_position)
                    pose_map[time_persent] = t_HC;
                    Vector3d Lie_ = SO3_::anti_hat(SO3_::log(start_rotation.transpose()*tmp_r));
                    rotation_map[time_persent] = Lie_;
                    //test_map[1] = start_acc_;
                    //test_map[0] = end_acc_;
                    //R3Spline计算参数与估计值
                    
                    spline.compute_spline_argument(start_pose,end_pose,start_velocity,end_velocity);
                    spline.acc_compute_spline_argument(start_pose, end_pose, start_velocity, end_velocity, start_acc_, end_acc_);

                  //  if(need_spline)
                   // {
                    //TicToc ceres_t;
                    //Spline_Least_Square_muti_ sq_spline(pose_map,start_pose,end_pose,start_velocity,end_velocity,end_time-start_time);
                    //sq_spline.compute_para_factor();
                    //轨迹飘飞，猜测可能是传入的test_map数据存在问题
                    //sq_spline.get_para(spline.C0,spline.C1,spline.C2,spline.C3);
                    velocity_map[time_persent] = current_V;
                    Vector3d vector_com = spline.C1+2*spline.C2+3*spline.C3;
                    //printf("IMU is:%f,%f,%f\n",current_V.x(),current_V.y(),current_V.z());
                    //printf("com V is:%f,%f,%f\n\n",vector_com.x(),vector_com.y(),vector_com.z());
                    
                    //    spline.compute_spline_argument(start_pose, end_pose, start_velocity, end_velocity);//*****************
                       // need_spline = false;

                        //ceres_t_sum += ceres_t.toc();
                        //ceres_op_time_sum += ceres_t_sum;
                       // printf("ceres op is:%f\n\n", ceres_op_time_sum);

                  //  }
                    

                    // 计算更新前的位姿
                    // 获取参数信息，acc_spline_para_指针指向一个18维数组地址
                    /* double *acc_spline_para_;
                    spline.get_acc_spline_para(acc_spline_para_); */

                    // 差值计算位姿
                    TicToc SO_tmp_t;
                    spline.compute_position(time_persent);
                    spline.acc_compute_position(time_persent);

                    // 计算 rotation 的spline 参数
                    Matrix3d delta_R = start_rotation.transpose() * end_rotation;
                    deltaR = delta_R;
                    Matrix3d M_start_gyc_ = SO3_::hat(start_gyc_);
                    Matrix3d M_end_gyc_ = SO3_::hat(end_gyc_);

                    Matrix3d alg_R = SO3_::log(delta_R);
                    alg_R_ = SO3_::anti_hat(alg_R);

                    // SO3Spline.compute_spline_para(delta_R,M_start_gyc_,M_end_gyc_,start_rotation);
                    SO3Spline.compute_spline_para(alg_R_, M_start_gyc_, M_end_gyc_, start_rotation);
                   // Eigen::Matrix3d spline_rotation = SO3Spline.compute_rotation(time_persent, start_rotation);
                    //  printf("SO3 time is :%f\n",SO_tmp_t.toc());


                    // Matrix3d SLERP_R = SO3Spline.SVDSpline(time_persent,start_time-(*keyframelist.begin())->time_stamp,end_time-(*keyframelist.begin())->time_stamp,start_rotation,end_rotation,start_gyc_,end_gyc_);

                    // SE3 LERP start
                    // ####################################################################################################
                    /*                   Sophus::SE3d com_curr_SEpose_ = SEspline.lerp_in_se(se_start_pose_,se_end_pose_,time_persent);
                                      Sophus::SE3d com_curr_SEpose_1 = SEspline.lerp_in_se(se_start_pose_,se_end_pose_,0.25);
                                      Sophus::SE3d com_curr_SEpose_2 = SEspline.lerp_in_se(se_start_pose_,se_end_pose_,0.75);

                                      Eigen::Vector3d SE_LERP_pose = com_curr_SEpose_.translation();
                                      Eigen::Matrix3d SE_LERP_rotation = com_curr_SEpose_.rotationMatrix();
                                      //Eigen::Vector3d SE_B_pose = test_Bslpine.translation();
                                      TicToc LERP_tmp_t;
                                      //double LERP_difference = (SE_LERP_pose-spline.compute_Liner_position(time_persent,start_pose,end_pose)).squaredNorm();
                                      double LERP_difference = (spline.position-tmp_t).squaredNorm();
                                      if(segment_index ==0)
                                          LERP_op_time_sum+=LERP_tmp_t.toc(); */
                    // ####################################################################################################
                    // SE3 SPline end

                    //***********************************************************************
                    // dejasulaju start
                    //***********************************************************************

                    // TicToc SE_tmp_t;
                    /*                     SEspline.Office_Phase(se_start_pose_,se_start_velocity,se_end_pose_,se_end_velocity,start_time,end_time);
                                        Sophus::SE3d spline_curr_SEpose_ =SEspline.Online_Phase(start_time,se_start_pose_,se_start_velocity,end_time,se_end_pose_,se_end_velocity,current_time);
                                       // if(segment_index ==0)
                                        //    SE_op_time_sum += SE_tmp_t.toc();
                                        Eigen::Vector3d SE_spline_pose = spline_curr_SEpose_.translation();
                                        printf("pose is : %f.%f.%f",SE_spline_pose.x(),SE_spline_pose.y(),SE_spline_pose.z());
                                        Eigen::Matrix3d SE_spline_rotation = spline_curr_SEpose_.rotationMatrix(); */

                    //***********************************************************************
                    // dejasulaju end
                    //***********************************************************************

                    //***********************************************************************
                    // cubic fuction start
                    //***********************************************************************

                    // 计算spline参数
                    TicToc SE_tmp_t;
                    SE3spline.computeSplineParameter(Lie_alg_delta, Lie_alg_start, Lie_alg_end);
                    // 计算delta
                    Matrix4d SE3_delta_Trans = SE3spline.compute_trans(time_persent);
                    // 计算真实数据
                    // Matrix4d SE3_Translation = start_T*SE3_delta_Trans;
                    // 输出信息
                    /*                          printf("translation is:\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n\n",
                                            SE3_delta_Trans(0,0),SE3_delta_Trans(0,1),SE3_delta_Trans(0,2),SE3_delta_Trans(0,3),
                                            SE3_delta_Trans(1,0),SE3_delta_Trans(1,1),SE3_delta_Trans(1,2),SE3_delta_Trans(1,3),
                                            SE3_delta_Trans(2,0),SE3_delta_Trans(2,1),SE3_delta_Trans(2,2),SE3_delta_Trans(2,3),
                                            SE3_delta_Trans(3,0),SE3_delta_Trans(3,1),SE3_delta_Trans(3,2),SE3_delta_Trans(3,3));   */
                    Matrix3d SE3SPline_rotation = start_rotation * SE3_::toRotationMatrix(SE3_delta_Trans);
                    Vector3d SE3Spline_pose = start_pose + SE3_::toMovementVector(SE3_delta_Trans);
                    // printf("SE time is :%f\n\n",SE_tmp_t.toc());
                    // Vector3d SE3Spline_pose = SE3SPline_rotation*SE3_::toRotationMatrix(start_T).transpose()*SE3_::toMovementVector(start_T)+SE3_::toMovementVector(SE3_delta_Trans);

                    //***********************************************************************
                    // cubic fuction end
                    //***********************************************************************

                    //***********************************************************************
                    // Bspline fuction start
                    //***********************************************************************

                    // Vector3d Bspline_pose = Bspline.Pose_Bspline(time_persent,pose1,pose2,pose3,pose4);
                    /*                     printf("time p is:%f\n",time_persent);
                                        printf("pose1  is : %f,%f,%f\n",pose1.x(),pose1.y(),pose1.z());
                                        printf("pose2  is : %f,%f,%f\n",pose2.x(),pose2.y(),pose2.z());
                                        printf("pose3  is : %f,%f,%f\n",pose3.x(),pose3.y(),pose3.z());
                                        printf("pose4  is : %f,%f,%f\n",pose4.x(),pose4.y(),pose4.z());
                                        printf("current  is : %f,%f,%f\n\n",Bspline_pose.x(),Bspline_pose.y(),Bspline_pose.z()); */
                    // printf("pose is:%f,%f,%f\n\n",Bspline_pose.x(),Bspline_pose.y(),Bspline_pose.z());

                    //***********************************************************************
                    // Bspline fuction end
                    //***********************************************************************

                    // double Spline_difference = (SE_spline_pose-spline.compute_Liner_position(time_persent,start_pose,end_pose)).squaredNorm();
                    double Spline_difference = (tmp_t - spline.position).squaredNorm();
                    // Vector3d c_V = spline.C3*time_persent*time_persent*3+time_persent*2*spline.C2+spline.C1;

                    // 确认差值方式
                    if (!std::isnan(SE3Spline_pose.x()) && !std::isnan(SE3Spline_pose.y()) && !std::isnan(SE3Spline_pose.z()))
                    // if(true)
                    {
                        
                        // tmp_t = SE_spline_pose;
                        // tmp_t =spline.acc_position;
                        // if(k>2&&abs(Spline_difference)<0.2)
                        // tmp_t = spline.compute_Liner_position(time_persent,start_pose,end_pose);
                        tmp_t = spline.position;
                        // tmp_t = t_HC;
                        // tmp_t = Bspline_pose;
                        //tmp_t = LERP_test_t;
                        // SE_count_A++;
                        // pubPointCloud_backend(LERP_test_t);

                        // if(abs(tmp_r.norm()-spline_rotation.norm())<0.0001)
                        //  {
                        // if(!tmp_r.hasNaN())

                        /*                         printf("tmp_r is:\n%f,%f,%f\n%f,%f,%f\n%f,%f,%f\n\n",
                                                tmp_r(0,0),tmp_r(0,1),tmp_r(0,2),tmp_r(1,0),tmp_r(1,1),tmp_r(1,2),tmp_r(2,0),tmp_r(2,1),tmp_r(2,2));
                                                printf("rotation is:\n%f,%f,%f\n%f,%f,%f\n%f,%f,%f\n\n",
                                                spline_rotation(0,0),spline_rotation(0,1),spline_rotation(0,2),spline_rotation(1,0),
                                                spline_rotation(1,1),spline_rotation(1,2),spline_rotation(2,0),spline_rotation(2,1),spline_rotation(2,2)); */
                        //tmp_r=spline_rotation;
                        /*                       Matrix3d test_r = tmp_r - start_rotation;
                                                printf("r is :\n %f.%f,%f,\n%f.%f,%f,\n%f.%f,%f,\n\n",
                                                test_r(0,0),test_r(0,1),test_r(0,2),
                                                test_r(1,0),test_r(1,1),test_r(1,2),
                                                test_r(2,0),test_r(2,1),test_r(2,2)
                                                ); */

                        // tmp_r = spline_rotation;
                        // tmp_r=SE3SPline_rotation;
                        if(!isnan(Quaterniond(final_R).w()))
                            tmp_r = final_R;
                        else
                            printf("rotation not an number!!!\n");
                        //  }
                        //  else{
                        // printf("e error is:%f\n\n",abs(tmp_r.norm()-spline_rotation.norm()));
                        //}

                        // Vector3d compute_V = spline.compute_velocity(time_persent);
                        // double factor_V = ((current_V.x()/compute_V.x())+(current_V.y()/compute_V.y())+(current_V.z()/compute_V.z()))/3;
                        // compute_V *= factor_V;
                        // printf("current velocity is:%f,%f.%f\n",current_V.x(),current_V.y(),current_V.z());
                        // printf("compute velocity is:%f,%f.%f\n\n",compute_V.x(),compute_V.y(),compute_V.z());
                        // tmp_r = SE_LERP_rotation;
                        // tmp_r = SE_spline_rotation;
                        // tmp_t = SE_LERP_pose;
                        // tmp_t = SE_B_pose;
                        // tmp_t = spline.compute_Liner_position(time_persent,start_pose,end_pose);
                    }
                    // else if(abs(LERP_difference)<0.002){
                    /* else if(LERP_difference<0.002){
                        //tmp_t = SE_LERP_pose;
                       //tmp_t = spline.compute_Liner_position(time_persent,start_pose,end_pose);
                       //tmp_t=LERP_test_t;
                        //tmp_t = spline.position;
                        //tmp_t =spline.acc_position;
                        SE_count_B++;
                    }*/
                    else
                    {
                        // if( (tmp_t - spline.compute_Liner_position(time_persent,start_pose,end_pose)).squaredNorm()< 2 ){
                        //  tmp_t = spline.compute_Liner_position(time_persent,start_pose,end_pose);
                        // tmp_t = SE_LERP_pose;
                        // tmp_t=LERP_test_t;
                        SE_count_A++;
                        tmp_t = spline.position;
                        // tmp_t = spline.acc_position;
                        // tmp_t = SE3Spline_pose;
                        //}
                    }

                    // tmp_t = test_spline;

                    // Spline in R3
                    // tmp_t = spline.position;
                    /*
                        Vector3d ab = end_pose - start_pose;
                        Vector3d ac = tmp_t - start_pose;
                        Vector3d M = ab.cross(ac);
                        double gongchang = M.norm()/ab.norm();
                    */

                    // LEPR in R3
                    // tmp_t = spline.compute_Liner_position(time_persent,start_pose,end_pose);
                    // printf("LERP_tmpt is:%f,%f,%f\n\n",tmp_t.x(),tmp_t.y(),tmp_t.z());
                    // Vector3d tmpt_ = tmp_r*tmp_t;
                    // tmp_t.z()+=1;
                    //(*it_)->getVioPose(tmp_t,tmp_r);
                    (*it_)->updatePose(tmp_t, tmp_r); //******************************

                    // spline.compute_spline_argument(no_start_pose,no_end_pose,start_velocity,end_velocity);

                    // spline.acc_compute_spline_argument(no_start_pose,no_end_pose,start_velocity,end_velocity,start_acc_test,end_acc_test);

                    /* Vector3d com_v = 5*spline.P5*time_persent*time_persent*time_persent*time_persent
                    +4*spline.P4*time_persent*time_persent*time_persent+3*spline.P3*time_persent*time_persent
                    +2*time_persent*spline.P2+spline.P1; */
                    // Vector3d com_a = 6*spline.C3*time_persent+2*spline.C2;
                    // printf("%f compute vo is:%f,%f,%f\n",time_persent,com_a.x(),com_a.y(),com_a.z());
                    // printf("%F current vo is:%f,%f,%f\n\n",time_persent,current_Acc.x(),current_Acc.y(),current_Acc.z());
                   // printf("inside map size is: %d\n\n",test_map.size());
                  // printf("loop segment is:%d\n",segment_index);
                   
                }
                else if (switch_lable)
                {
                    
                    TicToc polo_time;
                    if(Spline_para_buffer.empty()||segment_index>Spline_para_buffer.rbegin()->first){
                    start_pose = temp_pose_buffer[2 * segment_index];
                    end_pose = temp_pose_buffer[2 * segment_index + 1];
                    Eigen::Matrix3d start_rotation = temp_rotation_buffer[2 * segment_index];
                    Eigen::Matrix3d end_rotation = temp_rotation_buffer[2 * segment_index + 1];
                    Vector3d start_gyc_ = temp_gyc_buffer[2 * segment_index];
                    Vector3d end_gyc_ = temp_gyc_buffer[2 * segment_index + 1];

                   // printf("s_i is:%d\n",segment_index);
                   // printf("rebegin is:%d\n\n",Spline_para_buffer.rbegin()->first);
                    start_velocity = temp_velocity_buffer[2 * segment_index];
                   
                    //start_velocity = spline.C1 + 2*spline.C2 + 3*spline.C3;
                    //start_pose = spline.C0 + spline.C1 + spline.C2 + spline.C3;
                    //printf("\n\n!!!!!!!new v is:%f,%f,%f\n",start_velocity.x(),start_velocity.y(),start_velocity.z());
                    end_velocity = temp_velocity_buffer[2 * segment_index + 1]; 

                    
                    Spline_Least_Square_muti_ sq_spline(pose_map,start_pose,end_pose,start_velocity,end_velocity,temp_time_buffer[2 * segment_index+1]-temp_time_buffer[2 * segment_index]);

                    
                    sq_spline.compute_para_factor();
                    //printf("polo time cost is:%f\n\n",polo_time.toc());
                    
                    Vector3d start_rotation_lie = SO3_::anti_hat(SO3_::log(start_rotation));
                    Vector3d end_rotation_lie = SO3_::anti_hat(SO3_::log(end_rotation));
                   // rotation_map.clear();
                    SO3Spline_Least_Square_muti_ SO3Spline_LS_(rotation_map,start_rotation_lie,alg_R_,start_gyc_,end_gyc_,temp_time_buffer[2 * segment_index+1]-temp_time_buffer[2 * segment_index]);
                    SO3Spline_Least_Square_muti SO3Spline_LS(rotation_map,start_rotation_lie,alg_R_,start_gyc_,end_gyc_,temp_time_buffer[2 * segment_index+1]-temp_time_buffer[2 * segment_index]);
                    //Vector3d acc_gyc =(end_gyc_ -  start_gyc_)/(temp_time_buffer[2 * segment_index+1]-temp_time_buffer[2 * segment_index]);
                    //test_easy SO3test(start_gyc_,acc_gyc,deltaR);

                    Vector3d RO3_1;
                    Vector3d RO3_2;
                    Vector3d RO3_3;
                    Vector3d RO3_4;

                if(rotation_map.size()>1){
                        SO3Spline_LS_.compute_para_factor();
                        SO3Spline_LS_.get_para(RO3_4,RO3_3,RO3_2,RO3_1);
                        noV_sum++;
                    }
                    else
                   {
                        SO3Spline_LS.compute_para_factor();
                        SO3Spline_LS.get_para(RO3_4,RO3_3,RO3_2,RO3_1);
                        V_sum++;
                    }

                    //SO3test.get_para(RO3_3,RO3_2,RO3_1);

                    SO3Spline.updatePara(RO3_1,RO3_2,RO3_3);
                    vector<Vector3d> spline_rotation_para;

                    spline_rotation_para.push_back(RO3_1);
                    spline_rotation_para.push_back(RO3_2);
                    spline_rotation_para.push_back(RO3_3);
                    Spline_Rotation_para_buffer[segment_index] = spline_rotation_para;

                   // if(pose_map.size()!=0){
                        //spline.compute_spline_argument(start_pose, end_pose, start_velocity/10, end_velocity/10);
                        sq_spline.get_para(spline.C0,spline.C1,spline.C2,spline.C3);
                        vector<Vector3d> spline_para;
                        spline_para.push_back(spline.C1);
                        spline_para.push_back(spline.C2);
                        spline_para.push_back(spline.C3);
                        Spline_para_buffer[segment_index]=spline_para;//segment_index可以代表段落数，接下来需要考虑一下如何判断当此段已经有spline了就不要计算了，直接用已经有的数据就行
                        
                   // }
                    }
                    ceres_t_sum+=polo_time.toc();
                   

   
/*                  //****************************************************************
                    // ceres优化需要找18维数组指针，vector pose，vector timepersent
                    Eigen::Vector3d P0_, P1_, P2_, P3_, P4_, W1_, W2_, W3_;
                    spline.get_spline_para(P0_, P1_, P2_, P3_);
                    // double splinetest[18]={P0_.x(),P0_.y(),P0_.z(),P1_.x(),P1_.y(),P1_.z(),P2_.x(),P2_.y(),P2_.z(),P3_.x(),P3_.y(),P3_.z(),
                    //                       P4_.x(),P4_.y(),P4_.z(),P5_.x(),P5_.y(),P5_.z()};
                    double spline_pose_pa[12] = {P0_.x(), P0_.y(), P0_.z(), P1_.x(), P1_.y(), P1_.z(), P2_.x(), P2_.y(), P2_.z(), P3_.x(), P3_.y(), P3_.z()};
                    // SolveCeres solve_problem(end_Pose,spline_pose_pa,spline_pose_,spline_pose_time);
                    // solve_problem.Solve();
                    // 更新数据
                    //  spline.C0 = solve_problem.P0_op;
                    //   spline.C1 = solve_problem.P1_op;
                    //   spline.C2 = solve_problem.P2_op;
                    //  spline.C3 = solve_problem.P3_op;

                    //                     NomalSolveCeres solve_problem(spline_pose_pa,spline_pose_,spline_pose_time);
                    //                    solve_problem.Solve();
                    //                    //更新数据
                    //                    spline.C0 = solve_problem.P0_op;
                    //                    spline.C1 = solve_problem.P1_op;
                    //                    spline.C2 = solve_problem.P2_op;
                    //                    spline.C3 = solve_problem.P3_op; 

                    double spline_velocity_pa[12] = {spline.C0.x(), spline.C0.y(), spline.C0.z(),
                                                     spline.C1.x(), spline.C1.y(), spline.C1.z(), spline.C2.x(), spline.C2.y(), spline.C2.z(),
                                                     spline.C3.x(), spline.C3.y(), spline.C3.z()};
                    // VSolveCeres vsolve_problem(spline_velocity_pa,spline_velocity_,spline_pose_time);
                    // vsolve_problem.Solve();
                    // 更新数据
                    // spline.C0 = vsolve_problem.P0_op;
                    // spline.C1 = vsolve_problem.P1_op;
                    // spline.C2 = vsolve_problem.P2_op;
                    // spline.C3 = vsolve_problem.P3_op;

                    // printf("time is:%f,%f,%f\n",spline_pose_time[0]);
                    // printf("out spline P0 is:%f,%f.%f\n\n",splinetest[0],splinetest[1],splinetest[2]);
                    // solve_problem.get_para()[2]);
                    // printf("v0 befor is:%f,%f,%f\n",spline.C1.x(),spline.C1.y(),spline.C1.z());

                    VPSolveCeres vpsolve_problem(spline_pose_pa, spline_velocity_, spline_pose_, spline_pose_time);
                    // VPSolveCeres vpsolve_problem(spline_pose_pa,spline_acc_,spline_pose_,spline_pose_time);

                    vpsolve_problem.Solve();
                    // printf("origin p0 is:%f,%f,%f\n",spline.C1.x(),spline.C1.y(),spline.C1.z());
                    

                    spline.C0 = vpsolve_problem.P0_op;
                    spline.C1 = vpsolve_problem.P1_op;
                    spline.C2 = vpsolve_problem.P2_op;
                    spline.C3 = vpsolve_problem.P3_op;
                    // spline.P4 = solve_problem.P4_op;
                    // spline.P5 = solve_problem.P5_op;
                    // printf("v0 after is:%f,%f,%f\n\n",spline.C1.x(),spline.C1.y(),spline.C1.z());
                    spline.get_spline_para(P0_, P1_, P2_, P3_);

                    double _spline_pose_pa[12] = {P0_.x(), P0_.y(), P0_.z(), P1_.x(), P1_.y(), P1_.z(), P2_.x(), P2_.y(), P2_.z(), P3_.x(), P3_.y(), P3_.z()};

                    

                //    VVPSolveCeres vvpsolve_problem(_spline_pose_pa,spline_velocity_,spline_pose_,spline_pose_time);
                    // VPSolveCeres vpsolve_problem(spline_pose_pa,spline_acc_,spline_pose_,spline_pose_time);

                //    vvpsolve_problem.Solve();

                //    spline.C0 = vvpsolve_problem.P0_op;
                //    spline.C1 = vvpsolve_problem.P1_op;
                //    spline.C2 = vvpsolve_problem.P2_op;
                //    spline.C3 = vvpsolve_problem.P3_op;  

                    SO3Spline.get_para(W1_, W2_, W3_);
                    double spline_pose_pa_[18] = {P1_.x(), P1_.y(), P1_.z(), P2_.x(), P2_.y(), P2_.z(), P3_.x(), P3_.y(), P3_.z(),
                                                  W1_.x(), W1_.y(), W1_.z(), W2_.x(), W2_.y(), W2_.z(), W3_.x(), W3_.y(), W3_.z()};
                    VPWSolveCeres vpwsolve_problem(spline_pose_pa_, spline_velocity_, spline_pose_, spline_Rotaion_, spline_gyc_, spline_pose_time);
                    vpwsolve_problem.Solve();

                    SO3Spline.updatePara(vpwsolve_problem.W1_op, vpwsolve_problem.W2_op, vpwsolve_problem.W3_op); */

                    
                    int i = 0;
                    bool debug = true;
                    bool Debug = true;
                    // printf("size is:%d\n\n",spline_pose_time.size());
                    if (!spline_pose_time.empty())
                    {   
                        if(segment_index==0){
                            if(debug)
                                debug = false;
                            if(!debug)
                                Debug = false;
                        }
                        if(Debug&&pose_map.size()!=0)
                        {   
                            spline.C1 = Spline_para_buffer[segment_index][0];
                            spline.C2 = Spline_para_buffer[segment_index][1];
                            spline.C3 = Spline_para_buffer[segment_index][2];
                            spline.C0 = temp_pose_buffer[2 * segment_index];
                            //printf("%dtruestart pose is:%f\n\n",segment_index,C0.x(),C0.y(),C0.z());
                            SO3Spline.updatePara(Spline_Rotation_para_buffer[segment_index][0],Spline_Rotation_para_buffer[segment_index][1],Spline_Rotation_para_buffer[segment_index][2]);
                            //printf("P0 is:%f,%f,%f\n\n",spline.C0.x(),spline.C0.y(),spline.C0.z());
                        }
                        else if(pose_map.size()==0){
                            spline.compute_spline_argument(start_pose,end_pose,start_velocity,end_velocity);
                            printf("error! pose_map.size()==0 \n");
                        };

                        Vector3d Start_tmp_t;
                        Vector3d Start_tmp_v;
                        Matrix3d Start_tmp_r;
                        (*it_P)->getPose(Start_tmp_t, Start_tmp_r);
                        // printf("update time start is:%f\n\n",(*it_P)->time_stamp);
                       // printf("\nP0 is:%f,%f,%f\n\n",spline.P0.x(),spline.P0.y(),spline.P0.z());

                        for (it_P; i < spline_pose_time.size();)
                        {                             
                            Vector3d tmp_t;
                            Vector3d tmp_v;
                            Matrix3d tmp_r;
                            // double c_t = (*it_P)->time_stamp;
                            (*it_P)->getPose(tmp_t, tmp_r);
                            (*it_P)->getVelocity(tmp_v);

/*                             ofstream _loop_path_file("/home/wjb/output/omega.txt", ios::app);
                            _loop_path_file.setf(ios::fixed, ios::floatfield);
                            _loop_path_file.precision(5);
                            _loop_path_file <<SO3_::anti_hat(SO3_::log(tmp_r)).x()<<" ";
                            _loop_path_file <<SO3_::anti_hat(SO3_::log(tmp_r)).y()<<" ";
                            _loop_path_file <<SO3_::anti_hat(SO3_::log(tmp_r)).z()<<endl;
                            _loop_path_file.close(); */
                            
                            
                            //printf("time is:%f,and delta rotation is%f\n",spline_pose_time[i],alg_R_.squaredNorm());
                            // printf("time size is:%d,r size is: %d\n\n",spline_pose_time.size(),spline_Rotaion_.size());
/*                             if(rotation_map.size()>1){
                                tmp_r = SO3Spline.compute_rotation(spline_pose_time[i],Start_tmp_r);
                            } */
                                

                                
                            //Vector3d omega = SO3_::anti_hat(SO3_::log(tmp_r));
                            // tmp_r = spline_Rotaion_[i];
                            if(Debug){
                                TicToc polo_time_2;
                                spline.compute_position(spline_pose_time[i]);
                                total_ceres_t_sum+=polo_time_2.toc();
                                tmp_t = spline.position;
                                //tmp_t = spline_pose_[i];
                                //tmp_t = spline.correction_spline(start_pose,end_pose,spline_pose_time[i]);                       
                                }
                                
                            else{
                                printf("wait! not debug!");
                                tmp_t = spline_pose_[i];
                                }
                                

                            //Vector3d com_v = 3*spline.C3*spline_pose_time[i]*spline_pose_time[i]+2*spline.C2*spline_pose_time[i]+spline.C1;
                            //Vector3d com_a = 6*spline.C3*spline_pose_time[i]+2*spline.C2;
                            //(*it_P)->updateVelocity(com_v);
                            (*it_P)->updatePose(tmp_t,tmp_r);//***************************************
                            it_P++;
                            ++i;
                        }
                    }
                    // 估计是指针越界了
 
                    //***************************************************************************************
                    segment_index++;
                    switch_lable = false;
                    time_test_lable = 0;
                    spline_pose_.clear();
                    spline_pose_time.clear();
                    spline_velocity_.clear();
                    spline_Rotaion_.clear();
                    need_spline = true;
                    pose_map.clear();
                    velocity_map.clear();
                    rotation_map.clear();
                    // 记录时间
                    if(segment_index>test_sum)
                    {
                        test_sum = segment_index;
                    }
                }
                else
                {
                    pose_count++;
                }
                it_his = it_;
            }

            total_ceres_t_sum += ceres_t_sum;
           // printf("second loop time is:%f\n",SecondLoop.toc());
            printf("No add V sum is:%d,add V %d\n",noV_sum,V_sum);
            printf("polo time cost is:%f\n",total_ceres_t_sum);
            printf("lerp time cost is:%f\n\n",total_lerp_t_sum);
            // loop end
            // printf("algs cost time is: %f\n",test_algs.toc());
            //printf("se+r3 is:%d,se is:%d,R3 is :%d\n", SE_count_A, SE_count_B, R3_count);
            //printf("LERP time is :%f,R3 time is : %f.SE3 time is :%f\n\n", LERP_op_time_sum, spline_op_time_sum, SE_op_time_sum);
            // printf("pose num is :%d\n\n",pose_count);
            Vector3d cur_t, vio_t;
            Matrix3d cur_r, vio_r;
            cur_kf->getPose(cur_t, cur_r);
            cur_kf->getVioPose(vio_t, vio_r);
            m_drift.lock();
            yaw_drift = Utility::R2ypr(cur_r).x() - Utility::R2ypr(vio_r).x();
            r_drift = Utility::ypr2R(Vector3d(yaw_drift, 0, 0));
            t_drift = cur_t - r_drift * vio_t;
            m_drift.unlock();

            // cout << "t_drift " << t_drift.transpose() << endl;
            // cout << "r_drift " << Utility::R2ypr(r_drift).transpose() << endl;
            // cout << "yaw drift " << yaw_drift << endl;

            it++;

            for (; it != keyframelist.end(); it++)
            {
                Vector3d P;
                Matrix3d R;
                (*it)->getVioPose(P, R);
                P = r_drift * P + t_drift;
                R = r_drift * R;

                //不知道有啥用，目前注释掉了如果以后又有在取消注释
               // (*it)->updatePose(P, R);
            }

            m_keyframelist.unlock();
            //目前的想法是给updatePath函数传入一个list，里面包含了全部的需要补全的点信息
            updatePath();
            // if(path_change <=3);
            // part_optimize_lable=false;
            // ROS_INFO("solver costs: %fms", tmp_t.toc());
        }

        std::chrono::milliseconds dura(2000);
        std::this_thread::sleep_for(dura);
    }
}

void PoseGraph::updatePath()
{
    m_keyframelist.lock();
    p_v_SEspline.clear();
    point_cloud_SE.points.clear();
    list<KeyFrame *>::iterator it;
    list<KeyFrame *>::iterator it_add = additional_keyframelist.begin();
    for (int i = 1; i <= sequence_cnt; i++)
    {
        path[i].poses.clear();
    }
    base_path.poses.clear();
    posegraph_visualization->reset();

    if (SAVE_LOOP_PATH)
    {
        ofstream loop_path_file_tmp(VINS_RESULT_PATH, ios::out);
        loop_path_file_tmp.close();
    }

    bool path_change_ = false;

    //补充帧信息
/*     for(it_add = additional_keyframelist.begin();it_add != additional_keyframelist.end();it_add++){

        Vector3d P;
        Matrix3d R;
        (*it_add)->getPose(P, R);
        Quaterniond Q;
        Q = R;

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header.stamp = ros::Time((*it_add)->time_stamp);
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
        pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
        pose_stamped.pose.position.z = P.z();
        pose_stamped.pose.orientation.x = Q.x();
        pose_stamped.pose.orientation.y = Q.y();
        pose_stamped.pose.orientation.z = Q.z();
        pose_stamped.pose.orientation.w = Q.w();

        if (!(*it_add)->has_loop && (*it_add)->time_stamp >= time_stap_path)
        //if ((*it)->Spline_frame_ != "PoseGraph" && !(*it)->has_loop && (*it)->time_stamp >= time_stap_path)
        {
            path_change_ = true;
            path[path_change].poses.push_back(pose_stamped);
            path[path_change].header = pose_stamped.header;
            time_stap_path = (*it_add)->time_stamp;
        }
        else if (path_change_)
        {
            if (path_change <= 6)
                path_change++;
            path_change_ = false;
        }

        if ((*it_add)->sequence == 0)
        {
            base_path.poses.push_back(pose_stamped);
            base_path.header = pose_stamped.header;
        }
        else
        {

            path[(*it_add)->sequence].poses.push_back(pose_stamped);
            path[(*it_add)->sequence].header = pose_stamped.header;
        }

        // euroc
        if (SAVE_LOOP_PATH)
        {
            ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(9);
            loop_path_file << (*it_add)->time_stamp << " ";
            loop_path_file.precision(5);
            loop_path_file << P.x() << " "
                           << P.y() << " "
                           << P.z() << " "
                           << Q.x() << " "
                           << Q.y() << " "
                           << Q.z() << " "
                           << Q.w() << endl;
            loop_path_file.close();
        }

    }
    additional_keyframelist.clear(); */
    int add_index = 1;
    it_add = keyframelist.begin();
    bool switch_lable = false;
    int segment_index_update = 0;

    for (it = keyframelist.begin(); it != keyframelist.end(); it++)
    {   
        //判断段落切换
         if ((*it)->Spline_frame_ != "PoseGraph" && !(*it)->has_loop)
                {
                    switch_lable = true;
                }
            else if(switch_lable){
                segment_index_update++;
                switch_lable = false;
            }
        

        Vector3d P;
        Matrix3d R;
        Vector3d V;
        (*it)->getPose(P, R);
        (*it)->getVelocity(V);
        Quaterniond Q;
        Q = R;
        
        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header.stamp = ros::Time((*it)->time_stamp);
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
        pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
        pose_stamped.pose.position.z = P.z();
        pose_stamped.pose.orientation.x = Q.x();
        pose_stamped.pose.orientation.y = Q.y();
        pose_stamped.pose.orientation.z = Q.z();
        pose_stamped.pose.orientation.w = Q.w();

        Vector3d P_add;
        Matrix3d R_add;
        
        Quaterniond Q_add;
        Vector3d R_inter;
        geometry_msgs::PoseStamped pose_stamped_add;
        
         if(add_index)
        {   
            for( int i = 0;i < 3;i++ ){
            (*it_add)->getPose(P_add,R_add);

            Eigen::Quaterniond mamba = Eigen::Quaterniond(R_add.transpose()*R);
            //Eigen::Quaterniond Rotation_mamba_out = Eigen::Quaterniond(R_add)*Eigen::Quaterniond::Identity().slerp(1/4*(i+1),mamba);

            //Q_add = R_add;
            R_inter = P_add+(P-P_add)/4*(i+1);
           // printf("H is:%f,%f,%f\n",P_add.x(),P_add.y(),P_add.z());
           // printf("R is:%f,%f,%f\n",R_inter.x(),R_inter.y(),R_inter.z());
           // printf("T is:%f,%f,%f\n\n",P.x(),P.y(),P.z()); 

/*             pose_stamped_add.header.stamp = ros::Time((*it_add)->time_stamp+((*it)->time_stamp - (*it_add)->time_stamp)/4*(i+1));
            pose_stamped_add.header.frame_id = "world";
            pose_stamped_add.pose.position.x = R_inter.x();
            pose_stamped_add.pose.position.y = R_inter.y();
            pose_stamped_add.pose.position.z = R_inter.z();
            pose_stamped_add.pose.orientation.x = Rotation_mamba_out.x();
            pose_stamped_add.pose.orientation.y = Rotation_mamba_out.y();
            pose_stamped_add.pose.orientation.z = Rotation_mamba_out.z();
            pose_stamped_add.pose.orientation.w = Rotation_mamba_out.w();

            ofstream loop_path_file_(VINS_RESULT_PATH, ios::app);
            loop_path_file_.setf(ios::fixed, ios::floatfield);
            loop_path_file_.precision(9);
            loop_path_file_ << (*it_add)->time_stamp+((*it)->time_stamp - (*it_add)->time_stamp)/4*(i+1) << " ";
            loop_path_file_.precision(5);
            loop_path_file_ << R_inter.x() << " "
                            << R_inter.y() << " "
                            << R_inter.z() << " "
                            << Rotation_mamba_out.x() << " "
                            << Rotation_mamba_out.y() << " "
                            << Rotation_mamba_out.z() << " "
                            << Rotation_mamba_out.w() << endl;
            loop_path_file_.close();  */

            
        //}

           //  ofstream loop_path_file_(VINS_RESULT_PATH, ios::app);
           // loop_path_file_.setf(ios::fixed, ios::floatfield);
           // loop_path_file_.precision(9);
           // loop_path_file_ << pose_stamped_add.header.stamp << " ";
           // loop_path_file_.precision(5);
           // loop_path_file_ << R_inter.x() << " "
           //                << R_inter.y() << " "
           //                << R_inter.z() << " "
           //                << Q_add.x() << " "
            //               << Q_add.y() << " "
            //               << Q_add.z() << " "
            //               << Q_add.w() << endl;
            //loop_path_file_.close(); 
            if(switch_lable&&segment_index_update<test_sum){
                

                double start_time = temp_time_buffer[2 * segment_index_update];
                double end_time = temp_time_buffer[2 * segment_index_update + 1];
                Vector3d start_pose = temp_pose_buffer[2 * segment_index_update];
                double time_persent = ((*it)->time_stamp-start_time)/(end_time-start_time);
                double spline_time =(*it_add)->time_stamp+((*it)->time_stamp - (*it_add)->time_stamp)/4*(i+1)-start_time;
                double spline_time_rotation =((*it)->time_stamp - (*it_add)->time_stamp)/4*(i+1);

                Vector3d start_rotation = SO3_::anti_hat(SO3_::log(R_add));
                Vector3d alg_R_ = SO3_::anti_hat(SO3_::log(R_add.transpose()*R));

                Vector3d start_gyc_ ,start_acc;
                Vector3d end_gyc_ ,end_acc;

                (*it_add)->getIMU(start_acc,start_gyc_);
                (*it)->getIMU(start_acc,end_gyc_);

                SO3Spline_Least_Square_muti_test SO3Spline_test_test(start_rotation,alg_R_,start_gyc_,end_gyc_,(*it)->time_stamp-(*it_add)->time_stamp);

                Vector3d P0,P1,P2,P3;
                SO3Spline_test_test.compute_para_factor();
                SO3Spline_test_test.get_para(P3,P2,P1,P0);

                SO3Spline SO3Spline_mamba;
                SO3Spline_mamba.R0_ = R_add;
                //SO3Spline_mamba.updatePara(Spline_Rotation_para_buffer[segment_index_update][0],Spline_Rotation_para_buffer[segment_index_update][1],Spline_Rotation_para_buffer[segment_index_update][2]);
                SO3Spline_mamba.updatePara(P0,P1,P2);
                Matrix3d Rotation_mamba_out_ = SO3Spline_mamba.compute_rotation(spline_time_rotation,R_add);

                Eigen::Quaterniond Rotation_mamba_out = Eigen::Quaterniond(Rotation_mamba_out_);
                //double spline_time =1;
                //printf("%dspline time is:%f\n",segment_index_update,spline_time);
                //printf("%dstart pose is:%f\n",segment_index_update,start_pose.x(),start_pose.y(),start_pose.z());
                Vector3d C1 = Spline_para_buffer[segment_index_update][0];
                //printf("%dC1 is:%f\n\n",segment_index_update,C1.x(),C1.y(),C1.z());
                Vector3d C2 = Spline_para_buffer[segment_index_update][1];
                Vector3d C3 = Spline_para_buffer[segment_index_update][2];
                Vector3d spline_pose = start_pose+C1*spline_time+C2*spline_time*spline_time+C3*spline_time*spline_time*spline_time;
                //printf("error is:%f\n",(spline_pose-temp_pose_buffer[2 * segment_index_update]).squaredNorm());
                //printf("error is:%f\n",(spline_pose-R_inter).squaredNorm());

             //if((spline_pose-R_inter).squaredNorm()<0.0009){
/*                 ofstream loop_path_file_(VINS_RESULT_PATH, ios::app);
                loop_path_file_.setf(ios::fixed, ios::floatfield);
                loop_path_file_.precision(9);
                loop_path_file_ << (*it_add)->time_stamp+((*it)->time_stamp - (*it_add)->time_stamp)/4*(i+1) << " ";
                loop_path_file_.precision(5);
                loop_path_file_ << spline_pose.x() << " "
                                << spline_pose.y() << " "
                                << spline_pose.z() << " "
                                << Rotation_mamba_out.x() << " "
                                << Rotation_mamba_out.y() << " "
                                << Rotation_mamba_out.z() << " "
                                << Rotation_mamba_out.w() << endl;
                loop_path_file_.close(); */
           // }
            // else if(segment_index_update==31){
                //printf("P lerp is:%f,%f,%f\n",R_inter.x(),R_inter.y(),R_inter.z());
/*                 printf("%dspline time is:%f\n",segment_index_update,spline_time);
                printf("%dC0 is:%f,%f,%f\n",segment_index_update,start_pose.x(),start_pose.y(),start_pose.z());
                printf("%dC1 is:%f,%f,%f\n",segment_index_update,C1.x(),C1.y(),C1.z());
                printf("%dC2 is:%f,%f,%f\n",segment_index_update,C2.x(),C2.y(),C2.z());
                printf("%dC3 is:%f,%f,%f\n\n",segment_index_update,C3.x(),C3.y(),C3.z()); */
           // } 
            
                }
            }
            it_add=it;
        }
        

/*         Vector3d P_add;
        Matrix3d R_add;
        Quaterniond Q_add;
        geometry_msgs::PoseStamped pose_stamped_add;
        if(it_add!=additional_keyframelist.end()){

        (*it_add)->getPose(P_add,R_add);
        Q_add = R_add;

        pose_stamped_add.header.stamp = ros::Time((*it_add)->time_stamp);
        pose_stamped_add.header.frame_id = "world";
        pose_stamped_add.pose.position.x = P_add.x() + VISUALIZATION_SHIFT_X;
        pose_stamped_add.pose.position.y = P_add.y() + VISUALIZATION_SHIFT_Y;
        pose_stamped_add.pose.position.z = P_add.z();
        pose_stamped_add.pose.orientation.x = Q_add.x();
        pose_stamped_add.pose.orientation.y = Q_add.y();
        pose_stamped_add.pose.orientation.z = Q_add.z();
        pose_stamped_add.pose.orientation.w = Q_add.w();

        } */



        if (!(*it)->has_loop && (*it)->time_stamp >= time_stap_path)
        //if ((*it)->Spline_frame_ != "PoseGraph" && !(*it)->has_loop && (*it)->time_stamp >= time_stap_path)
        {
            path_change_ = true;
            path[path_change].poses.push_back(pose_stamped);
            path[path_change].header = pose_stamped.header;

/*             path[path_change].poses.push_back(pose_stamped_add);
            path[path_change].header = pose_stamped_add.header; */
            time_stap_path = (*it)->time_stamp;
        }
        else if (path_change_)
        {
            if (path_change <= 6)
                path_change++;
            path_change_ = false;
        }

        if ((*it)->sequence == 0)
        {
            base_path.poses.push_back(pose_stamped);
            base_path.header = pose_stamped.header;
        }
        else
        {   

/*             if(index)
            {   
               printf("inter R time is:%f\n",pose_stamped_add.header.stamp);
               printf("inter R is:%f,%f,%f\n\n",pose_stamped_add.pose.position.x,pose_stamped_add.pose.position.y,pose_stamped_add.pose.position.z);
                path[(*it)->sequence].poses.push_back(pose_stamped_add);
                path[(*it)->sequence].header = pose_stamped_add.header;
             } */
/*             printf("inter R time is:%f\n",pose_stamped.header.stamp);
            printf("R is:%f,%f,%f\n",pose_stamped.pose.position.x,pose_stamped.pose.position.y,pose_stamped.pose.position.z); */
            path[(*it)->sequence].poses.push_back(pose_stamped);
            path[(*it)->sequence].header = pose_stamped.header;



        }
        /*####################################################################################
                if (SAVE_LOOP_PATH)
                {
                    ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
                    loop_path_file.setf(ios::fixed, ios::floatfield);
                    loop_path_file.precision(0);
                    loop_path_file << (*it)->time_stamp * 1e9 << ",";
                    loop_path_file.precision(5);
                    loop_path_file  << P.x() << ","
                          << P.y() << ","
                          << P.z() << ","
                          << Q.w() << ","
                          << Q.x() << ","
                          << Q.y() << ","
                          << Q.z() << ","
                          << endl;
                    loop_path_file.close();
                }
         ######################################################################################       */
        // euroc
        if (SAVE_LOOP_PATH)
        {

            //输出轨迹值
            ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(9);
            loop_path_file << (*it)->time_stamp << " ";
            loop_path_file.precision(5);
            loop_path_file << P.x() << " "
                           << P.y() << " "
                           << P.z() << " "
                           << Q.x() << " "
                           << Q.y() << " "
                           << Q.z() << " "
                           << Q.w() << endl;
            loop_path_file.close();

/*             ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(9);
            loop_path_file << (*it)->time_stamp << " ";
            loop_path_file.precision(5);
            loop_path_file << P.x() << " "
                           << P.y() << " "
                           << P.z() << " "
                           << V.x() << " "
                           << V.y() << " "
                           << V.z() << endl;
            loop_path_file.close(); */

/*         if(it_add!=additional_keyframelist.end()){
            ofstream loop_path_file_(VINS_RESULT_PATH, ios::app);
            loop_path_file_.setf(ios::fixed, ios::floatfield);
            loop_path_file_.precision(9);
            loop_path_file_ << (*it_add)->time_stamp << " ";
            loop_path_file_.precision(5);
            loop_path_file_ << P_add.x() << " "
                           << P_add.y() << " "
                           << P_add.z() << " "
                           << Q_add.x() << " "
                           << Q_add.y() << " "
                           << Q_add.z() << " "
                           << Q_add.w() << endl;
            loop_path_file_.close();
        } */
        }

        /*####################################################################################
       //TUM
        if (SAVE_LOOP_PATH)
               {
                   ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
                   loop_path_file.setf(ios::fixed, ios::floatfield);cd
                   loop_path_file.precision(9);
                   loop_path_file << (*it)->time_stamp << " ";
                   loop_path_file.precision(5);
                   loop_path_file  << P.x() << " "
                                   << P.y() << " "
                                   << P.z() << " "
                                   << Q.x() << " "
                                   << Q.y() << " "
                                   << Q.z() << " "
                                   << Q.w() << endl;
                   loop_path_file.close();
               }
       ######################################################################################       */
        // draw local connection
        if (SHOW_S_EDGE)
        {
            list<KeyFrame *>::reverse_iterator rit = keyframelist.rbegin();
            list<KeyFrame *>::reverse_iterator lrit;
            for (; rit != keyframelist.rend(); rit++)
            {
                if ((*rit)->index == (*it)->index)
                {
                    lrit = rit;
                    lrit++;
                    for (int i = 0; i < 4; i++)
                    {
                        if (lrit == keyframelist.rend())
                            break;
                        if ((*lrit)->sequence == (*it)->sequence)
                        {
                            Vector3d conncected_P;
                            Matrix3d connected_R;
                            (*lrit)->getPose(conncected_P, connected_R);
                            posegraph_visualization->add_edge(P, conncected_P);
                        }
                        lrit++;
                    }
                    break;
                }
            }
        }
        if (SHOW_L_EDGE)
        {
            if ((*it)->has_loop && (*it)->sequence == sequence_cnt)
            {

                KeyFrame *connected_KF = getKeyFrame((*it)->loop_index);
                Vector3d connected_P;
                Matrix3d connected_R;
                connected_KF->getPose(connected_P, connected_R);
                //(*it)->getVioPose(P, R);
                (*it)->getPose(P, R);
                if ((*it)->sequence > 0)
                {
                    posegraph_visualization->add_loopedge(P, connected_P + Vector3d(VISUALIZATION_SHIFT_X, VISUALIZATION_SHIFT_Y, 0));
                }
            }
        }
        //删除下面这一行index的“//”可以用lerp补全部分威姿
        //index++;
    }
    publish();
    m_keyframelist.unlock();
}

void PoseGraph::part_updatePath(int frame_start, int frame_end)
{
    m_keyframelist.lock();

    for (int i = 1; i <= sequence_cnt; i++)
    {
        path[i].poses.erase(path[i].poses.begin(), path[i].poses.end() - 5);
    }
    base_path.poses.clear();
    posegraph_visualization->reset();

    if (SAVE_LOOP_PATH)
    {
        ofstream loop_path_file_tmp(VINS_RESULT_PATH, ios::out);
        loop_path_file_tmp.close();
    }

    for (int j = 0; j <= frame_end; j++)
    {
        KeyFrame *it = getKeyFrame(j);
        Vector3d P;
        Matrix3d R;
        it->getPose(P, R);
        Quaterniond Q;
        Q = R;
        // printf("path p: %f, %f, %f\n",  P.x(),  P.z(),  P.y() );

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header.stamp = ros::Time(it->time_stamp);
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.position.x = P.x() + VISUALIZATION_SHIFT_X;
        pose_stamped.pose.position.y = P.y() + VISUALIZATION_SHIFT_Y;
        pose_stamped.pose.position.z = P.z();
        pose_stamped.pose.orientation.x = Q.x();
        pose_stamped.pose.orientation.y = Q.y();
        pose_stamped.pose.orientation.z = Q.z();
        pose_stamped.pose.orientation.w = Q.w();
        if (it->sequence == 0)
        {
            base_path.poses.push_back(pose_stamped);
            base_path.header = pose_stamped.header;
        }
        else
        {
            path[it->sequence].poses.push_back(pose_stamped);
            path[it->sequence].header = pose_stamped.header;
        }

        if (SAVE_LOOP_PATH)
        {
            ofstream loop_path_file(VINS_RESULT_PATH, ios::app);
            loop_path_file.setf(ios::fixed, ios::floatfield);
            loop_path_file.precision(0);
            loop_path_file << it->time_stamp * 1e9 << ",";
            loop_path_file.precision(5);
            loop_path_file << P.x() << ","
                           << P.y() << ","
                           << P.z() << ","
                           << Q.w() << ","
                           << Q.x() << ","
                           << Q.y() << ","
                           << Q.z() << ","
                           << endl;
            loop_path_file.close();
        }
        // draw local connection
        /*if (SHOW_S_EDGE)
        {
            list<KeyFrame*>::reverse_iterator rit = keyframelist.rbegin();
            list<KeyFrame*>::reverse_iterator lrit;
            for (; rit != keyframelist.rend(); rit++)
            {
                if ((*rit)->index == (*it)->index)
                {
                    lrit = rit;
                    lrit++;
                    for (int i = 0; i < 4; i++)
                    {
                        if (lrit == keyframelist.rend())
                            break;
                        if((*lrit)->sequence == (*it)->sequence)
                        {
                            Vector3d conncected_P;
                            Matrix3d connected_R;
                            (*lrit)->getPose(conncected_P, connected_R);
                            posegraph_visualization->add_edge(P, conncected_P);
                        }
                        lrit++;
                    }
                    break;
                }
            }
        }*/
    }
    publish();
    m_keyframelist.unlock();
}

void PoseGraph::savePoseGraph()
{
    m_keyframelist.lock();
    TicToc tmp_t;
    FILE *pFile;
    printf("pose graph path: %s\n", POSE_GRAPH_SAVE_PATH.c_str());
    printf("pose graph saving... \n");
    string file_path = POSE_GRAPH_SAVE_PATH + "pose_graph.txt";
    pFile = fopen(file_path.c_str(), "w");
    // fprintf(pFile, "index time_stamp Tx Ty Tz Qw Qx Qy Qz loop_index loop_info\n");
    list<KeyFrame *>::iterator it;
    for (it = keyframelist.begin(); it != keyframelist.end(); it++)
    {
        std::string image_path, descriptor_path, brief_path, keypoints_path;
        if (DEBUG_IMAGE)
        {
            image_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_image.png";
            imwrite(image_path.c_str(), (*it)->image);
        }
        Quaterniond VIO_tmp_Q{(*it)->vio_R_w_i};
        Quaterniond PG_tmp_Q{(*it)->R_w_i};
        Vector3d VIO_tmp_T = (*it)->vio_T_w_i;
        Vector3d PG_tmp_T = (*it)->T_w_i;

        fprintf(pFile, " %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f %f %d\n", (*it)->index, (*it)->time_stamp,
                VIO_tmp_T.x(), VIO_tmp_T.y(), VIO_tmp_T.z(),
                PG_tmp_T.x(), PG_tmp_T.y(), PG_tmp_T.z(),
                VIO_tmp_Q.w(), VIO_tmp_Q.x(), VIO_tmp_Q.y(), VIO_tmp_Q.z(),
                PG_tmp_Q.w(), PG_tmp_Q.x(), PG_tmp_Q.y(), PG_tmp_Q.z(),
                (*it)->loop_index,
                (*it)->loop_info(0), (*it)->loop_info(1), (*it)->loop_info(2), (*it)->loop_info(3),
                (*it)->loop_info(4), (*it)->loop_info(5), (*it)->loop_info(6), (*it)->loop_info(7),
                (int)(*it)->keypoints.size());

        // write keypoints, brief_descriptors   vector<cv::KeyPoint> keypoints vector<BRIEF::bitset> brief_descriptors;
        assert((*it)->keypoints.size() == (*it)->brief_descriptors.size());
        brief_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_briefdes.dat";
        std::ofstream brief_file(brief_path, std::ios::binary);
        keypoints_path = POSE_GRAPH_SAVE_PATH + to_string((*it)->index) + "_keypoints.txt";
        FILE *keypoints_file;
        keypoints_file = fopen(keypoints_path.c_str(), "w");
        for (int i = 0; i < (int)(*it)->keypoints.size(); i++)
        {
            brief_file << (*it)->brief_descriptors[i] << endl;
            fprintf(keypoints_file, "%f %f %f %f\n", (*it)->keypoints[i].pt.x, (*it)->keypoints[i].pt.y,
                    (*it)->keypoints_norm[i].pt.x, (*it)->keypoints_norm[i].pt.y);
        }
        brief_file.close();
        fclose(keypoints_file);
    }
    fclose(pFile);

    printf("save pose graph time: %f s\n", tmp_t.toc() / 1000);
    m_keyframelist.unlock();
}
void PoseGraph::loadPoseGraph()
{
    TicToc tmp_t;
    FILE *pFile;
    string file_path = POSE_GRAPH_SAVE_PATH + "pose_graph.txt";
    printf("lode pose graph from: %s \n", file_path.c_str());
    printf("pose graph loading...\n");
    pFile = fopen(file_path.c_str(), "r");
    if (pFile == NULL)
    {
        printf("lode previous pose graph error: wrong previous pose graph path or no previous pose graph \n the system will start with new pose graph \n");
        return;
    }
    int index;
    double time_stamp;
    double VIO_Tx, VIO_Ty, VIO_Tz;
    double PG_Tx, PG_Ty, PG_Tz;
    double VIO_Qw, VIO_Qx, VIO_Qy, VIO_Qz;
    double PG_Qw, PG_Qx, PG_Qy, PG_Qz;
    double loop_info_0, loop_info_1, loop_info_2, loop_info_3;
    double loop_info_4, loop_info_5, loop_info_6, loop_info_7;
    int loop_index;
    int keypoints_num;
    Eigen::Matrix<double, 8, 1> loop_info;
    int cnt = 0;
    while (fscanf(pFile, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %d", &index, &time_stamp,
                  &VIO_Tx, &VIO_Ty, &VIO_Tz,
                  &PG_Tx, &PG_Ty, &PG_Tz,
                  &VIO_Qw, &VIO_Qx, &VIO_Qy, &VIO_Qz,
                  &PG_Qw, &PG_Qx, &PG_Qy, &PG_Qz,
                  &loop_index,
                  &loop_info_0, &loop_info_1, &loop_info_2, &loop_info_3,
                  &loop_info_4, &loop_info_5, &loop_info_6, &loop_info_7,
                  &keypoints_num) != EOF)
    {
        /*
        printf("I read: %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %d\n", index, time_stamp,
                                    VIO_Tx, VIO_Ty, VIO_Tz,
                                    PG_Tx, PG_Ty, PG_Tz,
                                    VIO_Qw, VIO_Qx, VIO_Qy, VIO_Qz,
                                    PG_Qw, PG_Qx, PG_Qy, PG_Qz,
                                    loop_index,
                                    loop_info_0, loop_info_1, loop_info_2, loop_info_3,
                                    loop_info_4, loop_info_5, loop_info_6, loop_info_7,
                                    keypoints_num);
        */
        cv::Mat image;
        std::string image_path, descriptor_path;
        if (DEBUG_IMAGE)
        {
            image_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_image.png";
            image = cv::imread(image_path.c_str(), 0);
        }

        Vector3d VIO_T(VIO_Tx, VIO_Ty, VIO_Tz);
        Vector3d PG_T(PG_Tx, PG_Ty, PG_Tz);
        Quaterniond VIO_Q;
        VIO_Q.w() = VIO_Qw;
        VIO_Q.x() = VIO_Qx;
        VIO_Q.y() = VIO_Qy;
        VIO_Q.z() = VIO_Qz;
        Quaterniond PG_Q;
        PG_Q.w() = PG_Qw;
        PG_Q.x() = PG_Qx;
        PG_Q.y() = PG_Qy;
        PG_Q.z() = PG_Qz;
        Matrix3d VIO_R, PG_R;
        VIO_R = VIO_Q.toRotationMatrix();
        PG_R = PG_Q.toRotationMatrix();
        Eigen::Matrix<double, 8, 1> loop_info;
        loop_info << loop_info_0, loop_info_1, loop_info_2, loop_info_3, loop_info_4, loop_info_5, loop_info_6, loop_info_7;

        if (loop_index != -1)
            if (earliest_loop_index > loop_index || earliest_loop_index == -1)
            {
                earliest_loop_index = loop_index;
                earliest_loop_index = -1;
            }

        // load keypoints, brief_descriptors
        string brief_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_briefdes.dat";
        std::ifstream brief_file(brief_path, std::ios::binary);
        string keypoints_path = POSE_GRAPH_SAVE_PATH + to_string(index) + "_keypoints.txt";
        FILE *keypoints_file;
        keypoints_file = fopen(keypoints_path.c_str(), "r");
        vector<cv::KeyPoint> keypoints;
        vector<cv::KeyPoint> keypoints_norm;
        vector<BRIEF::bitset> brief_descriptors;
        for (int i = 0; i < keypoints_num; i++)
        {
            BRIEF::bitset tmp_des;
            brief_file >> tmp_des;
            brief_descriptors.push_back(tmp_des);
            cv::KeyPoint tmp_keypoint;
            cv::KeyPoint tmp_keypoint_norm;
            double p_x, p_y, p_x_norm, p_y_norm;
            if (!fscanf(keypoints_file, "%lf %lf %lf %lf", &p_x, &p_y, &p_x_norm, &p_y_norm))
                printf(" fail to load pose graph \n");
            tmp_keypoint.pt.x = p_x;
            tmp_keypoint.pt.y = p_y;
            tmp_keypoint_norm.pt.x = p_x_norm;
            tmp_keypoint_norm.pt.y = p_y_norm;
            keypoints.push_back(tmp_keypoint);
            keypoints_norm.push_back(tmp_keypoint_norm);
        }
        brief_file.close();
        fclose(keypoints_file);

        KeyFrame *keyframe = new KeyFrame(time_stamp, index, VIO_T, VIO_R, PG_T, PG_R, image, loop_index, loop_info, keypoints, keypoints_norm, brief_descriptors);
        loadKeyFrame(keyframe, 0);
        if (cnt % 20 == 0)
        {
            publish();
        }
        cnt++;
    }
    fclose(pFile);
    printf("load pose graph time: %f s\n", tmp_t.toc() / 1000);
    base_sequence = 0;
}

void PoseGraph::publish()
{
    for (int i = 1; i <= sequence_cnt; i++)
    {
        // printf("sequence_cnt is :%d",sequence_cnt);
        // if (sequence_loop[i] == true || i == base_sequence)
        if (1 || i == base_sequence)
        {
            pub_pg_path.publish(path[i]);
            pub_path[i].publish(path[i]);
            for (int i = 2; i <= path_change; i++)
                pub_path[i].publish(path[i]);
            posegraph_visualization->publish_by(pub_pose_graph, path[sequence_cnt].header);
        }
    }
    base_path.header.frame_id = "world";
    pub_base_path.publish(base_path);
    // posegraph_visualization->publish_by(pub_pose_graph, path[sequence_cnt].header);
}

void PoseGraph::updateKeyFrameLoop(int index, Eigen::Matrix<double, 8, 1> &_loop_info)
{
    KeyFrame *kf = getKeyFrame(index);
    kf->updateLoop(_loop_info);
    if (abs(_loop_info(7)) < 30.0 && Vector3d(_loop_info(0), _loop_info(1), _loop_info(2)).norm() < 20.0)
    {
        if (FAST_RELOCALIZATION)
        {
            KeyFrame *old_kf = getKeyFrame(kf->loop_index);
            Vector3d w_P_old, w_P_cur, vio_P_cur;
            Matrix3d w_R_old, w_R_cur, vio_R_cur;
            old_kf->getPose(w_P_old, w_R_old);
            kf->getVioPose(vio_P_cur, vio_R_cur);

            Vector3d relative_t;
            Quaterniond relative_q;
            relative_t = kf->getLoopRelativeT();
            relative_q = (kf->getLoopRelativeQ()).toRotationMatrix();
            w_P_cur = w_R_old * relative_t + w_P_old;
            w_R_cur = w_R_old * relative_q;
            double shift_yaw;
            Matrix3d shift_r;
            Vector3d shift_t;
            shift_yaw = Utility::R2ypr(w_R_cur).x() - Utility::R2ypr(vio_R_cur).x();
            shift_r = Utility::ypr2R(Vector3d(shift_yaw, 0, 0));
            shift_t = w_P_cur - w_R_cur * vio_R_cur.transpose() * vio_P_cur;

            m_drift.lock();
            yaw_drift = shift_yaw;
            r_drift = shift_r;
            t_drift = shift_t;
            m_drift.unlock();
        }
    }
}