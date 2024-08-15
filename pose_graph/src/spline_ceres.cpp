#include <iostream>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include "so3.h"


using namespace std;
using namespace Eigen;
struct SplineCeres{
    public:
        SplineCeres(const double start_pose_x,const double start_pose_y,const double start_pose_z,
                const double pose_x,const double pose_y,const double pose_z,const double time_persent)
                : _start_pose_x(start_pose_x),_start_pose_y(start_pose_y),_start_pose_z(start_pose_z),
                    _pose_x(pose_x),_pose_y(pose_y),_pose_z(pose_z),_time_persent(time_persent){}

        template<typename T>
        bool operator()(
            const T *const spline_parameter,//18维->6个参数，每个参数三维
            T *residual
        )const{
            T com_pose[6];
            //compute_pose(spline_parameter,_time_persent,com_pose);
            //x-(a*t^5+b*t^4+c*t^3+d*t^2+e*t+f).x()
            //y-(a*t^5+b*t^4+c*t^3+d*t^2+e*t+f).y()
            //z-(a*t^5+b*t^4+c*t^3+d*t^2+e*t+f).z()
            /* com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
            spline_parameter[9]*_time_persent*_time_persent*_time_persent+
            spline_parameter[12]*_time_persent*_time_persent*_time_persent*_time_persent+
            spline_parameter[15]*_time_persent*_time_persent*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
            spline_parameter[10]*_time_persent*_time_persent*_time_persent+
            spline_parameter[13]*_time_persent*_time_persent*_time_persent*_time_persent+
            spline_parameter[16]*_time_persent*_time_persent*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
            spline_parameter[11]*_time_persent*_time_persent*_time_persent+
            spline_parameter[14]*_time_persent*_time_persent*_time_persent*_time_persent+
            spline_parameter[17]*_time_persent*_time_persent*_time_persent*_time_persent*_time_persent;
 */
            


            com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
                            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
                             spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
                            spline_parameter[11]*_time_persent*_time_persent*_time_persent; 

            com_pose[3] = spline_parameter[0]+spline_parameter[3]+spline_parameter[6]+
                            spline_parameter[9];

            com_pose[4] = spline_parameter[1]+spline_parameter[4]+spline_parameter[7]+
                            spline_parameter[10];

            com_pose[5] = spline_parameter[2]+spline_parameter[5]+spline_parameter[8]+
                            spline_parameter[11]; 

            //com_pose[3] = spline_parameter[0]+spline_parameter[3]+spline_parameter[6]+spline_parameter[9];

           // com_pose[4] = spline_parameter[1]+spline_parameter[4]+spline_parameter[7]+spline_parameter[10];

            //com_pose[5] = spline_parameter[2]+spline_parameter[5]+spline_parameter[8]+spline_parameter[11]; 

            residual[0] = _pose_x - com_pose[0];
            residual[1] = _pose_y - com_pose[1];
            residual[2] = _pose_z - com_pose[2];
            //residual[3] = _start_pose_x - com_pose[3];
            //residual[4] = _start_pose_y - com_pose[4];
            //residual[5] = _start_pose_z - com_pose[5];
           /*  residual[3]=_pose_x1- com_pose[3];
            residual[4]=_pose_y1- com_pose[4];
            residual[5]=_pose_z1- com_pose[5]; */
            return true;
        }

        static ceres::CostFunction *Create(const double _start_pose_x,const double _start_pose_y,const double _start_pose_z,const double _pose_x,const double _pose_y,const double _pose_z,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<SplineCeres,3,12>(new SplineCeres(_start_pose_x,_start_pose_y,_start_pose_z,_pose_x,_pose_y,_pose_z,_time_persent)));
        }

private:

        double _start_pose_x;
        double _start_pose_y;
        double _start_pose_z;
        double _pose_x;
        double _pose_y;
        double _pose_z;
        double _pose_x1;
        double _pose_y1;
        double _pose_z1;

        double _time_persent;
        
        vector<double>_data_y;

};

class SolveCeres{
    public:
    SolveCeres(Eigen::Vector3d end_pose,double *spline_parameter,vector<Eigen::Vector3d>spline_pose,const vector<double> time_persent):
                _end_pose(end_pose),spline_parameter_(spline_parameter),spline_pose_(spline_pose),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        for(int i=0;i<spline_pose_.size()&&i<time_persent_.size();i++){
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_pose = spline_pose_[i];
            double _time_persent = time_persent_[i];
            cost_funciton = SplineCeres::Create(_end_pose.x(),_end_pose.y(),_end_pose.z(),current_pose.x(),current_pose.y(),current_pose.z(),_time_persent);
            ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
            problem.AddResidualBlock(cost_funciton,loss_function,spline_para);
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options,&problem,&summary);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        //Eigen::Vector3d P0_op;
        P0_op.x()=spline_para[0];
        P0_op.y()=spline_para[1];
        P0_op.z()=spline_para[2];
        //Eigen::Vector3d P1_op;
        P1_op.x()=spline_para[3];
        P1_op.y()=spline_para[4];
        P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
        P2_op.x()=spline_para[6];
        P2_op.y()=spline_para[7];
        P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
        P3_op.x()=spline_para[9];
        P3_op.y()=spline_para[10];
        P3_op.z()=spline_para[11];
/*        // Eigen::Vector3d P4_op;
        P4_op.x()=spline_para[12];
        P4_op.y()=spline_para[13];
        P4_op.z()=spline_para[14];
        //Eigen::Vector3d P5_op;
        P5_op.x()=spline_para[15];
        P5_op.y()=spline_para[16];
        P5_op.z()=spline_para[17]; */
        //std::cout<<summary.FullReport()<<"\n";

    }

    double* get_para(){
        return spline_parameter_;
    }

    //private:
        Eigen::Vector3d _end_pose;
        Eigen::Vector3d P0_op;
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;
        Eigen::Vector3d P4_op;
        Eigen::Vector3d P5_op;
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_pose_;
        const vector<double> time_persent_;
        
};

struct VSplineCeres{
    public:
        VSplineCeres(const double velocity_x,const double velocity_y,const double velocity_z,const double time_persent)
         : _velocity_x(velocity_x),_velocity_y(velocity_y),_velocity_z(velocity_z),_time_persent(time_persent){}


        template<typename T>
        bool operator()(
            const T *const spline_parameter,//9维
            T *residual
        )const{
            //v=3a*t^2+2b*t+c
            //v0是真实值的速度，v1是估计值的速度
            //cos(angle)=a*b/|a|*|b| 
            //|a| = sqrt(a_x^2+a_y^2+a_z^2)
            //a*b = a_x*b_x+a_y*b_y+a_z*b_z
            T com_angle[1];
            T com_pose_x = spline_parameter[0]+T(2)*spline_parameter[3]*T(_time_persent)+
            T(3)*spline_parameter[6]*T(_time_persent)*T(_time_persent);

            T com_pose_y = spline_parameter[1]+T(2)*spline_parameter[4]*T(_time_persent)+
            T(3)*spline_parameter[7]*T(_time_persent)*T(_time_persent);

            T com_pose_z = spline_parameter[2]+T(2)*spline_parameter[5]*T(_time_persent)+
            T(3)*spline_parameter[8]*T(_time_persent)*T(_time_persent);

            T ab = T(_velocity_x)*T(com_pose_x)+T(_velocity_y)*T(com_pose_y)+T(_velocity_z)*T(com_pose_z);

            T a_norm = sqrt(T(_velocity_x)*T(_velocity_x)+T(_velocity_y)*T(_velocity_y)+T(_velocity_z)*T(_velocity_z));
            T b_norm = sqrt(T(com_pose_x)*T(com_pose_x)+T(com_pose_y)*T(com_pose_y)+T(com_pose_z)*T(com_pose_z));

            com_angle[0] = cos(ab/a_norm*b_norm);

            residual[0]=T(1) - com_angle[0];

            //printf("residual is %f,%f,%f\n\n",residual[0],residual[1],residual[2]);
            //printf("time persent is:%f\n\n",_time_persent);
            //residual[3]=_pose_x1- com_pose[3];
            //residual[4]=_pose_y1- com_pose[4];
            //residual[5]=_pose_z1- com_pose[5];
            return true;
        }

        static ceres::CostFunction *Create(const double _velocity_x,const double _velocity_y,const double _velocity_z,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<VSplineCeres,1,9>(new VSplineCeres(_velocity_x,_velocity_y,_velocity_z,_time_persent)));
        }

private:

        double _velocity_x;
        double _velocity_y;
        double _velocity_z;

        double _time_persent;
        
        //vector<double>_data_y;

};



class VSolveCeres{
    public:
    VSolveCeres(double *spline_parameter,vector<Eigen::Vector3d>spline_velocity,const vector<double> time_persent):
                spline_parameter_(spline_parameter),spline_velocity_(spline_velocity),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        //printf("spline_velocity_ size is:%d,time_persent_ size is:%d",spline_velocity_.size(),time_persent_.size());
        for(int i=0;i<spline_velocity_.size()&&i<time_persent_.size();i++){
/*             P1_op.x()=spline_para[3];
            P1_op.y()=spline_para[4];
            P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
            P2_op.x()=spline_para[6];
            P2_op.y()=spline_para[7];
            P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
            P3_op.x()=spline_para[9];
            P3_op.y()=spline_para[10];
            P3_op.z()=spline_para[11];
            Vector3d comp_V = 3*P3_op*time_persent_[i]*time_persent_[i]+2*time_persent_[i]*P2_op+P1_op; */
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_velocity = spline_velocity_[i];
            //printf("c_V is :%f.%f.%f\n",current_velocity.x(),current_velocity.y(),current_velocity.z());
            //printf("com_V is :%f.%f.%f\n\n",comp_V.x(),comp_V.y(),comp_V.z());
            double _time_persent = time_persent_[i];
            //if(_time_persent<=0.25&&_time_persent>=0.6)
            //{
                cost_funciton = VSplineCeres::Create(current_velocity.x(),current_velocity.y(),current_velocity.z(),_time_persent);
                ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
                problem.AddResidualBlock(cost_funciton,loss_function,spline_para);
            //}
            
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        options.max_num_iterations = 10;
        ceres::Solver::Summary summary;
        ceres::Solve(options,&problem,&summary);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        //Eigen::Vector3d P1_op;
        P1_op.x()=spline_para[3];
        P1_op.y()=spline_para[4];
        P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
        P2_op.x()=spline_para[6];
        P2_op.y()=spline_para[7];
        P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
        P3_op.x()=spline_para[9];
        P3_op.y()=spline_para[10];
        P3_op.z()=spline_para[11];
    }

    double* get_para(){
        return spline_parameter_;
    }

    //private:
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_velocity_;
        const vector<double> time_persent_;
        
};


struct VPSplineCeres{
    public:
        VPSplineCeres(const double pose_x,const double pose_y,const double pose_z,const double velocity_x,const double velocity_y,const double velocity_z,const double time_persent)
         : _pose_x(pose_x),_pose_y(pose_y),_pose_z(pose_z),_velocity_x(velocity_x),_velocity_y(velocity_y),_velocity_z(velocity_z),_time_persent(time_persent){}


        template<typename T>
        bool operator()(
            const T *const spline_parameter,//12维
            T *residual
        )const{
            //v=3a*t^2+2b*t+c
            T com_pose[4];
/*             com_pose[0] = spline_parameter[0]+T(2)*spline_parameter[3]*_time_persent+
            T(3)*spline_parameter[6]*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+T(2)*spline_parameter[4]*_time_persent+
            T(3)*spline_parameter[7]*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+T(2)*spline_parameter[5]*_time_persent+
            T(3)*spline_parameter[8]*_time_persent*_time_persent; */

            com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
            spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
            spline_parameter[11]*_time_persent*_time_persent*_time_persent;

            //velocity
            T com_pose_x = spline_parameter[3]+T(2)*spline_parameter[6]*T(_time_persent)+
            T(3)*spline_parameter[9]*T(_time_persent)*T(_time_persent);

            T com_pose_y = spline_parameter[4]+T(2)*spline_parameter[7]*T(_time_persent)+
            T(3)*spline_parameter[10]*T(_time_persent)*T(_time_persent);

            T com_pose_z = spline_parameter[5]+T(2)*spline_parameter[8]*T(_time_persent)+
            T(3)*spline_parameter[11]*T(_time_persent)*T(_time_persent);

/*             //acc
            T com_pose_x = T(2)*spline_parameter[3]+T(6)*spline_parameter[6]*T(_time_persent);
            T com_pose_y = T(2)*spline_parameter[4]+T(6)*spline_parameter[7]*T(_time_persent);
            T com_pose_z = T(2)*spline_parameter[5]+T(6)*spline_parameter[8]*T(_time_persent); */

            T ab = T(_velocity_x)*T(com_pose_x)+T(_velocity_y)*T(com_pose_y)+T(_velocity_z)*T(com_pose_z);

            T a_norm = sqrt(T(_velocity_x)*T(_velocity_x)+T(_velocity_y)*T(_velocity_y)+T(_velocity_z)*T(_velocity_z));
            T b_norm = sqrt(T(com_pose_x)*T(com_pose_x)+T(com_pose_y)*T(com_pose_y)+T(com_pose_z)*T(com_pose_z));

            com_pose[3] = cos(ab/a_norm*b_norm);


            residual[0]=_pose_x - com_pose[0];
            residual[1]=_pose_y - com_pose[1];
            residual[2]=_pose_z - com_pose[2];

/*             residual[0]=(_velocity_x - com_pose_x);
            residual[1]=(_velocity_y - com_pose_y);
            residual[2]=(_velocity_z - com_pose_z); */
           // printf("\ntrue v is:%f,%f.%f\n",_velocity_x,_velocity_y,_velocity_z);
           // printf("compute v is:%f,%f.%f\n",com_pose_x,com_pose_y,com_pose_z);

            //residual[3]=T(1) - com_pose[3];

            
            return true;
        }

        static ceres::CostFunction *Create(const double _pose_x,const double _pose_y,const double _pose_z,const double _velocity_x,const double _velocity_y,const double _velocity_z,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<VPSplineCeres,3,12>(new VPSplineCeres(_pose_x,_pose_y,_pose_z,_velocity_x,_velocity_y,_velocity_z,_time_persent)));
        }

private:

        double _velocity_x;
        double _velocity_y;
        double _velocity_z;
        double _pose_x;
        double _pose_y;
        double _pose_z;

        double _time_persent;
        
        //vector<double>_data_y;

};

class VPSolveCeres{
    public:
    VPSolveCeres(double *spline_parameter,vector<Eigen::Vector3d>spline_velocity,vector<Eigen::Vector3d>spline_pose,const vector<double> time_persent):
                spline_parameter_(spline_parameter),spline_velocity_(spline_velocity),spline_pose_(spline_pose),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        //printf("spline_velocity_ size is:%d,time_persent_ size is:%d",spline_velocity_.size(),time_persent_.size());
        for(int i=0;i<spline_velocity_.size()&&i<time_persent_.size();i++){
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_velocity = spline_velocity_[i];
            Eigen::Vector3d current_pose = spline_pose_[i];
            double _time_persent = time_persent_[i];
            cost_funciton = VPSplineCeres::Create(current_pose.x(),current_pose.y(),current_pose.z(),current_velocity.x(),current_velocity.y(),current_velocity.z(),_time_persent);
            ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
            problem.AddResidualBlock(cost_funciton,loss_function,spline_para);
            
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options,&problem,&summary);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        P0_op.x()=spline_para[0];
        P0_op.y()=spline_para[1];
        P0_op.z()=spline_para[2];
        //Eigen::Vector3d P1_op;
        P1_op.x()=spline_para[3];
        P1_op.y()=spline_para[4];
        P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
        P2_op.x()=spline_para[6];
        P2_op.y()=spline_para[7];
        P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
        P3_op.x()=spline_para[9];
        P3_op.y()=spline_para[10];
        P3_op.z()=spline_para[11];
    }

    double* get_para(){
        return spline_parameter_;
    }

    //private:
        Eigen::Vector3d P0_op;
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_velocity_;
        vector<Eigen::Vector3d>spline_pose_;
        const vector<double> time_persent_;
        
};

struct VVPSplineCeres{
    public:
        VVPSplineCeres(const double pose_x,const double pose_y,const double pose_z,const double velocity_x,const double velocity_y,const double velocity_z,const double time_persent)
         : _pose_x(pose_x),_pose_y(pose_y),_pose_z(pose_z),_velocity_x(velocity_x),_velocity_y(velocity_y),_velocity_z(velocity_z),_time_persent(time_persent){}


        template<typename T>
        bool operator()(
            const T *const spline_parameter,//12维
            T *residual
        )const{
            //v=3a*t^2+2b*t+c
            T com_pose[4];
/*             com_pose[0] = spline_parameter[0]+T(2)*spline_parameter[3]*_time_persent+
            T(3)*spline_parameter[6]*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+T(2)*spline_parameter[4]*_time_persent+
            T(3)*spline_parameter[7]*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+T(2)*spline_parameter[5]*_time_persent+
            T(3)*spline_parameter[8]*_time_persent*_time_persent; */

            com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
            spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
            spline_parameter[11]*_time_persent*_time_persent*_time_persent;

            //velocity
            T com_pose_x = spline_parameter[3]+T(2)*spline_parameter[6]*T(_time_persent)+
            T(3)*spline_parameter[9]*T(_time_persent)*T(_time_persent);

            T com_pose_y = spline_parameter[4]+T(2)*spline_parameter[7]*T(_time_persent)+
            T(3)*spline_parameter[10]*T(_time_persent)*T(_time_persent);

            T com_pose_z = spline_parameter[5]+T(2)*spline_parameter[8]*T(_time_persent)+
            T(3)*spline_parameter[11]*T(_time_persent)*T(_time_persent);

/*             //acc
            T com_pose_x = T(2)*spline_parameter[3]+T(6)*spline_parameter[6]*T(_time_persent);
            T com_pose_y = T(2)*spline_parameter[4]+T(6)*spline_parameter[7]*T(_time_persent);
            T com_pose_z = T(2)*spline_parameter[5]+T(6)*spline_parameter[8]*T(_time_persent); */

            T ab = T(_velocity_x)*T(com_pose_x)+T(_velocity_y)*T(com_pose_y)+T(_velocity_z)*T(com_pose_z);

            T a_norm = sqrt(T(_velocity_x)*T(_velocity_x)+T(_velocity_y)*T(_velocity_y)+T(_velocity_z)*T(_velocity_z));
            T b_norm = sqrt(T(com_pose_x)*T(com_pose_x)+T(com_pose_y)*T(com_pose_y)+T(com_pose_z)*T(com_pose_z));

            com_pose[3] = cos(ab/a_norm*b_norm);


/*             residual[0]=_pose_x - com_pose[0];
            residual[1]=_pose_y - com_pose[1];
            residual[2]=_pose_z - com_pose[2]; */

            residual[0]=(_velocity_x - com_pose_x);
            residual[1]=(_velocity_y - com_pose_y);
            residual[2]=(_velocity_z - com_pose_z);
            //residual[3]=();
           // printf("\ntrue v is:%f,%f.%f\n",_velocity_x,_velocity_y,_velocity_z);
           // printf("compute v is:%f,%f.%f\n",com_pose_x,com_pose_y,com_pose_z);

            //residual[3]=T(1) - com_pose[3];

            
            return true;
        }

        static ceres::CostFunction *Create(const double _pose_x,const double _pose_y,const double _pose_z,const double _velocity_x,const double _velocity_y,const double _velocity_z,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<VVPSplineCeres,3,12>(new VVPSplineCeres(_pose_x,_pose_y,_pose_z,_velocity_x,_velocity_y,_velocity_z,_time_persent)));
        }

private:

        double _velocity_x;
        double _velocity_y;
        double _velocity_z;
        double _pose_x;
        double _pose_y;
        double _pose_z;

        double _time_persent;
        
        //vector<double>_data_y;

};

class VVPSolveCeres{
    public:
    VVPSolveCeres(double *spline_parameter,vector<Eigen::Vector3d>spline_velocity,vector<Eigen::Vector3d>spline_pose,const vector<double> time_persent):
                spline_parameter_(spline_parameter),spline_velocity_(spline_velocity),spline_pose_(spline_pose),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        //printf("spline_velocity_ size is:%d,time_persent_ size is:%d",spline_velocity_.size(),time_persent_.size());
        for(int i=0;i<spline_velocity_.size()&&i<time_persent_.size();i++){
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_velocity = spline_velocity_[i];
            Eigen::Vector3d current_pose = spline_pose_[i];
            double _time_persent = time_persent_[i];
            cost_funciton = VVPSplineCeres::Create(current_pose.x(),current_pose.y(),current_pose.z(),current_velocity.x(),current_velocity.y(),current_velocity.z(),_time_persent);
            ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
            problem.AddResidualBlock(cost_funciton,loss_function,spline_para);
            
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options,&problem,&summary);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        P0_op.x()=spline_para[0];
        P0_op.y()=spline_para[1];
        P0_op.z()=spline_para[2];
        //Eigen::Vector3d P1_op;
        P1_op.x()=spline_para[3];
        P1_op.y()=spline_para[4];
        P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
        P2_op.x()=spline_para[6];
        P2_op.y()=spline_para[7];
        P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
        P3_op.x()=spline_para[9];
        P3_op.y()=spline_para[10];
        P3_op.z()=spline_para[11];
    }

    double* get_para(){
        return spline_parameter_;
    }

    //private:
        Eigen::Vector3d P0_op;
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_velocity_;
        vector<Eigen::Vector3d>spline_pose_;
        const vector<double> time_persent_;
        
};

struct NomalSplineCeres{
    public:
        NomalSplineCeres(const double pose_x,const double pose_y,const double pose_z,const double time_persent)
                : _pose_x(pose_x),_pose_y(pose_y),_pose_z(pose_z),_time_persent(time_persent){}

        template<typename T>
        bool operator()(
            const T *const spline_parameter,//18维->6个参数，每个参数三维
            T *residual
        )const{
            T com_pose[3];

            com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
                            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
                             spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
                            spline_parameter[11]*_time_persent*_time_persent*_time_persent; 

            //com_pose[3] = spline_parameter[0]+spline_parameter[3]+spline_parameter[6]+spline_parameter[9];

           // com_pose[4] = spline_parameter[1]+spline_parameter[4]+spline_parameter[7]+spline_parameter[10];

            //com_pose[5] = spline_parameter[2]+spline_parameter[5]+spline_parameter[8]+spline_parameter[11]; 

            residual[0] = _pose_x - com_pose[0];
            residual[1] = _pose_y - com_pose[1];
            residual[2] = _pose_z - com_pose[2];
           /*  residual[3]=_pose_x1- com_pose[3];
            residual[4]=_pose_y1- com_pose[4];
            residual[5]=_pose_z1- com_pose[5]; */
            return true;
        }

        static ceres::CostFunction *Create(const double _pose_x,const double _pose_y,const double _pose_z,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<NomalSplineCeres,3,12>(new NomalSplineCeres(_pose_x,_pose_y,_pose_z,_time_persent)));
        }

private:

        double _pose_x;
        double _pose_y;
        double _pose_z;

        double _time_persent;
        
        vector<double>_data_y;

};

class NomalSolveCeres{
    public:
    NomalSolveCeres(double *spline_parameter,vector<Eigen::Vector3d>spline_pose,const vector<double> time_persent):
                spline_parameter_(spline_parameter),spline_pose_(spline_pose),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        for(int i=0;i<spline_pose_.size()&&i<time_persent_.size();i++){
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_pose = spline_pose_[i];
            double _time_persent = time_persent_[i];
            cost_funciton = NomalSplineCeres::Create(current_pose.x(),current_pose.y(),current_pose.z(),_time_persent);
            ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
            problem.AddResidualBlock(cost_funciton,loss_function,spline_para);
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options,&problem,&summary);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        //Eigen::Vector3d P0_op;
        P0_op.x()=spline_para[0];
        P0_op.y()=spline_para[1];
        P0_op.z()=spline_para[2];
        //Eigen::Vector3d P1_op;
        P1_op.x()=spline_para[3];
        P1_op.y()=spline_para[4];
        P1_op.z()=spline_para[5];
        //Eigen::Vector3d P2_op;
        P2_op.x()=spline_para[6];
        P2_op.y()=spline_para[7];
        P2_op.z()=spline_para[8];
        //Eigen::Vector3d P3_op;
        P3_op.x()=spline_para[9];
        P3_op.y()=spline_para[10];
        P3_op.z()=spline_para[11];
/*        // Eigen::Vector3d P4_op;
        P4_op.x()=spline_para[12];
        P4_op.y()=spline_para[13];
        P4_op.z()=spline_para[14];
        //Eigen::Vector3d P5_op;
        P5_op.x()=spline_para[15];
        P5_op.y()=spline_para[16];
        P5_op.z()=spline_para[17]; */
        //std::cout<<summary.FullReport()<<"\n";

    }

    double* get_para(){
        return spline_parameter_;
    }

    //private:
        Eigen::Vector3d P0_op;
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_pose_;
        const vector<double> time_persent_;
        
};

struct VPWSplineCeres{
    public:
        VPWSplineCeres(const double pose_x,const double pose_y,const double pose_z,const double velocity_x,const double velocity_y,const double velocity_z,const Matrix3d rotation,const Vector3d Omega,const double time_persent)
         : _pose_x(pose_x),_pose_y(pose_y),_pose_z(pose_z),_velocity_x(velocity_x),_velocity_y(velocity_y),_velocity_z(velocity_z),Rotation(rotation),Omega_(Omega),_time_persent(time_persent){}


        template<typename T>
        bool operator()(
            const T *const spline_parameter,//18维，前9维是位移信息，后9维是旋转信息
            T *residual
        )const{
            //v=3a*t^2+2b*t+c
            T com_pose[3];
            T com_Rotation_lie[3];
            T com_omega[3];
            T com_Omega_[3];
           
           //提取出变化矩阵中的位移信息，并计算相关量
            com_pose[0] = spline_parameter[0]+T(2)*spline_parameter[3]*_time_persent+
            T(3)*spline_parameter[6]*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+T(2)*spline_parameter[4]*_time_persent+
            T(3)*spline_parameter[7]*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+T(2)*spline_parameter[5]*_time_persent+
            T(3)*spline_parameter[8]*_time_persent*_time_persent;


            com_Rotation_lie[0] = spline_parameter[15]*_time_persent+spline_parameter[12]*_time_persent*_time_persent+
            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_Rotation_lie[1] = spline_parameter[16]*_time_persent+spline_parameter[13]*_time_persent*_time_persent+
            spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_Rotation_lie[2] = spline_parameter[17]*_time_persent+spline_parameter[14]*_time_persent*_time_persent+
            spline_parameter[11]*_time_persent*_time_persent*_time_persent; 

            std::vector<T> com_rotation = compute_R(com_Rotation_lie[0],com_Rotation_lie[1],com_Rotation_lie[2]); 

  /*           com_pose[0] = spline_parameter[0]+spline_parameter[3]*_time_persent+spline_parameter[6]*_time_persent*_time_persent+
            spline_parameter[9]*_time_persent*_time_persent*_time_persent;

            com_pose[1] = spline_parameter[1]+spline_parameter[4]*_time_persent+spline_parameter[7]*_time_persent*_time_persent+
            spline_parameter[10]*_time_persent*_time_persent*_time_persent;

            com_pose[2] = spline_parameter[2]+spline_parameter[5]*_time_persent+spline_parameter[8]*_time_persent*_time_persent+
            spline_parameter[11]*_time_persent*_time_persent*_time_persent; */

/*             com_omega[0] = spline_parameter[12]*_time_persent+spline_parameter[15]*_time_persent*_time_persent+
            spline_parameter[18]*_time_persent*_time_persent*_time_persent;

            com_omega[1] = spline_parameter[13]*_time_persent+spline_parameter[16]*_time_persent*_time_persent+
            spline_parameter[19]*_time_persent*_time_persent*_time_persent;

            com_omega[2] = spline_parameter[14]*_time_persent+spline_parameter[17]*_time_persent*_time_persent+
            spline_parameter[20]*_time_persent*_time_persent*_time_persent; */

            //velocity
/*             T com_pose_x = spline_parameter[0]+T(2)*spline_parameter[3]*T(_time_persent)+
            T(3)*spline_parameter[6]*T(_time_persent)*T(_time_persent);

            T com_pose_y = spline_parameter[1]+T(2)*spline_parameter[4]*T(_time_persent)+
            T(3)*spline_parameter[7]*T(_time_persent)*T(_time_persent);

            T com_pose_z = spline_parameter[2]+T(2)*spline_parameter[5]*T(_time_persent)+
            T(3)*spline_parameter[8]*T(_time_persent)*T(_time_persent); */

            //omega
            

/*             com_omega[0] = spline_parameter[15]*_time_persent+spline_parameter[12]*T(_time_persent)*_time_persent+
            spline_parameter[9]*T(_time_persent)*T(_time_persent)*_time_persent;

            com_omega[1] = spline_parameter[16]*_time_persent+spline_parameter[13]*T(_time_persent)*_time_persent+
            spline_parameter[10]*T(_time_persent)*T(_time_persent)*_time_persent;

            com_omega[2] = spline_parameter[17]*_time_persent+spline_parameter[14]*T(_time_persent)*_time_persent+
            spline_parameter[11]*T(_time_persent)*T(_time_persent)*_time_persent; */

            com_Omega_[0] = spline_parameter[15]+T(2)*spline_parameter[12]*T(_time_persent)+
            T(3)*spline_parameter[9]*T(_time_persent)*T(_time_persent);

            com_Omega_[1] = spline_parameter[16]+T(2)*spline_parameter[13]*T(_time_persent)+
            T(3)*spline_parameter[10]*T(_time_persent)*T(_time_persent);

            com_Omega_[2] = spline_parameter[17]+T(2)*spline_parameter[14]*T(_time_persent)+
            T(3)*spline_parameter[11]*T(_time_persent)*T(_time_persent);

            //T test_cos = cos(com_omega[0]);

           // com_rotation = 

/*             T ab = T(_velocity_x)*T(com_pose_x)+T(_velocity_y)*T(com_pose_y)+T(_velocity_z)*T(com_pose_z);

            T a_norm = sqrt(T(_velocity_x)*T(_velocity_x)+T(_velocity_y)*T(_velocity_y)+T(_velocity_z)*T(_velocity_z));
            T b_norm = sqrt(T(com_pose_x)*T(com_pose_x)+T(com_pose_y)*T(com_pose_y)+T(com_pose_z)*T(com_pose_z));

            com_pose[3] = cos(ab/a_norm*b_norm); */


/*             residual[0]=_pose_x - com_pose[0];
            residual[1]=_pose_y - com_pose[1];
            residual[2]=_pose_z - com_pose[2];
 */
/*             residual[0]  = Rotation(0,0) - com_rotation[0];
            residual[1]  = Rotation(0,1) - com_rotation[1];
            residual[2]  = Rotation(0,2) - com_rotation[2];
            residual[3]  = Rotation(1,0) - com_rotation[3];
            residual[4]  = Rotation(1,1) - com_rotation[4];
            residual[5]  = Rotation(1,2) - com_rotation[5];
            residual[6]  = Rotation(2,0) - com_rotation[6];
            residual[7]  = Rotation(2,1) - com_rotation[7];
            residual[8]  = Rotation(2,2) - com_rotation[8];

            residual[1] *=T(0.1);
            residual[2] *=T(0.1);
            residual[3] *=T(0.1);
            residual[5] *=T(0.1);
            residual[6] *=T(0.1);
            residual[7] *=T(0.1); */
           // residual[3]=T(1) - com_pose[3];

           Vector3d _omega = SO3_::anti_hat(SO3_::log(Rotation));
/*            Matrix3d OR = SO3_::exp(_omega);
           Matrix3d test = OR.transpose()*Rotation;
           printf("test:\n,%f,%f,%f,\n%f,%f,%f,\n%f,%f,%f,\n\n",
            test(0,0),test(0,1),test(0,2),
            test(1,0),test(1,1),test(1,2),
            test(2,0),test(2,1),test(2,2)); */

           // residual[0]=( _omega.x() - com_omega[0] );
           // residual[1]=( _omega.y() - com_omega[1] );
           // residual[2]=( _omega.z() - com_omega[2] );

            residual[0]=( Omega_.x() - com_Omega_[0] );
            residual[1]=( Omega_.y() - com_Omega_[1] );
            residual[2]=( Omega_.z() - com_Omega_[2] );

            return true;
        }

        //根据Rotation的李代数的三维信息计算旋转旋转矩阵信息，用于计算残差方程
        template<typename T>
        static std::vector<T> compute_R(T para_1,T para_2,T para_3) {

            T Normalize = sqrt((para_1*para_1) + (para_2*para_2) + (para_3*para_3));
            T factor_A = sin(Normalize)/Normalize;
            T factor_B = (T(1)-cos(Normalize))/(Normalize*Normalize);

            std::vector<T> com_Rotation(9); // 创建一个指定大小的vector
            // 初始化数组或者填充数据到数组中
            com_Rotation[0] = -factor_B*para_3*para_3 - factor_B*para_2*para_2+ T(1);
            com_Rotation[1] = factor_B*para_1*para_2 - factor_A*para_3;
            com_Rotation[2] = factor_B*para_1*para_2 + factor_A*para_2;

            com_Rotation[3] = factor_B*para_1*para_2 + factor_A*para_3;
            com_Rotation[4] = -factor_B*para_3*para_3 - factor_B*para_1*para_1+ T(1);
            com_Rotation[5] = factor_B*para_2*para_3 - factor_A*para_1;

            com_Rotation[6] = factor_B*para_1*para_3 - factor_A*para_2;
            com_Rotation[7] = factor_B*para_2*para_3 + factor_A*para_1;
            com_Rotation[8] = -factor_B*para_2*para_2 - factor_B*para_1*para_1 + T(1);
            return com_Rotation;
        }


        static ceres::CostFunction *Create(const double _pose_x,const double _pose_y,const double _pose_z,const double _velocity_x,const double _velocity_y,const double _velocity_z,const Matrix3d Rotation,const Vector3d Omega_,const double _time_persent){
            return (new ceres::AutoDiffCostFunction<VPWSplineCeres,3,18>(new VPWSplineCeres(_pose_x,_pose_y,_pose_z,_velocity_x,_velocity_y,_velocity_z,Rotation,Omega_,_time_persent)));
        }

private:

        double _velocity_x;
        double _velocity_y;
        double _velocity_z;
        double _pose_x;
        double _pose_y;
        double _pose_z;
        double _omega_x;
        double _omega_y;
        double _omega_z;
        Matrix3d Rotation;
        Vector3d Omega_;

        double _time_persent;
        
        //vector<double>_data_y;

};

class VPWSolveCeres{
    public:
    VPWSolveCeres(double *spline_parameter,vector<Eigen::Vector3d>spline_velocity,vector<Eigen::Vector3d>spline_pose,vector<Eigen::Matrix3d>spline_rotation,vector<Eigen::Vector3d>spline_omega,const vector<double> time_persent):
                spline_parameter_(spline_parameter),spline_velocity_(spline_velocity),spline_pose_(spline_pose),spline_Rotation(spline_rotation),spline_omega_(spline_omega),time_persent_(time_persent){};

    void Solve(){
        vector<Eigen::Vector3d> op_pose_buffer;
        double *spline_para = this->spline_parameter_;

        //for(int i=0;i<12;i++)
            //printf("befor spline_para is:%f",spline_para[i]);

        ceres::Problem problem;
        //printf("spline_velocity_ size is:%d,time_persent_ size is:%d",spline_velocity_.size(),time_persent_.size());
        for(int i=0;i<spline_velocity_.size()&&i<time_persent_.size();i++){
            ceres::CostFunction *cost_funciton;
            Eigen::Vector3d current_velocity = spline_velocity_[i];
            Eigen::Vector3d current_pose = spline_pose_[i];
            Eigen::Vector3d current_omega = spline_omega_[i];

            Eigen::Matrix3d current_rotation = spline_Rotation[i];
            double _time_persent = time_persent_[i];
            cost_funciton = VPWSplineCeres::Create(current_pose.x(),current_pose.y(),current_pose.z(),current_velocity.x(),current_velocity.y(),current_velocity.z(),current_rotation,current_omega,_time_persent);
            ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
            problem.AddResidualBlock(cost_funciton,NULL,spline_para);
            
        }
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::LinearSolverType::SPARSE_SCHUR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
      //  printf("befor 0 is:%f\n\n",spline_para[0]);
        ceres::Solve(options,&problem,&summary);
       // printf("after 0 is:%f\n\n",spline_para[0]);
        //for(int i=0;i<12;i++)
            //printf("after spline_para is:%f",spline_para[i]);
        //输出完整优化信息
        
        P1_op.x()=spline_para[0];
        P1_op.y()=spline_para[1];
        P1_op.z()=spline_para[2];
        
        P2_op.x()=spline_para[3];
        P2_op.y()=spline_para[4];
        P2_op.z()=spline_para[5];
        
        P3_op.x()=spline_para[6];
        P3_op.y()=spline_para[7];
        P3_op.z()=spline_para[8];
        
        W1_op.x()=spline_para[9];
        W1_op.y()=spline_para[10];
        W1_op.z()=spline_para[11];
        
        W2_op.x()=spline_para[12];
        W2_op.y()=spline_para[13];
        W2_op.z()=spline_para[14];
        
        W3_op.x()=spline_para[15];
        W3_op.y()=spline_para[16];
        W3_op.z()=spline_para[17];

        
    }

    double* get_para(){
        return spline_parameter_;
    }



    //private:
       // Eigen::Vector3d P0_op;
        Eigen::Vector3d P1_op;
        Eigen::Vector3d P2_op;
        Eigen::Vector3d P3_op;

        Eigen::Vector3d W1_op;
        Eigen::Vector3d W2_op;
        Eigen::Vector3d W3_op;

        
        double *spline_parameter_;//spline parameter 0-17分别代表6个参数的三维，0-2是P0的三维，3-5是P1的三维....以此类推
        vector<Eigen::Vector3d>spline_velocity_;
        vector<Eigen::Vector3d>spline_pose_;
        vector<Eigen::Vector3d>spline_omega_;
        vector<Eigen::Matrix3d>spline_Rotation;
        const vector<double> time_persent_;
        
};

