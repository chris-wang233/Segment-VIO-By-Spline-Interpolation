#include <eigen3/Eigen/Dense>
#include <iostream>
#include <map>

class Spline_Least_Square{
    public:

        Spline_Least_Square(std::map<double,Eigen::Vector3d> time_pose_,Eigen::Vector3d start_pose_,
        Eigen::Vector3d end_pose_,Eigen::Vector3d start_velocity_,Eigen::Vector3d end_velocity_):
        time_pose(time_pose_),start_pose(start_pose_),end_pose(end_pose_),start_velocity(start_velocity_),
        end_velocity(end_velocity_),para_factor(4,6){};

        void compute_para(){
            //std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();
            P0=para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity+para_factor(0,4)*time_pose[0]+para_factor(0,5)*time_pose[1];
            P1=para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity+para_factor(1,4)*time_pose[0]+para_factor(1,5)*time_pose[1];
            P2=para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity+para_factor(2,4)*time_pose[0]+para_factor(2,5)*time_pose[1];
            P3=para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity+para_factor(3,4)*time_pose[0]+para_factor(3,5)*time_pose[1];
        }

        //利用最小二乘法计算矩阵系数，即(A^T*A)^(-1)*A^T
        void compute_para_factor(){
            std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();
            Eigen::MatrixXd A(6,4);
            A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3,
                0,0,2,6,
                0,0,2,0;
               // 1,it->first,it->first*it->first,it->first*it->first*it->first;

            Eigen::Matrix4d output_first = A.transpose()*A;
            para_factor = output_first.inverse()*A.transpose();
            compute_para();
        }

        //赋值
        void get_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_){
            P0_ = P0;
            P1_ = P1;
            P2_ = P2;
            P3_ = P3;

        }


    private:
        //map 键值为time_persent，内容为三维pose信息
        std::map<double,Eigen::Vector3d> time_pose;
        //端点信息
        Eigen::Vector3d start_pose;
        Eigen::Vector3d end_pose;
        Eigen::Vector3d start_velocity;
        Eigen::Vector3d end_velocity;
        Eigen::MatrixXd para_factor;
        Eigen::Vector3d P0;
        Eigen::Vector3d P1;
        Eigen::Vector3d P2;
        Eigen::Vector3d P3;
};

//origenal
class Spline_Least_Square_muti_{
    public:

        Spline_Least_Square_muti_(std::map<double,Eigen::Vector3d> time_pose_,Eigen::Vector3d start_pose_,
        Eigen::Vector3d end_pose_,Eigen::Vector3d start_velocity_,Eigen::Vector3d end_velocity_,double time):
        time_pose(time_pose_),start_pose(start_pose_),end_pose(end_pose_),start_velocity(start_velocity_),
        end_velocity(end_velocity_),para_factor(4,4+time_pose.size()),t(time){};

        void compute_para(){
            std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();
/*             P0=para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity+para_factor(0,4)*it->second;
            P1=para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity+para_factor(1,4)*it->second;
            P2=para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity+para_factor(2,4)*it->second;
            P3=para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity+para_factor(3,4)*it->second; */

            P0 =para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity;
            P1 =para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity;
            P2 =para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity;
            P3 =para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity;
            int i=1;
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
                P0 += para_factor(0,3+i)*it->second;
                P1 += para_factor(1,3+i)*it->second;
                P2 += para_factor(2,3+i)*it->second;
                P3 += para_factor(3,3+i)*it->second;
                i++;
               // printf("P0 is: %f,%f,%f\n",P0.x(),P0.y(),P0.z());
            }
            //printf("end loop\n\n");
        }

        //利用最小二乘法计算矩阵系数，即(A^T*A)^(-1)*A^T
        void compute_para_factor(){
            
            Eigen::MatrixXd A(4+time_pose.size(),4);

/*             A<<
                1,0,0,0,
                1,t,t^2,t^3,
                0,1,0,0,
                0,1,2*t,3*t^2; */

            A(0,0)=1;
            A(0,1)=A(0,2)=A(0,3)=0;

            //A(1,0)=A(1,1)=A(1,2)=A(1,3)=1;
            A(1,0)=1;
            A(1,1)=t;
            A(1,2)=t*t;
            A(1,3)=t*t*t;

            A(2,0)=A(2,2)=A(2,3)=0;
            A(2,1)=1;

            A(3,0)=0;
            A(3,1)=1;
            A(3,2)=2*t;
            A(3,3)=3*t*t;
           // printf("hey hey damm\n\n");
            int i =1;
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
            
                A(3+i,0)=1;
                A(3+i,1)=it->first;
                A(3+i,2)=it->first*it->first;
                A(3+i,3)=it->first*it->first*it->first;
                i++;
/*             A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3,
                1,it->first,it->first*it->first,it->first*it->first*it->first; */


            }
           Eigen::Matrix4d output_first = A.transpose()*A;
           para_factor = output_first.completeOrthogonalDecomposition().pseudoInverse()*A.transpose();
            compute_para();

        }

        //赋值
        void get_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_){
            P0_ = P0;
            P1_ = P1;
            P2_ = P2;
            P3_ = P3;

        }


    private:
        //map 键值为time_persent，内容为三维pose信息
        std::map<double,Eigen::Vector3d> time_pose;
        //端点信息
        Eigen::Vector3d start_pose;
        Eigen::Vector3d end_pose;
        Eigen::Vector3d start_velocity;
        Eigen::Vector3d end_velocity;
        Eigen::MatrixXd para_factor;
        Eigen::Vector3d P0;
        Eigen::Vector3d P1;
        Eigen::Vector3d P2;
        Eigen::Vector3d P3;
        double t;
};

//******************************************
class Spline_Least_Square_muti{
    public:

        Spline_Least_Square_muti(std::map<double,Eigen::Vector3d> time_pose_,std::map<double,Eigen::Vector3d> time_velocitry_,Eigen::Vector3d start_pose_,
        Eigen::Vector3d end_pose_,Eigen::Vector3d start_velocity_,Eigen::Vector3d end_velocity_,double time):
        time_pose(time_pose_),time_velocity(time_velocitry_),start_pose(start_pose_),end_pose(end_pose_),start_velocity(start_velocity_),
        end_velocity(end_velocity_),para_factor(4,4+time_pose.size()+time_velocity.size()),t(time){};

        void compute_para(){
            std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();
/*             P0=para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity+para_factor(0,4)*it->second;
            P1=para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity+para_factor(1,4)*it->second;
            P2=para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity+para_factor(2,4)*it->second;
            P3=para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity+para_factor(3,4)*it->second; */

            P0 =para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity;
            P1 =para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity;
            P2 =para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity;
            P3 =para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity;

            int i=0;
            std::map<double,Eigen::Vector3d>::iterator it_=time_velocity.begin();
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
                P0 += para_factor(0,4+2*i)*it->second;
                P1 += para_factor(1,4+2*i)*it->second;
                P2 += para_factor(2,4+2*i)*it->second;
                P3 += para_factor(3,4+2*i)*it->second;

                P0 += para_factor(0,4+2*i+1)*it_->second;
                P1 += para_factor(1,4+2*i+1)*it_->second;
                P2 += para_factor(2,4+2*i+1)*it_->second;
                P3 += para_factor(3,4+2*i+1)*it_->second;

                i++;
                it_++;
               //printf("P0 is: %f,%f,%f\n",P0.x(),P0.y(),P0.z());
            }
            //printf("end loop\n\n");
        }

        //利用最小二乘法计算矩阵系数，即(A^T*A)^(-1)*A^T
        void compute_para_factor(){
            
            Eigen::MatrixXd A(4+time_pose.size()+time_velocity.size(),4);

/*             A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3; */

            A(0,0)=1;
            A(0,1)=A(0,2)=A(0,3)=0;

            A(1,0)=1;
            A(1,1)=t;
            A(1,2)=t*t;
            A(1,3)=t*t*t;

            A(2,0)=A(2,2)=A(2,3)=0;
            A(2,1)=1;

            A(3,0)=0;
            A(3,1)=1;
            A(3,2)=2*t;
            A(3,3)=3*t*t;
           // printf("hey hey damm\n\n");

            int i =0;
            std::map<double,Eigen::Vector3d>::iterator it_=time_velocity.begin();
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
            
                A(4+2*i,0)=1;
                A(4+2*i,1)=it->first;
                A(4+2*i,2)=it->first*it->first;
                A(4+2*i,3)=it->first*it->first*it->first;

                A(4+2*i+1,0)=0;
                A(4+2*i+1,1)=1;
                A(4+2*i+1,2)=2*it_->first;
                A(4+2*i+1,3)=3*it_->first*it_->first;

                i++;
                it_++;
/*             A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3,
                1,it->first,it->first*it->first,it->first*it->first*it->first; */


            }
           // double lambda = 0.001; // 正则化项强度，根据实际情况选择合适的值
            //Eigen::Matrix4d output_first = A.transpose()*A;
            //Eigen::Matrix4d regularized_output = output_first + lambda * Eigen::Matrix4d::Identity();
            //Eigen::MatrixXd para_factor = regularized_output.inverse() * A.transpose();
           Eigen::Matrix4d output_first = A.transpose()*A;
           para_factor = output_first.completeOrthogonalDecomposition().pseudoInverse()*A.transpose();
            compute_para();

        }

        //赋值
        void get_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_){
            P0_ = P0;
            P1_ = P1;
            P2_ = P2;
            P3_ = P3;

        }


    private:
        //map 键值为time_persent，内容为三维pose信息
        std::map<double,Eigen::Vector3d> time_pose;
        std::map<double,Eigen::Vector3d> time_velocity;
        //端点信息
        Eigen::Vector3d start_pose;
        Eigen::Vector3d end_pose;
        Eigen::Vector3d start_velocity;
        Eigen::Vector3d end_velocity;
        Eigen::MatrixXd para_factor;
        Eigen::Vector3d P0;
        Eigen::Vector3d P1;
        Eigen::Vector3d P2;
        Eigen::Vector3d P3;
        double t;
};