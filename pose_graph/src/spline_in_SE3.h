#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include "sophus/se3.hpp"

class SESpline{
    public:
        //SE3 "LERP"
        Sophus::SE3d lerp_in_se(Sophus::SE3d start_pose,Sophus::SE3d end_pose,double delta_t){
            //Eigen::Matrix3d a;
            //Eigen::Vector3d b;
            //Sophus::SE3d start_pos(a,b);
            Eigen::Matrix<double,6,1> alg_start = start_pose.log();
            Eigen::Matrix<double,6,1> alg_end = end_pose.log();

            Eigen::Matrix<double,6,1> alg_com_start = alg_start*(1-delta_t);
            Eigen::Matrix<double,6,1> alg_com_end = alg_end*(delta_t);

            Sophus::SE3d compute_SE3 =Sophus::SE3d::exp(alg_com_start+alg_com_end);
            return compute_SE3;
        }

        //计算q(r)
        Sophus::SE3d Compute_Qr(double r , Sophus::SE3d SE_x , Sophus::SE3d SE_y){
            Eigen::Matrix<double,6,1> alg_start = SE_x.log();
            Eigen::Matrix<double,6,1> alg_end = SE_y.log();
            Sophus::SE3d compute_Qr = Sophus::SE3d::exp((SE_y.log() - SE_x.log())*r);

            return compute_Qr;
        }

        //SE3 "Spline"
        Sophus::SE3d function_c(double t, double r , Sophus::SE3d SE_x , Sophus::SE3d SE_y){
            Sophus::SE3d Qr = Compute_Qr(r,SE_x,SE_y);
            Eigen::Matrix<double,6,1> alg_x = SE_x.log();
            Eigen::Matrix<double,6,1> alg_y = SE_y.log();
            Eigen::Matrix<double,6,1> alg_qr = Qr.log();

            Eigen::Matrix<double,6,1> SE_vector = (alg_x-alg_qr)*(1-t) + t*(alg_y-alg_qr);
            Sophus::SE3d output = Sophus::SE3d::exp(alg_qr+SE_vector);
            return output;
        }
        
        //计算R()
        Sophus::SE3d function_R(Sophus::SE3d p,Eigen::Matrix<double,6,1> v){
            Eigen::Matrix<double,6,1> se_p = p.log();
            Sophus::SE3d output = Sophus::SE3d::exp(se_p+v);
            return output;
        }

        //offline 进程 计算当前分段的qi，wi+，wi+1-
        void Office_Phase(  Sophus::SE3d points_start,Eigen::Matrix<double,6,1>start_vector,Sophus::SE3d points_end,
                            Eigen::Matrix<double,6,1>end_vector,double start_time,double end_time){

            double hi = end_time - start_time;
            Sophus::SE3d pi_plus = function_R(points_start , start_vector*(hi/3));
            Sophus::SE3d pi_1_minus = function_R(points_end , end_vector*(-hi/3));
            Eigen::Matrix<double,6,1> alg_pi = pi_plus.log();
            Eigen::Matrix<double,6,1> alg_pi_1_minus = pi_1_minus.log();
            Eigen::Matrix<double,6,1> compute_pi_v = (alg_pi_1_minus - alg_pi)*0.5;

            qi = function_R(pi_plus,compute_pi_v);
            wi = alg_pi - qi.log();
            wi_1 = alg_pi_1_minus - qi.log();

        };

        Sophus::SE3d Online_Phase(  double start_time,Sophus::SE3d start_pose,Eigen::Matrix<double,6,1> start_velocity,
                                    double end_time,Sophus::SE3d end_pose,Eigen::Matrix<double,6,1> end_velocity,double current_time){
            double hi = end_time-start_time;
            //T is time's persentage
            double T = (current_time - start_time)/hi;
            
            Eigen::Matrix<double,6,1> beta_0_v = (current_time - start_time)/3*start_velocity;
            Sophus::SE3d beta_0 = function_R(start_pose,beta_0_v);

            Eigen::Matrix<double,6,1> beta_1_v = ((1-T) * wi) + (T * wi_1);
            Sophus::SE3d beta_1 = function_R(qi,beta_1_v);

            Eigen::Matrix<double,6,1> beta_2_v = -(current_time - start_time)/3*end_velocity;
            Sophus::SE3d beta_2 = function_R(end_pose,beta_2_v);

            Sophus::SE3d beta_01 = function_c(T,0,beta_0,beta_1);
            Sophus::SE3d beta_12 = function_c(T,1,beta_1,beta_2);

            Sophus::SE3d beta = function_c(T,T,beta_01,beta_12);
            return beta;
        }

        int num_difference = 0;
        double sum_difference = 0;
        double mean_difference = 0;

    private:
        Sophus::SE3d qi;
        Eigen::Matrix<double,6,1> wi;
        Eigen::Matrix<double,6,1> wi_1;


};
