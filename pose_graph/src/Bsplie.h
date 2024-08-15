#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include "sophus/se3.hpp"
#include "so3.h"

class Bspline{
    public:
      Sophus::SE3d spline_trajectory(double u,Sophus::SE3d T_i_before,Sophus::SE3d T_i,Sophus::SE3d T_i_1,Sophus::SE3d T_2){
            Eigen::Matrix4d C;
            C<< 6,0,0,0,
                5,3,-3,1,
                1,3,3,-2,
                0,0,0,1;
            C =C * 0.167; //初始化C
            Eigen::Vector4d vector_u ;
            vector_u<< 1,u,u*u,u*u*u;
            Eigen::Vector4d B_u = C*vector_u;
            Sophus::SE3d caculate_T_1 =Sophus::SE3d::exp((B_u.x()*(T_i_before.inverse()*T_i).log()));
            Sophus::SE3d caculate_T_2 =Sophus::SE3d::exp((B_u.y()*(T_i.inverse()*T_i_1).log()));
            Sophus::SE3d caculate_T_3 =Sophus::SE3d::exp((B_u.z()*(T_i_1.inverse()*T_2).log()));
            Sophus::SE3d T_u = T_i_before * caculate_T_1 *caculate_T_2*caculate_T_3;
      return T_u;

      }; 
};

class Bspline_{
    public:
      Eigen::Vector3d Pose_Bspline(double u,Eigen::Vector3d Pose_0,Eigen::Vector3d Pose_1,Eigen::Vector3d Pose_2,Eigen::Vector3d Pose_3){
            
            Eigen::MatrixXd C(3,4);
            C<<
                5,3,-3,1,
                1,3,3,-2,
                0,0,0,1;
            C = C*0.167;

            Eigen::Vector4d vt;
            vt<<
                  1,u,u*u,u*u*u;

            Eigen::Vector3d Ct = C*vt;

            Eigen::Vector3d Omega_0 = Pose_1 - Pose_0;
            Eigen::Vector3d Omega_1 = Pose_2 - Pose_1;
            Eigen::Vector3d Omega_2 = Pose_3 - Pose_2;
            


            //u就是timepersent ,p_i-1、p_i、p_i+1、p_i+2四个点0.167
            
            Eigen::Vector3d current_C = Pose_0 + (Ct.x()*Omega_0) + (Ct.y()*Omega_1)+(Ct.z()*Omega_2);

      return current_C;

      }; 

      Eigen::Matrix3d Rotation_Bspline(double u,Eigen::Matrix3d Rotation_0,Eigen::Matrix3d Rotation_1,Eigen::Matrix3d Rotation_2,Eigen::Matrix3d Rotation_3){

      Eigen::MatrixXd C(3,4);
            C<<
                5,3,-3,1,
                1,3,3,-2,
                0,0,0,1;
            C =C * 0.167;

            Eigen::Vector4d vt;
            vt<<
                  1,u,u*u,u*u*u;

            Eigen::Vector3d Ct = C*vt;

            Eigen::Vector3d omega0 = SO3_::anti_hat( SO3_::log(Rotation_0) );
            Eigen::Vector3d omega1 = SO3_::anti_hat( SO3_::log(Rotation_1) );
            Eigen::Vector3d omega2 = SO3_::anti_hat( SO3_::log(Rotation_2) );
            Eigen::Vector3d omega3 = SO3_::anti_hat( SO3_::log(Rotation_3) );

            Eigen::Vector3d delta_Omega_0 = SO3_::anti_hat( SO3_::log(Rotation_0.transpose()*Rotation_1) );
            Eigen::Vector3d delta_Omega_1 = SO3_::anti_hat( SO3_::log(Rotation_1.transpose()*Rotation_2) );
            Eigen::Vector3d delta_Omega_2 = SO3_::anti_hat( SO3_::log(Rotation_2.transpose()*Rotation_3) );
            
            Eigen::Vector3d test = (Ct.x()*delta_Omega_0) + (Ct.y()*delta_Omega_1) + (Ct.z()*delta_Omega_2);
            Eigen::Vector3d test_1 = (Ct.x()*delta_Omega_0);
            Eigen::Vector3d test_2 = (Ct.y()*delta_Omega_1);
            Eigen::Vector3d test_3 = (Ct.z()*delta_Omega_2);

            //printf("test1 is:%f,%f.%f\n",delta_Omega_0.x(),delta_Omega_0.y(),delta_Omega_0.z());
            //printf("test2 is:%f,%f.%f\n",delta_Omega_1.x(),delta_Omega_1.y(),delta_Omega_1.z());
            //printf("test3 is:%f,%f.%f\n\n",delta_Omega_2.x(),delta_Omega_2.y(),delta_Omega_2.z());
            Eigen::Matrix3d current_R;
            if(Ct.z()==0)
                   current_R = Rotation_0 * SO3_::exp(test_1)*SO3_::exp(test_2);
            else 
                   current_R = Rotation_0 * SO3_::exp(test_1)*SO3_::exp(test_2)*SO3_::exp(test_3);
            if(current_R.hasNaN()){
            Eigen::Matrix3d test = SO3_::exp(test_3);
            if(test.hasNaN()){
                  printf("%f\n\n",Ct.z());
            }
/*             printf(" THE rotation0 is :\n%f,%f,%f,\n%f,%f,%f,\n%f,%f,%f,\n\n",
            test_1.x(),test_1.y(),test_1.z(),
            test_2.x(),test_2.y(),test_2.z(),
            test_3.x(),test_3.y(),test_3.z()); */
            }
            return current_R;

      }
};