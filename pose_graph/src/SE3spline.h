#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include "sophus/se3.hpp"
#include "./SE3.h"

class SE3spline{
    private:
        Eigen::VectorXd para_P0;
        Eigen::VectorXd para_P1;
        Eigen::VectorXd para_P2;
    
    public:
        void computeSplineParameter(Eigen::VectorXd delta_T,Eigen::VectorXd LieAlg_H,
        Eigen::VectorXd LieAlg_T){
            //LieAlg前三维是位移，后三维是旋转
            //T(t)=T0*Exp(P2*t^3+P1*t^2+P0*t)
            //LieAlg_T_ = J_l(omega)*LieAlg_T
          //  delta_T(0,0)=delta_T(1,1)=delta_T(2,2)=1;
          //  delta_T(0,1)=delta_T(0,2)=delta_T(1,0)=delta_T(1,2)=delta_T(2,0)=delta_T(2,1)=0;
/*           printf("P0 is:%f,%f,%f,%f,\n%f,%f,%f,%f\n%f,%f,%f,%f,\n%f,%f,%f,%f\n\n",
          delta_T(0,0),delta_T(0,1),delta_T(0,2),delta_T(0,3),
          delta_T(1,0),delta_T(1,1),delta_T(1,2),delta_T(1,3),
          delta_T(2,0),delta_T(2,1),delta_T(2,2),delta_T(2,3),
          delta_T(3,0),delta_T(3,1),delta_T(3,2),delta_T(3,3)); */
            Eigen::VectorXd LieAlg_delta_T = delta_T;
           // LieAlg_delta_T(3,0) =  LieAlg_delta_T(4,0) =  LieAlg_delta_T(5,0) =0;
           // LieAlg_H(3,0) =  LieAlg_H(4,0) =  LieAlg_H(5,0) =0;
            //LieAlg_T(3,0) =  LieAlg_T(4,0) =  LieAlg_T(5,0) =0;
            //Eigen::MatrixXd left_jacobain = SE3_::SE3LeftJacobian(LieAlg_delta_T);
            Eigen::MatrixXd left_jacobain = SE3_::SE3RightJacobian(LieAlg_delta_T);
/*             printf("SE 3 jacobain is:\n%f,%f,%f,%f,%f,%f,\n%f,%f,%f,%f,%f,%f,\n%f,%f,%f,%f,%f,%f,\n%f,%f,%f,%f,%f,%f\n%f,%f,%f,%f,%f,%f\n%f,%f,%f,%f,%f,%f\n\n",
            left_jacobain(0,0),left_jacobain(0,1),left_jacobain(0,2),left_jacobain(0,3),left_jacobain(0,4),left_jacobain(0,5),
            left_jacobain(1,0),left_jacobain(1,1),left_jacobain(1,2),left_jacobain(1,3),left_jacobain(1,4),left_jacobain(1,5),
            left_jacobain(2,0),left_jacobain(2,1),left_jacobain(2,2),left_jacobain(2,3),left_jacobain(2,4),left_jacobain(2,5),
            left_jacobain(3,0),left_jacobain(3,1),left_jacobain(3,2),left_jacobain(3,3),left_jacobain(3,4),left_jacobain(3,5),
            left_jacobain(4,0),left_jacobain(4,1),left_jacobain(4,2),left_jacobain(4,3),left_jacobain(4,4),left_jacobain(4,5),
            left_jacobain(5,0),left_jacobain(5,1),left_jacobain(5,2),left_jacobain(5,3),left_jacobain(5,4),left_jacobain(5,5)); */
            //Eigen::VectorXd LieAlg_T_ = left_jacobain.inverse()*LieAlg_T;//left_jacobain.inverse() 应该是不对的
            Eigen::VectorXd LieAlg_T_ = LieAlg_T;

            
            para_P0 = LieAlg_H;
            para_P1 = 3*LieAlg_delta_T-2*LieAlg_H-LieAlg_T;
            para_P2 = -2*LieAlg_delta_T+LieAlg_H+LieAlg_T;
            // para_P1 = 3*LieAlg_delta_T-2*LieAlg_H-LieAlg_T_;
            // para_P2 = -2*LieAlg_delta_T+LieAlg_H+LieAlg_T_;
            // printf("P0 is:%f,%f,%f,%f,%f,%f\n\n",LieAlg_delta_T(0,0),LieAlg_delta_T(1,0),LieAlg_delta_T(2,0),LieAlg_delta_T(3,0),LieAlg_delta_T(4,0),LieAlg_delta_T(5,0));
        }

        Matrix4d compute_trans(double time_persent){
            //T(t)=T0*Exp(P2*t^3+P1*t^2+P0*t)
            Eigen::VectorXd lieAlg_trans = time_persent*time_persent*time_persent*para_P2+para_P1*time_persent*time_persent+para_P0*time_persent;
            Eigen::Matrix4d lie_trans = SE3_::hat(lieAlg_trans);
            Matrix4d Translation = SE3_::exp(lie_trans);
/*             printf("P0 is:%f,%f,%f,%f,%f,%f\n",para_P0(0,0),para_P0(1,0),para_P0(2,0),para_P0(3,0),para_P0(4,0),para_P0(5,0));
            printf("P1 is:%f,%f,%f,%f,%f,%f\n",para_P1(0,0),para_P1(1,0),para_P1(2,0),para_P1(3,0),para_P1(4,0),para_P1(5,0));
            printf("P2 is:%f,%f,%f,%f,%f,%f\n\n",para_P2(0,0),para_P2(1,0),para_P2(2,0),para_P2(3,0),para_P2(4,0),para_P2(5,0)); */
            
            return Translation;
        }

/*          Sophus::SE3::se3 computeInverseRightJacobins(Sophus::SE3::se3 LieAlg_T){
           Sophus::SE3::SE3d Lie_Group_T = Sophus::SE3::exp(LieAlg_T);
           Eigen::Matrix3d rotation_Matrix = Lie_Group_T.rotationMatrix();
           Sophus::SO3d rotation_Matrix_SO3 = rotation_Matrix;
           Sophus::SO3::so3 rotation_Matrix_so3 = rotation_Matrix_SO3.log();
           Sophus::SO3::so3 rotation_Jacobins_Right = Sophus::SO3d::Identity()-
           (1-cos(rotation_Matrix_so3.norm()))/(rotation_Matrix_so3.norm()*rotation_Matrix_so3.norm())*
        }  */

};