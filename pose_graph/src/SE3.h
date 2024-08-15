#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include "so3.h"

using namespace Eigen;

class SE3_{
    public:
  static Matrix4d setSE3(Vector3d pose,Matrix3d rotation)
    {
        Matrix4d T;
        T<< 
            rotation(0,0),rotation(0,1),rotation(0,2),pose.x(),
            rotation(1,0),rotation(1,1),rotation(1,2),pose.y(),
            rotation(2,0),rotation(2,1),rotation(2,2),pose.z(),
            0,0,0,1;

        return T;
    }

   static Matrix4d exp(Matrix4d hat_ksi)
    {
        Matrix3d so3 ;
        so3<<
                hat_ksi(0,0),hat_ksi(0,1),hat_ksi(0,2),
                hat_ksi(1,0),hat_ksi(1,1),hat_ksi(1,2),
                hat_ksi(2,0),hat_ksi(2,1),hat_ksi(2,2);
        Matrix3d Rotation = SO3_::exp(so3);
        Matrix4d T_out;
        Matrix3d Jacobi = SO3_::left_jacobian(so3);
        Vector3d trans_rou ;
        trans_rou<<
                    hat_ksi(0,3),hat_ksi(1,3),hat_ksi(2,3);
        Vector3d trans_t = Jacobi*trans_rou;
        T_out<< 
                Rotation(0,0),Rotation(0,1),Rotation(0,2),trans_rou.x(),
                Rotation(1,0),Rotation(1,1),Rotation(1,2),trans_rou.y(),
                Rotation(2,0),Rotation(2,1),Rotation(2,2),trans_rou.z(),
                0,0,0,1;
        return T_out;
    }   

    static VectorXd log(Matrix4d trasnform)
    {
        //李代数，前三维是“位移”，后三维是旋转
        VectorXd ksi(6);
        Vector3d displacement;
        Matrix3d rotation;
        displacement<<
                        trasnform(0,3),trasnform(1,3),trasnform(2,3);
        rotation<<
                    trasnform(0,0),trasnform(0,1),trasnform(0,2),
                    trasnform(1,0),trasnform(1,1),trasnform(1,2),
                    trasnform(2,0),trasnform(2,1),trasnform(2,2);
        Matrix3d LieAlg_R_Matrix = SO3_::log(rotation);

        Vector3d LieAlg_R_Vertex = SO3_::anti_hat(LieAlg_R_Matrix);
        Matrix3d SO3_Jacoabin_inverse = SO3_::left_jacobian_inverse(LieAlg_R_Matrix);
        Vector3d displacement_rou = SO3_Jacoabin_inverse * displacement;
         ksi<<
                displacement_rou.x(),displacement_rou.y(),displacement_rou.z(),
                LieAlg_R_Vertex.x(),LieAlg_R_Vertex.y(),LieAlg_R_Vertex.z(); 
/*          ksi<<
                displacement.x(),displacement.y(),displacement.z(),
                0,0,0;  */
        return ksi;
    }

    static Matrix4d hat(VectorXd ksi)
    {
        //李代数，前三维是“位移”，后三维是旋转
        Vector3d rotation_phi;
        rotation_phi << 
                        ksi(3,0),ksi(4,0),ksi(5,0);

        Matrix3d hat_SO3 = SO3_::hat(rotation_phi);
        Matrix4d ksi_out;
        ksi_out<<
                  hat_SO3(0,0),hat_SO3(0,1),hat_SO3(0,2),ksi(0,0),
                  hat_SO3(1,0),hat_SO3(1,1),hat_SO3(1,2),ksi(1,0),
                  hat_SO3(2,0),hat_SO3(2,1),hat_SO3(2,2),ksi(2,0),
                  0,0,0,0;
        return ksi_out;
    }

    static MatrixXd SE3LeftJacobian(VectorXd ksi){
        //6维李代数，前三维是“位移”，后三维是旋转
        Vector3d omega = {ksi(3,0),ksi(4,0),ksi(5,0)};
        Vector3d rou = {ksi(0,0),ksi(1,0),ksi(2,0)};
        //李代数的矩阵形式
        Matrix3d Matrix_omega = SO3_::hat(omega);
        //"位移"的反对称矩阵？？？
        Matrix3d Matrix_rou = SO3_::hat(rou);

        Matrix3d SO3_left_jacobain = SO3_::left_jacobian(Matrix_omega);
        //SE3的雅克比矩阵是9*9的
        MatrixXd SE3_left_jacobain(6,6);

        double factor_one = 0.5;
        double factor_two = (omega.norm()-sin(omega.norm()))/(omega.norm()*omega.norm()*omega.norm());
        double factor_three =(omega.norm()*omega.norm()+2*cos(omega.norm())-2)/(2*omega.norm()*omega.norm()*omega.norm()*omega.norm());
        double factor_four =(2*omega.norm()-3*sin(omega.norm())+omega.norm()*cos(omega.norm()))/(2*omega.norm()*omega.norm()*omega.norm()*omega.norm()*omega.norm());

        Matrix3d arguement_one = Matrix_rou;
        Matrix3d arguement_two = Matrix_omega*Matrix_rou + Matrix_rou*Matrix_omega + Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d arguement_three =Matrix_omega*Matrix_omega*Matrix_rou + Matrix_rou*Matrix_omega*Matrix_omega - 3*Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d arguement_four = Matrix_omega*Matrix_rou*Matrix_omega*Matrix_omega+Matrix_omega*Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d block_up_right = factor_one*arguement_one+factor_two*arguement_two+factor_three*arguement_three+factor_three*arguement_three;


        SE3_left_jacobain<<
                            SO3_left_jacobain(0,0),SO3_left_jacobain(0,1),SO3_left_jacobain(0,2),block_up_right(0,0),block_up_right(0,1),block_up_right(0,2),
                            SO3_left_jacobain(1,0),SO3_left_jacobain(1,1),SO3_left_jacobain(1,2),block_up_right(1,0),block_up_right(1,1),block_up_right(1,2),
                            SO3_left_jacobain(2,0),SO3_left_jacobain(2,1),SO3_left_jacobain(2,2),block_up_right(2,0),block_up_right(2,1),block_up_right(2,2),
                            0,0,0,SO3_left_jacobain(0,0),SO3_left_jacobain(0,1),SO3_left_jacobain(0,2),
                            0,0,0,SO3_left_jacobain(1,0),SO3_left_jacobain(1,1),SO3_left_jacobain(1,2),
                            0,0,0,SO3_left_jacobain(2,0),SO3_left_jacobain(2,1),SO3_left_jacobain(2,2); 

        return SE3_left_jacobain;

    }

    static MatrixXd SE3RightJacobian(VectorXd ksi){
        ksi = -1*ksi;
        //6维李代数，前三维是“位移”，后三维是旋转
        Vector3d omega = {ksi(3,0),ksi(4,0),ksi(5,0)};
        Vector3d rou = {ksi(0,0),ksi(1,0),ksi(2,0)};
        //李代数的矩阵形式
        Matrix3d Matrix_omega = SO3_::hat(omega);
        //"位移"的反对称矩阵？？？
        Matrix3d Matrix_rou = SO3_::hat(rou);

        Matrix3d SO3_left_jacobain = SO3_::left_jacobian(Matrix_omega);
        //SE3的雅克比矩阵是9*9的
        MatrixXd SE3_left_jacobain(6,6);

        double factor_one = 0.5;
        double factor_two = (omega.norm()-sin(omega.norm()))/(omega.norm()*omega.norm()*omega.norm());
        double factor_three =(omega.norm()*omega.norm()+2*cos(omega.norm())-2)/(2*omega.norm()*omega.norm()*omega.norm()*omega.norm());
        double factor_four =(2*omega.norm()-3*sin(omega.norm())+omega.norm()*cos(omega.norm()))/(2*omega.norm()*omega.norm()*omega.norm()*omega.norm()*omega.norm());

        Matrix3d arguement_one = Matrix_rou;
        Matrix3d arguement_two = Matrix_omega*Matrix_rou + Matrix_rou*Matrix_omega + Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d arguement_three =Matrix_omega*Matrix_omega*Matrix_rou + Matrix_rou*Matrix_omega*Matrix_omega - 3*Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d arguement_four = Matrix_omega*Matrix_rou*Matrix_omega*Matrix_omega+Matrix_omega*Matrix_omega*Matrix_rou*Matrix_omega; 
        Matrix3d block_up_right = factor_one*arguement_one+factor_two*arguement_two+factor_three*arguement_three+factor_three*arguement_three;


        SE3_left_jacobain<<
                            SO3_left_jacobain(0,0),SO3_left_jacobain(0,1),SO3_left_jacobain(0,2),block_up_right(0,0),block_up_right(0,1),block_up_right(0,2),
                            SO3_left_jacobain(1,0),SO3_left_jacobain(1,1),SO3_left_jacobain(1,2),block_up_right(1,0),block_up_right(1,1),block_up_right(1,2),
                            SO3_left_jacobain(2,0),SO3_left_jacobain(2,1),SO3_left_jacobain(2,2),block_up_right(2,0),block_up_right(2,1),block_up_right(2,2),
                            0,0,0,SO3_left_jacobain(0,0),SO3_left_jacobain(0,1),SO3_left_jacobain(0,2),
                            0,0,0,SO3_left_jacobain(1,0),SO3_left_jacobain(1,1),SO3_left_jacobain(1,2),
                            0,0,0,SO3_left_jacobain(2,0),SO3_left_jacobain(2,1),SO3_left_jacobain(2,2); 

        return SE3_left_jacobain;

    }

    static Eigen::Matrix3d toRotationMatrix(Eigen::Matrix4d Translation){
        Eigen::Matrix3d Rotation;
        Rotation<<
                    Translation(0,0),Translation(0,1),Translation(0,2),
                    Translation(1,0),Translation(1,1),Translation(1,2),
                    Translation(2,0),Translation(2,1),Translation(2,2);
        return Rotation;
    }

    static Eigen::Vector3d toMovementVector(Eigen::Matrix4d Translation){
        Eigen::Vector3d Movement;
        Movement<<
                    Translation(0,3),Translation(1,3),Translation(2,3);
        return Movement;
    }

};