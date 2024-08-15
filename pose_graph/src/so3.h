#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>

#pragma once

using namespace Eigen;

class SO3_{
    public:
   static Matrix3d hat(Eigen::Vector3d omega)
    {
        Matrix3d skewSymmetricMatrix;

        skewSymmetricMatrix<<
                                0 , -omega.z() , omega.y(),
                                omega.z() , 0 , -omega.x(),
                                -omega.y() , omega.x() , 0;
        return skewSymmetricMatrix;
    }

       static Vector3d anti_hat(Eigen::Matrix3d OMEGA)
    {
        Vector3d skewSymmetricMatrix;

        skewSymmetricMatrix.x() = OMEGA(2,1);
        skewSymmetricMatrix.y() = -OMEGA(2,0);
        skewSymmetricMatrix.z() = OMEGA(1,0);
        

        return skewSymmetricMatrix;
    }

   static Matrix3d exp(Vector3d lieAlgebra)
    {
        Matrix3d hat_Lie = hat(lieAlgebra);
        Matrix3d Identity = Matrix3d::Identity();
        Matrix3d exponential_1 = ( sin(lieAlgebra.norm())/lieAlgebra.norm() )*hat_Lie;
        Matrix3d exponential_2 = ( 1-cos(lieAlgebra.norm()) ) / ( lieAlgebra.norm()*lieAlgebra.norm() )*hat_Lie*hat_Lie;
        Matrix3d exponential = Identity+exponential_1+exponential_2;
        return exponential;
    }

    static Matrix3d exp(Matrix3d lieAlgebra)
    {
        Matrix3d hat_Lie= lieAlgebra;
        Vector3d man_lie = anti_hat(lieAlgebra);
        Matrix3d Identity = Matrix3d::Identity();
        Matrix3d exponential_1 = ( sin(man_lie.norm())/man_lie.norm() )*hat_Lie;
        Matrix3d exponential_2 = ( 1-cos(man_lie.norm()) ) / ( man_lie.norm()*man_lie.norm() )*hat_Lie*hat_Lie;
        Matrix3d exponential = Identity+exponential_1+exponential_2;
        return exponential;
    }

    static Matrix3d exp_first_order(Vector3d lieAlgebra)
    {
        Matrix3d hat_Lie= hat(lieAlgebra);
        Matrix3d Identity = Matrix3d::Identity();
        Matrix3d exponential = Identity+hat_Lie;
        return exponential;
    }

    static Matrix3d log(Matrix3d LieGroup){
        
        double theta = acos(( LieGroup.trace()-1 )/2);
        //printf("theta is %f\n\n",theta);
        Matrix3d result = ( theta/(2*sin(theta)) )*(LieGroup-LieGroup.transpose());
        return result;
    }

    static Matrix3d right_jacobian(Vector3d lieAlgebra){
        Matrix3d R = hat(lieAlgebra);
        Matrix3d I = Matrix3d::Identity();
        double factor_one = (1-cos(lieAlgebra.norm())) / (lieAlgebra.norm()*lieAlgebra.norm());
        double factor_two =  (lieAlgebra.norm()-sin(lieAlgebra.norm()))/(lieAlgebra.norm()*lieAlgebra.norm()*lieAlgebra.norm());
        Matrix3d result = I-(factor_one*R)+(factor_two*R*R);
        return result;
    }

    static Matrix3d right_jacobian(Matrix3d LieAlgebra){

        Matrix3d R = LieAlgebra;
        Matrix3d I = Matrix3d::Identity();
        Vector3d lieAlgebra = anti_hat(R);
        double factor_one = (1-cos(lieAlgebra.norm())) / (lieAlgebra.norm()*lieAlgebra.norm());
        double factor_two =  (lieAlgebra.norm()-sin(lieAlgebra.norm()))/(lieAlgebra.norm()*lieAlgebra.norm()*lieAlgebra.norm());
        Matrix3d result = I-(factor_one*R)+(factor_two*R*R);
        return result;
    }

    static Matrix3d right_jacobian_inverse(Matrix3d LieAlgebra){
        Matrix3d R = LieAlgebra;
        Vector3d lieAlgebra = anti_hat(R);
        Matrix3d I = Matrix3d::Identity();
        double factor_one = 0.5;
        double factor_two = 1/(lieAlgebra.norm()*lieAlgebra.norm())-( ( 1+cos(lieAlgebra.norm()) ) / ( 2*lieAlgebra.norm()*sin(lieAlgebra.norm()) ) );
        Matrix3d result = I+(factor_one*R)+(factor_two*R*R);
/*         printf("so3 jacobain is:\n%f,%f,%f,\n%f,%f,%f,\n%f,%f,%f\n\n",
            result(0,0),result(0,1),result(0,2),
            result(1,0),result(1,1),result(1,2),
            result(2,0),result(2,1),result(2,2)); */
        return result;
    }

    static Matrix3d left_jacobian(Matrix3d LieAlgebra){
        Matrix3d R = LieAlgebra;
        Vector3d lieAlgebra = anti_hat(R);
        Matrix3d I = Matrix3d::Identity();
        double factor_one = (1-cos(lieAlgebra.norm()))/(lieAlgebra.norm()*lieAlgebra.norm());
        double factor_two = (lieAlgebra.norm()-sin(lieAlgebra.norm()))/(lieAlgebra.norm()*lieAlgebra.norm()*lieAlgebra.norm());
        Matrix3d result = I+(factor_one*R)+(factor_two*R*R);
        return result;
    }

    static Matrix3d left_jacobian_inverse(Matrix3d LieAlgebra){
        Matrix3d R = LieAlgebra;
        Vector3d lieAlgebra = anti_hat(R);
        Matrix3d I = Matrix3d::Identity();
        double factor_one = 0.5;
        double factor_two = ( 1/( lieAlgebra.norm()*lieAlgebra.norm() ) )-
                            ( ( 1+cos(lieAlgebra.norm()) ) )/( 2*lieAlgebra.norm()*sin(lieAlgebra.norm()) );
        Matrix3d result = I-(factor_one*R)+(factor_two*R*R);
        return result;
    }

    static Matrix3d SLERP(double t,Matrix3d R1,Matrix3d R2)
    {
       
        Matrix3d deltaR = R1.transpose()*R2;
        Matrix3d r1 =t*log(deltaR);
        Matrix3d R3 = R1*exp(r1);
        return R3;
    }

    

};
