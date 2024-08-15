#include <iostream>
#include "eigen3/Eigen/Core"

class Spline{
    public:
    bool compute_spline_lable(Eigen::Vector3d velocity,Eigen::Vector3d Gy_velocity){
        double index = velocity.squaredNorm();
        double index_G = Gy_velocity.squaredNorm();
        double MAXMUM = 0.05;
        double MINIMUM = 0.005;
        //如果速度（变化过快）或过小（静止状态IMU会漂移）进行优化
        if (index_G>0.1||index >= MAXMUM ||index<=MINIMUM){
            return true;
        }
        else{
            return false;
        }

    };

    
    void compute_spline_argument(Eigen::Vector3d P0,Eigen::Vector3d P1,Eigen::Vector3d v0,Eigen::Vector3d v1){
        //计算P(t)=C3*t^3+C2*t^2+C1*t+C0中的自变量 C
        //根据自变量确定函数
        //先调用这个函数在计算position
        C3 = 2*P0+v0+v1-2*P1;
        C2 = -3*P0-2*v0-v1+3*P1;
        C1 = v0;
        C0 = P0;   
    }
    
    void acc_compute_spline_argument(Eigen::Vector3d P0,Eigen::Vector3d P1,Eigen::Vector3d v0,Eigen::Vector3d v1,Eigen::Vector3d a0,Eigen::Vector3d a1){
        //计算F(t)=P5*t^5+P4*t^4+P3*t^3+P2*t^2+P1*t+P0中的参数P
        //根据自变量确定函数
        //先调用这个函数在计算position
        //Eigen::Vector3d g;
        //g<<9.8,0,3;
        //a0-=g;
        //a1-=g;
        this->P0=P0;
        this->P1=v0;
        this->P2=0.5*a0;
        this->P3=-10*P0-6*v0-1.5*a0+10*P1-4*v1+0.5*a1;
        this->P4=15*P0+8*v0+1.5*a0-15*P1+7*v1-a1;
        this->P5=-6*P0-3*v0-0.5*a0+6*P1-3*v1+0.5*a1;
    }

    void compute_position(double t){
        position = (C3*t*t*t)+(C2*t*t)+(C1*t)+C0;

    }

    void acc_compute_position(double t){
        acc_position = P5*t*t*t*t*t+P4*t*t*t*t+P3*t*t*t+P2*t*t+P1*t+P0;
    }

    Eigen::Vector3d compute_velocity(double t){
        Eigen::Vector3d velocity = 3*C3*t*t+2*C2*t+C1;
        return velocity;

    }

    Eigen::Vector3d compute_Liner_position(double t_persent,Eigen::Vector3d start_pose,Eigen::Vector3d end_pose)
    {

        Eigen::Vector3d final_pose = t_persent*(end_pose - start_pose)+start_pose;
        return final_pose;
    }

/*     void get_acc_spline_para(double* acc_spline_para){
        double para[18] ={P0.x(),P0.y(),P0.z(),P1.x(),P1.y(),P1.z(),P2.x(),P2.y(),P2.z(),P3.x(),P3.y(),P3.z(),
        P4.x(),P4.y(),P4.z(),P5.x(),P5.y(),P5.z()};
        acc_spline_para = para;
        printf("acc_spline_para is:%f\n",acc_spline_para[0]);
        //return acc_spline_para[18];
    } */

    void get_acc_spline_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_,Eigen::Vector3d &P4_,Eigen::Vector3d &P5_){
            P0_ = P0;
            P1_ = P1;
            P2_ = P2;
            P3_ = P3;
            P4_ = P4;
            P5_ = P5;
    }

    void get_spline_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_){
        P0_=C0;
        P1_=C1;
        P2_=C2;
        P3_=C3;
    }

    Eigen::Vector3d correction_spline(Eigen::Vector3d P_start,Eigen::Vector3d P_end,double t){

        //2at^4+(3a+2b)t^3-(3a+2b)t^2+(a+b+c+P_A-P_B)t
        Eigen::Vector3d correction_pose = 2*C3*t*t*t*t+(3*C3+2*C2)*t*t*t-(3*C3+2*C2)*t*t+(C1+C2+C3+P_start-P_end)*t;
        return correction_pose;
    }


    Eigen::Vector3d C0,C1,C2,C3;
    Eigen::Vector3d P0,P1,P2,P3,P4,P5;
    Eigen::Vector3d A,B,C,D,E,F;
    double t0,t1,t2,t3,t4;
    Eigen::Vector3d position;
    Eigen::Vector3d acc_position;
};