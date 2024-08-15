由VINS-MONO (https://github.com/HKUST-Aerial-Robotics/VINS-Mono) 作为基础框架,主要修改了feature_tracker、pose_graph、vins_estimator这三个模块的代码。

修改后的代码实现了位姿图(pose graph)的分段优化:通过剔除部分误差较小的关键帧,在保证轨迹精度的同时有效的提升了后端优化速率。剔除的关键帧数据根据相应的时间信息，通过构建的样条函数（spline）插值计算得到。

其中，构建相机位移信息对应样条函数不仅参考了段落端点的位置信息与速度信息还参考了段中关键帧的相对位移信息，使段落对应的样条函数更加准确。
![image](https://github.com/user-attachments/assets/919eb807-1d59-4e62-ac0a-32bcb80002c2)

构建相机旋转信息对应的样条函数需要先将对应时刻的旋转矩阵转化为对应的李代数，即在李代数空间构建对应的三次样条函数
![image](https://github.com/user-attachments/assets/32ccb742-b3e2-4e36-9f2a-5daf927bde5b)

