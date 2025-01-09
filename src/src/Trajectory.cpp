#include <ros/ros.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "nav_msgs/Path.h"
#include <chrono>
// #include <qpOASES.hpp>

using namespace Eigen;
// 定义路径格式
typedef std::vector<Eigen::Vector2d> Path;

// Step 4: 自行实现轨迹生成类
class TrajectoryGenerator {
    // 轨迹规划的目标是根据A*算法给出的无碰撞路径，计算轨迹航点（控制点），从而生成一条以时间参数化的平滑轨迹，可以用于控制移动机器人跟踪
    // 本次作业中，我们要求生成一条分段多项式轨迹，即每段轨迹均为一个多项式函数
    // 你可以选择使用多项式、B样条、贝塞尔曲线、MINCO等多种轨迹基函数实现
    // 每段轨迹的连接处需要满足一定的连续性条件，如位置连续、速度连续、加速度连续等，这将构成轨迹优化的主要约束条件
    // 轨迹的初始和终止状态为到达指定位置，速度、加速度等状态均为0
    // 优化目标可以是最小化轨迹的加加速度（jerk）或加加加速度（snap），请自行选择合适的优化目标及多项式阶数
    // 本次作业对轨迹段的时间选择不做进一步要求，可以自行选择固定时间或自定义策略生成时间分配
    // 可以任意选用求解器，如qpOASES、OSQP-Eigen等，也可以自行实现闭式解法
public:
    TrajectoryGenerator() = default;

        // your code
    Path GeneratePolyTrajectory(Path path, VectorXd time, double duration, int minimum_order = 3) {
        Path trajectory;
        VectorXd PathX = VectorXd::Zero(path.size());
        VectorXd PathY = VectorXd::Zero(path.size());
        for (int i = 0; i < path.size(); i++) {
            PathX[i] = path[i][0];
            PathY[i] = path[i][1];
        }
        VectorXd poly_coef_x = PolyCoefficientsMinimum1D(PathX, time, minimum_order);
        VectorXd poly_coef_y = PolyCoefficientsMinimum1D(PathY, time, minimum_order);
        VectorXd generated_path_x = GeneratePolyPath1D(poly_coef_x, time, duration);
        VectorXd generated_path_y = GeneratePolyPath1D(poly_coef_y, time, duration);
        for (int i = 0; i < generated_path_x.size(); i++) {
            trajectory.push_back(Vector2d(generated_path_x[i], generated_path_y[i]));
        }
        return trajectory;
    }
    // 无约束贝塞尔曲线
    Path GenerateBezierTrajectory(Path control_points) {
        Path trajectory;
        int n = control_points.size() - 1;
        for (double t = 0; t <= 1 + 1e-6; t += 0.01) {
            Vector2d point = Vector2d::Zero();
            for (int i = 0; i <= n; i++) {
                point += control_points[i] * Combination(n, i) * pow(1 - t, n - i) * pow(t, i);
            }
            trajectory.push_back(point);
            std::cout << point.transpose() << std::endl;
        }
        return trajectory;
    }
    // 梯形时间规划，假设最大加速度和最大速度都为1
    VectorXd TimeAllocation(Path path) {
        VectorXd time = VectorXd::Ones(path.size() - 1);
        double max_velocity = 1.0, max_acceleration = 1.0;
        double t = max_velocity / max_acceleration;
        double distance_thereshold = 0.5 * max_velocity * t * 2;
        for (int i = 0; i < path.size() - 1; i++) {
            double distance = (path[i + 1] - path[i]).norm();
            if (distance < distance_thereshold) {
                time[i] = sqrt(distance / max_acceleration) * 2;
            } else {
                time[i] = 2 * t + (distance - distance_thereshold) / max_velocity;
            }
        }
        return time;
    }
private:
    int Factorial(int n) {
        int result = 1;
        for (int i = 1; i <= n; i++) {
            result *= i;
        }
        return result;
    }
    double Combination(int n, int k) {
        double result = 1.0;
        for (int i = 1; i <= k; i++) {
            result *= (n - i + 1) / static_cast<double>(i);
        }
        return result;
    }   
    
    double CalculatePolyValue(VectorXd poly_coef, double t) {
        double value = 0;
        for (int i = 0; i < poly_coef.size(); i++) {
            value += poly_coef[i] * pow(t, i);
        }
        return value;
    }
    // closed-form solution  
    // minimum_order: 最小阶数, 3阶时优化目标为最小化jerk，4阶时优化目标为最小化snap
    VectorXd PolyCoefficientsMinimum1D(VectorXd position, VectorXd time, int minimum_order = 3) {
        const int num_segments = position.size() - 1;
        const int poly_order = minimum_order * 2 - 1;
        const int num_poly_coeffs = poly_order + 1;
        MatrixXd M = GenerateM(num_segments, minimum_order, time);
        MatrixXd Ct = GenerateCt(num_segments, minimum_order);
        MatrixXd Q = GenerateQ(num_segments, minimum_order, time);
        MatrixXd R = Ct.transpose() * M.inverse().transpose() * Q * M.inverse() * Ct;
        int num_df = minimum_order * 2 + (num_segments - 1);
        int num_dp = (num_segments -1) * (minimum_order - 1);
        MatrixXd R_pp = R.bottomRightCorner(num_dp, num_dp);
        MatrixXd R_fp = R.topRightCorner(num_df, num_dp);
        VectorXd dF = VectorXd::Zero(num_df);
        VectorXd dP = VectorXd::Zero(num_dp);
        dF[0] = position[0];
        for (int i = 1; i <= num_segments; i++) {
            dF[minimum_order + i - 1] = position[i];
        }
        dP = -R_pp.inverse() * R_fp.transpose() * dF;
        VectorXd d = VectorXd::Zero(num_df + num_dp);
        d << dF, dP;
        VectorXd poly_coef = M.inverse() * Ct * d;
        return poly_coef;
    }

    // 生成映射矩阵
    MatrixXd GenerateM(int num_segments, int minimul_order, VectorXd time) {
        int num_poly_coeffs = minimul_order * 2;
        MatrixXd M = MatrixXd::Zero(num_segments * minimul_order * 2, num_segments * num_poly_coeffs);
        for (int i = 0; i < num_segments; i++) {
            for (int j = 0; j < minimul_order; j++) {
                for (int k = 0; k < num_poly_coeffs; k++) {
                    if (k < j){
                        continue;
                    }
                    M(i * minimul_order * 2 + j, i * num_poly_coeffs + k) = Factorial(k) / Factorial(k - j) * pow(0, k - j);
                    M(i * minimul_order * 2 + j + minimul_order, i * num_poly_coeffs + k) = Factorial(k) / Factorial(k - j) * pow(time[i], k - j);
                }
            }
        }
        return M;
    }
    // 生成置换矩阵
    MatrixXd GenerateCt(int num_segments, int minimul_order) {
        int ct_rows = num_segments * minimul_order * 2;
        int ct_cols = num_segments * minimul_order * 2 - (num_segments - 1) * minimul_order;
        int num_fixed = minimul_order * 2 + (num_segments -1);
        MatrixXd Ct = MatrixXd::Zero(ct_rows, ct_cols);
        for (int i = 0; i < ct_rows; i++) {
            if (i < minimul_order) {
                Ct(i, i) = 1;
            } 
            else if (i >= ct_rows - minimul_order) {
                Ct(i, num_fixed - minimul_order + (i - ct_rows + minimul_order)) = 1;
            }
            else if (((i % minimul_order) == 0) && (i / minimul_order % 2 == 1)) {
                Ct(i, minimul_order + (i / minimul_order)/2) = 1;
            }
            else if (((i % minimul_order) == 0) && (i / minimul_order % 2 == 0)) {
                Ct(i, minimul_order + (i / minimul_order)/2 - 1) = 1;
            }
            else if (((i % minimul_order) != 0) && (i / minimul_order % 2 == 1)) {
                Ct(i, num_fixed + (i % minimul_order - 1) + (i / minimul_order)/2 * (minimul_order - 1)) = 1;
            }
            else if (((i % minimul_order) != 0) && (i / minimul_order % 2 == 0)) {
                Ct(i, num_fixed + (i % minimul_order - 1) + (((i / minimul_order)/2 - 1) * (minimul_order - 1))) = 1;
            }
        }
        return Ct;
    }
    // 生成目标函数权重矩阵
    MatrixXd GenerateQ(int num_segments, int minimul_order, VectorXd time) {
        int num_poly_coeffs = minimul_order * 2;
        MatrixXd Q = MatrixXd::Zero(num_segments * num_poly_coeffs, num_segments * num_poly_coeffs);
        for (int i = 0; i < num_segments; i++) {
            for (int j = 0; j < num_poly_coeffs; j++) {
                for (int k = 0; k < num_poly_coeffs; k++) {
                    if (j < minimul_order || k < minimul_order) {
                        continue;
                    }
                    Q(i * num_poly_coeffs + j, i * num_poly_coeffs + k) = Factorial(j) / Factorial(j - minimul_order) * Factorial(k) / Factorial(k - minimul_order) * pow(time[i], j + k - 2 * minimul_order + 1) / (j + k - 2 * minimul_order + 1);
                }
            }
        }
        return Q;
    }
    VectorXd GeneratePolyPath1D(VectorXd poly_coef, VectorXd time, double duration) {
        VectorXd path;
        int num_poly_coeffs = poly_coef.size() / time.size();
        for (int i = 0; i < time.size(); i++) {
            double t = 0;
            // 1e-6避免浮点数精度问题
            while (t < time[i] - 1e-6) {
                path.conservativeResize(path.size() + 1);
                path[path.size() - 1] = CalculatePolyValue(poly_coef.segment(i * num_poly_coeffs, num_poly_coeffs), t);
                t += duration;
            }
            if (i == time.size() - 1) {
                path.conservativeResize(path.size() + 1);
                path[path.size() - 1] = CalculatePolyValue(poly_coef.segment(i * num_poly_coeffs, num_poly_coeffs), time[i]);
            }
        }
        return path;
    }
};



int main(int argc, char** argv) {
    // your code
    ros::init(argc, argv, "trajectory_generator");
    ros::NodeHandle nh;
    double duration, duration_step;
    int minimum_order;
    nh.param("trajectory_generator/duration_step", duration_step, 0.1);
    nh.param("trajectory_generator/minimum_order", minimum_order, 3);
    using clock = std::chrono::high_resolution_clock;
    TrajectoryGenerator trajectory_generator;
    // 订阅A*算法给出的无碰撞路径，然后调用GeneratePolyTrajectory函数生成轨迹
    // ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("path", 1);
    Path Astar_path;
    ros::Subscriber path_sub = nh.subscribe<nav_msgs::Path>("path", 1, [&Astar_path](const nav_msgs::Path::ConstPtr& msg) {
        for (const auto& pose : msg->poses) {
            Astar_path.push_back(Vector2d(pose.pose.position.x, pose.pose.position.y));
        }
    });
    // 发布生成的轨迹
    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1);
    ros::Rate rate(10);
    ros::Duration(5.0).sleep();
    ros::spinOnce();
    while (ros::ok()) {
        VectorXd time = trajectory_generator.TimeAllocation(Astar_path);
        // for (int i = 0; i < Astar_path.size(); i++) {
        //     std::cout << Astar_path[i].transpose() << std::endl;
        // }
        // 生成轨迹
        std::cout << "Start Calculate" << std::endl;
        auto start = clock::now();
        Path trajectory = trajectory_generator.GeneratePolyTrajectory(Astar_path, time, duration_step, minimum_order);
        // Path trajectory = trajectory_generator.GenerateBezierTrajectory(Astar_path);
        auto end = clock::now();
        std::cout << "End Calculate" << std::endl;
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
        // for (int i = 0; i < trajectory.size(); i++) {
        //     std::cout << trajectory[i].transpose() << std::endl;
        // }
        // 轨迹可视化
        if (trajectory.empty()) {
            continue;
        }
        nav_msgs::Path trajectory_msg;
        trajectory_msg.header.frame_id = "map";
        trajectory_msg.header.stamp = ros::Time::now();
        for (const auto& point : trajectory) {
            geometry_msgs::PoseStamped pose;
            pose.pose.position.x = point.x();
            pose.pose.position.y = point.y();
            pose.pose.position.z = 0.0; // 平面路径，z 设置为 0
            trajectory_msg.poses.push_back(pose);
        }
        path_pub.publish(trajectory_msg);
        rate.sleep();
        ros::Duration(1.0).sleep();
    }
    return 0;
}
