#include "../include/MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <Eigen/Core>

using CppAD::AD;

// 滚动时域长度
size_t N = 10;
// 控制周期
double dt = 0.1;  //tested with 0.3, 0.12, 0.1, 0.08

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
// start index for each state

size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
// size_t ey_start = psi_start + N;
// size_t epsi_start = ey_start + N;
size_t v_start = psi_start + N;
size_t omega_start = v_start + N - 1;

class FG_eval {
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd state_ref;
    FG_eval(Eigen::VectorXd state_ref) 
    { 
        this->state_ref = state_ref; 
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& vars) 
    {
        // fg的长度应等于约束条件数目+1（fg[0]是目标函数）
        assert(fg.size() == 3 * N + 1);
        assert(vars.size() == 3 * N + 2 * (N - 1));
        // 参考轨迹的初始状态
        double x_ref_start = state_ref[0];
        double y_ref_start = state_ref[1];
        double psi_ref_start = state_ref[2];
        double v_ref = state_ref[3];
        double omega_ref = state_ref[4];
        // std::cout << "i = 0: x_ref = " << x_ref_start << ", y_ref = " << y_ref_start << ", psi_ref = " << psi_ref_start << std::endl;
        size_t n_vars = vars.size();
        std::vector<double> vars_ref(n_vars);
        vars_ref[x_start] = state_ref[0];
        vars_ref[y_start] = state_ref[1];
        vars_ref[psi_start] = state_ref[2];
        vars_ref[v_start] = state_ref[3];
        vars_ref[omega_start] = state_ref[4];
        for(size_t i = 0; i < N-1; i++)
        {
            vars_ref[x_start+i+1] = vars_ref[x_start+i] + v_ref * CppAD::cos(vars_ref[psi_start+i]) * dt;
            vars_ref[y_start+i+1] = vars_ref[y_start+i] + v_ref * CppAD::sin(vars_ref[psi_start+i]) * dt;
            vars_ref[psi_start+i+1] = vars_ref[psi_start+i] + omega_ref * dt;
            // std::cout << "i = " << i+1 << ": x_ref = " << vars_ref[x_start+i+1] << ", y_ref = " << vars_ref[y_start+i+1] << ", psi_ref = " << vars_ref[psi_start+i+1] << std::endl;
        }
        for(size_t i = 0; i < N-2; i++)
        {
            vars_ref[v_start+i+1] = vars_ref[v_start+i];
            vars_ref[omega_start+i+1] = vars_ref[omega_start+i];
        }
        // fg[0]对应优化目标，fg[1]及其之后的元素对应约束条件
        // 添加初始条件约束
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        // fg[1 + ey_start] = vars[ey_start];
        // fg[1 + epsi_start] = vars[epsi_start];
        // 添加运动学约束
        for (size_t i = 0; i < N - 1; i++) 
        {
            // t+1时刻的状态变量
            AD<double> x1 = vars[x_start + i + 1];
            AD<double> y1 = vars[y_start + i + 1];
            AD<double> psi1 = vars[psi_start + i + 1];
            // AD<double> ey1 = vars[ey_start + i + 1];
            // AD<double> epsi1 = vars[epsi_start + i + 1];
            // t时刻的状态变量
            AD<double> x0 = vars[x_start + i];
            AD<double> y0 = vars[y_start + i];
            AD<double> psi0 = vars[psi_start + i];
            // AD<double> ey0 = vars[ey_start + i];
            // AD<double> epsi0 = vars[epsi_start + i];
            // 假设t时刻到t+1时刻控制变量保持不变
            AD<double> v0 = vars[v_start + i];
            AD<double> omega0 = vars[omega_start + i];

            // 暂时没有考虑时延
            // 离散方法用的是前向欧拉法，也可以考虑改成四阶龙格库塔法
            fg[2 + x_start + i]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[2 + y_start + i]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[2 + psi_start + i]  = psi1 - (psi0 + omega0 * dt);
            // TODO: 航迹误差和航向误差的计算方法还得再好好想清楚
            // fg[2 + ey_start + i]  = ey1 - ;
            // fg[2 + epsi_start + i] = epsi1 - ;
        }
        fg[0] = 0;
        // 跟踪误差尽可能小
        for (size_t i = 0; i < N; ++i) 
        {
            fg[0] += 200*CppAD::pow(vars[x_start+i]-vars_ref[x_start+i], 2);
            fg[0] += 200*CppAD::pow(vars[y_start+i]-vars_ref[y_start+i], 2);
            fg[0] += CppAD::pow(vars[psi_start+i]-vars_ref[psi_start+i], 2);
            // fg[0] += CppAD::pow(vars[ey_start + i], 2);
            // fg[0] += 200*CppAD::pow(vars[epsi_start + i], 2);
        }
        // 控制变量及其变化量尽可能小，最小能量控制
        for (size_t i = 0; i < N - 1; ++i) 
        {
            fg[0] += CppAD::pow(vars[v_start + i], 2);
            fg[0] += CppAD::pow(vars[omega_start + i], 2);
        }
        for (size_t i = 0; i < N - 2; ++i) 
        {
            fg[0] += CppAD::pow(vars[v_start + i + 1] - vars[v_start + i], 2);
            fg[0] += CppAD::pow(vars[omega_start + i +1] - vars[omega_start + i], 2);
        }
    }
};

MPC::MPC() {
    
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state_current, Eigen::VectorXd state_ref) {

    typedef CPPAD_TESTVECTOR(double) Dvector;
    // state是三维向量，依次是前轴中心（参考点）的横坐标，纵坐标，朝向（航迹误差，航向误差暂时不考虑）
    double x    = state_current[0];
    double y    = state_current[1];
    double psi  = state_current[2];
    // double ey   = state_current[3];
    // double epsi  = state_current[4];
    // 控制向量是二维向量，依次是线速度和角速度
    // std::cout << "x = " << x << ", y = " << y << ", psi = " << psi << std::endl;
    // 采用直接离散法（multi-shooting），把状态变量和控制变量都视为决策变量
    size_t n_vars = N * 3 + (N - 1) * 2;
    // 设置运动学约束的总数
    size_t n_constraints = N * 3;

    // 生成决策变量，这里把所有的状态变量和控制变量存储到同一个向量
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) 
    {
        vars[i] = 0;
    }
    // 对应0时刻的初始状态
    vars[x_start]    = x;
    vars[y_start]    = y;
    vars[psi_start]  = psi;
    // vars[ey_start]   = ey;
    // vars[epsi_start] = epsi;
    // 设置变量的上下界
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // 这里不考虑边界的约束，状态变量的上下界都可视为无穷
    for (size_t i = 0; i < v_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // 最大线速度取1.6 m/s
    for (size_t i = v_start; i < omega_start; i++) {
        vars_lowerbound[i] = -1.6;
        vars_upperbound[i] = 1.6;
    }

    // 最大角速度取0.4 rad/s
    for (size_t i = omega_start; i < n_vars; i++) {
        vars_lowerbound[i] = -0.4;
        vars_upperbound[i] = 0.4;
    }

    // 运动学约束都是等式约束
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }

    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    // constraints_lowerbound[ey_start] = ey;
    // constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    // constraints_upperbound[ey_start] = ey;
    // constraints_upperbound[epsi_start] = epsi;

    //构造优化问题
    FG_eval fg_eval(state_ref);

    // options for IPOPT solver
    std::string options;
    options += "Integer print_level  1\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time 0.1\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    std::cout << "solution.status == success: " << (solution.status  == CppAD::ipopt::solve_result<Dvector>::success) << std::endl;
    // Cost
    //auto cost = solution.obj_value;

    // Return the first actuator values
    double v = solution.x[v_start];
    double omega = solution.x[omega_start];
    vector<double> control_vec = {v, omega};
    // std::cout << "success\n";
    std::cout << "v = " << v << ", omega = " << omega << std::endl;
    return control_vec;
}
