#include <ros/ros.h>
#include <tf/transform_datatypes.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Float64.h>
#include <math.h>
#include "algorithm"
#include <iostream>
#include <fstream>
#include "../include/MPC.h"
#include <Eigen/Core>

int main(int argc, char **argv)
{
    MPC mpc;
    double x_current = 0;
    double y_current = 0;
    double psi_current = 0;
    double v_current = 0;
    double omega_current = 0;
    double x_ref = 2;
    double y_ref = 0;
    double psi_ref = 0;
    double v_ref = 0;
    double omega_ref = 0;
    Eigen::VectorXd state_current(5);
    Eigen::VectorXd state_ref(5);
    state_current << x_current, y_current, psi_current, v_current, omega_current;
    state_ref << x_ref, y_ref, psi_ref, v_ref, omega_ref;
    std::vector<double> control_vec = mpc.Solve(state_current, state_ref);
}
