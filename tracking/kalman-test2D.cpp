/**
 * Test for the KalmanFilter class with 2D motion with constant velocity
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "kalman.hh"

int main(int argc, char* argv[]) {

  int n = 4; // Number of states (x, y and velocities vx, vy)
  int m = 2; // Number of measurements (x and y in the 2D plane)

  double dt = 1.0; // Time step

  Eigen::MatrixXd A(n, n); // System dynamics matrix
  Eigen::MatrixXd C(m, n); // Output matrix
  Eigen::MatrixXd Q(n, n); // Process noise covariance
  Eigen::MatrixXd R(m, m); // Measurement noise covariance
  Eigen::MatrixXd P(n, n); // Estimate error covariance

  // Discrete LTI projectile motion, measuring position only
  A << 1, 0, dt, 0, 0, 1, 0, dt, 0, 0, 1, 0, 0, 0, 0, 1; 
  C << 1, 0, 0, 0, 0, 1, 0, 0;

  // Reasonable covariance matrices
  Q << 3., .0, .0, .0, .0, 3., .0, .0, .0, .0, 3., .0, .0, .0, .0, 3.;
  R << 2., 2., 2., 2.;
  P << 
    4., 4., 2., 2., 
    4., 4., 2., 2.,
    2., 2., 0.01, 2.,
    2., 2., 2., 0.01;

  std::cout << "A: \n" << A << std::endl;
  std::cout << "C: \n" << C << std::endl;
  std::cout << "Q: \n" << Q << std::endl;
  std::cout << "R: \n" << R << std::endl;
  std::cout << "P: \n" << P << std::endl;

  // Construct the filter
  KalmanFilter kf(dt, A, C, Q, R, P);

  // List of noisy position measurements (y)
  std::vector<double> xcoords = {
    1.04202710058, 2, 3, 5, 10
  };
  std::vector<double> ycoords = {
    1.1, 2.1, 2.8, 4.5, 8.9
  };

  // Best guess of initial states
  Eigen::VectorXd x0(n);
  x0 << xcoords[0], ycoords[0], 1, 1;
  kf.init(0.,x0);

  // Feed measurements into filter, output estimated states
  double t = 0;
  Eigen::VectorXd y(m);
  std::cout << "t = " << t << ", " << "x_hat[0]: " << kf.state().transpose() << std::endl;
  for(int i = 0; i < xcoords.size(); i++) {
    t += dt;
    y << xcoords[i], ycoords[i];
    kf.update(y);
    std::cout << "t = " << t << ", " << "y[" << i << "] = " << y.transpose()
        << ", x_hat[" << i << "] = " << kf.state().transpose() << std::endl;
  }

  return 0;
}
