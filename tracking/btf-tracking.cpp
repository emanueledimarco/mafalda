#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include <TFile.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

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

  auto myFile = TFile::Open("/Users/emanuele/Work/cygnus/data/btf-Dec2016/mafalda_output/MAFOutput_MPXNtuple_BTF_SF6_150Torr_gain1440_drift_field_086_Y_47_4600us_2MHz_trg_TOA.root");
  auto fChain = (TTree*)myFile->Get("BlobsFinder");

  if (!myFile || myFile->IsZombie()) {
    return 0;
  }

   Int_t           nBlobs;
   Float_t         geoCenter_x[1000];   //[nBlobs]
   Float_t         geoCenter_y[1000];   //[nBlobs]

   fChain->SetBranchAddress("nBlobs", &nBlobs);
   fChain->SetBranchAddress("geoCenter_x", geoCenter_x);
   fChain->SetBranchAddress("geoCenter_y", geoCenter_y);

   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntries();

   std::vector<float> xcoords,ycoords;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = fChain->LoadTree(jentry);
      if (ientry < 0) break;
      if (jentry !=2) continue;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::cout << "entry " << jentry << "nBlobs = " << nBlobs << std::endl;
      xcoords.clear(); ycoords.clear();
      for (int b = 0; b<nBlobs; ++b) {
        xcoords.push_back(geoCenter_x[b]);
        ycoords.push_back(geoCenter_y[b]);
      }
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
   }
  return 0;
}