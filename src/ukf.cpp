#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  n_aug_ = 7;
  n_x_ = 5;
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_,n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  lambda_ = 4;

  weights_ = VectorXd(15);

  weights_(0) = lambda_ / (lambda_ + 7);

  weights_.tail(2*n_aug_).fill( 0.5 / (lambda_ + n_aug_) );

  time_us_ = 0;
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

   Xsig_pred_ = MatrixXd(5, 15);
   is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) 
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if(!is_initialized_)
    {
       if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
        {
            std::cout<<"KF initialization with Radar"<<std::endl;
             x_ << cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0], 
                   sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0], 
                   abs(measurement_pack.raw_measurements_[2]), 
                    0.,
                    0.;
             time_us_ = measurement_pack.timestamp_;
            is_initialized_ = true;
         }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {   std::cout<<"KF initialization with Laser"<<std::endl;
           x_ << measurement_pack.raw_measurements_[0],
                measurement_pack.raw_measurements_[1],
               0.,
               0.,
               0.;
               time_us_ = measurement_pack.timestamp_;
               is_initialized_ = true; 
     }
       else
       {
            std::cout<<"KF initialization failed, sensor unknown"<<std::endl;
       }
      return; 
    }
    double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;  
    time_us_ = measurement_pack.timestamp_;
    
     std::cout << "delta_t = " << delta_t << std::endl;  
     while (delta_t > 0.1)
     {
        const double dt = 0.1;
        Prediction(dt);
        delta_t -= dt;
        }
     Prediction(delta_t);
    
     if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
     {
         std::cout<<"radar"<<std::endl;
     UpdateRadar(measurement_pack);
     } 
     else if (use_laser_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     UpdateLidar(measurement_pack);
      }
     else
     {
    
    std::cout << "Kalman Filter Update: sensor unkown or not used" << std::endl;
     }
      
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //set augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    //set augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    #if UKF_DEBUG
    std::cout << "P_aug = " << std::endl << P_aug << std::endl;
    #endif
    
    //calculate square root matrix
    MatrixXd A = P_aug.llt().matrixL();
    #if UKF_DEBUG
    std::cout << "A = " << std::endl << A << std::endl;
    #endif
    
    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    MatrixXd XsigDiff = sqrt(lambda_ + n_aug_) * A;
    for (int i=0; i<n_aug_; i++)
    {
        Xsig_aug.col(1+i)         = x_aug + XsigDiff.col(i);
        Xsig_aug.col(1+i+n_aug_)  = x_aug - XsigDiff.col(i);
    }

    //print augmented sigma points result
    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;


    /*****************************************************************************
    *  Predigt sigma points
    ****************************************************************************/
    //predict sigma points
    for (int i=0; i<2 * n_aug_ + 1; i++)
    {
        // for readability
        double px = Xsig_aug(0,i);
        double py = Xsig_aug(1,i);
        double v  = Xsig_aug(2,i);
        double psi = Xsig_aug(3,i);
        double psiDot = Xsig_aug(4,i);
        //psiDot = ((abs(psiDot)<(2.*M_PI)) ? psiDot : ((psiDot > 0) ? (2.*M_PI) : (-2.*M_PI))); //avoid to high yaw rates
        double nu_a = Xsig_aug(5,i);
        double nu_psiDDot = Xsig_aug(6,i);
        double px_pred, py_pred;
        //avoid division by zero
        if (fabs(psiDot) > 0.001) {
                px_pred = px + v/psiDot * ( sin (psi + psiDot*delta_t) - sin(psi));
                py_pred = py + v/psiDot * ( cos(psi) - cos(psi+psiDot*delta_t) );
        }
        else {
                px_pred = px + v*delta_t*cos(psi);
                py_pred = py + v*delta_t*sin(psi);
        }
        //add noise
        px_pred += 0.5*nu_a*delta_t*delta_t * cos(psi);
        py_pred += 0.5*nu_a*delta_t*delta_t * sin(psi);
            
        double v_pred = v + delta_t * nu_a;
        double psi_pred = psi + psiDot*delta_t + 0.5*delta_t*delta_t * nu_psiDDot;
        double psiDot_pred = psiDot + delta_t * nu_psiDDot;
        
        //write predicted sigma points into right column
        Xsig_pred_(0,i) = (((abs(px_pred)>0.0001) || (abs(py_pred)>0.0001)) ? px_pred : 0.0001); // avoid devision by zero later
        Xsig_pred_(1,i) = py_pred;
        Xsig_pred_(2,i) = v_pred;
        Xsig_pred_(3,i) = psi_pred;
        Xsig_pred_(4,i) = psiDot_pred;
    }  

    //print predicted sigma points
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;


    /*****************************************************************************
    *  Calculate predicted state mean and covariance
    ****************************************************************************/
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ += weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        //x_diff(3) = std::fmod (x_diff(3), (2.*M_PI));
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }

    //print predicted state mean and covariance
    #if UKF_DEBUG
        std::cout << "Predicted state" << std::endl;
        std::cout << x_ << std::endl;
        std::cout << "Predicted covariance matrix" << std::endl;
        std::cout << P_ << std::endl;
    #endif
                                                                                                                                                                                                                                                 
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    //set measurement dimension, lidar can measure px and py
  int n_z = 2;

  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      Zsig(0,i) = Xsig_pred_(0,i);  // px
      Zsig(1,i) = Xsig_pred_(1,i);  // py
  }
  
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  //and cross correlation matrix
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    //x_diff(3) = std::fmod (x_diff(3), (2.*M_PI));
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;

  //print predicted measurement mean and covariance matrix S
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

  /*****************************************************************************
  *  Update state with lidar measurement
  ****************************************************************************/
  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;
  //update state mean and covariance matrix
  VectorXd z_diff_mean =  measurement_pack.raw_measurements_ - z_pred;
  x_ += K * z_diff_mean;
  P_ -= K * S * K.transpose();

  //print updated state and covariance
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  //calculate the NIS
  NIS_laser_ = z_diff_mean.transpose() * S_inv * z_diff_mean;
  std::cout << "NIS lidar: " << std::endl << NIS_laser_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /*****************************************************************************
  *  Predict radar measurement
  ****************************************************************************/
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      // for readability
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v  = Xsig_pred_(2,i);
      double psi = Xsig_pred_(3,i);
      double psiDot = Xsig_pred_(4,i);

      Zsig(0,i) = sqrt(px*px + py*py);
      Zsig(1,i) = std::atan2(py,px);
      Zsig(2,i) = v*(px*cos(psi) + py*sin(psi)) / Zsig(0,i);
  }
  
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  //and cross correlation matrix
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    //x_diff(3) = std::fmod (x_diff(3), (2.*M_PI));
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    //z_diff(1) = std::fmod (z_diff(1), (2.*M_PI));
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

  //print predicted measurement mean and covariance matrix S
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

  /*****************************************************************************
  *  Update state with radar measurement
  ****************************************************************************/
  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;
  //update state mean and covariance matrix
  VectorXd z_diff_mean =  measurement_pack.raw_measurements_ - z_pred;
  x_ += K * z_diff_mean;
  P_ -= K * S * K.transpose();

  //print updated state and covariance
    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  //calculate the NIS
  NIS_radar_ = z_diff_mean.transpose() * S_inv * z_diff_mean;
  std::cout << "NIS radar: " << std::endl << NIS_radar_ << std::endl;
}
