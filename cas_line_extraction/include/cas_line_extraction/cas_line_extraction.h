// Line extraction algorithm from CAS Toolbox
// http://openslam.org
// https://svn.openslam.org/data/svn/cas-rnt
// Matlab code is licensed under the GPL
//
// converted to C++ by David Hodo - david.hodo@gmail.com
//
// MATLAB License below:
//
//%   Reference:
//%      K.O. Arras, "Feature-Based Robot Navigation in Known and Unknown
//%      Environments", Ph.D. dissertation, Nr. 2765, Swiss Federal Insti-
//%      tute of Technology Lausanne, Autonomous Systems Lab, June 2003.
//%
//% This file is part of the CAS Robot Navigation Toolbox.
//% Please refer to the license file for more infos.
//% Copyright (c) 2004 Kai Arras, ASL-EPFL
//% v.1.0, 1995, Kai Arras, IfR-ETHZ, diploma thesis
//% v.2.0-v.4.0, 1997, Kai Arras, ASL-EPFL: probabilistic version
//% v.4.1.1, 1.99/6.00/7.00, Kai Arras, ASL-EPFL
//% v.4.2, 05.12.02, Kai Arras, ASL-EPFL
//% v.4.3-4.4, Dec.2003, Kai Arras, CAS-KTH: minor adaptations for toolbox
#ifndef LINEEXTRACTOR_H
#define LINEEXTRACTOR_H

//#define EIGEN2_SUPPORT
//#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cas_line_extraction {

#define PI  3.14159265


class LineExtractor {
public:
	LineExtractor();
	~LineExtractor();

    void ExtractLines(Eigen::ArrayXd scan);
    void GenerateAngleVector(double angle_min, double angle_max, double angle_increment);
    void GenerateNoiseVector(double noise_std_deviation);
    void SetParameters(	unsigned int window_size, double threshold_fidelity,
        double fusion_alpha, double minimum_length);

private:
	//! Generates a range of angles for lidar data based on current parameters
	void FillLidarAngles();
	//! Generates a vector of repeated values of the lidar noise standard deviation
	void FillLidarNoiseVector();

    void FitLinePolar(Eigen::VectorXd range, Eigen::VectorXd theta, Eigen::VectorXd noise,
                     Eigen::VectorXd &parameters, Eigen::MatrixXd &covariance);


	////////////////////////////////////////////////////////
	// Line extraction parameters
	////////////////////////////////////////////////////////
	unsigned int window_size_;	//!< size of sliding window
	double threshold_fidelity_;  //!< threshold on compactness
	double fusion_alpha_;		//!< significance level for line fusion
	double minimum_length_;		//!< minimum length to accept a segment [m]
	double compensation_a_;		//!< compensation factor for angle [rad]
	double compesnation_r_;		//!< compensation factor for distance [m]


	////////////////////////////////////////////////////////
	// Lidar parameters
	////////////////////////////////////////////////////////
	double lidar_angle_min_;	//!< lidar minimum angle [rad]
	double lidar_angle_max_;	//!< lidar maximum angle [rad]
	double lidar_angle_increment_; //!< angular resolution [rad]
	double lidar_noise_std_dev_;	//!< standard deviation of lidar noise (1-sigma) [m]
    size_t number_of_scan_points_;      //!< number of points in one scan, calculated from angle parameters

    Eigen::VectorXd lidar_noise_vector_; //!< vector of noise std deviation [m]
    Eigen::VectorXd lidar_angles_; 	  //!< vector of angles corresponding to lidar measurements [rad]
};

}
#endif
