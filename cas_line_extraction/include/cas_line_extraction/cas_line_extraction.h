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

#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace cas_line_extraction {


class LineExtractor {
public:
	LineExtractor();
	~LineExtractor();

	void ExtractLines(VectorXd scan);

private:
	//! Generates a range of angles for lidar data based on current parameters
	void FillLidarAngles();
	//! Generates a vector of repeated values of the lidar noise standard deviation
	void FillLidarNoiseVector();


	////////////////////////////////////////////////////////
	// Line extraction parameters
	////////////////////////////////////////////////////////
	unsigned int window_size_;	//!< size of sliding window
	double threshold_fidelity_;  //!< threshold on compactness
	double fusion_alpha_;		//!< significance level for line fusion
	double minimum_length_;		//!< minimum length to accept a segment [m]
	double compensation_a_;		//!< compensation factor for angle [rad]
	double compesnation_r_;		//!< compensation factor for distance [m]
	bool cyclic_;				//!< are the scans cyclic?


	////////////////////////////////////////////////////////
	// Lidar parameters
	////////////////////////////////////////////////////////
	double lidar_angle_min_;	//!< lidar minimum angle [rad]
	double lidar_angle_max_;	//!< lidar maximum angle [rad]
	double lidar_angle_increment_; //!< angular resolution [rad]
	double lidar_noise_std_dev_;	//!< standard deviation of lidar noise (1-sigma) [m]

	VectorXd lidar_noise_vector_; //!< vector of noise std deviation [m]
	VectorXd lidar_angles_; 	  //!< vector of angles corresponding to lidar measurements [rad]
};

}
#endif
