#include <iostream>
#include <fstream>
#include "gtest/gtest.h"

// OMG this is so nasty...
#define private public
#define protected public

#include "cas_line_extraction/cas_line_extraction.h"
using namespace cas_line_extraction;

TEST(VectorGeneration, AngleVector) {
    LineExtractor test_extractor;
    double angle_min=-PI/2.;
    double angle_max=PI/2.;
    double angle_increment=1.*PI/180.;
    test_extractor.GenerateAngleVector(angle_min, angle_max, angle_increment);

    //std::cout << test_extractor.lidar_angles_ << std::endl;

    // test results
    ASSERT_EQ(181,test_extractor.lidar_angles_.rows());
}

TEST(VectorGeneration, NoiseVector) {
    LineExtractor test_extractor;
    double angle_min=-PI/2.;
    double angle_max=PI/2.;
    double angle_increment=1.*PI/180.;
    test_extractor.GenerateAngleVector(angle_min, angle_max, angle_increment);
    test_extractor.GenerateNoiseVector(0.05);
    //std::cout << test_extractor.lidar_noise_vector_ << std::endl;

    // test results
    ASSERT_EQ(181,test_extractor.lidar_noise_vector_.rows());
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_(0));
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_[4]);
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_[20]);
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_[55]);
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_[90]);
    ASSERT_EQ(0.05,test_extractor.lidar_noise_vector_[180]);
}

TEST(ExtractionTest, FitLines) {
    LineExtractor test_extractor;
    double angle_min=-PI/2.;
    double angle_max=PI/2.;
    double angle_increment=1.*PI/180.;
    VectorXd scan_window(9);
    VectorXd theta_window(9);
    VectorXd noise_window(9);
    VectorXd parameters;
    VectorXd covariance;

    scan_window << 4.0636,4.0687,4.09,4.0669,4.1029,4.1095,4.107,4.1361,4.1507;
    theta_window <<  -1.5708,-1.5533,-1.5359,-1.5184,-1.5010,-1.4835,-1.4661,-1.4486,-1.4312;
    noise_window <<   0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500;
    test_extractor.GenerateAngleVector(angle_min, angle_max, angle_increment);
    test_extractor.GenerateNoiseVector(0.05);
    test_extractor.SetParameters(9,0.1,23.0259,0.25);

    //test_extractor.FitLinePolar(scan_window,theta_window, noise_window, parameters, covariance);
    //ASSERT_EQ(181,test_extractor.lidar_angles_.rows());

    ASSERT_TRUE(true);
}


TEST(ExtractionTest, ExtractLines1) {
    LineExtractor test_extractor;
    double angle_min=-PI/2.;
    double angle_max=PI/2.;
    double angle_increment=1.*PI/180.;
    VectorXd scan(181);
    scan << 4.0636,4.0687,4.09,4.0669,4.1029,4.1095,4.107,4.1361,4.1507,4.1681,4.193,4.2068,4.2349,4.268,4.2822,4.3,4.3218,4.3766,4.3855,4.4327,4.4655,4.5004,4.5363,4.5674,4.6237,4.6601,4.7089,4.7584,4.8067,4.8545,4.9202,4.9977,5.0693,5.1161,5.1935,5.2658,5.323,5.4235,5.4989,5.583,5.6013,5.5287,5.4655,5.5577,5.6852,5.7988,5.9133,6.0516,6.1806,6.1226,6.0459,5.9839,5.9061,5.8442,5.7869,5.7278,5.6537,5.6182,5.5626,5.5211,4.2597,4.2038,4.1909,4.1555,4.1379,4.0995,4.0705,4.2715,4.5211,4.7882,5.1065,5.0876,5.0652,5.0604,5.0596,5.015,4.9925,4.9717,4.9741,4.9779,4.9619,4.9565,4.9573,4.9609,4.9242,7.4379,7.4353,7.4384,7.4244,7.4443,7.4788,7.4854,7.4936,7.5163,7.5142,7.5476,7.5691,5.0537,5.0801,4.8587,4.5529,4.2591,4.0223,3.799,3.5834,3.4358,3.2436,3.1165,2.9995,2.8493,2.7517,2.6507,2.5635,2.482,2.3939,2.3273,2.2759,2.2056,2.1474,2.0794,2.0144,1.9755,1.9321,1.8918,1.8477,1.7951,1.7513,1.7157,1.6908,1.6701,1.6509,1.6063,1.5898,1.5429,1.5444,1.5089,1.4987,1.4719,1.4394,1.4173,1.3861,1.391,1.3745,1.3566,1.355,1.3318,1.3255,1.3137,1.2968,1.2738,1.2776,1.2574,1.2375,1.2525,1.2586,1.218,1.2184,1.2233,1.1829,1.1787,1.1895,1.197,1.2153,1.2088,1.1826,1.1749,1.1592,1.144,1.158,1.1584,1.1555,1.1694,1.1526,1.1665,1.1575,1.1576,1.1355,1.1579,1.1573,1.1449,1.1571;
    test_extractor.GenerateAngleVector(angle_min, angle_max, angle_increment);
    test_extractor.GenerateNoiseVector(0.05);
    test_extractor.SetParameters(9,0.1,23.0259,0.25);

    //Eigen::Array3d v(M_PI, M_PI/2, M_PI/3);
    Eigen::Vector3d  v(M_PI, M_PI/2, M_PI/3);
    std::cout << v.array().cos() << std::endl;


    test_extractor.ExtractLines(scan);
    //std::cout << test_extractor.lidar_angles_ << std::endl;
    std::cout << scan.segment(0,9);
    // test results
    //ASSERT_EQ(181,test_extractor.lidar_angles_.rows());

    ASSERT_TRUE(true);
}


int main(int argc, char **argv) {
  try {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  } catch (std::exception &e) {
    std::cerr << "Unhandled Exception: " << e.what() << std::endl;
  }
  return 1;
}
