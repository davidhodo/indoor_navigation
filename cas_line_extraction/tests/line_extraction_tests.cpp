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

TEST(ExtractionTest, ExtractLines1) {
    LineExtractor test_extractor;
    double angle_min=-PI/2.;
    double angle_max=PI/2.;
    double angle_increment=1.*PI/180.;
    VectorXd scan(181);
    scan << -1.5708,-1.5533,-1.5359,-1.5184,-1.501,-1.4835,-1.4661,-1.4486,-1.4312,-1.4137,-1.3963,-1.3788,-1.3614,-1.3439,-1.3265,-1.309,-1.2915,-1.2741,-1.2566,-1.2392,-1.2217,-1.2043,-1.1868,-1.1694,-1.1519,-1.1345,-1.117,-1.0996,-1.0821,-1.0647,-1.0472,-1.0297,-1.0123,-0.99484,-0.97738,-0.95993,-0.94248,-0.92502,-0.90757,-0.89012,-0.87266,-0.85521,-0.83776,-0.8203,-0.80285,-0.7854,-0.76794,-0.75049,-0.73304,-0.71558,-0.69813,-0.68068,-0.66323,-0.64577,-0.62832,-0.61087,-0.59341,-0.57596,-0.55851,-0.54105,-0.5236,-0.50615,-0.48869,-0.47124,-0.45379,-0.43633,-0.41888,-0.40143,-0.38397,-0.36652,-0.34907,-0.33161,-0.31416,-0.29671,-0.27925,-0.2618,-0.24435,-0.22689,-0.20944,-0.19199,-0.17453,-0.15708,-0.13963,-0.12217,-0.10472,-0.087266,-0.069813,-0.05236,-0.034907,-0.017453,0,0.017453,0.034907,0.05236,0.069813,0.087266,0.10472,0.12217,0.13963,0.15708,0.17453,0.19199,0.20944,0.22689,0.24435,0.2618,0.27925,0.29671,0.31416,0.33161,0.34907,0.36652,0.38397,0.40143,0.41888,0.43633,0.45379,0.47124,0.48869,0.50615,0.5236,0.54105,0.55851,0.57596,0.59341,0.61087,0.62832,0.64577,0.66323,0.68068,0.69813,0.71558,0.73304,0.75049,0.76794,0.7854,0.80285,0.8203,0.83776,0.85521,0.87266,0.89012,0.90757,0.92502,0.94248,0.95993,0.97738,0.99484,1.0123,1.0297,1.0472,1.0647,1.0821,1.0996,1.117,1.1345,1.1519,1.1694,1.1868,1.2043,1.2217,1.2392,1.2566,1.2741,1.2915,1.309,1.3265,1.3439,1.3614,1.3788,1.3963,1.4137,1.4312,1.4486,1.4661,1.4835,1.501,1.5184,1.5359,1.5533,1.5708;
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
