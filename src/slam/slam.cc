//========================================================================
//  This software is free: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License Version 3,
//  as published by the Free Software Foundation.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  Version 3 in the file COPYING that came with this distribution.
//  If not, see <http://www.gnu.org/licenses/>.
//========================================================================
/*!
\file    slam.cc
\brief   SLAM Starter Code
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "shared/math/geometry.h"
#include "shared/math/math_util.h"
#include "shared/util/timer.h"

#include "slam.h"

#include "vector_map/vector_map.h"

using namespace math_util;
using Eigen::Affine2f;
using Eigen::Rotation2Df;
using Eigen::Translation2f;
using Eigen::Vector2f;
using Eigen::Vector2i;
using Eigen::VectorXf;
using Eigen::MatrixXf;
using std::cout;
using std::endl;
using std::string;
using std::swap;
using std::vector;
using vector_map::VectorMap;

namespace slam {

SLAM::SLAM() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false),
    prev_scans_(),
    prev_transforms_() {}

float ParticleFilter::_Distance(Vector2f p1, Vector2f p2) {
  return sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2));
}

void SLAM::GetPose(Eigen::Vector2f* loc, float* angle) const {
  // Return the latest pose estimate of the robot.
  *loc = Vector2f(0, 0);
  *angle = 0;
  //curr_pose.loc = *loc;
  //curr_pse.angle = *angle;
}

void SLAM::ObserveLaser(const vector<float>& ranges,
                        float range_min,
                        float range_max,
                        float angle_min,
                        float angle_max) {
  // A new laser scan has been observed. Decide whether to add it as a pose
  // for SLAM. If decided to add, align it to the scan from the last saved pose,
  // and save both the scan and the optimized pose.

  // Initialize starting transform matrix guess T0 with odometry
  Vector2f loc;
  float theta;

  // IMPORTANT: Assumed there was a variable called prev_loc and prev_angle
  GetPose(loc, theta);  // I presume this gives us the current coordinates of the robot
                        // If not, change

  Matrix3f T0;  
  float delta_x = loc.x() - prev_loc.x();
  float delta_y = loc.y() - prev_loc.y();
  float delta_theta = theta - prev_theta;
  T0  <<  cos(delta_theta), sin(delta_theta), -loc.x()*cos(delta_theta)-loc.y()*sin(delta_theta),
          -sin(delta_theta), cos(delta_theta), loc.x()*sin(delta_theta)-loc.y()*cos(delta_theta),
          0, 0, 1;

  // Initialize T and T_best
  Matrix3f T, T_best;

  // Initialize scan in current frame
  const float distance_base2lidar = 0.2; // From assignment 1
  float phi;

  float x_base2lidar = distance_base2lidar * cos(angle);
  float y_base2lidar = distance_base2lidar * sin(angle);

  vector<Vector2f> scan, T_scan;
  float increment = (angle_max - angle_min) / range.size();
  for (size_t i = 0; i < scan.size(); i++) {
    phi = angle + (angle_min + i * increment);
    scan(i) = Vector2f temp_loc(loc.x() + x_base2lidar + range(i)*cos(phi), loc.y() + y_base2lidar) + range(i)*sin(phi);
  }

  // Initialize randomized error generators
  // I'm unsure if this is the correct method, but this should sample a normal distribution 
  std::default_random_engine generator;

  float var = 0.1; // Just a guess
  std::normal_distribution<float> x_dist(prev_loc.x(),var);
  std::normal_distribution<float> y_dist(prev_loc.y(),var);
  std::normal_distribution<float> theta_dist(prev_theta,var);

  // Check odometry to see if pose is sufficiently divergent
  if ((_Distance(loc, prev_loc) > 0.5) || (abs(theta - prev_theta) > 30*2*M_PI/180)) {
    for (int N = 1; N < 15; N++) {
      if (N = 1) {
        T = T0;
      } else {
        // Adding 'error' to Trasformation Matrix  
        float delta_x = loc.x() - x_dist(generator);
        float delta_y = loc.y() - y_dist(generator);
        float delta_theta = theta - theta_dist(generator);
        
        T(0,0) = cos(delta_theta);
        T(0,1) = sin(delta_theta);
        T(0,2) = -delta_x*cos(delta_theta)-delta_y*sin(delta_theta);
        T(1,0) = -sin(delta_theta);
        T(1,1) = cos(delta_theta);
        T(1,2) = delta_x*sin(delta_theta) -delta_y*cos(delta_theta);
      }

      for (int M = 1; M <= length(scan), M++) {
        T_scan(M) = T * scan(i);
      }
      for (int M = 1; M <= length(T_scan), M++) {
        // Read pixel data
        int x_pixel = floor(raster.width()/2) + floor((prev_loc.x() - T_scan(M).x())/pixel_def);
        int y_pixel = floor(raster.height()/2) + floor((prev_loc.y() - T_scan(M).y())/pixel_def);
        prob += raster(x_pixel, y_pixel);
      }
      // Normalize probability
      prob = prob/scan.size();

      if (prob > prob_max) {
        prob_max = prob;
        T_best = T;
      }
    }
  }
}

void SLAM::ObserveOdometry(const Vector2f& odom_loc, const float odom_angle) {
  if (!odom_initialized_) {
    prev_odom_angle_ = odom_angle;
    prev_odom_loc_ = odom_loc;
    odom_initialized_ = true;
    return;
  }
  // Keep track of odometry to estimate how far the robot has moved between 
  // poses.
  // delta_x = odom_loc.x - prev_odom_loc.x;
  // delta_y = odom_loc.y - prev_odom_loc.y;
  // delta_angle = odom_angle - prev_odom_angle;
  // delta_dist = sqrt(pow((odom_loc.x - prev_odom_loc.x), 2) + pow((odom_loc.y - prev_odom_loc.y), 2));
  // prev_odom_loc = odom_loc;
  // prev_odom_angle = odom_angle;
}

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.
  return map;
}

}  // namespace slam
