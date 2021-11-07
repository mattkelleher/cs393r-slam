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
#include "visualization/CImg.h"

#include "slam.h"

#include "vector_map/vector_map.h"

using namespace math_util;
using Eigen::Affine2f;
using Eigen::Rotation2Df;
using Eigen::Translation2f;
using Eigen::Vector2f;
using Eigen::Vector2i;
using Eigen::VectorXf;
using Eigen::Matrix3f;
using std::cout;
using std::endl;
using std::string;
using std::swap;
using std::vector;
using vector_map::VectorMap;
using cimg_library::CImg;

namespace slam {

SLAM::SLAM() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false),
    add_pose_(true),
    prev_scans_(),
    prev_transforms_() {}

float SLAM::_Distance(Vector2f p1, Vector2f p2) {
  return sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2));
}

void SLAM::GetPose(Eigen::Vector2f* loc, float* angle) const {
  // Return the latest pose estimate of the robot.
  *loc = Vector2f(0, 0);  
  *angle = 0;
  //curr_pose.loc = *loc;
  //curr_pse.angle = *angle;
}


Vector2i SLAM::GetRasterIndex(Vector2f point) {
  return Vector2i(1,1);
}

void SLAM::MakeRaster(vector<Vector2f> pointCloud) {

}

void SLAM::ObserveLaser(const vector<float>& ranges,
                        float range_min,
                        float range_max,
                        float angle_min,
                        float angle_max) {
  // A new laser scan has been observed. Decide whether to add it as a pose
  // for SLAM. If decided to add, align it to the scan from the last saved pose,
  // and save both the scan and the optimized pose.
  if(add_pose_ == false) {
    return;
  }
  add_pose_ = false;
  
  // Convert LIDAR to point cloud scan in current frame
  const Vector2f laserLoc(0.2,0);
  float phi;

  vector<Vector2f> scan, T_scan;
  float increment = (angle_max - angle_min) / ranges.size();
  for (size_t i = 0; i < ranges.size(); i+=10) {
    phi = angle_min +  i * increment;
    Vector2f point(ranges[i] * cos(phi), ranges[i] * sin(phi));
    scan.push_back(point - laserLoc);
  } 
  prev_scans_.push_back(scan);
  if(prev_scans_.size() == 1) {
    Matrix3f hold; 
    hold << 0,0,0,0,0,0,0,0,0;
    prev_transforms_.push_back(hold);
    return;
  }
  
  // Declare candidate transfore matrix, T and best tranform matrix T_best
  Matrix3f T, T_best;
  T << 0,0,0,
       0,0,0,
       0,0,0;
  T_best << 0,0,0, 
            0,0,0,
            0,0,0; 

  //check CUBE 
  float prob_max = 0; 
  for(int i = -45; i <=45; i++) { // check potential theta's -45 degrees to 45
    float theta = M_PI / 180.0 * i; // angle in radians
    vector<Vector2f> t_scan;  // scan pointcloud transformed with only rotation
    //Define rotation matrix
    T(0,0) = cos(theta);
    T(0,1) = sin(theta);
    T(0,2) = 0;
    T(1,0) = -sin(theta);
    T(1,1) = cos(theta);
    T(1,2) = 0;
    // Rotate scan
    for (size_t m = 0; m < scan.size(); m++) {
      Vector2f point;
      point.x() = T(0,0) * scan[m].x() + T(0,1) * scan[m].y() + T(0,2);
      point.y() = T(1,0) * scan[m].x() + T(1,1) * scan[m].y() + T(1,2);
      t_scan.push_back(point);
    }
    //transform scan here (only once, we will then slide it around in inner loops)
    for(int x = -100; x < 101; x+=5) { //translation in x direction +- 1m (j is in cm)
      for(int y = -100; y < 101; y+=5) { // translation in y direction +-1m (k is in cm)
        float prob = 0;
        if (prob > prob_max) {
          prob_max = prob;
          T_best = T;
          T_best(0,2) = -x * cos(theta) - y * sin(theta);
          T_best(1,2) = x * sin(theta) - y * cos(theta);
        }
      }
    }
  }

 /* // Check odometry to see if pose is sufficiently divergent
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
      Vector2f point;
      point.x() = T(0,0) * scan[M].x() + T(0,1) * scan[M].y() + T(0,2);
      point.y() = T(1,0) * scan[M].x() + T(1,1) * scan[M].y() + T(1,2);
      T_scan.push_back(point)
    }
    for (int M = 1; M <= length(T_scan), M++) {
      // Read pixel data
      int x_pixel = std::floor(raster.width()/2) + std::floor((prev_loc.x() - T_scan[M].x())/pixel_def);
      int y_pixel = std::floor(raster.height()/2) + std::floor((prev_loc.y() - T_scan[M].y())/pixel_def); 
      prob += raster(x_pixel, y_pixel); //TODO what if x_pixel, y_pixel is off the raster?
    }
    // Normalize probability
    if (prob > prob_max) {
      prob_max = prob;
      T_best = T;
    }
*/
  
  prev_transforms_.push_back(T_best); 

  //TODO build raster from current scan (for next iteration of matching)
}

void SLAM::ObserveOdometry(const Vector2f& odom_loc, const float odom_angle) {
  if (!odom_initialized_) {
    prev_odom_angle_ = odom_angle;
    prev_odom_loc_ = odom_loc;
    odom_initialized_ = true;
    return;
  }
  curr_odom_loc_ = odom_loc;
  curr_odom_angle_ = odom_angle;
  // Keep track of odometry to estimate how far the robot has moved between 
  // poses.
  float delta_angle = odom_angle - prev_odom_angle_;
  float delta_dist = sqrt(pow((odom_loc.x() - prev_odom_loc_.x()), 2) + pow((odom_loc.y() - prev_odom_loc_.y()), 2));

  if((delta_angle > M_PI / 6.0) || delta_dist > 0.5) { //TODO what should thresholds be?
    add_pose_ = true;
    //prev_odom_loc_ = odom_loc;  
    //prev_odom_angle_ = odom_angle;
  }
}

vector<Vector2f> SLAM::GetMap() {
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.
  
  // Go through all previous scans except the first scan in reverse order and 
  // apply transformations 
  for(size_t i = prev_scans_.size() - 1; i > 0; i--) {
    // Update all points in current map
    for(size_t j = 0; j < map.size(); j++){
      map[j].x() = prev_transforms_[i](0,0) * map[j].x() + prev_transforms_[i](0,1) * map[j].y() + prev_transforms_[i](0,2);
      map[j].y() = prev_transforms_[i](1,0) * map[j].x() + prev_transforms_[i](1,1) * map[j].y() + prev_transforms_[i](1,2);
    }
    // Add new points from scan[i] into map
    for(auto p: prev_scans_[i]) {
      Vector2f point;
      point.x() = prev_transforms_[i](0,0) * p.x() + prev_transforms_[i](0,1) * p.y() + prev_transforms_[i](0,2);
      point.y() = prev_transforms_[i](1,0) * p.x() + prev_transforms_[i](1,1) * p.y() + prev_transforms_[i](1,2);
      
      map.push_back(point);
    }
  } 
  
  // Add points from first scan, don't need transformations
  for(auto p: prev_scans_[0]) {
    map.push_back(p);
  }
  
  return map;
}

}  // namespace slam
