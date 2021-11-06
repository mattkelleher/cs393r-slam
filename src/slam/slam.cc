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
    add_pose_(true),
    prev_scans_(),
    prev_transforms_() {}

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
  float delta_angle = odom_angle - prev_odom_angle_;
  float delta_dist = sqrt(pow((odom_loc.x() - prev_odom_loc_.x()), 2) + pow((odom_loc.y() - prev_odom_loc_.y()), 2));

  if((delta_angle > M_PI / 6.0) || delta_dist > 0.03) { //TODO what should thresholds be?
    add_pose_ = true;
    prev_odom_loc_ = odom_loc;  
    prev_odom_angle_ = odom_angle;
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
