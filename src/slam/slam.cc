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
using cimg_library::CImgDisplay;

namespace slam {

SLAM::SLAM() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false),
    add_pose_(true),
    prev_scans_(),
    prev_transforms_(),
    raster_(1600,1600,1,1,0) {}

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
  Vector2i index;
  index.x() = raster_.width() / 2 +  int(point.x() * 100) / 4;
  index.y() = raster_.height() / 2 - int(point.y()* 100) / 4;
  return index;  
}

void SLAM::MakeRaster(vector<Vector2f> pointCloud) {
  // We have a resolution of 4cm, total space is 62m (30m max range from center in all direction
  // + bonus 2m for translation (max should only be 1m so we should never get out of bounds indices.
  // As such our image only needs to be 6200/4 = 800px by 800px. Each pixel corresponds to a 
  // 4cm x 4cm area in the real world  
  std::cout << "Entering Make Raster" << std::endl;
  CImg<float> image(1600, 1600, 1, 1, 0);
  for(auto point : pointCloud) {
    Vector2i index = GetRasterIndex(point);
    if(index.x() < 0 || index.x() > (raster_.width() - 1) || index.y() < 0 || index.y() > (raster_.height() - 1)) {
      std::cout << "Encountered OOB index, this should not be happening" << std::endl;
      std::cout << "Point: (" << point.x() << ", " << point.y() << ")   Index[" << index.x() << "][" << index.y() << "]" << std::endl;
      continue;
    }
    float color = 0.5;
    cout << "@@@@@@@@@@@@@@@@@@ Drawing point (" << point.x() << ", " << point.y() << ") at: [" << index.x() << "][" << index.y() << "]" << endl;
    image.draw_circle(index.x(), index.y(), 5, &color); 
  }
  
  image.blur(4); //TODO blur over 10cm, not sure what value this should be
  raster_ = image; 
  raster_.save("raster.png");
  //CImgDisplay main_disp(raster_, "Raster Image");
  //while(!main_disp.is_close()) {
  //  main_disp.wait();
  //}  

}

void SLAM::ObserveLaser(const vector<float>& ranges,
                        float range_min,
                        float range_max,
                        float angle_min,
                        float angle_max) {
  // A new laser scan has been observed. Decide whether to add it as a pose
  // for SLAM. If decided to add, align it to the scan from the last saved pose,
  // and save both the scan and the optimized pose.
  std::cout << "Entering Observe Laser" << std::endl;
  if(add_pose_ == false) {
    return;
  }
  add_pose_ = false;
  
  // Convert LIDAR to point cloud scan in current frame
  const Vector2f laserLoc(0.2,0);
  float phi;

  vector<Vector2f> scan, T_scan;
  float increment = (angle_max - angle_min) / ranges.size();
  for (size_t i = 0; i < ranges.size(); i+=3) {
    phi = angle_min +  i * increment;
    if(ranges[i] > range_max) {
      continue; 
    }
    Vector2f point(ranges[i] * cos(phi), ranges[i] * sin(phi));
    scan.push_back(point - laserLoc);
  }
  cout << "Saving point cloud with: " << scan.size() << " points" << endl; 
  prev_scans_.push_back(scan);
  if(prev_scans_.size() == 1) {
    MakeRaster(scan);
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
    T(0,1) = -sin(theta);
    T(0,2) = 0;
    T(1,0) = sin(theta);
    T(1,1) = cos(theta);
    T(1,2) = 0;
    // Rotate scan
    // transform scan here (only once, we will then slide it around in inner loops)
    for (size_t m = 0; m < scan.size(); m++) {
      Vector2f point;
      point.x() = T(0,0) * scan[m].x() + T(0,1) * scan[m].y() + T(0,2);
      point.y() = T(1,0) * scan[m].x() + T(1,1) * scan[m].y() + T(1,2);
      t_scan.push_back(point);
    }
    for(int x = -100; x < 101; x+=4) { //translation in x direction +- 1m (j is in cm)
      for(int y = -100; y < 101; y+=4) { // translation in y direction +-1m (k is in cm)
        float prob = 0;
        for(auto p : t_scan) {
          Vector2i index = GetRasterIndex(p + Vector2f(x/100.0, y/100.0));
          prob += raster_(index.y(), index.x()); 
        }
        if (prob > prob_max) {
          prob_max = prob;
          T_best = T;
          T_best(0,2) = x/100.0 * cos(theta) - y/100.0 * sin(theta);
          T_best(1,2) = x/100.0 * sin(theta) + y/100.0 * cos(theta);
          cout << "########## " << prob << "/" << t_scan.size() << endl  << T_best << endl << "$$$$$$$$$$$$$$$" << endl;
        }
      }
    }
  }

/*      
      T(0,0) = cos(delta_theta);
      T(0,1) = -sin(delta_theta);
      T(0,2) = delta_x*cos(delta_theta)-delta_y*sin(delta_theta);
      T(1,0) = sin(delta_theta);
      T(1,1) = cos(delta_theta);
      T(1,2) = delta_x*sin(delta_theta) + delta_y*cos(delta_theta);
*/
  cout << "Best transformation matrix: " << endl << T_best << endl;  
  prev_transforms_.push_back(T_best); 
  MakeRaster(scan); 
}

void SLAM::ObserveOdometry(const Vector2f& odom_loc, const float odom_angle) {
  std::cout << "Entering Observe Odometry" << std::endl;
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
    cout << "Threshold exceeded, add pose !!" << endl;
    add_pose_ = true;
    prev_odom_loc_ = odom_loc;  
    prev_odom_angle_ = odom_angle;
  }
}

vector<Vector2f> SLAM::GetMap() {
  std::cout << "Entering Get Map" << std::endl;
  vector<Vector2f> map;
  // Reconstruct the map as a single aligned point cloud from all saved poses
  // and their respective scans.

  // No we dont have any scans so we cant build a map
  if(prev_scans_.size() == 0) {
    return map;
  } 
  // Go through all previous scans except the first scan in reverse order and 
  // apply transformations 
  cout << "Num Maps: " << prev_scans_.size() << endl;
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
