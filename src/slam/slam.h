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
\file    slam.h
\brief   SLAM Interface
\author  Joydeep Biswas, (C) 2018
*/
//========================================================================

#include <algorithm>
#include <vector>
#include <cmath>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "visualization/CImg.h"

#ifndef SRC_SLAM_H_
#define SRC_SLAM_H_

namespace slam {

class SLAM {
 public:
  // Default Constructor.
  SLAM();

  // Observe a new laser scan.
  void ObserveLaser(const std::vector<float>& ranges,
                    float range_min,
                    float range_max,
                    float angle_min,
                    float angle_max);

  // Distance between 2 2D points
  float _Distance(Eigen::Vector2f p1, Eigen::Vector2f p2);

  // Returns index corresponding to point in raster
  Eigen::Vector2i GetRasterIndex(Eigen::Vector2f point);

  // Create raster image from pointCloud
  void MakeRaster(std::vector<Eigen::Vector2f> pointCloud);
  
  // Observe new odometry-reported location.
  void ObserveOdometry(const Eigen::Vector2f& odom_loc,
                       const float odom_angle);

  // Get latest map.
  std::vector<Eigen::Vector2f> GetMap();

  // Get latest robot pose.
  void GetPose(Eigen::Vector2f* loc, float* angle) const;

 private:

  // Previous odometry-reported locations.
  Eigen::Vector2f prev_odom_loc_;
  float prev_odom_angle_;
  bool odom_initialized_;
  bool add_pose_;
  
  Eigen::Vector2f curr_odom_loc_;
  float curr_odom_angle_;
  // Previous scans and associated transforms
  std::vector<std::vector<Eigen::Vector2f>> prev_scans_;
  std::vector<Eigen::Matrix3f> prev_transforms_;

  // Raster image
  cimg_library::CImg<float> raster_;

};
}  // namespace slam

#endif   // SRC_SLAM_H_
