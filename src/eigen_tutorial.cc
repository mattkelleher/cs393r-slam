#include <stdio.h>
#include <iostream>
#include <cmath>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"

using std::cout;
using std::endl;
using std::sqrt;
using std::sin;
using std::cos;

using Eigen::Vector2f;
using Eigen::Matrix2f;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::Rotation2Df;

void DemoBasics() {
  cout << "Basic initialization" << endl;
  cout << "Initialize a 2D vector v1." << endl;
  Vector2f v1(1.0, 2.0);

  cout << "Read elements of the 3D vector:" << endl
       << "v1.x = " << v1.x() << endl
       << "v1.y = " << v1.y() << endl;

  cout << "Write 10 to the x coordinate of v1:" << endl;
  v1.x() = 10.0;
  cout << "v1.x = " << v1.x() << endl;

  cout << "Print the vector to stdout:\n" << v1 << endl;

  cout << "Initialize a 2x2 matrix m1." << endl;
  Matrix2f m1;
  m1  << 0, 2,
         3, 0;
  cout << "m1 = " << endl << m1 << endl;
  
  cout << "Initialize a 3D vector v2." << endl;
  VectorXf v2(3);
  v2 << 1, 2, 3;
  cout << "v2= " << endl << v2 << endl;
  
  cout << "Initialize a 3x3 matrix m2." << endl;
  MatrixXf m2(3,3);
  m2  << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
  cout << "m2 = " << endl << m2 << endl;
  
  cout << "Multiply 3x3 matrix by 3x1 vector" << endl;
  VectorXf v4(3);
  v4 = m2 * v2;
  cout << "Resulting vector:\n" << v4 << endl;

  cout << "Multiply matrix times vector." << endl;
  Vector2f v3 = m1 * v1;
  cout << "Resulting vector:\n" << v3 << endl;
}

void DemoRotations() {
  cout << "Rotations demonstration" << endl;
  // Note that pi/4 radians is 45 degrees.
  float angle1 = M_PI / 4.0;
  cout << "angle1 = " << angle1 << " radians = "
       << angle1 / M_PI * 180.0 << " degrees." << endl;

  cout << "Create a rotation" << endl;
  Rotation2Df r1(angle1);

  cout << "Apply that rotation to a vector" << endl;
  Vector2f v1(1.0, 0);
  Vector2f v2 = r1 * v1;
  cout << v2 << "\n";

  cout << "Convert rotation to a matrix" << endl;
  Matrix2f m1 = r1.toRotationMatrix();
  cout << "m1: \n" << m1 << "\n";
}

int main() {
  cout << "Basics: Vectors and Matrices.\n";
  DemoBasics();

  cout << "\n\n\nDifferent representations of rotation.\n";
  DemoRotations();

  return 0;
}
