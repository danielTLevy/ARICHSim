#ifndef UTILITY_INCLUDE
#define UTILITY_INCLUDE

#include "TVector3.h"
#include "TMatrixD.h"

static TMatrixD makeRotationMatrix(TVector3 dir) {
  // Rotation matrix to rotate (x,y,z)=(0,0,1) onto dir
  // https://math.stackexchange.com/questions/180418
  // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation

  TMatrixD rot(3, 3);
  // a = (0, 0, 1)
  // b = (dirx, diry, dirz)
  // v = a x b = (-diry, dirx, 0)
  TVector3 v(-dir.Y(), dir.X(), 0);
  double s = v.Mag(); // sin(a,b) = ||v|| 
  double c = dir.Z(); // cos(a,b) = a . b = (0,0,1) . (dirx, diry, dirz) = dirz

  // rot =   I   + v    + v^2                    * (1-c)/s^2
  TMatrixD V = TMatrixD(3, 3);

  double k = 1./(1.+c);
  rot(0,0) = 1.0        - (v[2]*v[2]+v[1]*v[1])*k;
  rot(0,1) =     - v[2] - (v[0]*v[1])          *k;
  rot(0,2) =     + v[1] + (v[0]*v[2])          *k;
  rot(1,0) =     + v[2] + (v[0]*v[1])          *k;
  rot(1,1) = 1.0        - (v[2]*v[2]+v[0]*v[0])*k; 
  rot(1,2) =     - v[0] - (v[2]*v[1])          *k;
  rot(2,0) =     - v[1] - (v[0]*v[2])          *k;
  rot(2,1) =     + v[0] + (v[1]*v[2])          *k;
  rot(2,2) = 1.0        - (v[1]*v[1]+v[0]*v[0])*k;
  return rot;
}

static TVector3 rotateVector(TMatrixD rot, double theta, double phi) {
  // Apply a rotation matrix to a theta and phi
  TVector3 dir;
  TVector3 dirRot;
  // Beam frame angle of photon
  dir.SetX(cos(phi)*sin(theta));
  dir.SetY(sin(phi)*sin(theta));
  dir.SetZ(cos(theta));

  // photon direction rotated onto particle direction
  //TVector3 dirRot = rot * dir;
  dirRot[0] = rot[0][0]*dir[0]+rot[0][1]*dir[1]+rot[0][2]*dir[2];
  dirRot[1] = rot[1][0]*dir[0]+rot[1][1]*dir[1]+rot[1][2]*dir[2];
  dirRot[2] = rot[2][0]*dir[0]+rot[2][1]*dir[1]+rot[2][2]*dir[2];
  return dirRot;
}

#endif