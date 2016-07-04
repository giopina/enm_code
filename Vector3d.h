/*
 * Vector3d.h
 *
 *  Created on: Jul 2, 2012
 *      Author: gpolles
 */

#ifndef POINT3D_H_
#define POINT3D_H_

#define SQ(i) (i)*(i)
#include <iostream>
#include <vector>
#include <cmath>

class Vector3d
{
  friend std::ostream& operator <<(std::ostream &out, const Vector3d& p);
  friend std::istream& operator >>(std::istream &in, Vector3d& p);


public:
  double X;
  double Y;
  double Z;
  public:
  Vector3d(){X=0;Y=0;Z=0;}
  Vector3d(double ax, double ay, double az) {X=ax;Y=ay;Z=az;}
  Vector3d(const Vector3d& a) {X = a.X ; Y = a.Y ; Z = a.Z;}

  double & getCoord(int mu){
    if(mu==0) return X;
    if(mu==1) return Y;
    if(mu==2) return Z;
  }



  bool operator ==(const Vector3d& a) { if(a.X == X && a.Y==Y && a.Z==Z) return true; return false;}
  inline Vector3d operator-(const Vector3d& a) const{
      return Vector3d(X-a.X,Y-a.Y,Z-a.Z);
  }
  inline Vector3d operator-() const{
        return Vector3d(-X,-Y,-Z);
  }
  inline Vector3d operator+(const Vector3d& a) const{
      return Vector3d(X+a.X,Y+a.Y,Z+a.Z);
  }
  inline Vector3d operator/(const Vector3d& a) const{
        return Vector3d(X/a.X,Y/a.Y,Z/a.Z);
  }
  inline Vector3d& operator-=(const Vector3d& a) {
      X-=a.X; Y-=a.Y; Z-=a.Z;
      return *this;
  }
  inline Vector3d& operator/=(const double a) {
      X/=a; Y/=a; Z/=a;
      return *this;
  }
  inline Vector3d& operator*=(const double a) {
    X*=a; Y*=a; Z*=a;
    return *this;
  }
  inline Vector3d& operator+=(const Vector3d& a) {
        X+=a.X; Y+=a.Y; Z+=a.Z;
        return *this;
  }
  Vector3d operator*(const double a) const{
        return Vector3d(a*X,a*Y,a*Z);
  }
  inline Vector3d operator/(const double a) const{
        return Vector3d(X/a,Y/a,Z/a);
  }

  double operator[](size_t i) {
    return *(&X+i);
  }
  double operator[](int i) {
    return *(&X+i);
  }


  //prodotti scalare e vettoriale (ho cambiato i nomi da dotProducts a dot)
  inline double dot(const Vector3d& p) const {return X*p.X + Y*p.Y + Z*p.Z;}
  inline Vector3d cross(const Vector3d& p) const {
    return Vector3d( Y*p.Z - Z*p.Y,
                    Z*p.X - X*p.Z,
                    X*p.Y - Y*p.X);
  }


  inline double distance(const Vector3d& p) const {return sqrt(SQ(p.X-X)+SQ(p.Y-Y)+SQ(p.Z-Z)); } //distanza fra due punti
  inline double distanceSQ(const Vector3d& p) const {return SQ(p.X-X)+SQ(p.Y-Y)+SQ(p.Z-Z); } //distanza quadra

  inline double normSQ() const {return X*X + Y*Y + Z*Z;} //norma/modulo quadra/o
  inline double norm() const {return sqrt(X*X + Y*Y + Z*Z);} //norma/modulo

  inline double angle(const Vector3d& p) const {
    double x = dot(p)/norm()/p.norm();
    return acos( x );}
};

//########################################################
//### funzioni utili che utilizzano la suddetta classe ###
//########################################################

inline double norm(Vector3d& v) {return sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z);}
inline double normSQ(Vector3d& v) {return (v.X*v.X + v.Y*v.Y + v.Z*v.Z);}

inline std::ostream& operator <<(std::ostream &out, const Vector3d& p){
      out << "(" << p.X << ", " <<  p.Y << ", " << p.Z << ")";
      return out;
}

inline std::istream& operator >>(std::istream &in, Vector3d& p){
        in >> p.X >>  p.Y >> p.Z;
        return in;
}

inline Vector3d operator *(const double x, const Vector3d& v) {
  return v*x;
}


inline double dot(const Vector3d& p,const Vector3d& q) {
  return q.X*p.X + q.Y*p.Y + q.Z*p.Z;
}

#endif /* POINT3D_H_ */
