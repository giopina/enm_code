/*
 * common.h
 *
 *  Created on: Jun 25, 2012
 *      Author: gpolles
 */

#ifndef COMMON_H_
#define COMMON_H_
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>


typedef double   real_t;
typedef int      int_t;
typedef unsigned int uint_t;

#define CLRLINE "\r                                                           \r"


inline void printProgressBar(long double stat, long double max, std::ostream& out, int lenght = 40){
  if (stat>max) stat = max;
  int nch = static_cast<int>(floor((stat/max)*lenght));
  out.setf(std::ios::fixed);
  out << "\r[" ;
  for (int i = 0; i < 40; ++i) {
    if (nch < i) out << " ";
    else out << "#";
  }
  out << "]  " << std::setw(6) << std::setprecision(2)
  << (stat/max)*100.0 << "%" << std::flush;
}


inline std::string str(double x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

#endif /* COMMON_H_ */
