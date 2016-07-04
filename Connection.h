/*
 * Connection.h
 *
 *  Created on: jun 2014
 *      Author: g. pinamonti
 */

#ifndef CONNECTION_H_
#define CONNECTION_H_

#include "io.h"
#include <iostream>

class Connection
{
  public:
  Connection(int i, int j) {
    _i1=i;
    _i2=j;
    //QUI DOVREI CONTROLLARE CHE SIANO POSITIVI E DIVERSI
  }
  
  inline int GetI1() {return _i1;}
  inline int GetI2() {return _i2;}
  
  inline bool SetI1(int i) { if((i>0)&&(i!=_i2)){ _i1=i; return true;} return false; /*CONTROLLA CHE SIANO POSITIVI E DIVERSI*/}
  inline bool SetI2(int i) { if((i>0)&&(i!=_i1)){ _i2=i; return true;} return false; /*CONTROLLA CHE SIANO POSITIVI E DIVERSI*/}
  
 private:
  int _i1;
  int _i2;
};


//########################################################
//### funzioni utili che utilizzano la suddetta classe ###
//########################################################

inline bool RandSwap(Connection &p,Connection &q) {
  double r=((double) rand() / ((double) RAND_MAX));
  if(r>0.5){
    int i1=p.GetI1();
    int j1=q.GetI1();
    if( (p.SetI1(j1)) &&
	(q.SetI1(i1)) )
      return true;
  }
  else{
    //  int j2=q.GetI2();
    int j1=q.GetI1();
    int i2=p.GetI2();
    if( (p.SetI2(j1)) &&
	(q.SetI1(i2)) )
      return true;
  }
  return false;
}

#endif /* CONNECTION_H_ */
