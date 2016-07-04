/*
 * Structure3d.h
 *
 *  Created on: Jul 3, 2012
 *      Author: gpolles
 */

// Modified
// by
// Giovanni Pinamonti
// 11/11/13

#ifndef STRUCTURE3D_H_
#define STRUCTURE3D_H_

#include <vector>

#include "Vector3d.h"
#include "Bead.h"
#include <cstdio>
#include "my_malloc.h"
#include "BeadList.h"


class Structure3d
{
  
 public:
  Structure3d();
  Structure3d(double cutoff);
  Structure3d(double cutoff, std::string fname);
  Structure3d(std::string fname);

  Structure3d(std::vector<Bead>);
  /* Structure3d(int,double*); */
  virtual ~Structure3d();
  
  
  std::string getTitle(){ return _title;}
  void setTitle(std::string title){_title=title;}
  void setTitle(const char* title){ std::string str_title(title); return setTitle(str_title);}

  int readFromPDBFile(FILE* fp);
  //  int readFromPDBFile(std::ifstream $fp);
  int readFromPDBFile(const char* fname);
  int readFromPDBFile(std::string fname) {return readFromPDBFile(fname.c_str());}

  int centersFromPDBFile(FILE* fp);
  int centersFromPDBFile(const char* fname);
  int centersFromPDBFile(std::string fname) {return centersFromPDBFile(fname.c_str());}

  int dumpPDBFile(std::ofstream &dumpfile);
  int dumpPDBFile(const char* fname);
  int dumpPDBFile(std::string fname) {return dumpPDBFile(fname.c_str());}
  
  void addBead(Bead );
  
  
  double    getDistance(size_t i, size_t j) const;
  double    getDistanceSQ(size_t i, size_t j) const;
  Vector3d  getDistanceVector(size_t i, size_t j);
  size_t    getSize() const {return _size;}
  size_t    getSizeRes() const {return _sizeRes;}
  void      updateNeighbors();
  void      updateCovBonded();
  
  std::vector<Bead>& getBeads() {
    return _beads;
  }
  Bead& getBead(size_t i){
    return _beads[i];
  }
  

  void setBeadList(vector<string> list)
  {
    _beadList=(BeadList(list));
  }

  void setBeadList(BeadList list)
  {
    _beadList=list;
  }


  bool isBond(int, int);



  double fit(Structure3d *ref_struc, double *w);

 private:
  std::vector<Bead> _beads;
  std::vector<Vector3d> _coordinates;
  std::vector<double>  _bFactors;

  std::string _title;

  size_t _sizeRes;
  size_t _size;
  double _cutoff;
  
  //  bool _isABead(char *);
  
  //  std::vector<string> _beadList;
  BeadList _beadList;

};

#endif /* STRUCTURE3D_H_ */
