/*
  Created by Giovanni Pinamonti
  PhD student @
  SISSA, Trieste, Italy
  November  2013
 */

#include "Structure3d.h"
#include "Matrix.h"
#include "Vector3d.h"
//#include "ListBead.h"
#include "BeadList.h"

class PrincipalComp
{
 private:
  int _N; //number of beads;
  int _size; //dimension of the matrix

  int _ntmax; //maximum number of time steps to read
  int _ntmin; //minimum time step to start the analisys

  int _NTOP; //N. of eigenvectors to print on file

  CovMatrix _covMat; //The correlation matrix
  IntMatrix _intMat; //The corresponding interaction matrix

  CovMatrix1d _forCovMat; //covariance matrix of the forces (w/ scalar prod)

  Structure3d _refStructure;
  Structure3d _mean;

  //  std::vector<std::string> _beadList;
  BeadList _beadList;

 public:
  PrincipalComp();
  PrincipalComp(int n);
  PrincipalComp(std::string fname);
  PrincipalComp(std::string fname, std::vector<std::string> beadlist);
  PrincipalComp(std::string fname, BeadList beadlist);
  virtual ~PrincipalComp();

  void setNTOP(int n){if(n<0) _NTOP=_size; else _NTOP=n;}
  void set_ntmax(int n){_ntmax=n;}
  void set_ntmin(int n){_ntmin=n;}

  double CovGetEigenval(int m) const {
    return _covMat.GetEigenval(m);}
  const std::vector<Vector3d>& CovGetEigenvec(int m) const { return _covMat.GetEigenvec(m);}

  CovMatrix & getCovMat(){return (_covMat);} //mica funzionano
  IntMatrix & getIntMat(){return (_intMat);} // ste due....

  void readTrajectory(const char* fname);
  void readTrajectory(std::string fname){return readTrajectory(fname.c_str());}
  //  void readTrajectory(const char* fname, std::vector<std::string> beadlist);
  void readForces(const char*);
  void readForces(std::string fname){return readForces(fname.c_str());}

  void Diagonalize();

  void Dump(const char* fname);
  void Dump(string fname){return Dump(fname.c_str());}
  void DumpCov(const char* fname);
  void DumpCov(string fname){return Dump(fname.c_str());}

  void dumpDistFluc(const char*name);
  void dumpMSF(const char*name);
  void dumpMSF_partial(const char*name);


  void dump_top_modes(int,double,const char*);

  void ComputeInt();

  void DumpInt(const char*fname);
  void DumpInt(string fname){return DumpInt(fname.c_str());}

/*   double Converge(double *, int); */

/*   void SetReference(double *); */

//###################################################################
/*   void dumpMean(char *name); */

  void Project(int n_comp, const char *fname){return Project(n_comp,fname,fname);}
  void Project(int n_comp, const char *iname, const char *oname );

  void readFromFile(const char *fname, int nmodes);
  void readFromFile(string fname, int nmodes){
    return readFromFile(fname.c_str(),nmodes); }

  void setBeadList(BeadList list){
    _beadList=list;}
  void setBeadList(std::vector<std::string> list){
    _beadList=(BeadList(list));}

};
