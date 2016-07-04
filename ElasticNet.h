/* This utterly awesome class         */
/* was created by Giovanni Pinamonti, */
/* PhD student at SISSA, Trieste.     */
/* November 20th 2013                 */


#ifndef ELASTICNET_H_
#define ELASTICNET_H_
#define TOL 0.000001

#include "Structure3d.h"
#include "Matrix.h"
#include <cstdio>
#include "Connection.h"
#include "my_malloc.h"

class ElasticNet
{

 public:
  ElasticNet();
  ElasticNet(Structure3d structure);
  ElasticNet(const char *filename); //TODO costruttore direttamente con nome PDB
  virtual ~ElasticNet();

  CovMatrix & getCovMatrix(){ return (_covMatrix);}
  IntMatrix & getIntMatrix(){ return (_intMatrix);}

  //    double CovGetEigenval(int m) const { return _covMatrix.GetEigenval(m);}

  //be careful when using this functions!!!
  double CovGetEigenval(int m) const { if(_intMatrix.GetEigenval(m)>TOL) return 1./_intMatrix.GetEigenval(m); else return 0.0; }
  //const std::vector<Vector3d>& CovGetEigenvec(int m) const { return _covMatrix.GetEigenvec(m);}
  const std::vector<Vector3d>& CovGetEigenvec(int m) const { return _intMatrix.GetEigenvec(m);}
  //  const std::vector<Vector3d>& CovGetEigenvec_fast(int m) const { return _intMatrix.GetEigenvec(m);}


  void readParameters(const char *filename);

  int readPDBFile(const char* fname);
  int readPDBFile(std::string fname){return readPDBFile(fname.c_str());}
  int readPDBFile(FILE *fp);

  void setFast(){_FAST=true;}
  void setFast(bool fast){_FAST=fast;}

  void setExp(){_EXP=true;}
  void setExp(bool fast){_EXP=fast;}

  int getN(){ return _N;}

  void setStructure(Structure3d structure);
  void setCutOff(double cutoff){_cutoff=cutoff;}
  void setRprot(double r){_R_prot=r;}
  void setNtopVectors(int n){_ntop_vectors=n;}

  void setBeadList(vector <string> list){_beadList=list; _structure.setBeadList(_beadList);}
  void setOutBeadList(vector<string> list){_outBeadList=list; _TRACCIAMENTO=true;} //ATTENZIONE: qui devo stare attento che non sia outbeadlist>beadlist

  Structure3d & getStructure(){return (_structure);}

  void constructContactMap();
  void constructIntMat();
  void computeCovar();
  int Solve();
 
  void Dump(const char *name);
  void Dump(std::string name){return Dump(name.c_str());}
  void dumpCovEigenvalues(const char*name){_intMatrix.dumpEigenvalues(name);}
  void dumpCovEigenvalues(std::string name){return dumpCovEigenvalues(name.c_str());}
  //  void dumpMSF(const char *name){_intMatrix.dumpMSF(name);}
  void dumpMSF(const char *name);
  void dumpMSF(std::string name){return dumpMSF(name.c_str());}
  void dumpDistFluc(const char *name);
  void dumpDistFluc(std::string name){return dumpDistFluc(name.c_str());}
  void dumpDistVecFluc(const char *name);
  void dumpDistVecFluc(std::string name){return dumpDistVecFluc(name.c_str());}
  void dump_top_modes(int, double, const char *name);
  void dump_top_modes(int ntop, double amp, std::string name){return dump_top_modes(ntop,amp,name.c_str());}

  double getDistFluc(int, int);

  void dumpNearestNeighbours(const char *name);
  void dumpNearestNeighbours(std::string name){return dumpNearestNeighbours(name.c_str());}

  int tracciamento();  //0 if success -1 if not

  void computeForCovMatrix();
  void dumpForCovMatrix(const char*);
  void dumpForCovMatrix(std::string name){return dumpForCovMatrix(name.c_str());}


 private:
  //  Matrix *_intMatrix;
  IntMatrix _intMatrix;
  CovMatrix _covMatrix; //ATTENZIONE: posso fare questa cosa? Lo dichiaro senza niente e poi quando lo chiamo gli assegno le cose?
  CovMatrix1d _forCovMatrix;

  double _cutoff;
  double _R_B;
  double _R_BB;
  double _R_PP;
  double _R_prot;
  double _R_N7;
  double _Kcov;
  double _KBB;
  double _KB;
  double _KPP;

  int _N;
  //  int _size;
  int _ntop_vectors;
  Structure3d _structure ; //, _out_struct; questa volendo la potrei fare. Ma alla fine una volta che traccio mi interessa solo la struttura "ridotta". O no?

  int **_cmap;

  vector<string> _beadList;

  vector<string> _outBeadList;
  vector<int> _outResList;

  bool _TRACCIAMENTO;

  bool _COM;

  //  bool _PROTEIN;

  bool _OUTRES; // true if you need to reduce the output to a set of residues

  bool _FAST;

  bool _DPF;

  bool _RANDOM;

  bool _EXP;

  bool _RARAND;

  bool _COVFOR;

  bool _DUMPMODES;
  bool _DUMPEIGENVALUES;
  bool _DUMPCOVMAT;
  bool _DUMPDISTFLUC;
  bool _DUMPREDCOVMAT;
  bool _DUMPNORMREDCOVMAT;
  bool _DUMPMSD;
 
  void constructRandMat(double **);
  void constructDemiRandMat();
  void constructExpMat(double **);


};




#endif /*ELASTICNET_H_*/
