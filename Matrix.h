/*
  Created by Giovanni Pinamonti
  PhD student @
  Scuola Internazionale Superiore di Studi Avanzati, Trieste, Italy
  November 21st 2013
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include "Vector3d.h"

class CovMatrix;
class IntMatrix;

class Matrix
{
 protected:
  int _N; //number of beads;
  int _size; //dimension of the matrix
  int _DIM;

  double **_matrix; //elements of the correaltion matrix
  // TODO Potrebbe servirmi salvare la matrice inversa per non ricalcolarla
  //bool _INV;  //  double **_invmat; 
  //eigenstuff for the spectral decomposition
  std::vector<double> _eigenvalues;
  std::vector<std::vector<Vector3d> > _eigenvectors;

  void _swapEigen();

 public:
  Matrix();
  Matrix(int);
  Matrix(int, double **);  
  virtual ~Matrix();
  Matrix(Matrix const &);

  Matrix& operator=(Matrix const &);

  void SetN(int);
  void SetMatrix(int, double**);
  double& GetElement(int i,int j){return _matrix[i][j];}

  int GetN(){return _N;}
  int GetSize(){return _size;}

  double GetEigenval(int i) const {return _eigenvalues.at(i);}
  const std::vector<Vector3d>& GetEigenvec(int i) const {return _eigenvectors[i];}


  void SetEigenvalues(std::vector<double> evals);
  void SetEigenvectors(std::vector<std::vector<Vector3d> > evec);

  void Compose(); // this function "composes" the matrix of which eigenvalues and eigenvectors are known
  void Compose(std::vector<double> evals, std::vector<std::vector<Vector3d> > evec);

  void Decompose();


  //  void GetInverse(Matrix &outMat);
  Matrix GetInverse();

  double GetTrace();
  double GetLogDet();

  void Times(double);
  void Plus(const Matrix &M);

  void Print();
  void dumpAll(const char *);//better not to use it
//###################################################################
  void dumpEigenvalues(const char *name);
  void dumpEigenvectors(const char *name);
  void dumpTopVectors(int i, const char *name);

  void dumpFullMatrix(const char *name);
  void dumpReducedMatrix(const char *name);
  void dumpNormalizedReducedMatrix(const char *name);

  void dumpFullInvMatrix(const char *name);
  void dumpReducedInvMatrix(const char *name);

  /*********************************/
  //di base stampa sempre le eigenvalues e gli eigenvectors nell'ordine inverso (che per la cov e' da maggiore a minore, per la INT viceversa)
  void dumpEigenvalues            (std::string name){return  dumpEigenvalues (name.c_str());}
  void dumpEigenvectors           (std::string name){return  dumpEigenvectors(name.c_str());}
  void dumpTopVectors      (int i, std::string name){return  dumpTopVectors(i,name.c_str());}
  void dumpFullMatrix             (std::string name){return  dumpFullMatrix  (name.c_str());}
  void dumpReducedMatrix          (std::string name){return dumpReducedMatrix(name.c_str());}
  void dumpNormalizedReducedMatrix(std::string name){return dumpNormalizedReducedMatrix (name.c_str());}
  void dumpFullInvMatrix          (std::string name){return dumpFullInvMatrix (name.c_str());}
  void dumpReducedInvMatrix       (std::string name){return dumpReducedInvMatrix (name.c_str());}

  //###################################################################

  friend Matrix dot_product(const Matrix &A,const Matrix &B);
  friend Matrix dot_product(const IntMatrix &A,const Matrix &B);
  friend Matrix dot_product(const Matrix &A,const IntMatrix &B);
  friend Matrix dot_product(const IntMatrix &A,const IntMatrix &B);

};


// NB: nelle CovMatrix gli eigenstuff sono salvati in ordine di eigenvalues crescente! (ma perche' poi?)


class IntMatrix : public Matrix {
 public:
  IntMatrix();
  IntMatrix(int i);
  IntMatrix(int i, double **d);


  friend Matrix dot_product(const IntMatrix &A,const Matrix &B);
  friend Matrix dot_product(const Matrix &A,const IntMatrix &B);
  friend Matrix dot_product(const IntMatrix &A,const IntMatrix &B);
  friend double RWSIP(const IntMatrix &A, const IntMatrix &B);
  friend double RWSIP(const IntMatrix &A, const CovMatrix &B);
  friend double RWSIP(const CovMatrix &A, const IntMatrix &B){return RWSIP(B,A);}
  friend double RMSIP(const IntMatrix &A, const IntMatrix &B, int ntop);
  friend double RMSIP(const IntMatrix &A, const CovMatrix &B, int ntop);
  friend double RMSIP(const CovMatrix &A, const IntMatrix &B, int ntop){return RMSIP(B,A,ntop);}

  friend double RWSIP_new(const IntMatrix &A, const IntMatrix &B);
  friend double RWSIP_new(const IntMatrix &A, const CovMatrix &B);
  friend double RWSIP_new(const CovMatrix &A, const IntMatrix &B){return RWSIP_new(B,A);}
  
  CovMatrix GetInverse();
  void Decompose();
  //  void Compose();
  void Compose(std::vector<double> evals, std::vector<std::vector<Vector3d> > evec);
  void dumpMSF(const char*);
  void dumpMSF(std::string name){return dumpMSF(name.c_str());}
  double GetMSF(int i);
  double GetMSF(int i, int mu, int nu);

  void dumpCovEigenvalues            (std::string name){return  dumpCovEigenvalues (name.c_str());}
  void dumpCovEigenvectors           (std::string name){return  dumpCovEigenvectors(name.c_str());}
  void dumpTopCovVectors      (int i, std::string name){return  dumpTopCovVectors(i,name.c_str());}
  void dumpCovEigenvalues(const char *name);
  void dumpCovEigenvectors(const char *name);
  void dumpTopCovVectors(int i, const char *name);

  
};


class CovMatrix : public Matrix {
 public:
   CovMatrix();
   CovMatrix(int i);
   CovMatrix(int i, double **d);
   
   IntMatrix GetInverse();
   void dumpMSF(const char*);
   void dumpMSF(std::string name){return dumpMSF(name.c_str());}
   double GetMSF(int i);
   double GetMSF(int i, int mu, int nu);

   void Reduce(int NTOP,const CovMatrix &D);

   friend double RWSIP(const CovMatrix &A, const CovMatrix &B);
   friend double RWSIP(const IntMatrix &A, const CovMatrix &B);
   friend double RWSIP_new(const CovMatrix &A, const CovMatrix &B);
   friend double RWSIP_new(const IntMatrix &A, const CovMatrix &B);
   friend double RMSIP(const CovMatrix &A, const CovMatrix &B, int ntop);
   friend double RMSIP(const IntMatrix &A, const CovMatrix &B, int ntop);
   //   friend double RWSIP(const CovMatrix &A, const IntMatrix &B){return RWSIP(B,A);}

   friend double ComputeBhatt(CovMatrix A,CovMatrix B);

};

 class CovMatrix1d : public CovMatrix {

  protected: 
   std::vector<std::vector<double> > _eigenvectors; 
    void _swapEigen();
  public:  
   //Costruttori 
   CovMatrix1d(); 
   CovMatrix1d(int); 
   CovMatrix1d(int, double **); 

   //func diverse perche' usano gli eigenvectors 
   const std::vector<double>& GetEigenvec(int i) const {return _eigenvectors[i];} 
   void SetEigenvectors(std::vector<std::vector<double> >); 
   void Compose(); 
   void Compose(std::vector<double>, std::vector<std::vector<double> >); 
   void Decompose(); 
   Matrix GetInverse();  //TODO 
   void dumpEigenvectors(const char *); 
   void dumpTopVectors(int ,const char *); 
   void dumpMatrix(const char *);
   void dumpReducedMatrix(const char *fname)
     //   {cout<<"*** WARNING: use CovMatrix1d->dumpMatrix() ***"<<endl; 
   {return dumpMatrix(fname);}
   void dumpFullMatrix(const char *fname)
   //   {cout<<"*** WARNING: use CovMatrix1d->dumpMatrix() ***"<<endl; 
   {return dumpMatrix(fname);}
   void dumpNormalizedReducedMatrix(const char *fname)
     //   {cout<<"*** WARNING: use CovMatrix1d->dumpMatrix() ***"<<endl; 
   {return dumpMatrix(fname);}
   void dumpMSF(const char*); 
   double GetMSF(int i); 
   //IntMatrix CovMatrix::GetInverse();//TODO: non so se mi serve sta cosa   
   void Reduce(int NTOP,const CovMatrix1d&D); 
   friend double RWSIP(const CovMatrix1d &A,const CovMatrix1d &B);
   friend double RWSIP_new(const CovMatrix1d &A,const CovMatrix1d &B);
   friend double RMSIP(const CovMatrix1d &A,const CovMatrix1d &B, int ntop); 
   friend double ComputeBhatt(CovMatrix1d A,CovMatrix1d B); 

};

#endif /* MATRIX_H_ */
