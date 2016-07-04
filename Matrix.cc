/*
  Created by Giovanni Pinamonti
  PhD student @
  Scuola Internazionale Superiore di Studi Avanzati, Trieste, Italy
  November 21st 2013
 */

//extern "C" {
//#include <mkl.h>
//}

#define TOL 0.000001

#include <fstream>
#include <iostream>
#include <cmath>
#include "omp.h"


#include "Matrix.h"
#include "my_malloc.h"
#include "io.h"
//#include "Vector3d.h"

#include "lapack_matrix_routines_wrapper_v3.cc"

using namespace std;

//Constructors for the class Matrix

Matrix:: Matrix():_DIM(3){ //forse sta cosa non si dovrebbe fare. Se non dichiaro la dimensione della matrice non la creo...
  _N=0;
  _size=0;
  cout<<"Constructing a null matrix"<<endl;

}

Matrix:: Matrix(int N):_DIM(3){
  //size of matrix
  _N=N;
  _size=_DIM*_N;
  cout<<"Constructing a zero matrix of size N = "<<N<<endl;
  //allocating memory and initializing arrays
  _matrix=d2t(_size,_size);  
  //options to print the output

}//END OF CONSTRUCTOR


Matrix:: Matrix(int N, double **matrix):_DIM(3){
  _N=N;
  _size=_DIM*_N;
  cout<<"Constructing a matrix of size N = "<<N<<endl;
  _matrix=d2t(_size,_size);

  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      _matrix[i][j]=matrix[i][j];
    }
  }
}

Matrix::~Matrix() //destructor
{
  if(_N>0) 
    free_d2t(_matrix);
  cout<<"Deleting matrix"<<endl;
  //delete [] _matrix[0];
  //  delete [] _matrix;
}

Matrix::Matrix(Matrix const &src):_DIM(3){   //copy constructor
  _N=src._N;
  _size=_DIM*_N;
  _matrix=d2t(_size,_size);
  cout<<"Copying matrix"<<endl;
  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      _matrix[i][j]=src._matrix[i][j];
    }
  }
  if(src._eigenvalues.size()>0){
    _eigenvalues=src._eigenvalues;
    _eigenvectors=src._eigenvectors;
  }
}



Matrix& Matrix::operator=(Matrix const&src){ //assignment operator
  if (this==&src) return *this;
  if (_N!=src._N){
    if(_N>0) free_d2t(_matrix);
    _N=src._N;
    _size=_DIM*_N;
    _matrix=d2t(_size,_size);
  }
  cout<<"Assigning matrix"<<endl;
  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      _matrix[i][j]=src._matrix[i][j];
    }
  }  
  if(src._eigenvalues.size()>0){
    _eigenvalues=src._eigenvalues;
    _eigenvectors=src._eigenvectors;
  }
  return *this;
}

//#############################################

void Matrix:: SetN(int N){
  //size of matrix
  if(_N==N) return;// AM I SURE THAT THIS IS SAFE??
  if(_size>0) free_d2t(_matrix);// AM I SURE THAT THIS IS SAFE??
  _N=N;
  _size=_DIM*_N;
  cout<<"Constructing a blank matrix of size N = "<<N<<endl;
  //allocating memory and initializing arrays
  _matrix=d2t(_size,_size);  
  return;
}


void Matrix:: SetMatrix(int N, double **matrix){
  if(_size>0) free_d2t(_matrix);// AM I SURE THAT THIS IS SAFE??
  _N=N;
  _size=_DIM*_N;
  cout<<"Matrix:: matrix of size N = "<<N<<endl;
  //  if(_matrix!=0) free_d2t(_matrix);// AM I SURE THAT THIS IS SAFE??
  _matrix=d2t(_size,_size);

  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      _matrix[i][j]=matrix[i][j];
    }
  }
}



 // std::vector<Vector3d> Matrix::GetEigenvec(int i){ 
 //   vector <Vector3d> eigenvec;
 //   for(int j=0;j<_N;++j){
 //     double x=_eigenvectors[i][_DIM*j+0];
 //     double y=_eigenvectors[i][_DIM*j+1];
 //     double z=_eigenvectors[i][_DIM*j+2];
 //     Vector3d ivec(x,y,z);
 //     eigenvec.push_back(ivec);
 //   }
 //   return eigenvec;}



void Matrix::SetEigenvalues(std::vector<double> eval){
  if(_size!=eval.size()){
    cerr<<endl<<"WARNING: the size of the matrix ("<<_size
	<<")does not correspond to the size of the eigenvalues("<<eval.size()<<")"<<endl;
  }
  if(_eigenvalues.size()>0) _eigenvalues.clear();

  try{
    _eigenvalues.reserve(eval.size());
  }catch(std::bad_alloc const&){
    cerr << endl << "NormalMode eigenvalues memory allocation failed!" << endl << flush;
    exit(1);
  }
  for(int i=0;i<eval.size();++i)
    _eigenvalues.push_back(eval[i]);
}


void Matrix::SetEigenvectors(std::vector<std::vector<Vector3d> > evec){ 
  int numModes=evec.size();
  int dim=evec.at(0).size();
  if(_size!=numModes){
    cerr<<endl<<"WARNING: the size of the matrix ("<<_size
	<<")does not correspond to the size of the eigenvectors("<<evec.size()<<")"<<endl;
  }
  if(_eigenvectors.size()>0) _eigenvectors.clear();
  
  try{
    _eigenvectors.resize(numModes);
    for (int i = 0; i < numModes; ++i) {
      _eigenvectors[i].reserve(dim);
      for (int j = 0; j < dim; ++j) {
        Vector3d p3d=evec.at(i).at(j);
        _eigenvectors[i].push_back(p3d);
      }
    }
    
  }catch(std::bad_alloc const&){
    cerr << endl << "NormalMode eigenvector memory allocation failed!" << endl << flush;
    exit(1);
  }  
}


void Matrix::Compose(){
  if(_size>0) free_d2t(_matrix);// AM I SURE THAT THIS IS SAFE??
  _size=_eigenvalues.size();
  if(_size<1){
      cout<<"Matrix.Compose --> ERROR: unvalid size of the matrix. Maybe you forgot to set the eigenvalues"<<endl;
      return;   }
  if(_eigenvectors.size()!=_size){
      cout<<"Matrix.Compose --> ERROR: unvalid size of the matrix. Maybe you forgot to set the eigenvectors"<<endl;
      return;   }
  _matrix=d2t(_size,_size);

  double** evec1=d2t(_size,_size);
  double** evec2=d2t(_size,_size);
  for(int alpha=0; alpha<_size; alpha++){
    for(int i=0; i<_N; i++){
      for(int mu=0;mu<_DIM;mu++){
	evec1[_DIM*i+mu][alpha]=_eigenvectors.at(alpha).at(i).getCoord(mu);
	evec2[alpha][_DIM*i+mu]=_eigenvectors.at(alpha).at(i).getCoord(mu)*_eigenvalues[alpha];
      }
    }
  }

  /* QUESTA PARTE E' LA VERSIONE VECCHIA E LENTA */  
// #pragma omp parallel for 
//   for(int i=0; i<_N; i++){
//     for(int j=i; j<_N; j++){
//       for(int k = 0;k< _size; k++){
//  	for(int mu=0;mu<_DIM;++mu)
//  	  for(int nu=mu;nu<_DIM;++nu)
//  	    _matrix[_DIM*i+mu][_DIM*j+nu] += _eigenvalues[k]*_eigenvectors.at(k).at(i).getCoord(mu)*_eigenvectors.at(k).at(j).getCoord(nu);
//       }//ENDDO K
//       for(int mu=0;mu<_DIM;++mu)
// 	for(int nu=mu+1;nu<_DIM;++nu)
// 	  _matrix[_DIM*i+nu][_DIM*j+mu]=_matrix[_DIM*i+mu][_DIM*j+nu];
//       if (i!=j){ /* symmetry */
//  	for(int mu=0;mu<_DIM;++mu){
//  	  for(int nu=0;nu<_DIM;++nu){
//  	    _matrix[_DIM*j+nu][_DIM*i+mu]=_matrix[_DIM*i+mu][_DIM*j+nu];
//  	    _matrix[j][i]=_matrix[i][j];
//  	}}}
//     }//ENDDO J
//   }//ENDDO I
  

  time_t  begin_time2=time(0);
  multiply_quad_matrix_lapack(evec2,evec1,_size,_matrix);
  cout << "tempo nel ciclo 2 = "<<difftime(time(0),begin_time2)<<endl;  
  free_d2t(evec1);
  free_d2t(evec2);

  // ofstream fout;
  // fout.open("cacca.tmp");
  // for(int i=0; i<_size; i++){
  //   for(int j=i; j<_size; j++){
  //     fout<<C[i][j]<<" "<<_matrix[i][j]<<endl;
  //   }
  // }
  // fout.close();


}//ENDFUNCTIONCOMPOSE

void Matrix::Compose(std::vector<double> eval, std::vector<std::vector<Vector3d> > evec){
  SetEigenvalues(eval);
  SetEigenvectors(evec);
  _swapEigen();//questo e' giusto solo per una covmatrix!!!
  //  cout<<"cacca"<<endl;
  return Compose();
  //  cout<<"cacca"<<endl;
}

//#############################################
 
void Matrix:: Decompose(){
  //  cout<<"The size of the matrix is "<<_size<<endl;

  if(_eigenvalues.size()>0){
    _eigenvalues.clear();
    _eigenvectors.clear();
    cout<<"deleting old eigenstuff"<<endl;
  }
  //  printf("allocating memory space for eigenvalues and eigenvectors\n");
  /* Use lapack routines: lapack_matrix_routines_wrapper_v3.c */
  
  double *temp_eigenvalues=d1t(_size);
  double **temp_eigenvectors=d2t(_size,_size);
  //  printf("taking the spectral decomposition...\n");
  //  cout<<_matrix[0][0]<<endl;
  decompose_symmetric_matrix_lapack(_matrix,_size,temp_eigenvalues,temp_eigenvectors);
  cout<<"done!"<<endl;

  try{
    _eigenvalues.reserve(_size);
  }catch(std::bad_alloc const&){
    cerr << endl << "eigenvalues memory allocation failed!" << endl << flush;
    exit(1);
  }  
  for(int i=0;i<_size;++i){
    _eigenvalues.push_back(temp_eigenvalues[i]);
  }  
  try{
    _eigenvectors.resize(_size);
    for (int i = 0; i < _size; ++i) {      
      _eigenvectors[i].reserve(_N);
      for (int j = 0; j < _N; ++j) {
	double x=temp_eigenvectors[i][_DIM*j+0];
	double y=temp_eigenvectors[i][_DIM*j+1];
	double z=temp_eigenvectors[i][_DIM*j+2];
	Vector3d p3d(x,y,z);
	_eigenvectors[i].push_back(p3d);	
      }      
    }    
  }catch(std::bad_alloc const&){
    cerr << endl<<"NormalMode eigenvector memory allocation failed!" << endl << flush;
    exit(1);
  }    
  // for(int i=0;i<_size;++i){
  //   std::vector<Vector3d> eigenvec;
  //   for(int j=0;j<_N;++j){
  //     double x=temp_eigenvectors[i][_DIM*j+0];
  //     double y=temp_eigenvectors[i][_DIM*j+1];
  //     double z=temp_eigenvectors[i][_DIM*j+2];
  //     Vector3d ivec(x,y,z);
  //     eigenvec.push_back(ivec);
  //   }
  //   _eigenvalues.push_back(temp_eigenvalues[i]);
  //   _eigenvectors.push_back(eigenvec);
  // }
  free_d1t(temp_eigenvalues);
  free_d2t(temp_eigenvectors);
}
//###########################################

void Matrix::dumpAll(const char *name){
  dumpEigenvalues(name);
  dumpEigenvectors(name);
  dumpFullMatrix(name);
  dumpReducedMatrix(name);
  dumpNormalizedReducedMatrix(name);
}//END FUNCTION DUMP



void Matrix:: Times(double a){
  /* This function multiplies the matrix for a certain scalar value.
     (be carefull. It modifies the matrix itself!!) */
  for(int i=0; i<_size; i++){
    for(int j=0; j<_size; j++){
      _matrix[i][j]*=a;
    }
  }
  if(_eigenvalues.size()>0){
    for(int i=0;i<_size;i++){
      _eigenvalues[i]*=a;
    }//enddo i
  }//endif
}//end function

void Matrix:: Plus(const Matrix &M){
  for(int i=0; i<_size; i++){
    for(int j=0; j<_size; j++){
      _matrix[i][j]+=M._matrix[i][j];
    }
  }
  if(_eigenvalues.size()>0){
    _eigenvalues.clear();
    _eigenvectors.clear();
    cout<<"Deleting eigenstuff: you need to compute them again"<<endl;
  }//endif
}//end function


//###########################################
void Matrix:: Print(){
  for(int i=0; i<_size; i++){
    for(int j=0; j<_size; j++){
      cout<<i<<" "<<j<<" "<<_matrix[i][j]<<endl;
    }//enddo
    cout<<endl;
  }//enddo
}//END FUNCTION PRINT


//TODO you have to decompose before invert. Can I modify this so that it check if it was done?
Matrix Matrix::GetInverse(){

  time_t begin_time2=time(0);
  Matrix tempMat(_N);
  vector<double>eval_tmp;
  vector<vector<Vector3d> >evec_tmp;
  int nzero=0;
  for(int k=0; k <_size; ++k){
    double eval= _eigenvalues[_size-k-1];
    if(fabs(eval)<TOL){++nzero; continue;} //skip zero eigenvals
    double tmp_val=1/eval; //evals of the inverse matrix;
    vector<Vector3d> tmp_vec=_eigenvectors[_size-k-1];
    eval_tmp.push_back(tmp_val);
    evec_tmp.push_back(tmp_vec);    
  }
  for(int k=0; k<nzero; ++k) {
    vector<Vector3d> tmp_vec=_eigenvectors[_size-k-1];
    eval_tmp.push_back(0.0);
    evec_tmp.push_back(tmp_vec);        
  } //put the zero eigenvalues at the end
  tempMat.Compose(eval_tmp,evec_tmp);
  cout << "GetInverse: tempo nel ciclo 2 = "<<difftime(time(0),begin_time2)<<endl;  
  return tempMat;
}

double Matrix::GetTrace(){
  double trace=0;
  for (int i=0; i<_size; ++i){
    trace+=_matrix[i][i];
  }
  return trace;
}

double Matrix::GetLogDet(){
  //  if(_eigenvalues.size()<1)
  double logdet=0;
  for (int i=0; i<_eigenvalues.size(); ++i){
    if(_eigenvalues[i]>TOL)
      logdet+=log(_eigenvalues.at(i));
  }
  return logdet;
}


//###################################################################
void Matrix::dumpEigenvalues(const char *name){
  char filename[200];
  FILE *fp;
  sprintf(filename,"%s_eigenvalues.dat",name);
  fp =open_file_w(filename);
  for(int i=0; i <_size; i++){
    fprintf(fp,"%4d %e\n",i,_eigenvalues.at(_size-i-1));
  }
  fclose(fp);
  printf("Eigenvalues written to file %s.\n",filename);
}
/******************************************/
void Matrix::dumpEigenvectors(const char *name){
  dumpTopVectors(-1,name);
}
/******************************************/
void Matrix::dumpTopVectors(int ntop,const char *name){
   char filename[200];
   FILE *fp;
   if (ntop<0) ntop=_size;
   for(int i=0; i < ntop; i++){
     if (i >= _size) break;
     sprintf(filename,"%s_eigenvector_%d.dat",name,i);
     fp =open_file_w(filename);
     for(int j=0; j < _N; j++){
       fprintf(fp,"%4d ",i);
       // for(int mu=0; mu < 3; mu++){
       //   fprintf(fp,"%lf ",_eigenvectors[_size-i-1][3*j+mu]);
       // }
       fprintf(fp,"%lf ",_eigenvectors.at(_size-i-1).at(j).X);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-i-1).at(j).Y);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-i-1).at(j).Z);

       fprintf(fp,"\n");
     }
     fclose(fp);
     printf("Eigenvector number %4d written to file %s\n",i,filename);
   }
}
/******************************************/
void Matrix::dumpFullMatrix(const char *name){
  char filename[200];
  //  FILE *fp;
  sprintf(filename,"%s_full_matrix.dat",name);
  ofstream fp;
  fp.open(filename);
  //  cout<<"cacca"<<endl;
  //  fp =open_file_w(filename);
  for (int mu=0; mu < _DIM; mu++){
    for (int nu=0; nu < _DIM; nu++){
      for(int i=0; i < _N; i++){
        for(int j=0; j < _N; j++){
	  //	  cout<<_DIM*i+mu<<" "<<_DIM*j+nu<<endl;
	  //	  cout<<_matrix[0][0]<<endl;
	  //	  cout<<_matrix[_DIM*i+mu][_DIM*j+nu]<<endl;
	  //          fprintf(fp,"%4d %4d %4d %4d %e\n",i,j,mu,nu,_matrix[_DIM*i+mu][_DIM*j+nu]);
	  fp<<i<<" "<<j<<" "<<mu<<" "<<nu<<" "<<_matrix[_DIM*i+mu][_DIM*j+nu]<<endl;
        }//enddo j
      }//enddo i
    }//enddo nu
  }//enddo mu
  //  fclose(fp);
  fp.close();
  printf("Full matrix written to file %s\n",filename);
}
/******************************************/
void Matrix::dumpReducedMatrix(const char *name){
  char filename[200];
  FILE *fp;
  sprintf(filename,"%s_reduced_matrix.dat",name);
  fp =open_file_w(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      double temp=0;
      for (int mu=0; mu < _DIM; mu++){
        temp += _matrix[_DIM*i+mu][_DIM*j+mu];
      }
      if(i==j) temp=0;
      fprintf(fp,"%4d %4d %e\n",i,j,temp);
    }
  }
  fclose(fp);
  printf("Reduced matrix written to file %s\n",filename);
}
/******************************************/
void Matrix::dumpNormalizedReducedMatrix(const char *name){
  char filename[200];
  FILE *fp;
  sprintf(filename,"%s_normalised_reduced_matrix.dat",name);
  fp =open_file_w(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      double temp=0.0;
      double norm1=0.0;
      double norm2=0.0;
      for (int mu=0; mu < _DIM; mu++){
        temp  += _matrix[_DIM*i+mu][_DIM*j+mu];
        norm1 += _matrix[_DIM*i+mu][_DIM*i+mu];
        norm2 += _matrix[_DIM*j+mu][_DIM*j+mu];
      }
      fprintf(fp,"%4d %4d %e\n",i,j,temp/(sqrt(norm1*norm2)));
    }
  }
  fclose(fp);
  printf("Normalised reduced matrix written to file %s\n",filename);
}
/******************************************/


//####################################################################
//####################################################################
//####################################################################

IntMatrix::IntMatrix(): Matrix() {}
IntMatrix::IntMatrix(int i): Matrix(i) {}
IntMatrix::IntMatrix(int i, double **d): Matrix(i,d) {}

void IntMatrix::Decompose(){
  Matrix::Decompose();
  _swapEigen();
}

void IntMatrix::Compose(std::vector<double> eval, std::vector<std::vector<Vector3d> > evec){
  SetEigenvalues(eval);
  SetEigenvectors(evec);
  //  _swapEigen();//questo e' giusto solo per una covmatrix!!!
  //  cout<<"cacca"<<endl;
  return Matrix::Compose();
  //  cout<<"cacca"<<endl;
}


void Matrix::_swapEigen(){
  /* this function revert the order of the eigenvalues and eigenvectors */
  double tempval;
  //  double *tempvec;
  std::vector<Vector3d> tempvec;
  
  for(int i=0; i<_size/2; ++i){
    // switch the eigenvalues
    tempval=_eigenvalues[i];
    _eigenvalues[i]=_eigenvalues[_size-1-i];
    _eigenvalues[_size-1-i]=tempval;
    //switch the eigenvectors
    tempvec=_eigenvectors[i];
    _eigenvectors[i]=_eigenvectors[_size-1-i];
    _eigenvectors[_size-1-i]=tempvec;
  }
}

CovMatrix IntMatrix::GetInverse(){
  time_t begin_time2=time(0);
  CovMatrix tempMat(_N);
  vector<double>eval_tmp;
  vector<vector<Vector3d> >evec_tmp;
  int nzero=0;
  for(int k=0; k <_size; ++k){
    double eval= _eigenvalues[_size-k-1];
    if(fabs(eval)<TOL){++nzero; continue;} //skip zero eigenvals
    double tmp_val=1/eval; //evals of the inverse matrix;
    vector<Vector3d> tmp_vec=_eigenvectors[_size-k-1];
    eval_tmp.push_back(tmp_val);
    evec_tmp.push_back(tmp_vec);    
  }
  for(int k=0; k<nzero; ++k) {
    vector<Vector3d> tmp_vec=_eigenvectors[_size-k-1];
    eval_tmp.push_back(0.0);
    evec_tmp.push_back(tmp_vec);        
  } //put the zero eigenvalues at the end
  //  cout<<"%@$%#@$%@#%@#size = "<<eval_tmp.size()<<" "<<_size<<endl;
  tempMat.Compose(eval_tmp,evec_tmp);
  cout << "GetInverse: tempo nel ciclo 2 = "<<difftime(time(0),begin_time2)<<endl;  
  return tempMat;
}


void IntMatrix::dumpMSF(const char*name){
  char filename[200];
  ofstream fout;
  sprintf(filename,"%s_mean_square_displ.dat",name);
  fout.open(filename);
  for(int i=0; i < _N; i++){
    double temp=GetMSF(i);
    //    fprintf(fp,"%4d %e\n",i,temp);
    fout<<i<<" "<<temp<<endl;
  }
  fout.close();
  printf("Beads' mean square displacement written to file %s\n",filename);
}

double IntMatrix::GetMSF(int i){
  double msf=GetMSF(i,0,0)+GetMSF(i,1,1)+GetMSF(i,2,2);
  return msf;
}

double IntMatrix::GetMSF(int i,int mu,int nu){
  double msf=0.0;
  if(_eigenvalues.size()<1){
    cerr<<"ERROR: you must compute eigenstuff in order to obtain the MSF from a interaction matrix!"<<endl;
    return 0.0;
  }
  for(int alpha=0; alpha<_size; ++alpha){
    double eval=_eigenvalues.at(alpha);
    if(fabs(eval)<TOL) continue;
    double temp1=_eigenvectors.at(alpha).at(i).getCoord(mu);
    double temp2=_eigenvectors.at(alpha).at(i).getCoord(nu);
    msf+=temp1*temp2/eval;
  }
  return msf;
}


void IntMatrix::dumpCovEigenvalues(const char *name){
  //  double tol=0.000001;
  char filename[200];
  FILE *fp;
  sprintf(filename,"%s_eigenvalues.dat",name);
  fp =open_file_w(filename);
  int nzero=0;
  for(int i=0; i <_size; ++i){
    double eval= _eigenvalues[_size-i-1];
    if(fabs(eval)<TOL){++nzero; continue;} //skip zero eigenvals
    double temp=1/eval; //evals of the inverse matrix;
    fprintf(fp,"%4d %e\n",i,temp);
  }
  for(int j=0; j<nzero; ++j) {fprintf(fp,"%4d 0.0\n",j);} //print the zero eigenvalues at the end
  fclose(fp);
  printf("Eigenvalues of the covariance matrix written to file %s.\n",filename);
}
/******************************************/
void IntMatrix::dumpCovEigenvectors(const char *name){
  dumpTopCovVectors(-1,name);
}
/******************************************/
void IntMatrix::dumpTopCovVectors(int ntop,const char *name){
   char filename[200];
   FILE *fp;
   if (ntop<0) ntop=_size;
   //   for(int i=6; i < ntop+6; i++){
   for(int i=0; i < ntop; i++){
     int ivec=i;
     if (i >= _size-6) ivec-=_size;  //devo stampare alla fine i 6 vettori nulli!
     sprintf(filename,"%s_eigenvector_%d.dat",name,i);
     fp =open_file_w(filename);
     for(int j=0; j < _N; j++){
       //       fprintf(fp,"%4d ",ivec);
       // for(int mu=0; mu < 3; mu++){
       //   //fprintf(fp,"%lf ",_eigenvectors[i][3*j+mu]);
       // 	 fprintf(fp,"%lf ",_eigenvectors[_size-i-1-6][3*j+mu]);
       // }
       fprintf(fp,"%4d ",ivec);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-1-ivec-6).at(j).X);
       fprintf(fp,"\n");
       fprintf(fp,"%4d ",ivec);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-1-ivec-6).at(j).Y);
       fprintf(fp,"\n");
       fprintf(fp,"%4d ",ivec);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-1-ivec-6).at(j).Z);
       fprintf(fp,"\n");
     }
     fclose(fp);
     printf("Eigenvector number %4d written to file %s\n",i,filename);
   }
}


//####################################################################
//####################################################################
//####################################################################

CovMatrix::CovMatrix():Matrix(){}
CovMatrix::CovMatrix(int i): Matrix(i){}
CovMatrix::CovMatrix(int i, double **d): Matrix(i,d){}

void CovMatrix::dumpMSF(const char*name){
  char filename[200];
  FILE *fp;
  sprintf(filename,"%s_mean_square_displ.dat",name);
  fp =open_file_w(filename);
  for(int i=0; i < _N; i++){
    double temp=GetMSF(i);
    fprintf(fp,"%4d %e\n",i,temp);
  }
  fclose(fp);
  printf("Beads' mean square displacement written to file %s\n",filename);
}


double CovMatrix::GetMSF(int i){
  double temp=GetMSF(i,0,0)+GetMSF(i,1,1)+GetMSF(i,2,2);
  return temp;}

double CovMatrix::GetMSF(int i, int mu, int nu){
  double temp= _matrix[_DIM*i+mu][_DIM*i+nu];
  return temp;}





IntMatrix CovMatrix::GetInverse(){
  IntMatrix tempMat(_N);
  time_t begin_time2=time(0);
  vector<double>eval_tmp;
  vector<vector<Vector3d> >evec_tmp;
  int nzero=0;
  for(int k=0; k <_size; ++k){
    double eval= _eigenvalues[k];
    if(fabs(eval)<TOL){++nzero; continue;} //skip zero eigenvals
    double tmp_val=1/eval; //evals of the inverse matrix;
    vector<Vector3d> tmp_vec=_eigenvectors[k];
    eval_tmp.push_back(tmp_val);
    evec_tmp.push_back(tmp_vec);    
  }
  for(int k=0; k<nzero; ++k) {
    vector<Vector3d> tmp_vec=_eigenvectors[k];
    //    eval_tmp.push_back(0.0);
    //    evec_tmp.push_back(tmp_vec);        
    eval_tmp.push_back(0.0);
    evec_tmp.push_back(tmp_vec);        
  } //put the zero eigenvalues at the end
  tempMat.Compose(eval_tmp,evec_tmp);
  cout << "GetInverse: tempo nel ciclo 2 = "<<difftime(time(0),begin_time2)<<endl;  
  return tempMat;
} 


void CovMatrix::Reduce(int NTOP, const CovMatrix&D){
  // reduces the matrix using the first NTOP eigenvectors of matrix D;
  // Matrix D has to be already decomposed;
  //  if D.GetEigenval
  double** temp_mat=d2t(NTOP,NTOP);

  for(int j=0;j<NTOP;++j){
    //vector<double> CQ_m_j;
    double * CQ_M_j=d1t(_size);
    for(int m=0;m<_size;++m){
      for(int l=0;l<_N;++l){
	CQ_M_j[m]+=_matrix[m][_DIM*l+0]*D._eigenvectors[_size-j-1].at(l).X;
	CQ_M_j[m]+=_matrix[m][_DIM*l+1]*D._eigenvectors[_size-j-1].at(l).Y;
	CQ_M_j[m]+=_matrix[m][_DIM*l+2]*D._eigenvectors[_size-j-1].at(l).Z;
      }//enddo l
    }//enddo m
    for(int i=0;i<NTOP;++i){
      for(int p=0;p<_N;++p){
	temp_mat[i][j]+=CQ_M_j[_DIM*p+0]*D._eigenvectors[_size-i-1].at(p).X;
	temp_mat[i][j]+=CQ_M_j[_DIM*p+1]*D._eigenvectors[_size-i-1].at(p).Y;
	temp_mat[i][j]+=CQ_M_j[_DIM*p+2]*D._eigenvectors[_size-i-1].at(p).Z;
      }//enddo i
    }//enddo p
  }//enddo j

  //Now I reset the matrix to its new value
  free_d2t(_matrix);
  _size=NTOP;
  _N=_size/_DIM;
  _matrix=d2t(NTOP,NTOP);
  if(_eigenvalues.size()>0){
    _eigenvectors.clear();
    _eigenvalues.clear();
  }
  for(int i=0;i<_size;++i){
    for(int j=0;j<_size;++j){
      _matrix[i][j]=temp_mat[i][j];
    }
  }

}


double RWSIP(const CovMatrix &A,const CovMatrix &B){
  bool NULL_EVAL=false;
  double RWSIP=0.0;
  cout<<"### RWSIP ###"<<endl;

  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  
  double *vec_a=d1t(nmodes*nmodes);
  double *vec_b=d1t(nmodes*nmodes);
  for (int n = 0; n < nmodes; ++n) {
    for (int m = 0; m < nmodes; ++m) {
      RWSIP+=A._matrix[n][m]*B._matrix[n][m];
      //      vec_a[nmodes*n+m]=A._matrix[n][m];
      //      vec_b[nmodes*n+m]=B._matrix[n][m];
    }
  }
  //  int inc=1;
  //  RWSIP=ddot(nmodes*nmodes,vec_a,inc,vec_b,inc);
  //  RWSIP=multiply_vector_lapack(vec_a,vec_b,nmodes*nmodes);
  double NORM=0.0;
  for (int m = 6; m < nmodes; ++m) {
    double eval_A=A._eigenvalues[m];
    double eval_B=B._eigenvalues[m];
    if(eval_A<TOL) NULL_EVAL=true;
    if(eval_B<TOL) NULL_EVAL=true;
    NORM+=eval_A*eval_B;
  }
  if(NULL_EVAL) return -999;
  RWSIP/=NORM;
  RWSIP=sqrt(RWSIP);
  return RWSIP; 
}

double RWSIP(const IntMatrix &A,const CovMatrix &B){
  //non andrebbe usato!!
  double RWSIP=0.0;
  bool NULL_EVAL=false;
  cout<<"### RWSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  
  double** evecA=d2t(nmodes,nmodes);
  double* evalA=d1t(nmodes);
  double** evecB=d2t(nmodes,nmodes);
  double* evalB=d1t(nmodes);
  for(int alpha=0; alpha<nmodes; alpha++){
    evalA[alpha]=A._eigenvalues[alpha];
    evalB[alpha]=B._eigenvalues[alpha];
    for(int i=0; i<n_beads; i++){
      for(int mu=0;mu<3;mu++){
	Vector3d dumvec=A._eigenvectors.at(alpha).at(i);
	evecA[alpha][3*i+mu]=dumvec.getCoord(mu);
	dumvec=B._eigenvectors.at(alpha).at(i);
	evecB[alpha][3*i+mu]=dumvec.getCoord(mu);
      }
    }
  }
   for (int m = 6; m < nmodes; ++m) {
     double eval_A=evalA[m-6];
     if(eval_A<TOL){
       NULL_EVAL=true;     
     }
     eval_A=1./eval_A;
     double temp_rwsip=0.0;
     for (int n = 6; n < nmodes; ++n) {
       double eval_B=evalB[n];
       double scal_prod=0;
       for(int i=0; i<nmodes; ++i){
 	scal_prod+=evecA[m-6][i]*evecB[n][i];
       }
       temp_rwsip+=eval_A*eval_B*scal_prod*scal_prod;
     }
     double eval_B_m=evalB[m];    
     if(eval_B_m<TOL){
       // #pragma omp critical
       NULL_EVAL=true;
     }
     RWSIP+=temp_rwsip;
     NORM+=eval_A*eval_B_m;
   }
  free_d2t(evecA);
  free_d2t(evecB);
  free_d1t(evalA);
  free_d1t(evalB);
   if(NULL_EVAL) return -999;
   RWSIP/=NORM;
   RWSIP=sqrt(RWSIP);
   return RWSIP; 
}

double RWSIP(const IntMatrix &A,const IntMatrix &B){
  bool NULL_EVAL=false;
  double RWSIP=0.0;
  cout<<"### RWSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
#pragma omp parallel for
  for (int m = 6; m < nmodes; ++m) {
    double eval_A=A._eigenvalues[m-6];
    if(eval_A<TOL) NULL_EVAL=true;
    eval_A=1./eval_A;
    vector<Vector3d> evec_A=A._eigenvectors[m-6];
    for (int n = 6; n < nmodes; ++n) {
      double eval_B=B._eigenvalues[n-6];
      if(eval_B<TOL) NULL_EVAL=true;
      eval_B=1./eval_B;
      vector<Vector3d> evec_B=B._eigenvectors[n-6];
      double scal_prod=0;
      for(int i=0; i<n_beads; ++i){
	scal_prod+=evec_A[i].dot(evec_B[i]);
      }
      double cacca=eval_A*eval_B;
      RWSIP+=cacca*scal_prod*scal_prod;
    }
    double eval_B_m=B._eigenvalues[m-6];
    NORM+=eval_A*eval_B_m;
  }
  if(NULL_EVAL) return -999;
  RWSIP/=NORM;
  RWSIP=sqrt(RWSIP);
  return RWSIP; 
}


double RMSIP(const CovMatrix &A,const CovMatrix &B,int ntop){
  double RMSIP=0.0;
  cout<<"### RMSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  if(ntop>nmodes) return -666;
  for (int m = 0; m < ntop; ++m) {
    double eval_A=A._eigenvalues[nmodes-1-m];
    if(eval_A<TOL) return -999;
    vector<Vector3d> evec_A=A._eigenvectors[nmodes-1-m];
    for (int n = 0; n < ntop; ++n){
      double eval_B=B._eigenvalues[nmodes-1-n];
          if(eval_B<TOL) return -999;
      vector<Vector3d> evec_B=B._eigenvectors[nmodes-1-n];
        double scal_prod=0;
        for(int i=0; i<n_beads; ++i){
          scal_prod+=evec_A[i].dot(evec_B[i]);
        }
        RMSIP+=scal_prod*scal_prod;
      }
    }
    RMSIP/=ntop;
    RMSIP=sqrt(RMSIP);
    return RMSIP; 
}

double RMSIP(const IntMatrix &A,const CovMatrix &B, int ntop){
  double RMSIP=0.0;
  cout<<"### RMSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  if(ntop>nmodes) return -666;
  for (int m = 0; m < ntop; ++m) {
    double eval_A=A._eigenvalues[nmodes-1-m-6];
    if(eval_A<TOL) return -999;
    vector<Vector3d> evec_A=A._eigenvectors[nmodes-1-m-6];
    for (int n = 0; n < ntop; ++n) {
      double eval_B=B._eigenvalues[nmodes-1-n];
      if(eval_B<TOL) return -999;
      vector<Vector3d> evec_B=B._eigenvectors[nmodes-1-n];
      double scal_prod=0;
      for(int i=0; i<n_beads; ++i){
	scal_prod+=evec_A[i].dot(evec_B[i]);
      }
      RMSIP+=scal_prod*scal_prod;
    }
  }
  RMSIP/=ntop;
  RMSIP=sqrt(RMSIP);
  return RMSIP; 
}

double RMSIP(const IntMatrix &A,const IntMatrix &B, int ntop){
  double RMSIP=0.0;
  cout<<"### RMSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  if(ntop>nmodes) return -666;
  for (int m = 0; m < ntop; ++m) {
    double eval_A=A._eigenvalues[nmodes-1-m-6];
    if(eval_A<TOL) return -999;
    vector<Vector3d> evec_A=A._eigenvectors[nmodes-1-m-6];
    for (int n = 0; n < ntop; ++n) {
      double eval_B=B._eigenvalues[nmodes-1-n-6];
      if(eval_B<TOL) return -999;
      vector<Vector3d> evec_B=B._eigenvectors[nmodes-1-n-6];
      double scal_prod=0;
      for(int i=0; i<n_beads; ++i){
	scal_prod+=evec_A[i].dot(evec_B[i]);
      }
      RMSIP+=scal_prod*scal_prod;
    }
  }
  RMSIP/=ntop;
  RMSIP=sqrt(RMSIP);
  return RMSIP; 
}


Matrix dot_product(const Matrix &A, const Matrix &B){
  if((A._N)!=(B._N)) {
    Matrix C(0); 
    return C;
  }
  Matrix C(A._N);
  for(int i=0;i<A._size;++i){
    for(int j=0;j<A._size;++j){
      for(int k=0;k<A._size;++k){
	C._matrix[i][j]+=A._matrix[i][k]*B._matrix[k][j];
      }
    }
  }
  return C;
}//endfunction dotproduct

Matrix dot_product(const IntMatrix &A, const Matrix &B){
  if((A._N)!=(B._N)) {
    Matrix C(0); 
    return C;
  }
  Matrix C(A._N);
  for(int i=0;i<A._size;++i){
    for(int j=0;j<A._size;++j){
      for(int k=0;k<A._size;++k){
	C._matrix[i][j]+=A._matrix[i][k]*B._matrix[k][j];
      }
    }
  }
  return C;
}//endfunction dotproduct

Matrix dot_product(const Matrix &A, const IntMatrix &B){
  if((A._N)!=(B._N)) {
    Matrix C(0); 
    return C;
  }
  Matrix C(A._N);
  for(int i=0;i<A._size;++i){
    for(int j=0;j<A._size;++j){
      for(int k=0;k<A._size;++k){
	C._matrix[i][j]+=A._matrix[i][k]*B._matrix[k][j];
      }
    }
  }
  return C;
}//endfunction dotproduct

Matrix dot_product(const IntMatrix &A, const IntMatrix &B){
  if((A._N)!=(B._N)) {
    Matrix C(0); 
    return C;
  }
  Matrix C(A._N);
  for(int i=0;i<A._size;++i){
    for(int j=0;j<A._size;++j){
      for(int k=0;k<A._size;++k){
	C._matrix[i][j]+=A._matrix[i][k]*B._matrix[k][j];
      }
    }
  }
  return C;
}//endfunction dotproduct


double ComputeBhatt(CovMatrix A,CovMatrix B){
  double Bhatt=0.0;
  cout<<"### Bhatt ###"<<endl;
  double norm_A=  A.GetTrace();
  double norm_B=  B.GetTrace();
  A.Times(1./norm_A);
  B.Times(1./norm_B);
  cout<<"creo D"<<endl;
  CovMatrix D=B;
  cout<<"creata D"<<endl;
  D.Plus(A);
  D.Times(0.5);
  D.Decompose(); // NB: tutto e' ordinato con gli eval crescenti! (ma perche'??!)
      
  int nmodes=D.GetSize();
  // calcolo il numero di autovalori per avere il *nu* per cento di RMSD
  double nu_frac=0.95; //fraction of RMSD to be reproduced by the eigenvalues
  double somma=0.0;
  for(int i=0; i<nmodes; i++){somma+=D.GetEigenval(nmodes-i-1); }
  double partsomma=0.0;
  int NTOP=0;
  for(int i=0;i<nmodes;i++){
    NTOP++;
    partsomma+=D.GetEigenval(nmodes-i-1);
    if(partsomma>nu_frac*somma) break;
  }
  cout<<"NTOP="<<NTOP<<endl;
  NTOP/=3;  NTOP*=3;
  cout<<"NTOP corretto="<<NTOP<<endl;
  
  cout<<"Reducing"<<endl;
  A.Reduce(NTOP,D);
  B.Reduce(NTOP,D);
  A.Decompose();
  B.Decompose();
  CovMatrix Dnew=B;
  Dnew.Plus(A);
  Dnew.Times(0.5);
  Dnew.Decompose();
      
  double log_det_A=0.0;
  double log_det_B=0.0;
  double log_det_D_new=0.0;
  for(int i=0;i<NTOP;++i){
    log_det_A+=log(A.GetEigenval(NTOP-i-1));
    log_det_B+=log(B.GetEigenval(NTOP-i-1));
    log_det_D_new+=log(Dnew.GetEigenval(NTOP-i-1));
  }
  double D_B=(log_det_D_new - log_det_A/2. - log_det_B/2.)/(2*NTOP);
  
  Bhatt=exp(-D_B);
  return Bhatt; 
}
