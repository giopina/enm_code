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

#include "Matrix.h"
#include "my_malloc.h"
#include "io.h"
//#include "Vector3d.h"

#include "lapack_matrix_routines_wrapper_v3.cc"

using namespace std;

//Constructors for the class CovMatrix1d
CovMatrix1d:: CovMatrix1d():CovMatrix(){ //forse sta cosa non si dovrebbe fare. Se non dichiaro la dimensione della matrice non la creo...
  _DIM=1;
}//END OF CONSTRUCTOR

CovMatrix1d:: CovMatrix1d(int N):CovMatrix(){
  //size of matrix
  _DIM=1;
  _N=N;
  _size=_DIM*_N;
  _matrix=d2t(_size,_size);  
}//END OF CONSTRUCTOR


CovMatrix1d:: CovMatrix1d(int N, double **matrix):CovMatrix(){
  _DIM=1;
  _N=N;
  _size=_DIM*_N;
  _matrix=d2t(_size,_size);
  for(int i=0; i<_size; ++i){
    for(int j=0; j<_size; ++j){
      _matrix[i][j]=matrix[i][j];
    }//enddo
  }//enddo
}//END OF CONSTRUCTOR


//#############################################

void CovMatrix1d::SetEigenvectors(std::vector<std::vector<double> > evec){ 
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
        double p1d=evec.at(i).at(j);
        _eigenvectors[i].push_back(p1d);
      }
    }
    
  }catch(std::bad_alloc const&){
    cerr << endl << "NormalMode eigenvector memory allocation failed!" << endl << flush;
    exit(1);
  }  
}


void CovMatrix1d::Compose(){
  if(_size>0) free_d2t(_matrix);// AM I SURE THAT THIS IS SAFE??
  _size=_eigenvalues.size();
  if(_size<1){
      cout<<"Matrix.Compose --> ERROR: unvalid size of the matrix. Maybe you forgot to set the eigenvalues"<<endl;
      return;   }
  if(_eigenvectors.size()!=_size){
      cout<<"Matrix.Compose --> ERROR: unvalid size of the matrix. Maybe you forgot to set the eigenvectors"<<endl;
      return;   }
  _matrix=d2t(_size,_size);
  // for(int i=0; i<_size; i++){
  //   for(int j=i; j<_size; j++){
  for(int i=0; i<_N; i++){
    for(int j=i; j<_N; j++){
      /* check se deve essere >= k 0 no */
      for(int k = 0;k< _size; k++){
	// imnotdoingthis?? /* Run the sum backwards to minimize roundoff errors in the sum of inverse eigenvalues */
	//        if (fabs(_eigenvalues[k]) < tol) continue;
	
        _matrix[i][j] += _eigenvalues[k]*_eigenvectors.at(k).at(i)*_eigenvectors.at(k).at(j);

      }
      if (i!=j){ /* symmetry */
	_matrix[j][i]=_matrix[i][j];	
      }//endif
    }//enddo j
  }//endoi
}//end function

void CovMatrix1d::Compose(std::vector<double> eval, std::vector<std::vector<double> > evec){
  SetEigenvalues(eval);
  SetEigenvectors(evec);
  _swapEigen();
  return Compose();
}

//#############################################
 
void CovMatrix1d:: Decompose(){

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
  //  cout<<"done!"<<endl;

  for(int i=0;i<_size;++i){
    std::vector<double> eigenvec;
    for(int j=0;j<_N;++j){
      double ivec=temp_eigenvectors[i][j];
      eigenvec.push_back(ivec);
    }
    _eigenvalues.push_back(temp_eigenvalues[i]);
    _eigenvectors.push_back(eigenvec);
  }
  free_d1t(temp_eigenvalues);
  free_d2t(temp_eigenvectors);
}

//TODO you have to decompose before invert. Can I modify this so that it check if it was done?
Matrix CovMatrix1d::GetInverse(){
  //TODO: sta cosa non so se mi serve in fondo...
  Matrix tempMat(0);
  return tempMat;
  // double **invmat=d2t(_size,_size);
  // double *temp_eigenvalues=d1t(_size);
  // double **temp_eigenvectors=d2t(_size,_size);
  // for(int i=0;i<_size;++i){
  //   temp_eigenvalues[i]=_eigenvalues.at(i);
  //   for(int j=0;j<_N;++j){
  //     temp_eigenvectors[i][_DIM*j+0]=_eigenvectors.at(i).at(j).X;
  //     temp_eigenvectors[i][_DIM*j+1]=_eigenvectors.at(i).at(j).Y;
  //     temp_eigenvectors[i][_DIM*j+2]=_eigenvectors.at(i).at(j).Z;
  //   }
  // }    
  // spectral_singular_inversion(invmat,_size,temp_eigenvalues,temp_eigenvectors,TOL);
  // free_d1t(temp_eigenvalues);
  // free_d2t(temp_eigenvectors);
  // //  cout<<"caccacca"<<endl;
  // //  cout<< invmat[0][0]<<endl;
  // Matrix tempMat(_N,invmat);
  // free_d2t(invmat);
  // //Matrix schifo=tempMat;
  // //outMat=schifo;
  // return tempMat;
}

void CovMatrix1d::dumpEigenvectors(const char *name){
  dumpTopVectors(-1,name);
}
/******************************************/
void CovMatrix1d::dumpTopVectors(int ntop,const char *name){
   char filename[200];
   FILE *fp;
   if (ntop<0) ntop=_size;
   for(int i=0; i < ntop; i++){
     if (i >= _size) break;
     sprintf(filename,"%s_eigenvector_%d.dat",name,i);
     fp =open_file_w(filename);
     for(int j=0; j < _N; j++){
       fprintf(fp,"%4d ",i);
       fprintf(fp,"%lf ",_eigenvectors.at(_size-i-1).at(j));
       fprintf(fp,"\n");
     }
     fclose(fp);
     printf("Eigenvector number %4d written to file %s\n",i,filename);
   }
}


void CovMatrix1d::dumpMatrix(const char *name){
  char filename[200];
  //  FILE *fp;
  sprintf(filename,"%s_matrix.dat",name);
  ofstream fp;
  fp.open(filename);
  for(int i=0; i < _N; i++){
    for(int j=0; j < _N; j++){
      fp<<i<<" "<<j<<" "<<_matrix[i][j]<<endl;
    }//enddo j
  }//enddo i
  fp.close();
  printf("Matrix written to file %s\n",filename);
}
/******************************************/


/******************************************/

void CovMatrix1d::_swapEigen(){
  /* this function revert the order of the eigenvalues and eigenvectors */
  double tempval;
  //  double *tempvec;
  std::vector<double> tempvec;
  
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



//####################################################################
//####################################################################
//####################################################################

void CovMatrix1d::dumpMSF(const char*name){
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


double CovMatrix1d::GetMSF(int i){
  double temp= _matrix[i][i];
  return temp;}


//TODO: non so se mi serve sta cosa
//IntMatrix CovMatrix::GetInverse(){
  // double **invmat=d2t(_size,_size);
  // double *temp_eigenvalues=d1t(_size);
  // double **temp_eigenvectors=d2t(_size,_size);
  // for(int i=0;i<_size;++i){
  //   temp_eigenvalues[i]=_eigenvalues.at(i);
  //   for(int j=0;j<_N;++j){
  //     temp_eigenvectors[i][_DIM*j+0]=_eigenvectors.at(i).at(j).X;
  //     temp_eigenvectors[i][_DIM*j+1]=_eigenvectors.at(i).at(j).Y;
  //     temp_eigenvectors[i][_DIM*j+2]=_eigenvectors.at(i).at(j).Z;
  //   }
  // }    
  // spectral_singular_inversion(invmat,_size,temp_eigenvalues,temp_eigenvectors,TOL);
  // free_d1t(temp_eigenvalues);
  // free_d2t(temp_eigenvectors);
  // IntMatrix tempMat(_N,invmat);
  // free_d2t(invmat);
  // return tempMat;
//}


void CovMatrix1d::Reduce(int NTOP, const CovMatrix1d&D){
  // reduces the matrix using the first NTOP eigenvectors of matrix D;
  // Matrix D has to be already decomposed;
  //  if D.GetEigenval
  double** temp_mat=d2t(NTOP,NTOP);

  for(int j=0;j<NTOP;++j){
    //vector<double> CQ_m_j;
    double * CQ_M_j=d1t(_size);
    for(int m=0;m<_size;++m){
      for(int l=0;l<_N;++l){
	CQ_M_j[m]+=_matrix[m][l]*D._eigenvectors[_size-j-1].at(l);
      }//enddo l
    }//enddo m
    for(int i=0;i<NTOP;++i){
      for(int p=0;p<_N;++p){
	temp_mat[i][j]+=CQ_M_j[p]*D._eigenvectors[_size-i-1].at(p);
      }//enddo i
    }//enddo p
  }//enddo j
  //Now I reset the matrix to its new value
  free_d2t(_matrix);
  _size=NTOP;
  _N=_size;
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


double RWSIP(const CovMatrix1d &A,const CovMatrix1d &B){
  double RWSIP=0.0;
  cout<<"### RWSIP ###"<<endl;
  double NORM=0.0;
  int nmodes=A._size;
  int n_beads=nmodes/A._DIM;
  for (int m = 6; m < nmodes; ++m) {
    double eval_A=A._eigenvalues[m];
    vector<double> evec_A=A._eigenvectors[m];
    for (int n = 6; n < nmodes; ++n) {
      double eval_B=B._eigenvalues[n];
      vector<double> evec_B=B._eigenvectors[n];
        double scal_prod=0;
        for(int i=0; i<n_beads; ++i){
          scal_prod+=evec_A[i]*evec_B[i];
        }
        double cacca=eval_A*eval_B;
        RWSIP+=cacca*scal_prod*scal_prod;
      }
      double eval_B_m=B._eigenvalues[m];
      NORM+=eval_A*eval_B_m;
    }
    RWSIP/=NORM;
    RWSIP=sqrt(RWSIP);
    return RWSIP; 
}
