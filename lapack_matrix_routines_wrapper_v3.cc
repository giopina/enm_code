extern "C" {
#include <mkl.h>
}


//#include <fstream>
#include <stdio.h>
#include <math.h>                                                          
//#include <malloc.h>                                                        
//#include <string.h>   
#include <stdlib.h>                                                        


/* // copy this in you .h
extern "C" void dgetrf(int *,int *, double *, int *, int *, int *);
extern "C" void dgetri(int *,double *,int *,int *, double *, int*,int *);
extern "C" void dsyevd(char *,char *,int *, double *, int *,double *,double *, int *, int *,int *, int * );
*/


static int invert_symmetric_matrix_lapack( double **m, int dim_m, double **inv_m){
  
/*
 M is a NON-SINGULAR symmetric matrix of size  (dim_m * dim_m).
 The program calculates the inverse and returns it in inv_m

 written by CM,  20/07/2008
 tested by CM,  21/07/2008
*/

  int c, i,j,k,l,lwork, info;
  double *array_m;
  double *work, tmp1, tmp2, worksize;
  int *iwork;
  int *isupps;

  //  cout<<"cacca2"<<endl;  
  array_m = d1t(dim_m*dim_m);
  //  cout<<"cacca3"<<endl;    

  printf("copio l'array\n");
  c=0;
  for(i=0; i< dim_m; i++){
    for(j=0; j< dim_m; j++){
      //      cout<<i<<j<<endl;
      //      cout<<m[i][j]<<endl;
      /* copy matrix in array */
      array_m[c]=m[i][j];
      c++;
    }
  }

  isupps = i1t(1+dim_m); /* OK for ATLAS and MKL */

  printf("dgetrf\n");
  dgetrf(&dim_m,&dim_m,array_m,&dim_m,isupps,&info);

  if ( info!=0 ) {
    //    fprintf(stderr,"ERROR in dgetrf. Bailing out.\n");
    cout<<"ERROR in dgetrf. Bailing out. Info: "<<info<<endl;
    //    fprintf(stderr,"\tInfo: %d\n", info);
    //    exit(1);
    return -1;
  }


  lwork=-1; /*workspace query */

  printf("dgetri 1\n");
  //lapack function
  dgetri(&dim_m,array_m,&dim_m,isupps,&worksize,&lwork,&info);
  if ( info!=0 ) {
    cout<<"ERROR in workspace query for dgetri. Bailing out."<<endl;
    cout<<"Info:"<<info<<endl;
    //    exit(1);
    return -1;
  }

  lwork=(int) worksize+1;
  printf("dgetri info: %d  ...allocating workspace for work of size: %d (%lf)\n",info,lwork,worksize+1);
  work = d1t(lwork);
  printf("dgetri 2\n");
  //lapack function
  dgetri(&dim_m,array_m,&dim_m,isupps,work,&lwork,&info);
  if ( info!=0 ) {
    cout<<"ERROR in dgetri. Bailing out.Info: "<<info<<endl;
    // exit(1);
    return -1;
  }


  printf("ricopio\n");
  c=0;
  for(i=0; i< dim_m; i++){
    for(j=0; j< dim_m; j++){
      /* copy array in inverse matrix */
      inv_m[i][j]=array_m[c];
      c++;
    }
  }
  printf("...ricopiato!\n");
  
  free_d1t(array_m);
  free_d1t(work);
  free_i1t(isupps);
  
  printf("finito!\n");
  return(0); /* Success */
}

static int multiply_quad_matrix_lapack( double **m_a,double **m_b,int dim_m, double**m_c){
/*
non so cosa cazzo sto facendo
*/
  double *A,*B,*C;
  char TRANS = 'N';
  double ALPHA = 1.0;
  double BETA = 0.0;
  int k=0;
  A = d1t(dim_m*dim_m);
  B = d1t(dim_m*dim_m);
  C = d1t(dim_m*dim_m);
  printf("copio l'array\n");
  k=0;
  for(int i=0; i< dim_m; i++){
    for(int j=0; j< dim_m; j++){
      /* copy matrix in array */
      A[k]=m_a[i][j];
      B[k]=m_b[i][j];
      k++;
    }
  }

  dgemm(&TRANS, &TRANS, &dim_m, &dim_m, &dim_m, &ALPHA, A, &dim_m, B, &dim_m, &BETA, C, &dim_m);
  
  printf("ricopio\n");
  k=0;
  for(int i=0; i< dim_m; i++){
    for(int j=0; j< dim_m; j++){
      /* copy array in inverse matrix */
      m_c[i][j]=C[k];
      k++;
    }
  }
  printf("...ricopiato!\n");
  free_d1t(A);
  free_d1t(B);
  free_d1t(C);  
  printf("finito!\n");
  return(0); /* Success */
}


static double multiply_vector_lapack( double *v_a,double *v_b,int dim_v){
/*
non so cosa cazzo sto facendo
*/
  int inc=1;
  double temp=ddot(&dim_v,v_a,&inc,v_b,&inc);
  return temp;
}



static void decompose_symmetric_matrix_lapack( double **m, int dim_m, double *eigenvalues_m, double **eigenvectors){
  //  cout<<"cacca"<<endl;
  /* M is a SINGULAR (or NON-SINGULAR) symmetric matrix of size dim_m * dim_m.
     The program calculates the eigenvalues, eigenvectors, the pseudoinverse.
     returns the number of zero eigenvalues (i.e. eigenvalues that are less than "tol" in modulus.
  */

  int c,i,j,k,l,lwork, info;
  double *array_m;
  double *work, worksize;
  int liwork, *iwork, iworksize, nzeros;
  char c1, c2;
  //  cout<<"cacca"<<endl;
  array_m = d1t(dim_m*dim_m);
  //  cout<<"cacca"<<endl;
  c=0;
  for(i=0; i< dim_m; i++){
    for(j=0; j< dim_m; j++){
      /* copy matrix in array */
      //      cout<<i<<j<<endl;
      //      cout<<m[i][j]<<endl;
      array_m[c]=m[i][j];
      c++;
    }
  }
  //  printf("...computing size of workspace\n");
  c1='v';
  c2='u';
  lwork=liwork=-1; /* workspace query */

  //### lapack function
  //  printf("dimensione %d\n",dim_m);
  dsyevd( &c1, &c2, &dim_m, array_m, &dim_m, eigenvalues_m,&worksize, &lwork, &iworksize, &liwork, &info );
  if (info!=0) {printf("Error in dsyevd. Workspace query failed. Info value: %d\n",info); exit(1);}
  lwork=(int) (worksize+1);
  liwork=iworksize;
  
  c1='v';
  c2='u';
  //  printf("dsyevd info: %d ...allocating workspace for work of size: %d (%lf) and iwork (size: %d )\n",info,lwork,worksize+1,liwork);
  work=d1t(lwork);
  iwork=i1t(liwork);
  //  printf("...decomposing matrix\n");
  dsyevd( &c1, &c2, &dim_m, array_m, &dim_m, eigenvalues_m,work, &lwork, iwork, &liwork, &info );
  //  printf("done\n");

  if ( info!=0 ) {
    fprintf(stderr,"ERROR in dsyevd.\n");
    fprintf(stderr,"\tInfo: %d\n", info);
    exit(1);
  }
  
  /* eigenvalues, in ascending order, are in eigenvalues_m */
  /* eigenvectors, in ascending order, are in array. */
  
  /* For a given eigenvector we should have x1,y1,z1, x2,y2,z2, etc... 
  */
  for( c=0, j=0; j< dim_m; j++)
    for(l=0; l< dim_m; l++, c++ ) {
      eigenvectors[j][l] = array_m[c];
    }

  free_d1t(work);
  free_i1t(iwork);
  free_d1t(array_m);
  
}


static void spectral_singular_inversion(double **inv_m, int dim_m, double *eigenvalues_m, double **eigenvectors, double tol){

int i, j, k;

  for(i=0; i< dim_m; i++){
    for(j=i; j< dim_m; j++){
      /* copy matrix in array */
      inv_m[i][j]=0;
      /* check se deve essere >= k 0 no */
      for( k = 0;k< dim_m; k++){
        /* Run the sum backwards to minimize roundoff errors in the sum of inverse eigenvalues */
        if (fabs(eigenvalues_m[k]) < tol) continue;
        inv_m[i][j] += 1.0/eigenvalues_m[k]*eigenvectors[k][i]*eigenvectors[k][j];
      }
      if (i!=j) inv_m[j][i]=inv_m[i][j]; /* symmetry */
    }
  }

}
