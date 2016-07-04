
double scal_d (double *a, double *b, int dim)
{

  int i;
  double temp;

  temp = 0.0;
  for (i = 0; i < dim; i++)
    {
      temp += a[i] * b[i];
    }
  return (temp);
}
/*******************************/

double norm_d (double *a, int dim)
{

  return (sqrt (scal_d (a, a, dim)));
}
/*******************************/

double dist_d (double *a, double *b, int dim)
{

  int i;
  double temp;

  temp = 0.0;
  for (i = 0; i < dim; i++)
    {
      temp += (a[i] - b[i]) * (a[i] - b[i]);
    }

  temp = sqrt (temp);
  return (temp);
}

/*******************************/
void my_invert_symmetric_matrix(double **m,double **invm, int size, double *eigenvalues, double **eigenvectors){

  int p,q,r,nrot,i,j,l, nsweeps;
  double theta, t, c, s, tau, Sum, a, avg_entry;
  double *rowp,*rowq;
  
  rowp = d1t(size);
  rowq = d1t(size);
  
  for(i=0; i < size; i++){      
    eigenvectors[i][i]=1.0;
    for(j=i+1; j < size; j++){    
      eigenvectors[i][j]=0.0;
      eigenvectors[j][i]=0.0;
    }
  }
  
  nrot =0;
  for(nsweeps=0;;nsweeps++){
    if (nsweeps > 50){
      printf("Error. Reached maximum number of sweeps for Jacobi convergence. Exiting...\n");
      exit(1);
      
    }   

    Sum=0.0;
    avg_entry=0.0;
    for(i=0; i < size; i++){    
      for(j=i+1; j < size; j++){    
        Sum +=  2.0 * m[i][j]*m[i][j];      
        avg_entry += 2.0*fabs(m[i][j]);
      }
    }
    avg_entry = avg_entry/(size*(size-1.0));
    
    for(p=0; p < size; p++){
      for(q=p+1; q < size; q++){
        /* For the first sweeps only perform rotations to kill the 
           off-diagonal matrix elements which are larger than average
        */

        if ((nsweeps <3) && (fabs(m[p][q]) < avg_entry)) continue;
        
        theta = (m[q][q]-m[p][p])/(2.0*m[p][q]);
        
        if (fabs(theta) > 1.0e+10) t = 1.0/(2*theta);
        else t = 1.0/(fabs(theta) + sqrt (theta*theta+1.0));
        if (theta < 0.0) t = -t;
        
        c = 1.0/sqrt(t*t+1);
        s = t*c;
        tau = s/(1.0+c);
        
        /* Apply the Jacobi rotation to get the new M matrix */
        for(r=0; r < size; r++){
          rowp[r]=m[r][p];
          rowq[r]=m[r][q];
        }
        
        m[p][p] = rowp[p]- t*rowp[q];
        m[q][q] = rowq[q]+ t*rowq[p];
        
        for(r=0; r < size; r++){
          if (r==p) continue;
          if (r==q) continue;
          m[r][p]= rowp[r] - s*(rowq[r] + tau*rowp[r]);
          m[p][r]= m[r][p];
          m[r][q]= rowq[r] + s*(rowp[r] - tau*rowq[r]);
          m[q][r]= m[r][q];
        }
        m[p][q]=0.0;
        m[q][p]=0.0;
        Sum = Sum - 2*rowp[q]*rowp[q];
        
        
        /* Combine the various rotations to get the matrix of eigenvectors
         */
        
        for(r=0; r < size; r++){
          rowp[r]=eigenvectors[p][r];
          rowq[r]=eigenvectors[q][r];
        }
        
        for(r=0; r < size; r++){
          eigenvectors[p][r] = rowp[r] - s*(rowq[r] + tau*rowp[r]);
          eigenvectors[q][r] = rowq[r] + s*(rowp[r] - tau*rowq[r]);
        }
        
        if (Sum < 1.0e-12) goto end;       
        nrot++;  
      }
      
    }
  }
 end:
  

  /* the diagonal elements of the rotated matrix contain the
     eigenvalues */
  
  for(r=0; r < size; r++){
    eigenvalues[r] = m[r][r];
  }
  
  
  
  /* Construct the inverse matrix by doing the spectral decomposition
     and eliminating the null-eigenvalue space */
  
  for (j=0;j<size;j++) {
    for (i=0;i<size;i++) {
      invm[j][i]=0.0;
      for (l=0;l<size;l++) {
        if (fabs(eigenvalues[l])>0.0000001) invm[j][i]+= 1.0/eigenvalues[l]*eigenvectors[l][j]*eigenvectors[l][i];
      }
    }
  }        
  
  printf("End matrix inversion. Sorting eigenvalues... \n");
  /* Now let's sort the eigenvalues (and eigenvectors) in descending
     order */
  
  
  for(j=1; j < size; j++){
    
    i = j-1;
    a = eigenvalues[j];
    for(r=0; r < size; r++) rowp[r]=eigenvectors[j][r];
    
    while ((i >=0) && (eigenvalues[i] < a)){
      eigenvalues[i+1] = eigenvalues[i];
      for(r=0; r < size; r++){
        eigenvectors[i+1][r]=eigenvectors[i][r];
      }
      i--;
    }
    eigenvalues[i+1]=a;
    for(r=0; r < size; r++){
      eigenvectors[i+1][r]= rowp[r];
    }
  }

 free_d1t(rowp);
 free_d1t(rowq);

}

