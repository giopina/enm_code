/* Written by C. Micheletti, michelet@sissa.it  */
/* Last revision Oct 2010 */

#include <fstream>
#include <iostream>
#include <math.h>

#include "kabsch.h"
#include "io.h"
#include "my_malloc.h"
#include "vectop.cc"
#include "myjacobi.cc"

/* ============================================ */
double rmsd_without_alignement(double str1[][3], double str2[][3],int length){
  /* ritorna l'rmsd tra str1 e str2. Non fa nessuna sovrapposizione ottimale */
  double  rmsd, msd, v[3];
  int i,j;

  msd=0.0;
  for(i=0 ; i < length; i++){
    for(j=0; j < 3; j++) {
      v[j]= str1[i][j]-str2[i][j];
      msd+=v[j]*v[j];
    }
  }
  rmsd=sqrt(msd/length);
  return(rmsd);
}

/* ============================================ */

double drms(double str1[][3], double str2[][3],int length){
/* ritorna l'rms delle distanze tra le N*(N-1) coppie di punti delle due strutture */

        double drms, d1,d2;
        int i,j;


        drms=0.0;
        for(i=0 ; i < length; i++){
           for(j=i+1 ; j < length; j++){
              d1 = dist_d(str1[i],str1[j],3);                     
              d2 = dist_d(str2[i],str2[j],3);                     
              drms +=  (d1-d2)*(d1-d2);
           }                       
        }               
        drms = drms*2/(length*(length-1));
        drms = sqrt(drms);                
        return(drms);
}

/* ============================================ */

double  best_rmsd(double str1[][3], double str2[][3],int length)
{
/* ritorna il migliore rmsd (Kabsch) tra str1 e str2. Le coordinate di str1 e str2 vengono preservate  */

        double str3[2000][3], cm1[3],cm2[3], rmsd;
        double u[3][3];

        if (length > 2000) {
                printf("aumenta la size vettore in best_rmsd\n");
                exit(1);
        }       
        rmsd=optimal_alignment(str1,str2,str3,length,u,cm1,cm2);
        return(rmsd);
}

/* ============================================ */
double  kabsch(double str1[][3], double str2[][3], double str3[][3],int length)
{
/* ritorna il migliore rmsd tra str1 e str2 senza allineare niente */

        double rmsd;
        double u[3][3], cm1[3],cm2[3];

        rmsd=optimal_alignment(str1,str2,str3,length, u, cm1,cm2);
        return(rmsd);

}

/* ============================================ */

double  optimal_alignment(double str1[][3], double str2[][3], double str3[][3],int length, double u[][3], double *cm1, double *cm2)   {

/* ritorna il migliore rmsd tra str1 e str2 e mette in str3 la
prima struttura allineata sulla seconda, in u la matrice di rotazione */

  void myjacobi (double a[][3], int n, double *d, double v[][3], int *nrot);

  int i, j, k, sign[3], order[3], nrot;
  double e, e0;
  double r[3][3], rt[3][3], temp, **x, **y;
  double a[3][3], eval[3], evec[3][3];
  double eigenvalues[3], eigenvectors[3][3], b[3][3];


  x = d2t(length,3);
  y = d2t(length,3);

  // zero_vec_d (cm1, 3);
  // zero_vec_d (cm2, 3);
  cm1=d1t(3);
  cm2=d1t(3);


  for (i = 0; i < length; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  cm1[j] += str1[i][j] / length;
	  cm2[j] += str2[i][j] / length;

	}
    }
/*    printf("cm1 %lf %lf %lf\n",cm1[0],cm1[1],cm1[2]);
    printf("cm2 %lf %lf %lf\n",cm2[0],cm2[1],cm2[2]);
*/
  for (i = 0; i < length; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  x[i][j] = str1[i][j] - cm1[j];
	  y[i][j] = str2[i][j] - cm2[j];
	}
    }
    
  /* Mettiamo un pelo di noise per evitare allineamenti perfetti */
    
  e0 = 0.0;
  for (i = 0; i < length; i++)
    {
      e0 += 0.5 * norm_d (x[i], 3) * norm_d (x[i], 3);
      e0 += 0.5 * norm_d (y[i], 3) * norm_d (y[i], 3);
    }

  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  r[i][j] = 0.0;
	  for (k = 0; k < length; k++)
	    {
	      r[i][j] += y[k][i] * x[k][j];
	    }
	  rt[j][i] = r[i][j];
	}
    }



  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  a[i][j] = 0;
	  for (k = 0; k < 3; k++)
	    {
	      a[i][j] += rt[i][k] * r[k][j];
	    }
	}
    }


  myjacobi (a, 3, eval, evec, &nrot);
/* aggiungiamo delle piccole quantita' per rimuovere degenerazioni */
  eval[0]+=0.0000000000000001;
  eval[1]+=0.00000000000000001;
  eval[2]+=0.00000000000000002;
        
  if ((eval[0] < eval[1]) && (eval[0] < eval[2]))
    {
      order[0] = 1;
      order[1] = 2;
      order[2] = 0;
    }

  if ((eval[1] < eval[0]) && (eval[1] < eval[2]))
    {
      order[0] = 0;
      order[1] = 2;
      order[2] = 1;
    }

  if ((eval[2] < eval[0]) && (eval[2] < eval[1]))
    {
      order[0] = 0;
      order[1] = 1;
      order[2] = 2;
    }


  for (i = 0; i < 3; i++)
    {
      eigenvalues[i] = eval[order[i]];
      for (j = 0; j < 3; j++)
	{
	  eigenvectors[i][j] = evec[j][order[i]];
	}
    }


  normalize_d (eigenvectors[0], 3);
  normalize_d (eigenvectors[1], 3);
  vecprod_d (eigenvectors[0], eigenvectors[1], eigenvectors[2]);
  normalize_d (eigenvectors[2], 3);

  for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  b[i][j] = 0;
	  for (k = 0; k < 3; k++)
	    {
	      b[i][j] += r[j][k] * eigenvectors[i][k];
	    }
	}
      normalize_d (b[i], 3);
    }


  vecprod_d (b[0], b[1], b[2]);
  normalize_d (b[2], 3);


  temp = 0.0;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  temp += b[2][i] * r[i][j] * eigenvectors[2][j];
	}
    }
  sign[2] = +1;
  if (temp < 0)
  sign[2] = -1;
  sign[0]=sign[1]=1;
        
  e = e0 - sqrt (eigenvalues[0]) - sqrt (eigenvalues[1]) - sign[2] * sqrt (eigenvalues[2]);



/********************/
  e =2.0 * e / length;
  if (e <0.0) {	
  if (fabs(e) < 1.0e-3) { printf("Warning. In Kabsch alignment found slightly negative value of e (%e). Roundoff error? I will set it equal to zero.\n",e); e=0.0;}/* occasionally, when dealing with two practically identical configurations
				  the value of e may be slightly negative due to the small offsets and roundoff errors. 
				  In this case we set it equal to zero. */
  else {fprintf(stderr,"ERROR. In Kabsch alignment found negative value of e: (%lf)\n",e); exit(1);} 
  }
  e = sqrt (e);
/********************/


  for(i=0; i < 3; i++){
    for(j=0; j < 3; j++){
      u[i][j]=0.0;
        for(k=0; k < 3; k++){
          u[i][j]+= b[k][i]*eigenvectors[k][j];

        }
/*       printf("%10.4lf ",u[i][j]); */

    }
/*        printf("\n"); */
  }

/* allinea la 1 sulla 2 - compreso shift del centro di massa */
  for (i = 0; i < length; i++){
      for (j = 0; j < 3; j++){
        str3[i][j]=cm2[j]; 
        for (k = 0; k < 3; k++){
           str3[i][j]+=u[j][k]*x[i][k];
        }
      }
/*        printf("%8.3lf %8.3lf %8.3lf\n",str3[i][0],str3[i][1],str3[i][2]);
        printf("AAAA %8.3lf %8.3lf %8.3lf\n",x[i][0],x[i][1],x[i][2]);
*/
  }
	
  free_d2t(x);
  free_d2t(y);      
  return (e);

}



/* ============================================ */

double  best_weighted_rmsd(double str1[][3], double str2[][3],double weight[],int length)
{
/* ritorna il migliore rmsd tra str1 e str2 senza allineare niente */

        double str3[2000][3], cm1[3],cm2[3], rmsd;
        double u[3][3];

        rmsd=optimal_weighted_alignment(str1,str2,str3,weight,length,u,cm1,cm2);
        return(rmsd);
}

/* ============================================ */
double  weighted_kabsch(double str1[][3], double str2[][3], double str3[][3],double weight[],int length)
{
/* fa un allineamento ottimale pesato tra str1 e str2 e mette in str3 la prima struttura allineata sulla seconda. Ritorna l'rmsd pesato */

        double rmsd;
        double u[3][3], cm1[3],cm2[3];

        rmsd=optimal_weighted_alignment(str1,str2,str3,weight,length, u, cm1,cm2);
        return(rmsd);

}

/* ============================================ */


double  optimal_weighted_alignment(double str1[][3], double str2[][3], double str3[][3],double *weight,int length, double u[][3], double *cm1, double *cm2)   {

/*fa un allineamento ottimale pesato tra str1 e str2 e mette in str3 la
prima struttura allineata sulla seconda, in u la matrice di rotazione. Ritorna l'rmsd pesato */

/* Corretta il 11/6/2010. Non veniva considerato il peso per calcolare il centro di massa, cioe' l'origine attorno a cui ruotare le due strutture */ 


  void myjacobi (double a[][3], int n, double *d, double v[][3], int *nrot);

  int i, j, k, sign[3], order[3], nrot;
  double e, e0;
  double r[3][3], rt[3][3], temp, **x, **y;
  double a[3][3], eval[3], evec[3][3];
  double eigenvalues[3], eigenvectors[3][3], b[3][3];
  double tot_weight;

  x = d2t(length,3);
  y = d2t(length,3);

  //  zero_vec_d (cm1, 3);
  //  zero_vec_d (cm2, 3);
  cm1=d1t(3);
  cm2=d1t(3);

  tot_weight=0.0;
  for (i = 0; i < length; i++)     tot_weight+=weight[i];

  for (i = 0; i < length; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  cm1[j] += weight[i]*str1[i][j]/tot_weight;
	  cm2[j] += weight[i]*str2[i][j]/tot_weight;

	}
  
    }
  
 
/*    printf("cm1 %lf %lf %lf\n",cm1[0],cm1[1],cm1[2]);
    printf("cm2 %lf %lf %lf\n",cm2[0],cm2[1],cm2[2]);
*/
  for (i = 0; i < length; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  x[i][j] = str1[i][j] - cm1[j];
	  y[i][j] = str2[i][j] - cm2[j];
	}
    }


  e0 = 0.0;
  for (i = 0; i < length; i++)
    {
      e0 += 0.5 * weight[i]*norm_d (x[i], 3) * norm_d (x[i], 3);
      e0 += 0.5 * weight[i]*norm_d (y[i], 3) * norm_d (y[i], 3);
    }

  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  r[i][j] = 0.0;
	  for (k = 0; k < length; k++)
	    {
	      r[i][j] += weight[k]*y[k][i] * x[k][j];
	    }
	  rt[j][i] = r[i][j];
	}
    }



  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  a[i][j] = 0;
	  for (k = 0; k < 3; k++)
	    {
	      a[i][j] += rt[i][k] * r[k][j];
	    }
	}
    }


  myjacobi (a, 3, eval, evec, &nrot);
/* aggiungiamo delle piccole quantita' per rimuovere degenerazioni */
  eval[0]+=0.0000000000000001;
  eval[1]+=0.00000000000000001;
  eval[2]+=0.00000000000000002;
        
  if ((eval[0] < eval[1]) && (eval[0] < eval[2]))
    {
      order[0] = 1;
      order[1] = 2;
      order[2] = 0;
    }

  if ((eval[1] < eval[0]) && (eval[1] < eval[2]))
    {
      order[0] = 0;
      order[1] = 2;
      order[2] = 1;
    }

  if ((eval[2] < eval[0]) && (eval[2] < eval[1]))
    {
      order[0] = 0;
      order[1] = 1;
      order[2] = 2;
    }


  for (i = 0; i < 3; i++)
    {
      eigenvalues[i] = eval[order[i]];
      for (j = 0; j < 3; j++)
	{
	  eigenvectors[i][j] = evec[j][order[i]];
	}
    }


  normalize_d (eigenvectors[0], 3);
  normalize_d (eigenvectors[1], 3);
  vecprod_d (eigenvectors[0], eigenvectors[1], eigenvectors[2]);
  normalize_d (eigenvectors[2], 3);

  for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  b[i][j] = 0;
	  for (k = 0; k < 3; k++)
	    {
	      b[i][j] += r[j][k] * eigenvectors[i][k];
	    }
	}
      normalize_d (b[i], 3);
    }


  vecprod_d (b[0], b[1], b[2]);
  normalize_d (b[2], 3);


  temp = 0.0;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  temp += b[2][i] * r[i][j] * eigenvectors[2][j];
	}
    }
  sign[2] = +1;
  if (temp < 0)
  sign[2] = -1;
  sign[0]=sign[1]=1;
        
  e = e0 - sqrt (eigenvalues[0]) - sqrt (eigenvalues[1]) - sign[2] * sqrt (eigenvalues[2]);



/********************/
  e =2.0 * e / tot_weight;
  if (e <0.0) {	
  if (fabs(e) < 1.0e-3) { printf("Warning. In Kabsch alignment found slightly negative value of e (%e). Roundoff error? I will set it equal to zero.\n",e); e=0.0;}/* occasionally, when dealing with two practically identical configurations
				  the value of e may be slightly negative due to the small offsets and roundoff errors. 
				  In this case we set it equal to zero. */
  else {fprintf(stderr,"ERROR. In Kabsch alignment found negative value of e: (%lf)\n",e); exit(1);} 
  }
  e = sqrt (e);
/********************/



  for(i=0; i < 3; i++){
    for(j=0; j < 3; j++){
      u[i][j]=0.0;
        for(k=0; k < 3; k++){
          u[i][j]+= b[k][i]*eigenvectors[k][j];
        }
/*        printf("%10.4lf ",u[i][j]); */

    }
/*        printf("\n"); */
  }

/* allinea la 1 sulla 2 - compreso shift del centro di massa */
  for (i = 0; i < length; i++){
      for (j = 0; j < 3; j++){
        str3[i][j]=cm2[j]; 
        for (k = 0; k < 3; k++){
           str3[i][j]+=u[j][k]*x[i][k];
        }
      }
/*        printf("%8.3lf %8.3lf %8.3lf\n",str3[i][0],str3[i][1],str3[i][2]);
        printf("AAAA %8.3lf %8.3lf %8.3lf\n",x[i][0],x[i][1],x[i][2]);
*/
  }
	
        free_d2t(x);        
        free_d2t(y);
  return (e);

}



