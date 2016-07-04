/* Written by C. Micheletti, michelet@sissa.it  */
/* Last revision June 2010 */


double  best_rmsd(double str1[][3], double str2[][3],int length);
/* migliore rmsd di allineamento tra 2 strutture */
double  kabsch(double str1[][3], double str2[][3], double str3[][3],int length);
/* da il migliore rmsd e mette in str3 l'allineamento della str1 sulla str2*/
double  optimal_alignment(double str1[][3], double str2[][3], double str3[][3],int length, double u[][3], double *cm1, double *cm2);
/* da il migliore rmsd e mette in str3 l'allineamento della str1 sulla str2,in u la matrice di rotazione e le traslazioni 
del centro di massa da fare per portare la 1 sulla 2 */

double  best_weighted_rmsd(double str1[][3], double str2[][3],double weight[],int length);
/* migliore rmsd pesato di allineamento tra 2 strutture */
double  weighted_kabsch(double str1[][3], double str2[][3], double str3[][3],double weight[],int length);
/* da il migliore rmsd pesato e mette in str3 l'allineamento della str1 sulla str2*/
double  optimal_weighted_alignment(double str1[][3], double str2[][3], double str3[][3],double weight[],int length, double u[][3], double *cm1, double *cm2);
/* da il migliore rmsd pesato e mette in str3 l'allineamento della str1 sulla str2,in u la matrice di rotazione e le traslazioni 
del centro di massa da fare per portare la 1 sulla 2 */
