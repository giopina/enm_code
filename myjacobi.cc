/* Written by C. Micheletti, michelet@sissa.it  */
/* Last revision January 2007 */

/* ============================================ */

void myjacobi (double a[][3], int n, double *d, double v[][3], int *nrot)
{
  
/* modificata 24/1/2007 per eliminare controlli di uguaglianza tra floats */

  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  double b[3], z[3];
  
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
  
  for (ip = 0; ip <= (n - 1); ip++)
    {
      for (iq = 0; iq <= (n - 1); iq++)
	v[ip][iq] = 0.0;
      v[ip][ip] = 1.0;
    }
  for (ip = 0; ip <= (n - 1); ip++)
    {
      b[ip] = d[ip] = a[ip][ip];
      z[ip] = 0.0;
    }
  *nrot = 0;
  for (i = 1; i <= 500; i++)
    {
      sm = 0.0;
      for (ip = 0; ip <= n - 2; ip++)
	{
	  for (iq = ip + 1; iq <= (n - 1); iq++)
	    sm += fabs (a[ip][iq]);
	}
      if (sm == 0.0)
	{
	  return;
	}
      if (i < 4)
	tresh = 0.2 * sm / (n * n);
      else
	tresh = 0.0;
      for (ip = 0; ip <= n - 2; ip++)
	{
	  for (iq = ip + 1; iq <= (n - 1); iq++)
	    {
	      g = 100.0 * fabs (a[ip][iq]);
	      if (i > 4 && (fabs ((fabs (d[ip]) + g) - fabs (d[ip])) < 1.0e-6)
		  && (fabs( (fabs(d[iq]) + g)- fabs (d[iq])) < 1.0e-6))
		a[ip][iq] = 0.0;
	      else if (fabs (a[ip][iq]) > tresh)
		{
		  h = d[iq] - d[ip];
		  if (  fabs((fabs (h) + g) - fabs (h)) < 1.0e-6)
		    t = (a[ip][iq]) / h;
		  else
		    {
		      theta = 0.5 * h / (a[ip][iq]);
		      t = 1.0 / (fabs (theta) + sqrt (1.0 + theta * theta));
		      if (theta < 0.0)
			t = -t;
		    }
		  c = 1.0 / sqrt (1 + t * t);
		  s = t * c;
		  tau = s / (1.0 + c);
		  h = t * a[ip][iq];
		  z[ip] -= h;
		  z[iq] += h;
		  d[ip] -= h;
		  d[iq] += h;
		  a[ip][iq] = 0.0;
		  for (j = 0; j <= ip - 1; j++)
		    {
		      ROTATE (a, j, ip, j, iq)
		    }
		  for (j = ip + 1; j <= iq - 1; j++)
		    {
		      ROTATE (a, ip, j, j, iq)
		    }
		  for (j = iq + 1; j <= (n - 1); j++)
		    {
		      ROTATE (a, ip, j, iq, j)
		    }
		  for (j = 0; j <= (n - 1); j++)
		    {
		      ROTATE (v, j, ip, j, iq)
		    }
		  ++(*nrot);
		}
	    }
	}
      for (ip = 0; ip <= (n - 1); ip++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }
  printf ("Too many iterations in routine JACOBI %lf",sm);
  /*  exit (1); */
#undef ROTATE
}

