/*****        program geneigen        *****/
/*    tridiagonalization and Wilkinson    */
/*  shift approach for symmetric matrices */
/*   t.r.chandrupatla and a.d.belegundu   */
/******************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   double *s,*gm,*d,*b;
   double omega,freq,c;
   int *nord,nq,nbw,i,j,k,n,ii,jn,ifl,iter;
   double pi = 3.14159;
   char dummy[81], file1[81], file2[81];
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr1 = fopen(file1, "r");
   fptr2 = fopen(file2, "w");
   fgets(dummy,80,fptr1);
   fgets(dummy,80,fptr1);
/* --- read in number of equations and bandwidth --- */
   fscanf(fptr1,"%d %d\n", &nq, &nbw);
/* -----         memory allocation         ----- */
   s = (double *) calloc(nq*nq, sizeof(double));
   gm = (double *) calloc(nq*nq, sizeof(double));
   d = (double *) calloc(nq, sizeof(double));
   b = (double *) calloc(nq*nq, sizeof(double));
   nord = (int *) calloc(nq, sizeof(int));
/* --------------------------------------------- */
   for (i = 0; i < nq; i++) {
	   nord[i] = i;
	  }
/* ----- banded stiffness matrix into s(nq,nq) ----- */
   fgets(dummy,80,fptr1);
   for (i = 1; i <= nq; i++) {
	 for (jn = 1; jn <= nbw; jn++) {
         fscanf(fptr1, "%lg\n", &c);
	    j = i + jn - 1;
        if (j <= nq) {
           s[nq*(i-1)+j-1] = c;
	       s[nq*(j-1)+i-1] = c;
           }
         }
      }
/* ----- banded mass matrix into gm(nq,nq) ----- */
   fgets(dummy,80,fptr1);
   for (i = 1; i <= nq; i++) {
	 for (jn = 1; jn <= nbw; jn++) {
         fscanf(fptr1, "%lg\n", &c);
	    j = i + jn - 1;
        if (j <= nq) {
           gm[nq*(i-1)+j-1] = c;
	       gm[nq*(j-1)+i-1] = c;
           }
	 }
	   }
     fclose(fptr1);
   /* ----- Cholesky Factorization of Mass Matrix */
     cholesky(gm, nq);
   /* ----- Update of Stiffness Matrix - Standard Form Ax=(lambda)x */
     updatestiff(s, gm, nq);
   /* ----- Tri-diagonalize D() Diagonal, B() Sub-diagonal */
   /* ----- S() has the rotation matrix */ 
    tridiag(s, d, b, nq);   
   /* ----- Find Eigenvalues and Eigenvectors */
     eigentd(d, b, s, nq, iter);
   /* ----- Determine Eigenvectors */
     eigenvec(s, gm, nq);
   /* ----- Ascending Order of Eigenvalues */
     ascendeigen(d, nord, nq);
   /* -----   results   ----- */
   printf("eigenvalues and eigenvectors for data in file %s\n)", file1);
   fprintf(fptr2, "eigenvalues and eigenvectors for data in file %s\n", file1);
     for (i = 0; i < nq; i++) {
        ii = nord[i];
	printf( "eigenvalue number  %d\n ", i+1);
	fprintf(fptr2, "eigenvalue number  %d\n ", i+1);
        c = d[ii];
	    omega = sqrt((float) c);
	    freq = .5 * omega / pi;
	printf("eigenvalue = %11.4e",  c);
	fprintf(fptr2, "eigenvalue = %11.4e",  c);
        printf("   omega = %10.3e", omega);
        fprintf(fptr2, "   omega = %10.3e", omega);
	printf("   freq = %10.3e hz\n", freq);
	fprintf(fptr2, "   freq = %10.3e hz\n", freq);
	printf( "eigenvector \n");
	fprintf(fptr2, "eigenvector \n");
	ifl = 0;
	for (j = 0; j < nq; j++) {
	   printf( "%10.3e ", s[nq*j+ii]);
	   fprintf(fptr2, "%10.3e ", s[nq*j+ii]);
           ifl = ifl + 1;
           if (ifl > 5)
              ifl = 0;
           if (ifl == 0) {
             printf( "\n");
	         fprintf(fptr2, "\n");
	         }
           }
        printf("\n");
	    fprintf(fptr2, "\n");
        }
        fclose(fptr2);
     return(0);
}

cholesky(a, n)
    int n;
    double *a;
    {
      int i, j, k;

     /* ----- L into lower left triangle of a */

     for ( k = 0; k < n; k++ ) {
	a[k*n + k] = sqrt(a[k*n + k]);

	for ( i = k + 1; i < n; i++) {
           a[i*n + k] = a[i*n + k]/a[k*n + k];
	 }
	                              
	for ( j = k + 1; j < n; j++ ) {

	   for ( i = j; i < n; i++) {
	      a[i*n + j] = a[i*n + j] - a[i*n + k] * a[j*n + k];
	   }
	 }
     }
    return(0);
  }

updatestiff(a, b, n)
  int n;
  double *a, *b;
  {
     int i,j,k;
     /* --- forward substitution i  (invL)*a */
     for ( j = 0; j < n; j++) {
	a[j] = a[j] / b[0];
	for ( i = 1; i < n; i++) {
	   for ( k = 0; k < i; k++) {
              a[i*n + j] = a[i*n + j] - b[i*n + k] * a[k*n + j];
           }
           a[i*n + j] = a[i*n + j] / b[i*n + i];
	}
      }

     /* --- forward substitution ii   (invl)*a*(invl)' */
     for ( j = 0; j < n; j++) {
        a[j*n] = a[j*n] / b[0];
        for ( i = 1; i < n; i++ ) {
           for ( k = 0; k < i; k++ ) {
              a[j*n + i] = a[j*n + i] - b[i*n + k] * a[j*n + k];
           }
           a[j*n + i] = a[j*n + i] / b[i*n + i];
        }
      }
    return(0);
  }

 tridiag(a, d, b, n)
    int n;
    double *a, *d, *b;
   {
    int i,j,k,ia, ii, ij, ik;
    double aa,ww,bet, c1;
    for ( i = 0; i < n - 2; i++ ) {
       aa = 0;

       for ( j = i + 1; j < n; j++ ) {
          aa = aa + a[j*n + i] * a[j*n + i];
       }
       aa = sqrt(aa);
       ww = 2 * aa * (aa + fabs(a[(i + 1)*n + i]));
       ww = sqrt(ww);
       ia = 1;
       if (a[(i + 1)*n + i] < 0)
          ia = -1;
       /* ----- diagonal and next to diagonal term */
       d[i] = a[i*n + i];
       b[i] = -ia * aa;
       /* ----- unit vector w() in column i from row i+1 to n */
       for ( j = i + 1; j < n; j++ ) {
	  a[j*n + i] = a[j*n + i] / ww;
       }
       a[(i + 1)*n + i] = a[(i + 1)*n + i] + ia * aa / ww;
       /* ----- w'a in row i from col i+1 to n */
       bet = 0;

       for ( j = i + 1; j < n; j++ ) {
          a[i*n + j] = 0;
          for ( k = i + 1; k < n; k++ ) { 
              a[i*n + j] = a[i*n + j] + a[k*n + i] * a[k*n + j];
          }
          bet = bet + a[i*n + j] * a[j*n + i];
       }
       /* ----- modified a() */
       for ( j = i + 1; j < n; j++ ) {
          for ( k = i + 1; k < n; k++ ) {
	     c1 = - 2 * a[j*n + i] * a[i*n + k];
	     c1 = c1 - 2 * a[i*n + j] * a[k*n + i];
	     c1 = c1 + 4 * bet * a[k*n + i] * a[j*n + i];
	     a[j*n + k] = a[j*n + k] + c1;
          }
       }
    }
    d[n - 2] = a[(n - 2)*n + n - 2];
    b[n - 2] = a[(n - 2)*n + n - 1];
    d[n - 1] = a[(n - 1)*n +n - 1];
    a[(n - 2)*n + n - 2] = 1;
    a[(n - 1)*n +n - 1] = 1;
    a[(n - 2)*n + n - 1] = 0;
    a[(n - 1)*n +n - 2] = 0;

    /* ----- now create the q matrix in a() */
     for (i = 0; i < n - 2; i++) {
       ii = n - i - 3;
       a[ii*n + ii] = 1;
       for (j = 1; j <= i + 2; j++) {
          ij = ii + j;
          c1 = 0;
          for (k = 1; k <= i + 2; k++) {
             ik = ii + k;
	     c1 = c1 + a[ik*n + ij] * a[ik*n + ii];
          }
	  for (k = 1; k <= i + 2; k++) {
             ik = ii + k;
	     a[ik*n + ij] = a[ik*n + ij] - 2 * c1 * a[ik*n + ii];
          }
       }
       for (j = 1; j <= i + 2; j ++) {
          ij = ii + j;
	  a[ii*n + ij] = 0;
	  a[ij*n + ii] = 0;
       }
    }
    return(0);
 }

 eigentd(a, b, q, n, iter)
  int n,iter;
  double *a, *b, *q;

  {
     int m, i, k;
     double dd,bb,p,pp,bot,x,a1,b1,a2,sn,cs;
     iter = 0;
     m = n;
  do {
     iter = iter + 1;
     dd = 0.5 * (a[m - 2] - a[m - 1]);
     bb = b[m - 2] * b[m - 2];
     p = 1;
     if ( dd < 1 )
        p = -1;
     bot = dd + p * sqrt(dd * dd + bb);
     p = a[0] - a[m - 1] + bb / bot;
     x = b[0];
     for ( i = 0; i < m - 1; i++ ) {
	pp = sqrt(p * p + x * x);
        cs = -p / pp;
        sn = x / pp;
        if (i > 0)
           b[i - 1] = cs * b[i - 1] - sn * x;
        a1 = a[i];
        a2 = a[i + 1];
        b1 = b[i];
        a[i] = a1 * cs * cs - 2 * b1 * cs * sn + a2 * sn * sn;
        b[i] = (a1 - a2) * cs * sn + b1 * (cs * cs - sn * sn);
        a[i + 1] = a1 * sn * sn + 2 * b1 * cs * sn + a2 * cs * cs;
        /* ----- update q()  */
        for ( k = 0; k < n; k++) {
           a1 = q[k*n + i];
           a2 = q[k*n + i + 1];
           q[k*n + i] = cs * a1 - sn * a2;
           q[k*n + i + 1] = sn * a1 + cs * a2;
        }
        if (i == m - 2)
            break;
        x = -b[i + 1] * sn;
        b[i + 1] = b[i + 1] * cs;
        p = b[i];
     }
    while (m > 1 && fabs(b[m - 2]) < 0.000001) {
           m = m - 1;
      }
  } while (m > 1);
  m = m - 1;
  return(0);
 }

 eigenvec(a, b, n)
  int n;
  double *a, *b;
  {
    int i,j,k;
     /* --- backsubstitution --- */
     for ( j = 0; j < n; j++ ) {
        a[(n - 1)*n + j] = a[(n - 1)*n + j] / b[(n - 1)*n + n - 1];
        for ( i = n - 2; i >= 0; i--) {
           for ( k = n - 1; k >= i + 1; k-- ) {
              a[i*n + j] = a[i*n + j] - b[k*n + i] * a[k*n + j];
           }
           a[i*n + j] = a[i*n + j] / b[i*n + i];
        }
     }
    return(0);
  }

ascendeigen(d, nord, n)
  int *nord, n;
  double *d;
  {
     int i, j, i1, j1, ii, ij;
     double c1;
     /* --- Initialization of nord[] --- */
     for ( i = 0; i < n; i++ ) {
      nord[i] = i;
      }
     /* ----- ascending order of eigenvalues */
     for ( i = 0; i < n - 1; i++ ) {
        ii = nord[i];
        i1 = ii;
        c1 = d[ii];
        j1 = i;
        for ( j = i + 1; j < n; j++) {
           ij = nord[j];
           if (c1 > d[ij]) {
             c1 = d[ij];
             i1 = ij;
             j1 = j;
           }
        }
        nord[i] = i1;
        nord[j1] = ii;
     }
    return(0);
  }
