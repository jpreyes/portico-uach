/*****       program invitr        *****/
/*      inverse iteration method       */
/*  for eigenvalues and eigenvectors   */
/*        searching in subspace        */
/*         for banded matrices         */
/* t.r.chandrupatla and a.d.belegundu  */
/***************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   float *s,*gm,*ev1,*ev2,*evt,*evs,*evc,*evl,*st,c;
   float cv,el1,el2,c1,c2,omega,freq;
   int nev=0,nq,nbw,i,j,ich,n,nv,iter,itmax=50,k,ka,kz,l;
   int k1,l1,ia,iz,i1,ifl;
   float pi = 3.14159, tol = 0.000001, sh = 0;
   char dummy[81], title[81], file1[81], file2[81];
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
   printf("<enter 0 for default tolerance of 1e-6>\n");
   scanf("%f", &c);
   if (c != 0)
      tol = c;
  printf("number of eigenvalues desired\n");
  scanf("%d", &nev);
/* -----         memory allocation         ----- */
   s = (float *) calloc(nbw*nq, sizeof(float));
   gm = (float *) calloc(nbw*nq, sizeof(float));
   ev1 = (float *) calloc(nq, sizeof(float));
   ev2 = (float *) calloc(nq, sizeof(float));
   evt = (float *) calloc(nq, sizeof(float));
   evs = (float *) calloc(nq, sizeof(float));
   st = (float *) calloc(nq, sizeof(float));
   evc = (float *) calloc(nev*nq, sizeof(float));
   evl = (float *) calloc(nev, sizeof(float));   
/* --------------------------------------------- */
   /* ----- read in stiffness matrix ----- */
   fgets(dummy,80,fptr1);
     for (i = 0; i < nq; i++) {
        for (j = 0; j < nbw; j++) {
           fscanf(fptr1, "%f\n", &c);
           s[nbw*i+j] = c;
           }
        }
   /* ----- read in mass matrix ----- */
   fgets(dummy,80,fptr1);
     for (i = 0; i < nq; i++) {
        for (j = 0; j < nbw; j++) {
           fscanf(fptr1, "%f\n", &c);
           gm[nbw*i+j] = c;
           }
        }
   /* ----- read in starting vector ----- */
   fgets(dummy,80,fptr1);
     for (i = 0; i < nq; i++) {
           fscanf(fptr1, "%f\n", &c);
           st[i] = c;
        }
     fclose(fptr1);
           printf( "default shift value  %f\n", sh);
	   printf( "enter desired shift value \n");
           scanf("%f", &c);
           if (c != 0) {
              sh = c;
             for(i = 0; i < nq; i++) {
                for(j = 0; j < nbw; j++) {
		   s[nbw*i+j] = s[nbw*i+j] - sh * gm[nbw*i+j];
		 }
               }
             }
   printf("eigenvalues and eigenvectors for data in file %s\n", file1);
   fprintf(fptr2, "eigenvalues and eigenvectors for data in file %s\n", file1);
   bansol2(s,nq,nbw);        /*<----stiffness to upper triangle */
     for (nv = 1; nv <= nev; nv++) {
        /* --- starting value for eigenvector --- */
        for (i = 0; i < nq; i++) {
           ev1[i] = st[i];
           }
        el2 = 0;
        iter = 0;
	do  {
	   el1 = el2;
	   iter = iter + 1;
	   if (iter > itmax) {
	      printf( "no convergence for %d  iterations\n", iter);
	      exit(1);
	     }
	   if (nv > 1) {
	      /*----  starting vector orthogonal to ----*/
	      /*----       evaluated vectors ----*/
	      for (i = 1; i < nv; i++) {
		 cv = 0;
		 for (k = 1; k <= nq; k++) {
		    ka = k - nbw + 1;
		    kz = k + nbw - 1;
		    if (ka < 1)
		       ka = 1;
		    if (kz > nq)
		       kz = nq;
		    for (l = ka; l <= kz; l++) {
		       if (l < k) {
			  k1 = l - 1;
			  l1 = k - l;
			  } else {
			  k1 = k - 1;
			  l1 = l - k;
			  }
		       cv = cv + evs[k-1] * gm[nbw*k1+l1] * evc[nev*(l-1)+i-1];
		       }
		    }
		 for (k = 0; k < nq; k++) {
		    ev1[k] = ev1[k] - cv * evc[nev*k+i-1];
		    }
		 }
	      }
	   for (i = 1; i <= nq; i++) {
	      ia = i - nbw + 1;
		  iz = i + nbw - 1;
	      evt[i-1] = 0;
	      if (ia < 1)
		 ia = 1;
	      if (iz > nq)
		 iz = nq;
	      for (k = ia; k <= iz; k++) {
		 if (k < i) {
		    i1 = k - 1;
			k1 = i - k;
			}
		 else {
		    i1 = i - 1;
			k1 = k - i;
		    }
		 evt[i-1] = evt[i-1] + gm[nbw*i1+k1] * ev1[k-1];
		}
	      ev2[i-1] = evt[i-1];
	     }
      /* --- reduce right side and solve --- */
	   rhsolve(s,nq,nbw,ev2);
	   c1 = 0;
	       c2 = 0;
           for (i = 0; i < nq; i++) {
              c1 = c1 + ev2[i] * evt[i];
              }
           for (i = 1; i <= nq; i++) {
              ia = i - nbw + 1;
	          iz = i + nbw - 1;
    	     evt[i-1] = 0;
              if (ia < 1)
                 ia = 1;
              if (iz > nq)
                 iz = nq;
              for (k = ia; k <= iz; k++) {
                 if (k < i) {
                    i1 = k - 1;
	                k1 = i - k;
	                }
                 else {
                    i1 = i - 1;
	                k1 = k - i;
                   }
                 evt[i-1] = evt[i-1] + gm[nbw*i1+k1] * ev2[k-1];
               }
             }
           for (i = 0; i < nq; i++) {
              c2 = c2 + ev2[i] * evt[i];
             }
           el2 = c1 / c2;
           c2 = sqrt((double) c2);
           for (i = 0; i < nq; i++) {
              ev1[i] = ev2[i] / c2;
              evs[i] = ev1[i];
              }
	 } while (fabs(el2 - el1) / fabs(el2) > tol);
	 
        for (i = 0; i < nq; i++) {
	   evc[nev*i + nv - 1] = ev1[i];
           }
        printf("eigenvalue number %d", nv);
        fprintf(fptr2, "eigenvalue number %d", nv);
        printf("     iteration number %d\n", iter);
        fprintf(fptr2, "     iteration number %d\n", iter);
        el2 = el2 + sh;
	evl[nv-1] = el2;
	printf("eigenvalue = %10.4e ", el2);
	fprintf(fptr2, "eigenvalue = %10.4e ", el2);
        omega = sqrt((double) el2);
	    freq = .5 * omega / pi;
        printf("   omega = %10.3e ", omega);
        fprintf(fptr2, "   omega = %10.3e ", omega);
	printf("   freq = %10.3e hz\n", freq);
	fprintf(fptr2, "   freq = %10.3e hz\n", freq);
	printf("eigenvector \n");
	fprintf(fptr2, "eigenvector \n");
        ifl = 0;
	for (i = 0; i < nq; i++) {
	   printf("%10.3e ", evc[nev*i + nv - 1]);
	   fprintf(fptr2, "%10.3e ", evc[nev*i + nv - 1]);
	   ifl = ifl + 1;
	   if (ifl > 5)
	      ifl = 0;
           if (ifl == 0) {
              printf("\n");
	          fprintf(fptr2, "\n");
	         }
          }
          printf("\n");
          fprintf(fptr2, "\n");          
       }
     fclose(fptr2);
     return(0);
}     
bansol2(s,nq,nbw)
   int nq,nbw;
   float *s;
{
   int k,i,nk,i1,j,j1;
   float c1;
/* --- gauss elimination ldu approach (for symmetric banded matrices) --- */
     /* ----- reduction to upper triangular form ----- */
     for (k = 1; k < nq; k++) {
        nk = nq - k + 1;
        if (nk > nbw)
           nk = nbw;
        for (i = 2; i <= nk; i++) {
	   c1 = s[nbw*(k-1)+i-1] / s[nbw*(k-1)];
           i1 = k + i - 2;
           for (j = i; j <= nk; j++) {
              j1 = j - i;
              s[nbw*i1+j1] = s[nbw*i1+j1] - c1 * s[nbw*(k-1)+j-1];
              }
           }
        }
     return(0);
}
rhsolve(s,nq,nbw,ev2)
 int nq,nbw;
 float *s,*ev2;
  {
   int k,nk,i,ni,i1;
   float c1;
     /* ----- reduction of the right hand side ----- */
     for (k = 1; k < nq; k++) {
        nk = nq - k + 1;
        if (nk > nbw)
          nk = nbw;
        for (i = 2; i <= nk; i++) {
	       i1 = k + i - 1;
           c1 = 1 / s[nbw*(k-1)];
           ev2[i1-1] = ev2[i1-1] - c1 * s[nbw*(k-1)+i-1] * ev2[k-1];
          }
       }
     /* ----- back substitution ----- */
     ev2[nq-1] = ev2[nq-1] / s[nbw*(nq-1)];
     for (i = nq - 1; i >= 1; i--) {
	    c1 = 1 / s[nbw*(i - 1)];
        ni = nq - i + 1;
        if (ni > nbw)
           ni = nbw;
        ev2[i-1] = c1 * ev2[i-1];
        for (k = 2; k <= ni; k++) {
           ev2[i-1] = ev2[i-1] - c1 * s[nbw*(i-1)+k-1] * ev2[i + k - 2];
           }
       }
     return(0);
     }

