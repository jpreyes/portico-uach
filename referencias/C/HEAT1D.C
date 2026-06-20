/****************************************/
/*           program  heat1d           */
/*  t.r.chandrupatla and a.d.belegundu  */
/****************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
main()
{  
   FILE *fptr;
   int ne,nbc,nq,*nb,k,i1,i2,nn,nbw,i,j,n,*bc;
   float *x,*f,*tc,*v,*h,*s,c,ell,ekl,cnst;
   char dummy[81], title[81], file1[81], file2[81];
   puts("\nInput file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr = fopen(file1, "r");
   fgets(dummy,80,fptr);
   fgets(title,80,fptr);
   fgets(dummy,80,fptr);
   fscanf(fptr,"%d %d %d %d %d %d\n", &ne, &nbc, &nq);
   nn = ne + 1;
   nbw = 2;    /* --- bandwidth for nodes in natural order --- */
/* ------------- memory allocation ------------- */   
   x = (float *) calloc(nn, sizeof(float));
   f = (float *) calloc(nn, sizeof(float));
   tc = (float *) calloc(ne, sizeof(float));
   v = (float *) calloc(nbc, sizeof(float));
   h = (float *) calloc(nbc, sizeof(float));
   s = (float *) calloc(nbw*nn, sizeof(float));
   nb = (int *) calloc(nbc, sizeof(int));
   bc = (int *) calloc(nbc, sizeof(int));
/* ---------------------------------------------- */
/* --- thermal conductivity ---*/
   fgets(dummy, 80, fptr);
   for (i = 0; i < ne; i++) {
      fscanf(fptr, "%d %f\n", &k, &c);
      tc[k-1] = c;
   }
 /* ---  coordinates --- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nn; i++) {
      fscanf(fptr, "%d %f\n",&k, &c);
      x[k-1] = c;
      }
/* ---  boundary conditions --- */
     fgets(dummy,80,fptr);
     for (i = 0; i < nbc; i++) {
	 fscanf(fptr, "%d  %s\n", &k, dummy);
	 nb[i] = k;
	 bc[i] = 0;
printf("%d  %s\n",k,dummy);
	  if (strcmp(dummy,"temp")==0 || strcmp(dummy,"TEMP")==0){
	    bc[i] = 1;
	    fscanf(fptr,"%f\n", &c);
	    v[i] = c;
	    }
	  if (strcmp(dummy,"flux")==0 || strcmp(dummy,"FLUX")==0){
	    bc[i] = 2;
	    fscanf(fptr,"%f\n", &c);
	    v[i] = c;
	    }
	  if (strcmp(dummy,"conv")==0 || strcmp(dummy,"CONV")==0){
	    bc[i] = 3;
	    fscanf(fptr,"%f\n", &c);
	    h[i] = c;
	    fscanf(fptr,"%f\n", &c);
	    v[i] = c;
	    }
	}
     /* --- calculate and input nodal heat source vector --- */
     for (i = 0; i < nn; i++) {
	    f[i] = 0;
	    }	    
     fgets(dummy,80,fptr);
     if (nq > 0) {
       for (i = 0; i < nq; i++) {
           fscanf(fptr, "%d %f\n",&k,&c);
           f[k-1] = c;
	       }
       }
     fclose(fptr);
   /* --- stiffness matrix --- */
	 for (i = 0; i < nn; i++) {
	    for (j = 0; j < nbw; j++) {
	       s[nbw*i+j] = 0;
	       }
	    }
     for (i = 0; i < ne; i++) {
	  i1 = i;
	  i2 = i + 1;
	  ell = fabs(x[i2] - x[i1]);
	  ekl = tc[i] / ell;
	  s[nbw*i1] = s[nbw*i1] + ekl;
	  s[nbw*i2] = s[nbw*i2] + ekl;
	  s[nbw*i1+1] = s[nbw*i1+1] - ekl;
	 }
     /* --- account for b.c.'s --- */
     cnst = 0;
     for (i = 0; i < nn; i++) {
        if (s[nbw*i] > cnst)
           cnst = s[nbw*i];
        }
     cnst = cnst * 10000;
     for (i = 0; i < nbc; i++) {
	n = nb[i]-1;
	/* --- temperature --- */
	if (bc[i] == 1) {
	   s[nbw*n] = s[nbw*n] + cnst;
	   f[n] = f[n] + cnst * v[i];
	   }
	/* --- heat flux --- */
	if (bc[i] == 2)
	   f[n] = f[n] - v[i];
	/* --- convection --- */
	if (bc[i] == 3) {
	  s[nbw*n] = s[nbw*n] + h[i];
	  f[n] = f[n] + h[i] * v[i];
	  }
	}
   bansol(s,f,nn,nbw);
     /* --- f contains the solution. 'rhs' is over-written --- */
     fptr = fopen(file2, "w");
     printf("\n%s\n", title);
     fprintf(fptr, "\n%s\n", title); 
     printf("node#  temperature\n");
     fprintf(fptr,"node#  temperature\n");
     for (i = 0; i < nn; i++) {
	     printf("%3d    %8.4f\n", i+1, f[i]);
	     fprintf(fptr,"%3d    %8.4f\n", i+1, f[i]);
	    }
	 fclose(fptr);
	 printf("results are in file %s\n",file2);
     return(0);
}
/* ----- band solver ----- */
bansol(s,f,nq,nbw)
  int nq, nbw;
  float *s, *f;
{
 int n1,k,nk,i,i1,j,j1,kk;
  float c1;
  /* ----- band solver ----- */
  n1 = nq - 1;
  /* --- forward elimination --- */
  for (k = 1; k <= n1; k++) {
     nk = nq - k + 1;
     if (nk > nbw)
	nk = nbw;
     for (i = 2; i <= nk; i++) {
       c1 = s[nbw*(k-1)+i-1] / s[nbw*(k-1)];
       i1 = k + i - 1;
       for (j = i; j <= nk; j++) {
	j1 = j - i + 1;
	s[nbw*(i1-1)+j1-1] = s[nbw*(i1-1)+j1-1] - c1 * s[nbw*(k-1)+j-1];
	}
       f[i1-1] = f[i1-1] - c1 * f[k-1];
       }
     }
  /* --- back-substitution --- */
  f[nq-1] = f[nq-1] / s[nbw*(nq-1)];
  for (kk = 1; kk <= n1;kk++) {
     k = nq - kk;
     c1 = 1 / s[nbw*(k-1)];
     f[k-1] = c1 * f[k-1];
     nk = nq - k + 1;
     if (nk > nbw)
       nk = nbw;
       for (j = 2; j <= nk; j++) {
	 f[k-1] = f[k-1] - c1 * s[nbw*(k-1)+j-1] * f[k + j - 2];
	}
     }
    return(0);
}