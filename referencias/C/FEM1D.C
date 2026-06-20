/*==============================*/
/*            fem1d.c           */
/*   CHANDRUPATLA & BELEGUNDU   */
/*          (C) 2001            */
/* =============================*/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr;
   int n,i,j,k,m,i1,i2,m1;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc;
   int *noc, *nu, *mat, *mpc;
   float *x, *area, *pm, *u, *tempr, *stiff, *force, *beta;
   float c, x21, eal, tld, cnst, eps, stress, reaction;
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr = fopen(file1, "r");
   fgets(dummy,80,fptr);
   fgets(title,80,fptr);
   fgets(dummy,80,fptr);
   fscanf(fptr,"%d %d %d %d %d %d\n", &nn, &ne, &nm, &ndim, &nen, &ndn);
   fgets(dummy, 80, fptr);
   fscanf(fptr,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
   npr = 2;    /* Material Properties E, Alpha */
/* ----- memory allocation ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   area = (float *) calloc(ne, sizeof(float));
   force = (float *) calloc(nn*ndn, sizeof(float));
   tempr = (float *) calloc(ne, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
/* ----- coordinates ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nn; i++) {
      fscanf(fptr, "%d %f\n",&n, &c);
      x[n-1] = c;
      }
/* ----- connectivity etc ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < ne; i++) {
      fscanf(fptr,"%d", &n);
      for (j = 0; j < nen; j++) {
         fscanf(fptr,"%d", &k);
         noc[(n-1)*nen+j]=k;
      }
         fscanf(fptr,"%d", &k);
         mat[n-1] = k;
         fscanf(fptr,"%f",&c);
         area[n-1] = c;
         fscanf(fptr,"%f\n",&c);
         tempr[n-1] = c;
   }
/* ----- boundary conditions ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nd; i++) {
      fscanf(fptr, "%d %f\n", &k, &c);
      nu[i] = k;
      u[i] = c;
   }
/* ----- component loads ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nl; i++) {
      fscanf(fptr, "%d %f\n", &k, &c);
      force[k-1] = c;
   }
/* ----- material properties ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nm; i++){
      fscanf(fptr, "%d", &k);
      for (j = 0; j < npr; j++) {
         fscanf(fptr, "%f\n", &c);
	 pm[(k-1)*npr+j] = c;
      }
   }
/* ----- multipoint constraints ----- */
   if (nmpc > 0) 
      { fgets(dummy,80,fptr);
        for(j=0;j<nmpc;j++){
           fscanf(fptr,"%f",&c);
	   beta[3*j]=c;
	   fscanf(fptr,"%d",&k);
	   mpc[2*j]=k;
           fscanf(fptr,"%f",&c);
	   beta[3*j+1]=c;
           fscanf(fptr,"%d",&k);
           mpc[2*j+1]=k;
           fscanf(fptr,"%f",&c);
	   beta[3*j+2]=c;
	   }
	}
   fclose (fptr);
/* ----- bandwidth evaluation ----- */
   nbw = 0;
   for (i = 0; i < ne; i++) {
      n = abs(noc[nen*i] - noc[nen*i+1])+1;
      if (nbw < n)
         nbw = n;
   }
   for (i = 0; i < nmpc; i++) {
      n = abs(mpc[2*i] - mpc[2*i+1])+1;
      if (nbw < n)
         nbw = n;
   }
   fptr = fopen(file2, "w");
   printf("\n%s\n", title);
   fprintf(fptr, "\n%s\n", title);
   printf("bandwidth = %d\n",nbw);
   fprintf(fptr, "bandwidth = %d\n",nbw);
/* ----- allocate memory for stiffness ----- */
   stiff = (float *) calloc(nn*nbw, sizeof(float));
/* ----- assemble the stiffness matrix ----- */
   for (i = 0; i < ne; i++) {
      i1 = noc[nen*i];
      i2 = noc[nen*i+1];
      m1 = mat[i];
      x21 = x[i2-1] - x[i1-1];
      eal = pm[npr*(m1-1)] * area[i] / fabs(x21);
      c = 0;
      if (npr > 1)
	 c = pm[npr*(m1-1)+1];
      tld = eal*c*tempr[i]*x21;
      stiff[(i1-1)*nbw] = stiff[(i1-1)*nbw] + eal;
      stiff[(i2-1)*nbw] = stiff[(i2-1)*nbw] + eal;
      n = i1;
      if (i2 < i1)
         n = i2;
         m = abs(i2-i1);
      stiff[(n-1)*nbw+m] = stiff[(n-1)*nbw+m] - eal;
   /* --- temperature forces --- */
      force[i1-1] = force[i1-1] - tld;
      force[i2-1] = force[i2-1] + tld;
   }
/* ----- decide penalty parameter cnst ----- */
   cnst = 0;
   for (i = 0;i < nn; i++){
       if (cnst < stiff[i*nbw])
          cnst = stiff[i*nbw];
       }
   cnst = 10000 * cnst;
/* ----- modify for displacement boundary conditions ----- */
   for (i = 0; i < nd; i++) {
      k = nu[i];
      stiff[(k-1)*nbw] = stiff[(k-1)*nbw] + cnst;
      force[k-1] = force[k-1] + cnst * u[i];
   }
/* ----- modify for multipoint constraints ----- */
   for (i = 0; i < nmpc; i++){
       i1 = mpc[2*i];
       i2 = mpc[2*i+1];
       stiff[(i1-1)*nbw] = stiff[(i1-1)*nbw] + cnst*beta[3*i]*beta[3*i];
       stiff[(i2-1)*nbw] = stiff[(i2-1)*nbw] + cnst*beta[3*i+1]*beta[3*i+1];
       n=i1;
       if (n > i2)
	  n = i2;
       m = abs(i2-i1);
       stiff[(n-1)*nbw+m] = stiff[(n-1)*nbw+m]+cnst*beta[3*i]*beta[3*i+1];
       force[i1-1] = force[i1-1] + cnst*beta[3*i]*beta[3*i+2];
       force[i2-1] = force[i2-1] + cnst*beta[3*i+1]*beta[3*i+2];
       }
/* ----- solution of equations using band solver ----- */
   bansol(stiff,force,nn,nbw);
/* ----- printing displacements ----- */
   printf("node#  displacement\n");
   fprintf(fptr, "node#  displacement\n");
   for (i = 0; i < nn; i++) {
      printf(" %d     %e\n",i+1,force[i]);
      fprintf(fptr, " %d     %e\n",i+1,force[i]);
   }
/* ----- stress calculation ----- */
   printf("elem#  stress\n");
   fprintf(fptr, "elem#  stress\n");
   for (i = 0; i < ne; i++) {
      i1 = noc[nen*i];
      i2 = noc[nen*i+1];
      m1 = mat[i];
      x21 = x[i2-1] - x[i1-1];
      eps = (force[i2-1] - force[i1-1]) / x21;
      stress = pm[npr*(m1-1)] * (eps - pm[npr*(m1-1)+1]*tempr[i]);
      printf(" %d     %e\n", i+1, stress);
      fprintf(fptr, " %d     %e\n", i+1, stress);
   }
/* ----- reaction calculation ----- */
   printf("node#  reaction\n");
   fprintf(fptr, "node#  reaction\n");
   for (i = 0; i < nd; i++) {
      k = nu[i];
      reaction = cnst * (u[i] - force[k-1]);
      printf(" %d     %e\n", k, reaction);
      fprintf(fptr, " %d     %e\n", k, reaction);
   }
   fclose (fptr);
   printf("\n output is in file %s \n",file2);
   return 0;
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