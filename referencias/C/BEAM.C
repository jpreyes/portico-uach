/*==============================*/
/*            beam.c            */
/*   chandrupatla & belegundu   */
/*          (C) 2001            */
/* =============================*/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr;
   int n,nq,i,j,k,m,i1,i2,ii,nrt,it,nr,jj,nct,jt,nc;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc;
   int *noc, *nu, *mat, *mpc;
   float *x, *smi, *pm, *u, *s, *f, *beta, se[4][4];
   float c, el, eil, cnst, reaction;
   printf("\n");
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
   npr = 1;  /* Material property E */
   /* ----- total dof is nq ----- */
   nq = ndn*nn;
/* ----- memory allocation ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   smi = (float *) calloc(ne, sizeof(float));
   f = (float *) calloc(nn*ndn, sizeof(float));
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
	 fscanf(fptr,"%f\n",&c);
         smi[n-1] = c;
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
      f[k-1] = c;
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
      n = ndn * (abs(noc[nen*i] - noc[nen*i+1])+1);
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
   s = (float *) calloc(nq*nbw, sizeof(float));
/* ----- assemble the stiffness matrix ----- */
   for (n = 0; n < ne; n++) {
      printf("forming stiffness matrix of element... %d\n", n+1);
	 i1 = noc[nen*n];
         i2 = noc[nen*n+1];
         m = mat[n];
         el = fabs(x[i1-1] - x[i2-1]);
	 eil = pm[m-1]* smi[n] / (el*el*el);
         se[0][0] = 12 * eil;
         se[0][1] = eil * 6 * el;
         se[0][2] = -12 * eil;
         se[0][3] = eil * 6 * el;
            se[1][0] = se[0][1];
            se[1][1] = eil * 4 * el * el;
            se[1][2] = -eil * 6 * el;
            se[1][3] = eil * 2 * el * el;
         se[2][0] = se[0][2];
         se[2][1] = se[1][2];
         se[2][2] = eil * 12;
         se[2][3] = -eil * 6 * el;
            se[3][0] = se[0][3];
            se[3][1] = se[1][3];
	    se[3][2] = se[2][3];
            se[3][3] = eil * 4 * el * el;
	printf (".... placing in global locations\n");
	for (ii = 0; ii < nen; ii++) {
	   nrt = ndn * (noc[nen*n + ii] - 1);
	   for (it = 0; it < ndn; it++) {
	      nr = nrt + it;
	      i = ndn * ii + it;
	      for (jj = 0; jj < nen; jj++) {
		 nct = ndn * (noc[nen*n+jj] - 1);
		 for (jt = 0; jt < ndn; jt++) {
		    j = ndn * jj + jt;
		    nc = nct + jt - nr;
		    if (nc >= 0)
		       s[nbw*nr+nc] = s[nbw*nr+nc] + se[i][j];
		    }
		  }
	     }
	  }
   }
/* ----- decide penalty parameter cnst ----- */
   cnst = 0;
   for (i = 0;i < nq; i++){
       if (cnst < s[i*nbw])
          cnst = s[i*nbw];
       }
   cnst = 10000 * cnst;
/* ----- modify for displacement boundary conditions ----- */
   for (i = 0; i < nd; i++) {
      k = nu[i];
      s[(k-1)*nbw] = s[(k-1)*nbw] + cnst;
      f[k-1] = f[k-1] + cnst * u[i];
   }
/* ----- modify for multipoint constraints ----- */
   for (i = 0; i < nmpc; i++){
       i1 = mpc[2*i];
       i2 = mpc[2*i+1];
       s[(i1-1)*nbw] = s[(i1-1)*nbw] + cnst*beta[3*i]*beta[3*i];
       s[(i2-1)*nbw] = s[(i2-1)*nbw] + cnst*beta[3*i+1]*beta[3*i+1];
       n=i1;
       if (n > i2)
	  n = i2;
       m = abs(i2-i1);
       s[(n-1)*nbw+m] = s[(n-1)*nbw+m]+cnst*beta[3*i]*beta[3*i+1];
       f[i1-1] = f[i1-1] + cnst*beta[3*i]*beta[3*i+2];
       f[i2-1] = f[i2-1] + cnst*beta[3*i+1]*beta[3*i+2];
       }
/* ----- solution of equations using band solver ----- */
   bansol(s,f,nq,nbw);
/* ----- printing displacements ----- */
   printf("node#   displ.       rotation\n");
   fprintf(fptr, "node#   displ.       rotation\n");
   for (i = 0; i < nn; i++) {
      printf(" %3d  %11.4e  %11.4e\n",i+1,f[2*i],f[2*i+1]);
      fprintf(fptr, " %3d  %11.4e  %11.4e\n",i+1,f[2*i],f[2*i+1]);
   }
/* ----- reaction calculation ----- */
   printf("node#   reaction\n");
   fprintf(fptr, "node#   reaction\n");
   for (i = 0; i < nd; i++) {
      k = nu[i];
      reaction = cnst * (u[i] - f[k-1]);
      printf(" %3d   %11.4e\n", k, reaction);
      fprintf(fptr, " %3d   %11.4e\n", k, reaction);
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