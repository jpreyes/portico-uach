/***************************************/
/*         program  truss2d            */
/* t.r.chandrupatla and a.d.belegundu  */
/***************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   int n,i,j,k,m,i1,i2,ii,jj,m1,nmin,nmax,nrt,nct,it,jt;
   int nr,nc,j1,j2,k1,k2;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc;
   int *noc, *nu, *mat, *mpc;
   float *x, *area, *pm, *u, *tempr, *s, *f, *beta;
   float c, al, e, tld, cnst, reaction, stress, x21, y21;
   float se[4][4], tl[4], sn, cs, ee0, el, eal;
/*-------------------------------------------------------*/
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr1 = fopen(file1, "r");
   fgets(dummy,80,fptr1);
   fgets(title,80,fptr1);
   fgets(dummy,80,fptr1);
   fscanf(fptr1,"%d %d %d %d %d %d\n", &nn, &ne, &nm, &ndim, &nen, &ndn);
   fgets(dummy, 80, fptr1);
   fscanf(fptr1,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
   npr = 2;  /* Material Properties  E, Alpha
/* ----- memory allocation ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   area = (float *) calloc(ne, sizeof(float));
   f = (float *) calloc(nn*ndn, sizeof(float));
   tempr = (float *) calloc(ne, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
/* ----- total dof is  nq ----- */
     nq = ndn * nn;
/* ===============  read data  ==================== */
/* ----- coordinates ----- */
     fgets(dummy,80,fptr1);
     for (i = 0; i < nn; i++){
        fscanf(fptr1, "%d", &n);
        for (j = 0; j < ndim; j++){
		fscanf(fptr1, "%f\n", &c);
                x[ndim*(n-1)+j] = c;
            }
         }
/* ----- connectivity, material, thickness, temp-change ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < ne; i++) {
       fscanf(fptr1,"%d", &n);
       for (j = 0; j < nen; j++) {
           fscanf(fptr1,"%d", &k);
           noc[(n-1)*nen+j]=k;
       }
       fscanf(fptr1,"%d", &k);
       mat[n-1] = k;
       fscanf(fptr1,"%f\n",&c);
       area[n-1] = c;
       fscanf(fptr1,"%f\n",&c);
       tempr[n-1] = c;
   }
/* ----- displacement bc  ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nd; i++) {
      fscanf(fptr1, "%d %f\n", &k, &c);
      nu[i] = k;
      u[i] = c;
   }
/* ----- component loads ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nl; i++) {
      fscanf(fptr1, "%d %f\n", &k, &c);
      f[k-1] = c;
   }
/* ----- material properties ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nm; i++){
      fscanf(fptr1, "%d", &k);
      for (j = 0; j < npr; j++) {
         fscanf(fptr1, "%f\n", &c);
	 pm[(k-1)*npr+j] = c;
      }
   }
/* ----- multipoint constraints ----- */
   if (nmpc > 0) 
      { fgets(dummy,80,fptr1);
	for(j=0;j<nmpc;j++){
	   fscanf(fptr1,"%f",&c);
	   beta[3*j]=c;
	   fscanf(fptr1,"%d",&k);
	   mpc[2*j]=k;
           fscanf(fptr1,"%f",&c);
	   beta[3*j+1]=c;
           fscanf(fptr1,"%d",&k);
           mpc[2*j+1]=k;
           fscanf(fptr1,"%f",&c);
	   beta[3*j+2]=c;
	   }
	}
   fclose (fptr1);
/* ----- bandwidth nbw from connectivity noc() and mpc ----- */
   nbw = 0;
   for (i = 0; i < ne; i++) {
	n = ndn * (abs(noc[nen*i] - noc[nen*i+1]) + 1);
        if (nbw < n)
           nbw = n;
   }    
   for (i = 0; i < nmpc; i++) {
        n = abs(mpc[2*i] - mpc[2*i+1]) + 1;
        if (nbw < n)
           nbw = n;
        }
     printf ("the bandwidth is %d\n", nbw);
/* ----- allocate memory for stiffness ----- */
   s = (float *) calloc(nq*nbw, sizeof(float));
/* ----- global stiffness matrix -----*/
   for (n = 0; n < ne; n++) {
	 printf("forming stiffness matrix of element %d\n", n+1);
     /* --- element stiffness --- */
     
        i1 = noc[nen*n] - 1;
        i2 = noc[nen*n+1] - 1;
        m = mat[n] - 1;
        x21 = x[ndim*i2] - x[ndim*i1];
        y21 = x[ndim*i2+1] - x[ndim*i1+1];
        el = sqrt((double) (x21 * x21 + y21 * y21));
        eal = pm[npr*m] * area[n] / el;
        cs = x21 / el;
        sn = y21 / el;
    /* ----------- element stiffness matrix se() ----------- */
        se[0][0] = cs * cs * eal;
        se[0][1] = cs * sn * eal;
	se[0][2] = -cs * cs * eal;
	se[0][3] = -cs * sn * eal;
	se[1][0] = se[0][1];
	se[1][1] = sn * sn * eal;
	se[1][2] = -cs * sn * eal;
	se[1][3] = -sn * sn * eal;
	se[2][0] = se[0][2];
	se[2][1] = se[1][2];
	se[2][2] = cs * cs * eal;
	se[2][3] = cs * sn * eal;
	se[3][0] = se[0][3];
	se[3][1] = se[1][3];
	se[3][2] = se[2][3];
	se[3][3] = sn * sn * eal;
     /* --- temperature load vector --- */
        ee0 = pm[npr*m+1] * tempr[n] * pm[npr*m] * area[n];
        tl[0] = -ee0 * cs;
        tl[1] = -ee0 * sn;
        tl[2] = ee0 * cs;
        tl[3] = ee0 * sn;
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
	      f[nr] = f[nr] + tl[i];
	      }
	   }
     }
/* ----- decide penalty parameter cnst ----- */
     cnst = 0.;
     for (i = 0; i < nq; i++) {
	 if (cnst < s[i*nbw])
	    cnst = s[i*nbw];
	 }
     cnst = cnst * 10000.;
/* ----- modify for displacement boundary conditions ----- */
   for (i = 0; i < nd; i++) {
      k = nu[i];
      s[(k-1)*nbw] = s[(k-1)*nbw] + cnst;
      f[k-1] = f[k-1] + cnst * u[i];
   }
/* ----- modify for multipoint constraints ----- */
   for (i = 0; i < nmpc; i++){
       i1 = mpc[2*i]-1;
       i2 = mpc[2*i+1]-1;
       s[i1*nbw] = s[i1*nbw] + cnst*beta[3*i]*beta[3*i];
       s[i2*nbw] = s[i2*nbw] + cnst*beta[3*i+1]*beta[3*i+1];
       n=i1;
       if (n > i2)
       n = i2;
       m = abs(i2-i1);
       s[n*nbw+m] = s[n*nbw+m]+cnst*beta[3*i]*beta[3*i+1];
       f[i1] = f[i1] + cnst*beta[3*i]*beta[3*i+2];
       f[i2] = f[i2] + cnst*beta[3*i+1]*beta[3*i+2];
       }
/* ----- solution of equations using band solver ----- */
       bansol(s,f,nq,nbw);
/* ----- printing displacements ----- */
   fptr1 = fopen(file2, "w");
   printf("\n%s\n", title);
   fprintf(fptr1, "\n%s\n", title);
   fprintf(fptr1, "bandwidth = %d\n",nbw);
   fprintf(fptr1, "node#     x-displ      y-displ\n");
   printf ("node#     x-displ      y-displ\n");
   for (i = 0; i < nn; i++) {
	 printf(" %4d  %11.4e  %11.4e\n",i+1,f[2*i],f[2*i+1]);
	 fprintf(fptr1," %4d  %11.4e  %11.4e\n",i+1,f[2*i],f[2*i+1]);
	 }
/* ----- reaction calculation ----- */
   printf("node#    reaction\n");
   fprintf(fptr1, "node#     reaction\n");
   for (i = 0; i < nd; i++) {
      k = nu[i];
      reaction = cnst * (u[i] - f[k-1]);
      printf(" %4d  %11.4e\n", k, reaction);
      fprintf(fptr1, " %4d  %11.4e\n", k, reaction);
   }
/* -----  stress calculations ----- */
     printf("elem#    stress\n");
     fprintf(fptr1, "elem#    stress\n");
     for (i = 0; i < ne; i++) {
	i1 = noc[nen*i] - 1;
	i2 = noc[nen*i+1] - 1;
	m = mat[i] - 1;
	x21 = x[ndim*i2] - x[ndim*i1];
	y21 = x[ndim*i2+1] - x[ndim*i1+1];
	el = sqrt((double) (x21 * x21 + y21 * y21));
	cs = x21 / el;
	sn = y21 / el;
	j1 = 2 * i1;
	j2 = j1 + 1;
	k1 = 2 * i2;
	k2 = k1 + 1;
        c = (f[k1] - f[j1]) * cs + (f[k2] - f[j2]) * sn;
        stress = pm[npr*m] * (c / el - pm[npr*m+1] * tempr[i]);
	printf(" %4d  %11.4e\n", i+1, stress);
        fprintf(fptr1, " %4d  %11.4e\n", i, stress);
        }
     fclose(fptr1);
     printf( "complete results are in file %s\n", file2);
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