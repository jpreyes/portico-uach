/***************************************/
/*          program  cstkm             */
/*    stiffness and mass matrices      */
/*   2-d  constant strain triangle     */
/* t.r.chandrupatla and a.d.belegundu  */
/***************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr;
   int n,i,j,k,m,ii,jj,m1,nmin,nmax,nrt,nct,it,jt,i1,i2;
   int nr,nc;
   char dummy[81], title[81], file1[81], file2[81], file3[81];
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nbw,nch,nmpc,lc,ipl;
   int *noc, *nu, *mat, *mpc;
   float *x, *thick, *pm, *u, *s, *gm, *beta;
   float c, dj, al, pnu, cnst;
   float se[6][6], em[6][6];
/*-------------------------------------------------------*/
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   printf("\n");
   printf("  1) plane stress analysis\n");
   printf("  2) plane strain analysis\n");
   printf("     choose <1 or 2>  ");
   scanf("%d", &lc);
   if (lc < 1 || lc > 2)
      lc = 1;
   fptr = fopen(file1, "r");
   fgets(dummy,80,fptr);
   fgets(title,80,fptr);
   fgets(dummy,80,fptr);
   fscanf(fptr,"%d %d %d %d %d %d\n", &nn, &ne, &nm, &ndim, &nen, &ndn);
   fgets(dummy, 80, fptr);
   fscanf(fptr,"%d %d %d \n", &nd, &nl, &nmpc);
   npr = 4;    /* Material Properties E, Nu, Alpha, Density(Rho) */
/* -----         memory allocation          ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   thick = (float *) calloc(ne, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
/* ----------------------------------------------- *
/* ----- total dof is  nq ----- */
     nq = ndn * nn;
/* ===============  read data  ==================== */
/* ----- coordinates ----- */
     fgets(dummy,80,fptr);
     for (i = 0; i < nn; i++){
        fscanf(fptr, "%d", &n);
        for (j = 0; j < ndim; j++){
		fscanf(fptr, "%f\n", &c);
                x[ndim*(n-1)+j] = c;
            }
         }
/* ----- connectivity, material, thickness, temp-change ----- */
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
       thick[n-1] = c;
       fscanf(fptr,"%f\n",&c);
   }
/* ----- displacement bc  ----- */
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
      /* dummy read */
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
/* ----- bandwidth nbw from connectivity noc() and mpc ----- */
   nbw = 0;
   for (i = 0; i < ne; i++) {
        nmin = noc[nen*i];
        nmax = nmin;
        for (j = 1; j < 3;j++) {
            n =noc[nen*i+j];
            if (nmin > n)
               nmin = n;
            if (nmax < n)
               nmax = n;
            }  
        n= ndn * (nmax - nmin + 1);
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
   gm = (float *) calloc(nq*nbw, sizeof(float));
/* ----- global stiffness matrix -----*/
     for (n = 0; n < ne; n++) {
	printf("forming stiffness matrix of element %d\n", n+1);
	elkm(lc,n,mat,pm,npr,x,noc,thick,se,em);
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
		    if (nc >= 0) {
		       s[nbw*nr+nc] = s[nbw*nr+nc] + se[i][j];
		       gm[nbw*nr+nc] = gm[nbw*nr+nc] + em[i][j];
		       }
		    }
		  }
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
       }
     /* -----  additional springs and lumped masses  ----- */
	 printf("spring supports  < dof# = 0 exits this mode >\n");
	 do {
	printf( "   dof#  ");
	scanf("%d", &k);
	if (k > 0) {
	   printf("   spring const ");
	   scanf("%f", &c);
	   s[nbw*(k-1)] = s[nbw*(k-1)] + c;
	   }
	} while (k > 0);
	 printf("lumped masses  < dof# = 0 exits this mode >\n");
	 do {
	printf( "   dof#  ");
	scanf("%d", &k);
	if (k > 0) {
	   printf("   lumped mass ");
           scanf("%f", &c);
           gm[nbw*(k-1)] = gm[nbw*(k-1)] + c;
           }
        } while (k > 0);
/* --- print banded stiffness and mass matrices in output file --- */
       fptr = fopen(file2, "w");
       fprintf(fptr, "stiffness and mass for data in file %s\n", file1);
       fprintf(fptr,"%d  %d\n", nq, nbw);
       fprintf(fptr, "banded stiffness matrix\n");
       for (i = 0; i < nq; i++) {
          for (j = 0; j < nbw; j++) {
	     fprintf(fptr, "%g ", s[nbw*i+j]);
	     }
	      fprintf(fptr, "\n");
	     }
      fprintf(fptr, "banded mass matrix\n");
       for (i = 0; i < nq; i++) {
          for (j = 0; j < nbw; j++) {
	     fprintf(fptr, "%g ", gm[nbw*i+j]);
	     }
	      fprintf(fptr, "\n");
	     }
      fprintf(fptr, "starting vector for inverse iteration\n");
       for (i = 0; i < nq; i++) {
           fprintf(fptr, "1 ");
           }
	 fprintf(fptr, "\n");
   fclose (fptr);
   printf("global stiffness and mass matrices are in file  %s\n", file2);
   printf("run invitr or jacobi or geneigen program to get \n");
   printf("eigenvalues and eigenvectors\n");
   return(0);
}
/* ----- element stiffness and mass matrices ----- */
elkm(lc,n,mat,pm,npr,x,noc,thick,se,em)
    int lc,n,npr;
    int *mat,*noc;
    float *x,*pm,*thick;
    float se[][6],em[][6];
{
    int i,j,k,m,i1,i2,i3;
    float c,e,c1,c2,c3,dj,pnu,rho,cm;
    float b[3][6],db[3][6],d[3][3];
    float x1,x2,x3,y1,y2,y3,x21,x32,x13,y12,y23,y31;
/* ----- d(), b() and db() matrices ----- */
     /* --- first the d-matrix --- */
     m = mat[n]-1;
     e = pm[npr*m];
     pnu = pm[npr*m+1];
     /* --- d() matrix --- */
     if (lc == 1) {
        /* --- plane stress --- */
        c1 = e / (1 - pnu * pnu);
        c2 = c1 * pnu;
        }
     else {
        /* --- plane strain --- */
        c = e / ((1 + pnu) * (1 - 2 * pnu));
        c1 = c * (1 - pnu);
        c2 = c * pnu;
        }
     c3 = .5 * e / (1 + pnu);
     d[0][0] = c1;
     d[0][1] = c2;
     d[0][2] = 0;
     d[1][0] = c2;
     d[1][1] = c1;
     d[1][2] = 0;
     d[2][0] = 0;
     d[2][1] = 0;
     d[2][2] = c3;
/* --- strain-displacement matrix b() --- */
     i1 = noc[3*n]-1;
     i2 = noc[3*n+1]-1;
     i3 = noc[3*n+2]-1;
     x1 = x[2*i1];
     y1 = x[2*i1+1];
     x2 = x[2*i2];
     y2 = x[2*i2+1];
     x3 = x[2*i3];
     y3 = x[2*i3+1];
     x21 = x2 - x1;
     x32 = x3 - x2;
     x13 = x1 - x3;
     y12 = y1 - y2;
     y23 = y2 - y3;
     y31 = y3 - y1;
     dj = x13 * y23 - x32 * y31; /* dj is determinant of jacobian */
/* --- definition of b() matrix --- */
     b[0][0] = y23 / dj;
     b[1][0] = 0.;
     b[2][0] = x32 / dj;
     b[0][1] = 0.;
     b[1][1] = x32 / dj;
     b[2][1] = y23 / dj;
     b[0][2] = y31 / dj;
     b[1][2] = 0.;
     b[2][2] = x13 / dj;
     b[0][3] = 0;
     b[1][3] = x13 / dj;
     b[2][3] = y31 / dj;
     b[0][4] = y12 / dj;
     b[1][4] = 0.;
     b[2][4] = x21 / dj;
     b[0][5] = 0.;
     b[1][5] = x21 / dj;
     b[2][5] = y12 / dj;
/* --- db matrix db = d*b ---*/
     for (i = 0; i < 3; i++) {
         for (j = 0; j < 6; j++) {
             c = 0.;
             for (k = 0; k < 3; k++) {
		 c = c + d[i][k] * b[k][j];
		 }
	     db[i][j] = c;
	     }
	 }
     /* --- element stiffness --- */
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                c = 0.;
                for (k = 0; k < 3; k++) {
                    c = c + .5 * fabs(dj) * b[k][i] * db[k][j] * thick[n];
                    }
		se[i][j] = c;
		}
	    }
     /* ----- element mass  em[][] ----- */
     rho = pm[npr*m+3];
     cm = rho * thick[n] * .5 * fabs(dj) / 12;
     for (i = 0; i < 6; i++) {
	    for (j = 0; j < 6; j++) {
	       em[i][j] = 0;
	       }
	    }
     /* --- non-zero elements of mass matrix are defined --- */
     em[0][0] = 2 * cm;
	 em[0][2] = cm;
	 em[0][4] = cm;
     em[1][1] = 2 * cm;
	 em[1][3] = cm;
	 em[1][5] = cm;
     em[2][0] = cm;
	 em[2][2] = 2 * cm;
	 em[2][4] = cm;
     em[3][1] = cm;
	 em[3][3] = 2 * cm;
	 em[3][5] = cm;
     em[4][0] = cm;
	 em[4][2] = cm;
	 em[4][4] = 2 * cm;
     em[5][1] = cm;
	 em[5][3] = cm;
	 em[5][5] = 2 * cm;
    return(0);
}
