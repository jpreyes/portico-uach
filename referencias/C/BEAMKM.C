/*==============================*/
/*          beamkm .c           */
/*   shaft vibration analysis   */
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
   float *x, *smi, *pm, *u, *s, *beta, *gm, *area;
   float c, el, eil, cnst, se[4][4], em[4][4];
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
   fscanf(fptr,"%d %d %d n", &nd, &nl, &nmpc);
   npr = 2;   /* Material Properties E, Density(Rho) */
   /* ----- total dof is nq ----- */
   nq = ndn*nn;
/* ----- memory allocation ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   smi = (float *) calloc(ne, sizeof(float));
   area = (float *) calloc(ne, sizeof(float));
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
         smi[n-1] = c;
     fscanf(fptr,"%f\n",&c);
         area[n-1] = c;
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
      /* this is dummy read */
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
/* ----- allocate memory for stiffness and global mass ----- */
   s = (float *) calloc(nq*nbw, sizeof(float));
   gm = (float *) calloc(nq*nbw, sizeof(float));
/* ----- assemble the stiffness matrix ----- */
   for (n = 0; n < ne; n++) {
      printf("forming stiffness and mass matrices of element... %d\n", n+1);
      elkm(n,noc,mat,pm,x,smi,area,se,em);
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
       fprintf(fptr, "num. of dof    bandwidth");
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
elkm(n,noc,mat,pm,x,smi,area,se,em)
int n,*noc,*mat;
float *pm,*x,*smi,*area,se[][4],em[][4];
{
  int i1,i2,m;
  float el,eil,rho,c1;
     i1 = noc[2*n];
     i2 = noc[2*n+1];
     m = mat[n];
     el = fabs(x[i1-1] - x[i2-1]);
     eil = pm[2*(m-1)]* smi[n] / (el*el*el);
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
      /* --- element mass --- */
      rho = pm[2*(m-1)+1];
      c1 = rho * area[n] * el / 420;
      em[0][0] = 156 * c1;
      em[0][1] = 22 * el * c1;
      em[0][2] = 54 * c1;
      em[0][3] = -13 * el * c1;
	 em[1][0] = em[0][1];
	 em[1][1] = 4 * el * el * c1;
	 em[1][2] = 13 * el * c1;
	 em[1][3] = -3 * el * el * c1;
      em[2][0] = em[0][2];
      em[2][1] = em[1][2];
      em[2][2] = 156 * c1;
      em[2][3] = -22 * el * c1;
         em[3][0] = em[0][3];
         em[3][1] = em[1][3];
         em[3][2] = em[2][3];
         em[3][3] = 4 * el * el * c1;
  return(0);
}