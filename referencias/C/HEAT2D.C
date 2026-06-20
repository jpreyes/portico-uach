/***************************************/
/*          program  heat2d            */
/*   heat 2-d  with 3-noded triangles  */
/* t.r.chandrupatla and a.d.belegundu  */
/***************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   int i,j,k,m,n,ii,jj,nconv,nhf,nbw,n1,n2,nmax,nmin,i1,i2,i3,ii1,ii2;
   char dummy[121], title[81], file1[81], file2[81], file3[81];
   int nn,ne,nm,ndim,nen,ndn,nd,nl,npr,nmpc,*noc,*mat,*nu,ipl;
   float *x,*pm,*f,*u,*ehs,*s,v,h,tinf,bt[2][3],detj;
   float elen,x32,x13,x21,y23,y31,y12,area,sum,cnst,eq,c,qx,qy;
/*-------------------------------------------------------*/
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   printf("\n");
   fptr1 = fopen(file1, "r");
   fgets(dummy,80,fptr1);
   fgets(title,80,fptr1);
   fgets(dummy,80,fptr1);
   fscanf(fptr1,"%d %d %d %d %d %d\n", &nn, &ne, &nm, &ndim, &nen, &ndn);
   fgets(dummy, 80, fptr1);
   fscanf(fptr1,"%d %d %d %d %d\n", &nd, &nmpc);
   npr = 1;    /* Material property Thermal Conductivity */
   nmpc = 0;
   /* ---  nd = no. of specified temperatures --- */
   /* ---  nl = no. of nodal heat sources --- */
   /* note!! npr = 1 (thermal conductivity) and nmpc = 0 for this program */
   /* ehs(i) = element heat source, i = 1,...,ne */
/* -----------    memory allocation    ------------ */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   f = (float *) calloc(nn*ndn, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   ehs = (float *) calloc(ne, sizeof(float));
/* ------------------------------------------------ */ 
   printf("\n\n    PLOT CHOICE\n");
   printf("  1) no plot data\n");
   printf("  2) create data file containing nodal temperatures\n");
   printf("     choose <1 or 2> ");
   scanf("%d%*c", &ipl);
   if(ipl < 1 || ipl > 2)
     ipl = 1;         /* --- default is no data ---*/
   if(ipl > 1){
      printf("Output file name < dr:fn.ext >:\n");
      gets(file3);
      }
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
/* ----- connectivity, material#, element heat source ----- */
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
       ehs[n-1] = c;
   }
/* ----- temperature bc  ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nd; i++) {
      fscanf(fptr1, "%d %f\n", &k, &c);
      nu[i] = k;
      u[i] = c;
   }
/* ----- nodal heat sources ----- */
   for (i = 0; i < nn; i++) {
      f[i] = 0;
      }
   fgets(dummy,80,fptr1);
   for (i = 0; i < nl; i++) {
      fscanf(fptr1, "%d %f\n", &k, &c);
      f[k-1] = c;
   }	 
/* ----- thermal conductivity of material ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nm; i++){
      fscanf(fptr1, "%d", &k);
      for (j = 0; j < npr; j++) {
         fscanf(fptr1, "%f\n", &c);
	 pm[(k-1)*npr+j] = c;
      }
   }
/* ----- bandwidth nbw from connectivity noc() ----- */
   nbw = 0;
   for (i = 0; i < ne; i++) {
        nmin = noc[nen*i];
        nmax = nmin;
	for (j = 0; j < 3;j++) {
	    n = noc[nen*i+j];
	    if (nmin > n)
	       nmin = n;
	    if (nmax < n)
	       nmax = n;
	    }
	n = ndn * (nmax - nmin + 1);
	if (nbw < n)
           nbw = n;
   }    
   printf ("the bandwidth is %d\n", nbw);
/* ----- allocate memory for stiffness ----- */
   s = (float *) calloc(nn*nbw, sizeof(float));
  /* --- flux data if exist --- */
   fgets(dummy,120,fptr1);
   fscanf(fptr1, "%d\n", &nhf);
  if (nhf > 0) {
     for (i = 0; i < nhf; i++) {
	 fscanf(fptr1, "%d %d %f\n", &n1,&n2,&v);
         n1 = n1 -1;
         n2 = n2 -1;
         elen = (x[2*n1] - x[2*n2])*(x[2*n1] - x[2*n2]);
         elen = elen + (x[2*n1+1] - x[2*n2+1])*(x[2*n1+1] - x[2*n2+1]);
         elen = sqrt((double) elen);
         f[n1] = f[n1] - elen * v / 2;
         f[n2] = f[n2] - elen * v / 2;
         }
     }
  /* --- convection data if exist --- */
   fgets(dummy,80,fptr1);
   fscanf(fptr1, "%d\n", &nconv);
   if (nconv > 0) {
      for (i = 0; i < nconv; i++) {
	 fscanf(fptr1, "%d %d %f %f\n", &n1, &n2, &h, &tinf);
	 n1 = n1 -1;
	 n2 = n2 -1;
	 elen = (x[2*n1] - x[2*n2])*(x[2*n1] - x[2*n2]);
	 elen = elen + (x[2*n1+1] - x[2*n2+1])*(x[2*n1+1] - x[2*n2+1]);
	 elen = sqrt((double) elen);
	 f[n1] = f[n1] + elen * h * tinf / 2;
	 f[n2] = f[n2] + elen * h * tinf / 2;
	 s[nbw*n1] = s[nbw*n1] + h * elen / 3;
	 s[nbw*n2] = s[nbw*n2] + h * elen / 3;
	 if (n1 >= n2) {
	    k = n1;
	    n1 = n2;
	    n2 = k;
           }
         s[nbw*n1+n2-n1] = s[nbw*n1+n2-n1] + h * elen / 6;
         }
      }
   fclose(fptr1);
     /* --- element heat sources and conductivity --- */
     for (i = 0; i < ne; i++) {
         i1 = noc[nen*i]-1;
	 i2 = noc[nen*i+1]-1;
	 i3 = noc[nen*i+2]-1;
         x32 = x[2*i3] - x[2*i2];
	 x13 = x[2*i1] - x[2*i3];
	 x21 = x[2*i2] - x[2*i1];
         y23 = x[2*i2+1] - x[2*i3+1];
	 y31 = x[2*i3+1] - x[2*i1+1];
	 y12 = x[2*i1+1] - x[2*i2+1];
	 detj = x13 * y23 - x32 * y31;
	 area = .5 * fabs(detj);
      /* --- element heat sources if exist --- */
         eq = ehs[i];
	 if (fabs(eq) > 0) {
            c = eq * area / 3;
            f[i1] = f[i1] + c;
	    f[i2] = f[i2] + c;
	    f[i3] = f[i3] + c;
	 }
	 /* --- element stiffness and placing in global loc. --- */
	 bt[0][0] = y23 / detj;
	 bt[0][1] = y31 / detj;
	 bt[0][2] = y12 / detj;
	 bt[1][0] = x32 / detj;
	 bt[1][1] = x13 / detj;
	 bt[1][2] = x21 / detj;
	 for (ii = 0; ii < 3; ii++) {
	     ii1 = noc[nen*i+ii]-1;
	     for (jj = 0; jj < 3; jj++) {
		 ii2 = noc[nen*i+jj]-1;
		 if (ii1 <= ii2) {
		    sum = 0;
		    for (j = 0; j < 2; j++) {
		       sum = sum + bt[j][ii] * bt[j][jj];
			   }
		    n = nbw*ii1+ii2-ii1;
		    m = mat[i] - 1;
		    s[n] = s[n]+sum*area*pm[npr*m];
		   }

		}
	     }
	}
     if (nd > 0) {
       /* --- modify for temp. boundary conditions --- */
       sum = 0;
       for (i = 0; i < nn; i++) {
          sum = sum + s[nbw*i]/nn;
	      }
       cnst = sum * 1000000;
       for (i = 0; i < nd; i++) {
           n = nu[i]-1;
	   s[nbw*n] = s[nbw*n] + cnst;
           f[n] = f[n] + cnst * u[i];
	       }
     }
/* --- equation solving --- */
     bansol(s,f,nn,nbw);
   fptr1 = fopen(file2, "w");
   fprintf(fptr1, "%s\n", title);
   printf("%s\n", title);
   fprintf(fptr1, "node no.  temperature\n");
   printf("node no.  temperature\n");
   for (i = 0; i < nn; i++) {
     fprintf(fptr1, "%4d      %8.4f\n", i+1, f[i]);
     printf("%4d      %8.4f\n", i+1, f[i]);
     }
     if (ipl == 2 ) {
       fptr2 = fopen(file3, "w");
       fprintf( fptr2,  "nodal temperatures\n");
       for (i = 0; i < nn; i++) {
	      fprintf(fptr2, " %8.4f\n", f[i]);
	      }
	   fclose(fptr2);
       printf("\n");
	   printf ("nodal temperature data in file  %s \n", file3);
       printf ("run contour1 or contour2 to plot isotherms\n");
       }
   fprintf(fptr1, "conduction heat flow per unit area in each element \n");
   fprintf(fptr1, "element#   qx= -k*dt/dx  qy= -k*dt/dy \n");
     for (i = 0; i < ne; i++) {
         i1 = noc[nen*i]-1;
	 i2 = noc[nen*i+1]-1;
	 i3 = noc[nen*i+2]-1;
	 x32 = x[2*i3] - x[2*i2];
	 x13 = x[2*i1] - x[2*i3];
	 x21 = x[2*i2] - x[2*i1];
	 y23 = x[2*i2+1] - x[2*i3+1];
	 y31 = x[2*i3+1] - x[2*i1+1];
	 y12 = x[2*i1+1] - x[2*i2+1];
	 detj = x13 * y23 - x32 * y31;
	 bt[0][0] = y23 / detj;
	 bt[0][1] = y31 / detj;
	 bt[0][2] = y12 / detj;
	 bt[1][0] = x32 / detj;
	 bt[1][1] = x13 / detj;
	 bt[1][2] = x21 / detj;
     qx = bt[0][0] * f[i1] + bt[0][1] * f[i2] + bt[0][2] * f[i3];
     qx = -qx * pm[npr*(mat[i]-1)];
     qy = bt[1][0] * f[i1] + bt[1][1] * f[i2] + bt[1][2] * f[i3];
     qy = -qy * pm[npr*(mat[i]-1)];
     fprintf(fptr1,"%4d       %8.4f     %8.4f\n", i+1, qx, qy);
     }
     fclose(fptr1);
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