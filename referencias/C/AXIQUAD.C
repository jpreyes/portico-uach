/********     program axiquad       **********/
/* axisymmetric stress analysis using 4-node */
/*  quadrilateral elements with temperature  */
/*    t.r.chandrupatla and a.d.belegundu     */
/*********************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   int n,i,j,k,m,i1,i2,i3,ii,jj,m1,nmin,nmax,nrt,nct,it,jt;
   int ip,nr,nc,in;
   char dummy[81], title[81], file1[81], file2[81], file3[81];
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc,lc,ipl;
   int *noc, *nu, *mat, *mpc;
   float *x, *thick, *pm, *u, *tempr, *s, *f, *beta;
   float xi,eta,th,q[8],c1,c2,sv,rad;
   float c, dj, al, pnu, tld, cnst, reaction, s1, s2, s3, ang, r;
   float b[4][8],d[4][4],db[4][8],se[8][8],str[4],tl[8],xni[4][2];
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
   fscanf(fptr1,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
   npr = 3;   /* Material properties E, Nu, Alpha */
/* ----- memory allocation ----- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   f = (float *) calloc(nn*ndn, sizeof(float));
   tempr = (float *) calloc(ne, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
   printf("\n\n    PLOT CHOICE\n");
   printf("  1) no plot data\n");
   printf("  2) create data file for von mises stress\n");
   printf("     choose <1 or 2> ");
   scanf("%d%*c", &ipl);
   if(ipl < 1 || ipl > 2)
     ipl = 1;         /* --- default is no data ---*/
   if(ipl > 1){
      printf("Output file name for plot data < dr:fn.ext >:\n");
      gets(file3);
      }
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
/* ----- connectivity, material, temp-change ----- */
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
/* ----- global stiffness matrix -----*/
    /* ----- corner nodes and integrationpoints ----- */
    integ(xni);
    for (n = 0; n < ne; n++) {
	printf("forming stiffness matrix of element %d\n", n+1);
	dmatrix(n,pm,mat,npr,&al,d);
	/* --- element stiffness --- */
	elstif(n,se,tl,xni,d,tempr,x,al,noc);
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
     fprintf(fptr1, "node#     r-displ      z-displ\n");
     printf ("node#     r-displ      z-displ\n");
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
     if (ipl > 1){
        fptr2 = fopen(file3, "w");
	  fprintf(fptr2, "von Mises stress (element) for data in file  %s\n", file1);
        }
/* -----  stress calculations ----- */
     fprintf (fptr1, "elem#    von mises stresses at 4 integration points\n");
   /* ----- stresses at integration points ----- */

     for (n = 0; n < ne; n++) {
        fprintf (fptr1, "%4d ", n+1);
        for (ip = 0; ip < 4; ip++) {
	   xi = xni[ip][0];
	   eta = xni[ip][1];
	   dmatrix(n,pm,mat,npr,&al,d);
	   dbmat(n,x,noc,d,b,db,&dj,&rad,xi,eta);
	   /* --- stress evaluation --- */
	   for (i = 0; i < nen; i++) {
	      in = ndn * (noc[nen*n+i] - 1);
	      ii = ndn * i;
	      for (j = 0; j < ndn; j++) {
		 q[ii + j] = f[in + j];
		 }
	      }
	   c1 = al * tempr[n];
           for (i = 0; i < 4; i++) {
              c = 0;
              for (k = 0; k < 8; k++) {
                  c = c + db[i][k] * q[k];
                  }
                 str[i] = c - c1 * (d[i][0] + d[i][1] + d[i][3]);
              }
        /* --- von mises stress at integration point --- */
	   c1 = str[0] + str[1] + str[3];
           c2 = str[0] * str[1] + str[1] * str[3] + str[3] * str[0];
           c2 = c2 - str[2] * str[2];
           sv = sqrt((double)(c1 * c1 - 3 * c2));
	   fprintf(fptr1, " %10.4e ", sv);
	   if (ipl == 2)
	      fprintf(fptr2, " %f ", sv);
       }   
	   fprintf(fptr1, "\n");
	   if (ipl > 1)
	      fprintf(fptr2, "\n");
     }
     fclose(fptr1);
     printf("complete results are in file %s\n", file2);
     printf("view using a text processor\n");
     if (ipl > 1) {
        fclose(fptr2);
        printf("element stress data in file %s\n", file3);
        printf("run bestfitq and then contourA or contourB to plot stresses\n");
        }
     return(0);
   }
integ(xni)
  float xni[][2];
  {
    float c;
   /* ----- integration points xni() ----- */
     c = .57735026919;
     xni[0][0] = -c;
     xni[0][1] = -c;
     xni[1][0] = c;
     xni[1][1] = -c;
     xni[2][0] = c;
     xni[2][1] = c;
     xni[3][0] = -c;
     xni[3][1] = c;
     return(0);
  }

dmatrix(n,pm,mat,npr,al1,d)
   int n,npr,*mat;
   float *pm,*al1,d[][4];
  {
     int m;
     float e,c,c1,c2,c3,pnu,al;
   /* -----  d() matrix  ----- */
     /* --- material properties --- */
     m = mat[n]-1;
     e = pm[npr*m];
     pnu= pm[npr*m+1];
     al = pm[npr*m+2];
     *al1 = al;
     /* --- d() matrix --- */
     c1 = e * (1 - pnu) / ((1 + pnu) * (1 - 2 * pnu));
     c2 = pnu / (1 - pnu);
     d[0][0] = c1;
     d[0][1] = c1*c2;
     d[0][2] = 0;
     d[0][3] = c1*c2;
     d[1][0] = c1*c2;
     d[1][1] = c1;
     d[1][2] = 0;
     d[1][3] = c1*c2;
     d[2][0] = 0;
     d[2][1] = 0;
     d[2][2] = .5*e/(1+pnu);
     d[2][3] = 0;
     d[3][0] = c1*c2;
     d[3][1] = c1*c2;
     d[3][2] = 0;
     d[3][3] = c1;
     return(0);
  }

elstif(n,se,tl,xni,d,tempr,x,al,noc)
    int n,*noc;
    float al;
    float *x,*tempr,d[][4],tl[8],se[][8],xni[][2];
  {
    int i,j,k,ip;
    float dte,c,xi,eta,dj,rad,b[4][8],db[4][8];
    float pi = 3.14159;
    /* -----  element stiffness and temperature load  ----- */
     for (i = 0; i < 8; i++) {
	 for (j = 0; j < 8; j++) {
	     se[i][j] = 0.;
	     }
	 tl[i] = 0.;
	 }
     dte = tempr[n];
     /* --- weight factor is one --- */
     /* --- loop on integration points --- */
     for (ip = 0; ip < 4; ip++) {
        /* ---  get db matrix at integration point ip --- */
        xi = xni[ip][0];
	eta = xni[ip][1];
	dbmat(n,x,noc,d,b,db,&dj,&rad,xi,eta);
	/* --- element stiffness matrix  se --- */
	for (i = 0; i < 8; i++) {
           for (j = 0; j < 8; j++) {
              c = 0;
              for (k = 0; k < 4; k++) {
                 c = c + 2 * pi * rad * b[k][i] * db[k][j] * dj;
                 }
              se[i][j] = se[i][j] + c;
              }
	   }
	/* --- determine temperature load tl --- */
	c = 2 * pi * rad * dj * al * dte;
        for (i = 0; i < 8; i++) {
           tl[i] = tl[i] + c * (db[0][i] + db[1][i]);
           }
	}
     return(0);
   }
dbmat(n,x,noc,d,b,db,dj1,rad1,xi,eta)
  float *x,*dj1,*rad1,xi,eta;
  float d[][4],b[][8],db[][8];
  int n,*noc;
  {
   int n1,n2,n3,n4,i,j,k;
   float x1,y1,x2,y2,x3,y3,x4,y4,tj11,tj12,tj21,tj22,dj,c;
   float rad,sh1,sh2,sh3,sh4;
   float a[3][4],g[4][8];
   /* -----  db()  matrix  ----- */
     /* --- nodal coordinates --- */     
     n1 = noc[4*n];
     n2 = noc[4*n+1];
     n3 = noc[4*n+2];
     n4 = noc[4*n+3];
     x1 = x[2*(n1-1)];
     y1 = x[2*(n1-1)+1];
     x2 = x[2*(n2-1)];
     y2 = x[2*(n2-1)+1];
     x3 = x[2*(n3-1)];
     y3 = x[2*(n3-1)+1];
     x4 = x[2*(n4-1)];
     y4 = x[2*(n4-1)+1];
     /* --- formation of jacobian  tj --- */
     tj11 = ((1 - eta) * (x2 - x1) + (1 + eta) * (x3 - x4)) / 4;
     tj12 = ((1 - eta) * (y2 - y1) + (1 + eta) * (y3 - y4)) / 4;
     tj21 = ((1 - xi) * (x4 - x1) + (1 + xi) * (x3 - x2)) / 4;
     tj22 = ((1 - xi) * (y4 - y1) + (1 + xi) * (y3 - y2)) / 4;
     /* --- determinant of the jacobian --- */
     dj = tj11 * tj22 - tj12 * tj21;
     *dj1 = dj;
     /* --- a[3,4] matrix relates strains to --- */
     /* --- local derivatives of u --- */
     a[0][0] = tj22 / dj;
     a[1][0] = 0;
     a[2][0] = -tj21 / dj;
     a[0][1] = -tj12 / dj;
     a[1][1] = 0;
     a[2][1] = tj11 / dj;
     a[0][2] = 0;
     a[1][2] = -tj21 / dj;
     a[2][2] = tj22 / dj;
     a[0][3] = 0;
     a[1][3] = tj11 / dj;
     a[2][3] = -tj12 / dj;
     /* --- g[4,8] matrix relates local derivatives of u --- */
     /* --- to local nodal displacements q[8] --- */
     for (i = 0; i < 4; i++) {
	     for (j = 0; j < 8; j++) {
             g[i][j] = 0;
	         }
	     }
     g[0][0] = -(1 - eta) / 4;
     g[1][0] = -(1 - xi) / 4;
     g[2][1] = -(1 - eta) / 4;
     g[3][1] = -(1 - xi) / 4;
     g[0][2] = (1 - eta) / 4;
     g[1][2] = -(1 + xi) / 4;
     g[2][3] = (1 - eta) / 4;
     g[3][3] = -(1 + xi) / 4;
     g[0][4] = (1 + eta) / 4;
     g[1][4] = (1 + xi) / 4;
     g[2][5] = (1 + eta) / 4;
     g[3][5] = (1 + xi) / 4;
     g[0][6] = -(1 + eta) / 4;
     g[1][6] = (1 - xi) / 4;
     g[2][7] = -(1 + eta) / 4;
     g[3][7] = (1 - xi) / 4;
     /* --- Shape Function Values --- */
     sh1 = .25 * (1 - xi) * (1 - eta);
     sh2 = .25 * (1 + xi) * (1 - eta);
     sh3 = .25 * (1 + xi) * (1 + eta);
     sh4 = .25 * (1 - xi) * (1 + eta);
     rad = sh1 * x1 + sh2 * x2 + sh3 * x3 + sh4 * x4;
     *rad1 = rad;
     /* --- b[4,8] matrix relates strains to q --- */
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 8; j++) {
           c = 0;
           for (k = 0; k < 4; k++) {
               c = c + a[i][k] * g[k][j];
               }
           b[i][j] = c;
           }
        }
     b[3][0] = sh1 / rad; 
     b[3][1] = 0;
     b[3][2] = sh2 / rad; 
     b[3][3] = 0;
     b[3][4] = sh3 / rad;
     b[3][5] = 0;
     b[3][6] = sh4 / rad; 
     b[3][7] = 0;
     /* --- db[3,8] matrix relates stresses to q[8] --- */
     for (i = 0; i < 4; i++) {
        for (j = 0; j < 8; j++) {
           c = 0;
           for (k = 0; k < 4; k++) {
               c = c + d[i][k] * b[k][j];
               }
	   db[i][j] = c;
	   }
	}
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