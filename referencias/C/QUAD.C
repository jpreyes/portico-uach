/********       program quad        **********/
/*     2-d stress analysis using 4-node      */
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
   float xi,eta,th,q[8],c1,sv;
   float c, dj, al, pnu, tld, cnst, reaction, s1, s2, s3, ang, r;
   float b[3][8],d[3][3],db[3][8],se[8][8],str[3],tl[8],xni[4][2];
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
   thick = (float *) calloc(ne, sizeof(float));
   f = (float *) calloc(nn*ndn, sizeof(float));
   tempr = (float *) calloc(ne, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
   printf("\n\n    PLOT CHOICE\n");
   printf("  1) no plot data\n");
   printf("  2) create data file for in-plane shear stress\n");
   printf("  3) create data file for von mises stress\n");
   printf("     choose <1 or 2 or 3> ");
   scanf("%d%*c", &ipl);
   if(ipl < 1 || ipl > 3)
     ipl = 1;         /* --- default is no data ---*/
   if(ipl > 1){
      printf("Output file name < dr:fn.ext >:\n");
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
       fscanf(fptr1,"%f",&c);
       thick[n-1] = c;
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
	dmatrix(n,pm,mat,npr,&pnu,&al,lc,d);
	/* --- element stiffness --- */
	elstif(n,lc,se,tl,xni,d,thick,tempr,x,al,pnu,noc);
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
     if (lc == 1)
         fprintf(fptr1, "plane stress analysis\n");
     if (lc == 2)
         fprintf(fptr1, "plane strain analysis\n");
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
     if (ipl > 1){
        fptr2 = fopen(file3, "w");
        if (ipl == 2)
	    fprintf(fptr2, "max. in-plane Shear Stress");
	    if (ipl == 3)
	       fprintf(fptr2, "von Mises stress");
	    fprintf(fptr2, "(element) for data in file  %s\n", file1);
        }
/* -----  stress calculations ----- */
     fprintf (fptr1, "elem#    von mises stresses at 4 integration points\n");
   /* ----- stresses at integration points ----- */
     for (n = 0; n < ne; n++) {
        fprintf (fptr1, "%4d ", n+1);
        for (ip = 0; ip < 4; ip++) {
	   xi = xni[ip][0];
	   eta = xni[ip][1];
	   dmatrix(n,pm,mat,npr,&pnu,&al,lc,d);
	   dbmat(n,x,noc,thick,&th,d,b,db,&dj,xi,eta);
	   /* --- stress evaluation --- */
	   for (i = 0; i < nen; i++) {
	      in = ndn * (noc[nen*n+i] - 1);
	      ii = ndn * i;
	      for (j = 0; j < ndn; j++) {
		 q[ii + j] = f[in + j];
		 }
	      }
	   c1 = al * tempr[n];
	   if (lc == 2)
              c1 = c1 * (1 + pnu);
           for (i = 0; i < 3; i++) {
              c = 0;
              for (k = 0; k < 8; k++) {
                  c = c + db[i][k] * q[k];
                  }
                 str[i] = c - c1 * (d[i][0] + d[i][1]);
              }
        /* --- von mises stress at integration point --- */
           c = 0;
	   if (lc == 2)
	      c = pnu * (str[0] + str[1]);
	   c1 = (str[0] - str[1]) * (str[0] - str[1]);
	   c1 = c1 + (str[1] - c) * (str[1] - c);
	   c1 = c1 + (c - str[0]) * (c - str[0]);
	   sv = sqrt((double)(.5 * c1 + 3 * str[2] * str[2]));
	   fprintf(fptr1, " %10.4e ", sv);
        /* --- maximum shear stress r --- */
           c = .25 * (str[0]-str[1])*(str[0]-str[1]);
           c = c + str[2]*str[2];
           r = sqrt((double) c);
	   if (ipl == 2)
              fprintf(fptr2," %f ", r);
	   if (ipl == 3)
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
	printf("run bestfit and then contourA or contourB to plot stresses\n");
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

dmatrix(n,pm,mat,npr,pnu1,al1,lc,d)
   int lc,n,npr,*mat;
   float *pm,*pnu1,*al1,d[][3];
  {
     int m;
     float e,c,c1,c2,c3,pnu,al;
   /* -----  d() matrix  ----- */
     /* --- material properties --- */
     m = mat[n]-1;
     e = pm[npr*m];
     pnu= pm[npr*m+1];
     al = pm[npr*m+2];
     *pnu1 = pnu;
     *al1 = al;
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
     return(0);
  }

elstif(n,lc,se,tl,xni,d,thick,tempr,x,al,pnu,noc)
    int n,lc,*noc;
    float al,pnu;
    float *x,*tempr,*thick,d[][3],tl[8],se[][8],xni[][2];
  {
    int i,j,k,ip;
    float dte,c,xi,eta,th,dj,b[3][8],db[3][8];
    /* -----  element stiffness and temperature load  ----- */
     for (i = 0; i < 8;i++) {
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
	dbmat(n,x,noc,thick,&th,d,b,db,&dj,xi,eta);
	/* --- element stiffness matrix  se --- */
	for (i = 0; i < 8; i++) {
           for (j = 0; j < 8; j++) {
              c = 0;
              for (k = 0; k < 3; k++) {
                 c = c + b[k][i] * db[k][j] * dj * th;
                 }
              se[i][j] = se[i][j] + c;
              }
	   }
	/* --- determine temperature load tl --- */
	c = al * dte;
	if (lc == 2)
	    c = (1 + pnu) * c;
        for (i = 0; i < 8; i++) {
           tl[i] = tl[i] + th * dj * c * (db[0][i] + db[1][i]);
           }
	}
     return(0);
   }
dbmat(n,x,noc,thick,th1,d,b,db,dj1,xi,eta)
  float *x,*dj1,*thick,*th1,xi,eta;
  float d[][3],b[][8],db[][8];
  int n,*noc;
  {
   int n1,n2,n3,n4,i,j,k;
   float x1,y1,x2,y2,x3,y3,x4,y4,tj11,tj12,tj21,tj22,dj,c;
   float th,a[3][4],g[4][8];
   /* -----  db()  matrix  ----- */
     /* --- nodal coordinates --- */     
     th = thick[n];
     *th1 = th;
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
     /* --- b[3,8] matrix relates strains to q --- */
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 8; j++) {
           c = 0;
           for (k = 0; k < 4; k++) {
               c = c + a[i][k] * g[k][j];
               }
           b[i][j] = c;
           }
        }
     /* --- db[3,8] matrix relates stresses to q[8] --- */
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 8; j++) {
           c = 0;
           for (k = 0; k < 3; k++) {
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