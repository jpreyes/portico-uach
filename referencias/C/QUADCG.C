/********     program quadcg        **********/
/*     2-d stress analysis using 4-node      */
/*  quadrilateral elements with temperature  */
/*         conjugate gradient approach       */
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
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nmpc,lc,ipl;
   int *noc, *nu, *mat, *mpc;
   double *x, *thick, *pm, *u, *tempr, *s, *f, *beta;
   double *ad, *dd, *gg, *q;
   double xi,eta,th,bt,qe[8],c1,sv;
   double c, dj, al, pnu, tld, cnst, reaction, s1, s2, s3, ang, r;
   double b[3][8],d[3][3],db[3][8],str[3],tl[8],xni[4][2];
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
   npr = 3;  /* Material properties E, Nu, Alpha */
/* ----- memory allocation ----- */
   x = (double *) calloc(nn*ndim, sizeof(double));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (double *) calloc(nd, sizeof(double));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   thick = (double *) calloc(ne, sizeof(double));
   f = (double *) calloc(nn*ndn, sizeof(double));
   ad = (double *) calloc(nn*ndn, sizeof(double));
   dd = (double *) calloc(nn*ndn, sizeof(double));
   gg = (double *) calloc(nn*ndn, sizeof(double));
   q = (double *) calloc(nn*ndn, sizeof(double));
   tempr = (double *) calloc(ne, sizeof(double));
   pm = (double *) calloc(nm*npr, sizeof(double));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (double *) calloc(3*nmpc, sizeof(double));
/* ----- allocate memory for stiffness ----- */
   s = (double *) calloc(ne*8*8, sizeof(double));
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
		fscanf(fptr1, "%lf\n", &c);
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
       fscanf(fptr1,"%lf",&c);
       thick[n-1] = c;
       fscanf(fptr1,"%lf\n",&c);
       tempr[n-1] = c;
   }
/* ----- displacement bc  ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nd; i++) {
      fscanf(fptr1, "%d %lf\n", &k, &c);
      nu[i] = k;
      u[i] = c;
   }
/* ----- component loads ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nl; i++) {
      fscanf(fptr1, "%d %lf\n", &k, &c);
      f[k-1] = c;
   }
/* ----- material properties ----- */
   fgets(dummy,80,fptr1);
   for (i = 0; i < nm; i++){
      fscanf(fptr1, "%d", &k);
      for (j = 0; j < npr; j++) {
         fscanf(fptr1, "%lf\n", &c);
	 pm[(k-1)*npr+j] = c;
      }
   }
/* ----- multipoint constraints ----- */
   if (nmpc > 0) 
      { fgets(dummy,80,fptr1);
	for(j=0;j<nmpc;j++){
	   fscanf(fptr1,"%lf",&c);
	   beta[3*j]=c;
	   fscanf(fptr1,"%d",&k);
	   mpc[2*j]=k;
           fscanf(fptr1,"%lf",&c);
	   beta[3*j+1]=c;
           fscanf(fptr1,"%d",&k);
           mpc[2*j+1]=k;
	   fscanf(fptr1,"%lf",&c);
	   beta[3*j+2]=c;
	   }
	}
   fclose (fptr1);
/* ----- global stiffness matrix -----*/
/* ----- corner nodes and integrationpoints ----- */
    integ(xni);
    for (n = 0; n < ne; n++) {
	printf("forming stiffness matrix of element %d\n", n+1);
	dmatrix(n,pm,mat,npr,&pnu,&al,lc,d);

	/* --- element stiffness --- */
	elstif(n,lc,s,tl,xni,d,thick,tempr,x,al,pnu,noc);
	printf (".... placing in global locations\n");
	for (ii = 0; ii < nen; ii++) {
	   nrt = ndn * (noc[nen*n + ii] - 1);
	   for (it = 0; it < ndn; it++) {
	      nr = nrt + it;
	      i = ndn * ii + it;
	      f[nr] = f[nr] + tl[i];
	      gg[nr] = gg[nr] + s[64*n + 8*i + i];
	      }
	   }
    }

    /* ----- decide penalty parameter cnst ----- */
/* ----- GG() diagonal stiffness summation   */
     cnst = 0.;
     for (i = 0; i < nq; i++) {
	 if (cnst < gg[i])
	    cnst = gg[i];
	 }
     cnst = cnst * 10000.;
/* ----- Modify right hand side F() for Boundary Conditions ----- */
/* ----- Displacement BC */
   for (i = 0; i < nd; i++) {
      k = nu[i];
      f[k-1] = f[k-1] + cnst * u[i];
   }
/* ----- modify for multipoint constraints ----- */
   for (i = 0; i < nmpc; i++){
       i1 = mpc[2*i]-1;
       i2 = mpc[2*i+1]-1;
       f[i1] = f[i1] + cnst*beta[3*i]*beta[3*i+2];
       f[i2] = f[i2] + cnst*beta[3*i+1]*beta[3*i+2];
       }

/* ----- solution of equations using conjgate gradient method ----- */
   cgsolve(s,f,q,ad,dd,gg,noc,cnst,nu,mpc,nmpc,beta,nd,nq,ne);
/* ----- printing displacements ----- */
   fptr1 = fopen(file2, "w");
   printf("\n%s\n", title);
   fprintf(fptr1, "\n%s\n", title);
     if (lc == 1)
         fprintf(fptr1, "plane stress analysis\n");
     if (lc == 2)
         fprintf(fptr1, "plane strain analysis\n");
     fprintf(fptr1, "node#     x-displ      y-displ\n");
     printf ("node#     x-displ      y-displ\n");
     for (i = 0; i < nn; i++) {
	 printf(" %4d  %11.4e  %11.4e\n",i+1,q[2*i],q[2*i+1]);
	 fprintf(fptr1," %4d  %11.4e  %11.4e\n",i+1,q[2*i],q[2*i+1]);
	 }
/* ----- reaction calculation ----- */
   printf("node#    reaction\n");
   fprintf(fptr1, "node#     reaction\n");
   for (i = 0; i < nd; i++) {
      k = nu[i];
      reaction = cnst * (u[i] - q[k-1]);
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
		 qe[ii + j] = q[in + j];
		 }
	      }
	   c1 = al * tempr[n];
	   if (lc == 2)
              c1 = c1 * (1 + pnu);

           for (i = 0; i < 3; i++) {
              c = 0;
              for (k = 0; k < 8; k++) {
                  c = c + db[i][k] * qe[k];
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
	printf("run bestfitq and then contourA or contourB to plot stresses\n");
	}

     return(0);
   }

integ(xni)
  double xni[][2];
  {
    double c;
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
   double *pm,*pnu1,*al1,d[][3];
  {
     int m;
     double e,c,c1,c2,c3,pnu,al;
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

  elstif(n,lc,s,tl,xni,d,thick,tempr,x,al,pnu,noc)
    int n,lc,*noc;
    double al,pnu;
    double *x,*tempr,*thick,d[][3],tl[8],s[][8][8],xni[][2];
  {
    int i,j,k,ip;
    double dte,c,xi,eta,th,dj,b[3][8],db[3][8];
    /* -----  element stiffness and temperature load  ----- */
     for (i = 0; i < 8;i++) {
	 for (j = 0; j < 8; j++) {
	     s[n][i][j] = 0.;
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
              s[n][i][j] = s[n][i][j] + c;
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
  double *x,*dj1,*thick,*th1,xi,eta;
  double d[][3],b[][8],db[][8];
  int n,*noc;
  {
   int n1,n2,n3,n4,i,j,k;
   double x1,y1,x2,y2,x3,y3,x4,y4,tj11,tj12,tj21,tj22,dj,c;
   double th,a[3][4],g[4][8];
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

/* ----- cgsolve ----- */
cgsolve(s,f,q,ad,dd,gg,noc,cnst,nu,mpc,nmpc,beta,nd,nq,ne)
  int nq, *nu, nmpc, nd, ne, *noc, mpc[][2];
  double s[][8][8], *f, *q, *ad, *dd, *gg, cnst, beta[][3];
{
 int i,j,n,ii,jj,i1,j1,i2,il,jl,ig,jg,igy,ilt;
 int igt,jgt,jlt,iter=0;
  double gg1,gg2,dad,c,al,bta;
  /* ----- Conjugate Gradient Method ----- */
  for (i = 0; i < nq; i++) {
      gg[i] = -f[i];
      dd[i] = f[i];
      q[i] = 0.;
      gg1 = gg1 + gg[i] * gg[i];
     }
  /* --- iteration loop --- */
  do {
     iter = iter + 1;
     /* =====  element loop  ===== */
     for (n = 0; n < ne; n++) {
	for (i = 0; i < 4; i++) {

           igt = 2 * (noc[4*n + i] - 1);
	   ilt = 2 * i;
           for (ii = 0; ii < 2; ii++){
              ig = igt + ii;
              il = ilt + ii;
              for (j = 0; j < 4; j++) {
                 jgt = 2 * (noc[4*n +j] - 1);
		 jlt = 2 * j;
                 for (jj = 0; jj < 2; jj++){
                    jg = jgt + jj;
                    jl = jlt + jj;
                    ad[ig] = ad[ig] + s[n][il][jl] * dd[jg];
                 }
              }
           }
         }
      }
     /* --- displacement bc --- */
     for (i = 0; i < nd; i++) {
        n = nu[i] - 1;
	ad[n] = ad[n] + cnst * dd[n];
      }
     /* --- multi-point constraints --- */
     for (i = 0; i < nmpc; i++) {
        i1 = mpc[i][1];
        i2 = mpc[i][2];
	c = beta[i][1] * dd[i1] + beta[i][2] * dd[i2];
        ad[i1] = ad[i1] + cnst * beta[i][1] * c;
        ad[i2] = ad[i2] + cnst * beta[i][2] * c;
      }
        dad = 0.;
        for (i = 0; i < nq; i++) {
           dad = dad + dd[i] * ad[i];
        }
        al = gg1 / dad;
        gg2 = 0.;
        for (i = 0; i < nq; i++) {
           gg[i] = gg[i] + al * ad[i];
           q[i] = q[i] + al * dd[i];
           gg2 = gg2 + gg[i] * gg[i];
	 }
	if (gg2 > 0.00000001) {
          bta = gg2 / gg1;
          gg1 = gg2;
          for (i = 0; i < nq; i++) {
           dd[i] = -gg[i] + bta * dd[i];
          }
          for (i = 0; i < nq; i++) {
            ad[i] = 0.;
          }
        }
   } while (gg2 > 0.00000001);
    return(0);
}