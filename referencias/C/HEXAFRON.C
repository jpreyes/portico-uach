/*****        program hexafron         *****/
/*   3-d stress analysis using  8-node     */
/*    isoparametric hexahedral element     */
/*          using frontal solver           */
/*   t.r.chandrupatla and a.d.belegundu    */
/*******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
struct data
  {
  int variable;
  double coefft;
  };
  FILE *fptr;
main()
{
   FILE *fptr1;
   long int icount;
   int n,i,j,k,m,nfron,ntogo,ndcnt,in,ii;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nmpc,ibl;
   int mtn,mtn1,ip;
   int *noc, *nu, *mat, *mpc, *isbl, *iebl, *indx;
   double *x, *pm, *u, *tempr, *s, *f, *beta;
   double c,dj,al,tld,cnst,reaction,s1,s2,s3,pi,cal,siv1,siv2,vm;
   double se[24][24],xi[3][8],xni[3][8],d[6][6];
   double b[6][24],db[6][24],qt[24],str[6];
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
   npr = 3;    /* Material properties E, Nu, Alpha */
   nq = nn * ndn;
   fgets(dummy, 80, fptr1);
   fscanf(fptr1,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
/* ----- memory allocation ----- */
   x = (double *) calloc(nn*ndim, sizeof(double));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (double *) calloc(nd, sizeof(double));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   f = (double *) calloc(nn*ndn, sizeof(double));
   tempr = (double *) calloc(ne, sizeof(double));
   pm = (double *) calloc(nm*npr, sizeof(double));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (double *) calloc(3*nmpc, sizeof(double));
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
   prefront(nn,ne,nen,ndn,nq,nmpc,noc,mpc,&ibl);
/* ------------------------------------------------- */
   isbl = (int *) calloc(ibl+1, sizeof(int));
   indx = (int *) calloc(ibl+1, sizeof(int));
   s = (double *) calloc(ibl*ibl, sizeof(double));
/* ------------------------------------------------- */
     nfron = 0;
	 ntogo = 0;
	 ndcnt = 0;
     for (i = 1; i <= ibl; i++) {
	    indx[i] = i;
	    }
     icount = 0;
/* =====  frontal assembly & eliminaton etc.  ===== */
/* ----- corner nodes and integration points */
     integ(xi,xni);
/* ----- open scratch file for writing ----- */
     if ((fptr = fopen("scratch.dat","wb"))==NULL)
       {printf("Can't open file scratch.dat"); exit(1);}
/* ----- element loop ----- */
     mtn1 = 0;
     for (n = 0; n < ne; n++) {
        printf("... forming stiffness matrix of element %d\n", n+1);
        mtn = mat[n];
	    if (mtn != mtn1)
	       dmat(mtn,&al,pm,d);
	       mtn1 = mtn;
	elstif(n,se,qt,xi,xni,d,b,db,tempr,x,al,noc);
        if (n == 0) {
           cnst = 0;
	   for (i = 0; i < 24; i++) {
	           cnst = cnst + se[i][i];
	           }
           cnst = 1e+11 * cnst;
           mpcfron(indx,isbl,mpc,nmpc,&nfron,s,f,ibl,beta,cnst);
	  }
     /* ----- add temperature load to the force vector f[] ----- */
       for (i = 0; i < 8; i++) {
	  ii = 3*(abs(noc[nen*n])-1);
	  for (j = 0; j < 3; j++) {
	     f[ii+j] = f[ii+j] + qt[3*i+j];
	     }
	  }
 /* frontal assembly  and forward elimination */
front(n,noc,nen,ndn,nd,&icount,indx,isbl,ibl,s,f,&nfron,&ntogo,&ndcnt,se,nu,cnst,u);
      }
fclose(fptr);
/* ----- assembly and reduction are complete */
/* ----- now backsubstitute */
     if ((fptr = fopen("scratch.dat","rb"))==NULL)
       {printf("Can't open file scratch.dat"); exit(1);}
     backsub(icount,f);
     fclose(fptr);
/* ----- printing displacements ----- */
   fptr1 = fopen(file2, "w");
   printf("\n%s\n", title);
   fprintf(fptr1, "\n%s\n", title);
 fprintf(fptr1, "node#     x-displ      y-displ      z-displ\n");
 printf ("node#     x-displ      y-displ      z-displ\n");
for (i = 0; i < nn; i++) {
 printf(" %4d  %11.4e  %11.4e  %11.4e\n",i+1,f[3*i],f[3*i+1],f[3*i+2]);
 fprintf(fptr1," %4d  %11.4e  %11.4e  %11.4e\n",i+1,f[3*i],f[3*i+1],f[3*i+2]);
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
     mtn1 = 0;
     for (n = 0; n < ne; n++) {
     fprintf(fptr1, "von mises stress at 8 int. pts. in elem# %d\n", n+1);
        mtn = mat[n];
	  if (mtn != mtn1)
	     dmat(mtn,&al,pm,d);
	     mtn1 = mtn;
      cal = al * tempr[n];
        for (ip = 0; ip < 8; ip++) {
           /* --- von mises stress at integration points */
	   dbmat(n,ip,x,noc,d,b,db,&dj,xi,xni);
      /* --- element nodal displacements stored in qt() */
           for (i = 0; i < 8; i++) {
              in = 3 * (abs(noc[nen*n+i]) - 1);
	      ii = 3 * i;
              for (j = 0; j < 3; j++) {
                 qt[ii + j] = f[in + j];
                }
             }
           /* --- stress calculation str = db * q */
           for (i = 0; i < 6; i++) {
              str[i] = 0;
              for (j = 0; j < 24; j++) {
                 str[i] = str[i] + db[i][j] * qt[j];
                 }
              str[i] = str[i] - cal * (d[i][1] + d[i][2] + d[i][3]);
              }
           /* --- calculation of von mises stress at ip */
           siv1 = str[0] + str[1] + str[2];
           siv2 = str[0] * str[1] + str[1] * str[2] + str[2] * str[0];
           siv2 = siv2 - str[3]*str[3] - str[4]*str[4] - str[5]*str[5];
           vm = sqrt((double) siv1 * siv1 - 3 * siv2);
           if (ip == 4)
              fprintf(fptr1,"\n");
           fprintf(fptr1, "   %11.4e", vm);
          }
        fprintf(fptr1,"\n");
       }
     printf("the results are saved in the file %s\n", file2);
     fclose(fptr1);
     return(0);
   }
/* end of main function */

integ(xi,xni)
double xi[][8], xni[][8];
{
  int i;
  double c;
/* ------- integration points xni() -------- */
     c = .57735026919;
     xi[0][0] = -1;
     xi[1][0] = -1;
     xi[2][0] = -1;
     xi[0][1] = 1;
     xi[1][1] = -1;
     xi[2][1] = -1;
     xi[0][2] = 1;
     xi[1][2] = 1;
     xi[2][2] = -1;
     xi[0][3] = -1;
     xi[1][3] = 1;
     xi[2][3] = -1;
     xi[0][4] = -1;
     xi[1][4] = -1;
     xi[2][4] = 1;
     xi[0][5] = 1;
     xi[1][5] = -1;
     xi[2][5] = 1;
     xi[0][6] = 1;
     xi[1][6] = 1;
     xi[2][6] = 1;
     xi[0][7] = -1;
     xi[1][7] = 1;
     xi[2][7] = 1;
     for (i = 0; i < 8; i++) {
        xni[0][i] = c * xi[0][i];
        xni[1][i] = c * xi[1][i];
        xni[2][i] = c * xi[2][i];
       }
     return(0);
 }
dmat(mtn,al,pm,d)
   int mtn;
   double *al,*pm,d[][6];
  {
   double e,pnu,c1,c2;
   int i,j;
   /* --- d() matrix relating stresses to strains */
     e = pm[3*(mtn-1)];
	 pnu = pm[3*(mtn-1)+1];
	 *al = pm[3*(mtn-1)+2];
     c1 = e / ((1 + pnu) * (1 - 2 * pnu));
     c2 = .5 * e / (1 + pnu);
     for (i = 0; i < 6; i++) {
	     for (j = 0; j < 6; j++) {
	        d[i][j] = 0;
	       }
	    }
     d[0][0] = c1 * (1 - pnu);
     d[0][1] = c1 * pnu;
     d[0][2] = d[0][1];
     d[1][0] = d[0][1];
     d[1][1] = d[0][0];
     d[1][2] = d[0][1];
     d[2][0] = d[0][2];
     d[2][1] = d[1][2];
     d[2][2] = d[0][0];
     d[3][3] = c2;
     d[4][4] = c2;
     d[5][5] = c2;
     return(0);
 }
elstif(n,se,qt,xi,xni,d,b,db,tempr,x,al,noc)
   int n,*noc;
   double al,*x,*tempr,d[][6],qt[],xi[][8],xni[][8];
   double se[][24],b[][24],db[][24];
{
   int i,j,k,ip;
   double dj,dte,c,dsum;
/* --------  element stiffness and temperature load  ----- */
     for (i = 0; i < 24; i++) {
	    for (j = 0; j < 24; j++) {
           se[i][j] = 0;
	       }
	     qt[i] = 0;
	     }
     dte = tempr[n];
     /* --- weight factor is one */
     /* --- loop on integration points */
     for (ip = 0; ip < 8; ip++) {
        printf("integration point = %d\n", ip);
        dbmat(n,ip,x,noc,d,b,db,&dj,xi,xni);
        /* --- element stiffness matrix  se */
        for (i = 0; i < 24; i++) {
           for (j = 0; j < 24; j++) {
              for (k = 0; k < 6; k++) {
                 se[i][j] = se[i][j] + b[k][i] * db[k][j] * dj;
                }
             }
          }
        /* --- determine temperature load qt() */
        c = al * dte;
        for (i = 0; i < 24; i++) {
	   dsum = db[0][i] + db[1][i] + db[2][i];
	   qt[i] = qt[i] + c * fabs(dj) * dsum / 6;
           }
     }
     return(0);
 }
dbmat(n,ip,x,noc,d,b,db,dj,xi,xni)
  int n,*noc,ip;
  double d[][6],b[][24],db[][24],xi[][8],xni[][8],*x,*dj;
{
 int i,j,k,kn,ir,ic;
 double aj[3][3],tj[3][3],gn[3][8],h[9][24],g[6][9];
 double dj1,dj2,dj3,c;
/* -------  db()  matrix  ------ */
/* --- gradient of shape functions - the gn() matrix */
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 8; j++) {
           c = 1;
           for (k = 0; k < 3; k++) {
              if (k != i)
                 c = c * (1 + xi[k][j] * xni[k][ip]);
              }
           gn[i][j] = .125 * xi[i][j] * c;
           }
        }
     /* --- formation of jacobian  tj */
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
           tj[i][j] = 0;
           for (k = 0; k < 8; k++) {
              kn = abs(noc[8*n+k]);
              tj[i][j] = tj[i][j] + gn[i][k] * x[3*(kn-1)+j];
              }
           }
         }
     /* --- determinant of the jacobian */
     dj1 = tj[0][0] * (tj[1][1] * tj[2][2] - tj[2][1] * tj[1][2]);
     dj2 = tj[0][1] * (tj[1][2] * tj[2][0] - tj[2][2] * tj[1][0]);
     dj3 = tj[0][2] * (tj[1][0] * tj[2][1] - tj[2][0] * tj[1][1]);
     *dj = dj1 + dj2 + dj3;
     /* --- inverse of the jacobian aj() */
     aj[0][0] = (tj[1][1] * tj[2][2] - tj[1][2] * tj[2][1]) / *dj;
     aj[0][1] = (tj[2][1] * tj[0][2] - tj[2][2] * tj[0][1]) / *dj;
     aj[0][2] = (tj[0][1] * tj[1][2] - tj[0][2] * tj[1][1]) / *dj;
     aj[1][0] = (tj[1][2] * tj[2][0] - tj[1][0] * tj[2][2]) / *dj;
     aj[1][1] = (tj[0][0] * tj[2][2] - tj[0][2] * tj[2][0]) / *dj;
     aj[1][2] = (tj[0][2] * tj[1][0] - tj[0][0] * tj[1][2]) / *dj;
     aj[2][0] = (tj[1][0] * tj[2][1] - tj[1][1] * tj[2][0]) / *dj;
     aj[2][1] = (tj[0][1] * tj[2][0] - tj[0][0] * tj[2][1]) / *dj;
     aj[2][2] = (tj[0][0] * tj[1][1] - tj[0][1] * tj[1][0]) / *dj;
     /* --- h() matrix relates local derivatives of  u  to local */
     /*     displacements  q */
     for (i = 0; i < 9; i++) {
        for (j = 0; j < 24; j++) {
           h[i][j] = 0;
           }
        }
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
           ir = 3 * i  + j ;
           for (k = 0; k < 8; k++) {
              ic = 3 * k + i;
              h[ir][ic] = gn[j][k];
              }
           }
        }
     /* --- g() matrix relates strains to local derivatives of  u */
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 9; j++) {
           g[i][j] = 0;
          }
        }
     g[0][0] = aj[0][0];
     g[0][1] = aj[0][1];
     g[0][2] = aj[0][2];
     g[1][3] = aj[1][0];
     g[1][4] = aj[1][1];
     g[1][5] = aj[1][2];
     g[2][6] = aj[2][0];
     g[2][7] = aj[2][1];
     g[2][8] = aj[2][2];
     g[3][3] = aj[2][0];
     g[3][4] = aj[2][1];
     g[3][5] = aj[2][2];
     g[3][6] = aj[1][0];
     g[3][7] = aj[1][1];
     g[3][8] = aj[1][2];
     g[4][0] = aj[2][0];
     g[4][1] = aj[2][1];
     g[4][2] = aj[2][2];
     g[4][6] = aj[0][0];
     g[4][7] = aj[0][1];
     g[4][8] = aj[0][2];
     g[5][0] = aj[1][0];
     g[5][1] = aj[1][1];
     g[5][2] = aj[1][2];
     g[5][3] = aj[0][0];
     g[5][4] = aj[0][1];
     g[5][5] = aj[0][2];
     /* --- b() matrix relates strains to  q */
     for (i = 0; i < 6; i++) {
        for (j = 0; j <  24; j++) {
            b[i][j] = 0;
            for (k = 0; k < 9; k++) {
              b[i][j] = b[i][j] + g[i][k] * h[k][j];
             }
           }
        }
     /* --- db() matrix relates stresses to  q */
     for (i = 0; i <  6; i++) {
        for (j = 0; j <  24; j++ ) {
           db[i][j] = 0;
           for (k = 0; k < 6; k++) {
              db[i][j] = db[i][j] + d[i][k] * b[k][j];
             }
           }
        }
     return(0);
}
prefront(nn,ne,nen,ndn,nq,nmpc,noc,mpc,ibl)
  int nn,ne,nen,ndn,nq,nmpc,*ibl;
  int *noc,*mpc;
 {
   int i,j,k,n,ifron,*ide,i1,ia,ineg;
        /* ----- mark last appearance of node / make it negative in noc() */
        /*  last appearance is first appearance for reverse element order */
        for (i = 1; i <= nn; i++) {
           for (j = ne-1; j >= 0; j--) {
              for (k = 0; k < nen; k++) {
                 if (i == noc[nen*j+k])
                    goto label1;
              }
           }
label1:
           noc[nen*j+k] = -i;
        }
    /* ===== block size determination */
     ide = (int *) calloc(nq+1, sizeof(int));
     for (i = 1; i <= nq; i++) {
	    ide[i] = 0;
	    }
     for (i = 0; i < nmpc; i++) {
	    for (j = 0; j < 2; j++) {
	       ide[mpc[2*i+j]] = 1;
	    }
	 }
    ifron = 0;
	for (i = 1; i <= nq; i++) {
	    ifron = ifron + ide[i];
	   }
    *ibl = ifron;
    for (n = 0; n < ne; n++) {
	   ineg = 0;
	   for (i = 0; i < nen; i++) {
	      i1 = noc[nen*n+i];
	          ia = ndn * (abs(i1) - 1);
              for (j = 0; j < ndn; j++) {
                 ia = ia + 1;
                 if (ide[ia] == 0) {
                    ifron = ifron + 1;
	                ide[ia] = 1;
                    }
               }
              if (i1 < 0)
                 ineg = ineg + 1;
           }
           if (*ibl < ifron)
              *ibl = ifron;
           ifron = ifron - ndn * ineg;
        }
        printf( "block size = %d\n", *ibl);
        return(0);
 }
mpcfron(indx,isbl,mpc,nmpc,nfron,s,f,ibl,beta,cnst)
   int *indx,*isbl,*mpc,nmpc,*nfron,ibl;
   double *s,*f,*beta,cnst;
   {
    int i,j,i1,j1,ifl,k,k1,i2;
/* ----- modifications for multipoint constraints by penalty method */
    for (i = 0; i < nmpc; i++) {
       i1 = mpc[2*i];
           ifl = 0;
           for (j = 1; j <= *nfron; j++) {
              j1 = indx[j];
              if (i1 == isbl[j1]) {
                 ifl = 1;
	             break;
                }
              }
	   if (ifl == 0) {
	      *nfron = *nfron + 1;
		  j1 = indx[*nfron];
		  isbl[j1] = i1;
	      }
	   i2 = mpc[2*i+1];
	   ifl = 0;
	   for (k = 1; k <= *nfron; k++) {
	      k1 = indx[k];
	      if (k1 == isbl[k1]) {
		 ifl = 1;
		     break;
		 }
	      }
	   if (ifl == 0) {
	      *nfron = *nfron + 1;
		  k1 = indx[*nfron];
		  isbl[k1] = i2;
	      }
	   /* ----- stiffness modification */
	   j1 = j1 - 1;
	   k1 = k1 - 1;
	   s[ibl*j1+j1] = s[ibl*j1+j1] + cnst * beta[3*i] * beta[3*i];
	   s[ibl*k1+k1] = s[ibl*k1+k1] + cnst * beta[3*i+1] * beta[3*i+1];
	   s[ibl*j1+k1] = s[ibl*j1+k1] + cnst * beta[3*i] * beta[3*i+1];
	   s[ibl*k1+j1] = s[ibl*j1+k1];
	   /* ----- force modification */
	   f[i1-1] = f[i1-1] + cnst * beta[3*i+2] * beta[3*i];
	   f[i2-1] = f[i2-1] + cnst * beta[3*i+2] * beta[3*i+1];
     }
    return(0);
   }


front(n,noc,nen,ndn,nd,icount,indx,isbl,ibl,s,f,nfron,ntogo,ndcnt,se,nu,cnst,u)
  int n,nen,ndn,nd,*noc,*indx,*isbl,ibl,*nfron,*ntogo,*ndcnt,*nu,*icount;
  double *s,*f,se[][24],cnst,*u;
 {
   struct data record;
   double pivot,c;
   int i,i1,ia,is1,idj,idf,ie1,ifl,ii,ix,j,j1;
   int itemp,ipv,ipg,ig,iebl[24],ntg1,iba;
   /* ----- frontal method assembly and elimination ----- */
   /* -------------  assembly of element n  ------------- */
        for (i = 0; i < nen; i++) {
           i1 = noc[nen*n+i];
	       ia = abs(i1);
	       is1 = 1;
	       if(i1 < 0)
	         is1 = -1;
           idf = ndn * (ia - 1);
	       ie1 = ndn * i;
           for (j = 0; j < ndn; j++) {
              idf = idf + 1;
	          ie1 = ie1 + 1;
	          ifl = 0;
              if (*nfron > *ntogo) {
                 for (ii = *ntogo+1; ii <= *nfron; ii++) {
                    ix = indx[ii];
		    if (idf == isbl[ix]) {
                       ifl = 1;
	                   break;
                       }
                   }
                }
              if (ifl == 0) {
		*nfron = *nfron + 1;
		    ii = *nfron;
	            ix = indx[ii];
                }
              isbl[ix] = idf;
	      iebl[ie1] = ix;
	      if (is1 == -1) {
                 *ntogo = *ntogo + 1;
                 itemp = indx[*ntogo];
                 indx[*ntogo] = indx[ii];
                 indx[ii] = itemp;
               }
           }
	}
	for (i = 0; i < 24; i++) {
           i1 = iebl[i+1]-1;
	   for (j = 0; j < 24; j++) {
              j1 = iebl[j+1]-1;
              s[ibl*i1+j1] = s[ibl*i1+j1] + se[i][j];
             }
          }
/* ------------------------------------------------------------------ */
     if (*ndcnt < nd) {
/* -----  modification for displacement bcs / penalty approach  ----- */
        for (i = 1; i <= *ntogo; i++) {
	   i1 = indx[i];
	   ig = isbl[i1];
	      for (j = 0; j < nd; j++) {
		 if (ig == nu[j]) {
		    i1 = i1 - 1;
		    s[ibl*i1+i1] = s[ibl*i1+i1] + cnst;
		    f[ig-1] = f[ig-1] + cnst * u[j];
		    *ndcnt = *ndcnt + 1;      /* counter for check */
		    break;
		   }
		}
	    }
       }
/* ------------   elimination of completed variables   --------------- */
        ntg1 = *ntogo;
        for (ii = 1; ii <= ntg1; ii++) {
           ipv = indx[1];
	       ipg = isbl[ipv];
           pivot = s[(ibl+1)*(ipv-1)];
        /* -----  write separator "0" and pivot value to disk  ----- */
     *icount = *icount + 1;
     record.variable = 0;
     record.coefft = pivot;
     fwrite(&record, sizeof(record),1,fptr);
       s[(ibl+1)*(ipv-1)] = 0;
           for (i = 2; i <= *nfron; i++) {
              i1 = indx[i];
	          ig = isbl[i1];
              if (s[ibl*(i1-1)+ipv-1] != 0) {
                  c = s[ibl*(i1-1)+ipv-1] / pivot;
	              s[ibl*(i1-1)+ipv-1] = 0;
                  for (j = 2; j <= *nfron; j++) {
                     j1 = indx[j];
                     if (s[ibl*(ipv-1)+j1-1] != 0)
         s[ibl*(i1-1)+j1-1] = s[ibl*(i1-1)+j1-1] - c * s[ibl*(ipv-1)+j1-1];
                     }
                  f[ig-1] = f[ig-1] - c * f[ipg-1];
                 }
              }
           for (j = 2; j <= *nfron; j++) {
        /* -----  write variable# and reduced coeff/pivot to disk  ----- */
              j1 = indx[j];
              if (s[ibl*(ipv-1)+j1-1] != 0) {
                 *icount = *icount + 1;
	             iba = isbl[j1];
	    record.variable = iba;
	    record.coefft = s[ibl*(ipv-1)+j1-1]/pivot;
	    fwrite(&record, sizeof(record),1,fptr);
	    s[ibl*(ipv-1)+j1-1] = 0;
                }
            }
           *icount = *icount + 1;
        /* -----  write eliminated variable# and rhs/pivot to disk  ----- */
     record.variable = ipg;
     record.coefft = f[ipg-1]/pivot;
     fwrite(&record, sizeof(record),1,fptr);
     f[ipg-1] = 0;
	/* ----- (*ntogo) into (1); (*nfron) into (*ntogo) */
	/* ----- ipv into (*nfron) and reduce front & *ntogo sizes by 1 */
	   if (*ntogo > 1)
	      indx[1] = indx[*ntogo];
	   indx[*ntogo] = indx[*nfron];
	       indx[*nfron] = ipv;
           *nfron = *nfron - 1;
	       *ntogo = *ntogo - 1;
        }
        return(0);
 }

 backsub(icount,f)
   long int icount;
   double *f;
   {
   struct data record;
   long int offset;
   int n1,n2;
     /* ===== backsubstitution */;
    while (icount > 0) {
	offset = (icount-1) * sizeof(record);
	fseek(fptr,offset,0);
	fread(&record,sizeof(record),1,fptr);
	icount = icount -1;
	n1 = record.variable;
	f[n1-1] = record.coefft;
     while (icount > 0) {
       offset = (icount-1) * sizeof(record);
       fseek(fptr,offset,0);
       fread(&record,sizeof(record),1,fptr);
       icount = icount -1;
       n2 = record.variable;
	if (n2 == 0)
	   break;
	f[n1-1] = f[n1-1] - record.coefft * f[n2-1];
	  }
	}
	return(0);
 }
