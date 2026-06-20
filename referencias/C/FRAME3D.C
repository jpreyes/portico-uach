/*======================================*/
/*               frame3d.c              */
/*     3-d  frame analysis by fem       */
/*   t.r.chandrupatla & a.d.belegundu   */
/* =====================================*/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr;
   int n,nq,i,j,k,m,i1,i2,ii,nrt,it,nr,jj,nct,jt,nc,nnt;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nm,nd,nl,nen,ndn,ndim,npr,nbw,nch,nmpc,nnref;
   int *noc,*nu,*mat,*mpc;
   float *x,*arin,*pm,*u,*s,*f,*beta,*udl,se[12][12],ed[12];
   float edp[12],alambda[12][12],sep[12][12],ef[12];
   float c,el,eil,cnst,reaction;
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr = fopen(file1, "r");
   fgets(dummy,80,fptr);
   fgets(title,80,fptr);
   fgets(dummy,80,fptr);
   fscanf(fptr,"%d %d %d %d %d %d %d\n", &nn,&ne,&nm,&ndim,&nen,&ndn,&nnref);
   fgets(dummy, 80, fptr);
   fscanf(fptr,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
   npr = 2;   /* Material properties E, G (shear modulus) */
   nch = 4;   /* Area, Iy,  Iz,  J  for each element */
   /* ----- total dof is nq ----- */
   nq = ndn*nn;
   nnt = nn + nnref;
/* --------------- memory allocation ---------------- */
   x = (float *) calloc(nnt*ndim, sizeof(float));
   noc = (int *) calloc((nen+1)*ne, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   arin = (float *) calloc(nch*ne, sizeof(float));
   f = (float *) calloc(nn*ndn, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
   udl = (float *) calloc(2*ne, sizeof(float));
/* ------------------ read data --------------------- */
/* ----- coordinates ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nnt; i++) {
      fscanf(fptr, "%d", &n);
      for (j = 0; j < ndim; j++) {
	  fscanf(fptr, "%f\n",&c);
          x[ndim*(n-1)+j] = c;
          }
      }
/* --- el# n1 n2 ref_pt  mat#  area  Iy  Iz  J  UDLy'  UDLz' --- */
   fgets(dummy,80,fptr);
   for (i = 0; i < ne; i++) {
      fscanf(fptr,"%d", &n);
      for (j = 0; j <= nen; j++) {
         fscanf(fptr,"%d", &k);
	 noc[(nen+1)*(n-1)+j]=k;
	  }
	 fscanf(fptr,"%d", &k);
	 mat[n-1] = k;
      for (j = 0; j < nch; j++) {
	     fscanf(fptr,"%f\n",&c);
         arin[nch*(n-1)+j] = c;
         }
         fscanf(fptr,"%f\n",&c);
	 udl[2*(n-1)] = c;
         fscanf(fptr,"%f\n",&c);
	 udl[2*(n-1)+1] = c;
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
      n = ndn * (abs(noc[3*i] - noc[3*i+1])+1);
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
/* ----- global stiffness matrix ----- */
   for (n = 0; n < ne; n++) {
      printf("forming stiffness matrix of element... %d\n", n+1);
    elstif(2,se,alambda,sep,&el,x,noc,mat,pm,arin,nch,n);
	printf (".... placing in global locations\n");
	for (ii = 0; ii < nen; ii++) {
	   nrt = ndn * (noc[3*n + ii] - 1);
	   for (it = 0; it < ndn; it++) {
	      nr = nrt + it;
	      i = ndn * ii + it;
	      for (jj = 0; jj < nen; jj++) {
		 nct = ndn * (noc[3*n+jj] - 1);
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
   /* ----- loads due to uniformly distributed load on element ----- */
     for (n = 0; n < ne; n++) {
	if (fabs(udl[2*n]) > 0 || fabs(udl[2*n+1]) > 0) {
	   elstif(1,se,alambda,sep,&el,x,noc,mat,pm,arin,nch,n);
	   ed[0] = 0;
	   ed[3] = 0;
	   ed[6] = 0;
	   ed[9] = 0;
	   ed[1] = udl[2*n] * el / 2;
	   ed[7] = ed[1];
	   ed[5] = udl[2*n] * el * el / 12;
	   ed[11] = -ed[5];
	   ed[2] = udl[2*n+1] * el * el / 12;
	   ed[8] = ed[2];
	   ed[4] = -udl[2*n+1] * el * el / 12;
	   ed[10] = -ed[4];
	   for (i = 0; i < 12; i++) {
	      edp[i] = 0;
	      for (k = 0; k < 12; k++) {
		 edp[i] = edp[i] + alambda[k][i] * ed[k];
		 }
	       }
	   i1 = 6*(noc[3*n]-1);
	   i2 = 6*(noc[3*n+1]-1);
	    for (i = 0; i < 6; i++) {
	      f[i1 + i] = f[i1 + i] + edp[i];
	      f[i2 + i] = f[i2 + i] + edp[i + 6];
	      }
           }
	}
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
   printf("node#    x-displ.    y-displ.    z-displ.");
   printf("   xrot.       y-rot.     z-rot.\n");
   fprintf(fptr, "node#    x-displ.    y-displ.   z-displ.");
   fprintf(fptr, "   xrot.       y-rot.     z-rot.\n");
   for (i = 0; i < nn; i++) {
    i1 = 6*i;
    printf(" %3d  %10.3e  %10.3e  %10.3e" ,i+1,f[i1],f[i1+1],f[i1+2]);
    printf(" %10.3e  %10.3e  %10.3e\n", f[i1+3],f[i1+4],f[i1+5]);
    fprintf(fptr," %3d  %10.3e  %10.3e  %11.4e",i+1,f[i1],f[i1+1],f[i1+2]);
    fprintf(fptr," %10.3e  %10.3e  %11.4e\n", f[i1+3],f[i1+4],f[i1+5]);
   }
/* ----- reaction calculation ----- */
   printf(" dof#   reaction\n");
   fprintf(fptr, "dof#    reaction\n");
   for (i = 0; i < nd; i++) {
      k = nu[i];
      reaction = cnst * (u[i] - f[k-1]);
      printf(" %3d   %11.4e\n", k, reaction);
      fprintf(fptr, " %3d   %11.4e\n", k, reaction);
   }
/* ---- member end-actions ----- */
     fprintf(fptr, " member end-forces\n");
     for (n = 0; n < ne; n++) {
       elstif(1,se,alambda,sep,&el,x,noc,mat,pm,arin,nch,n);
       i1 = 6*(noc[3*n]-1);
	   i2 = 6*(noc[3*n+1]-1);
       for (i = 0; i < 6; i++) {
	       ed[i] = f[i1 + i];
	       ed[i + 6] = f[i2 + i];
           }
       for (i = 0; i < 12; i++) {
           edp[i] = 0;
           for (k = 0; k < 12; k++) {
              edp[i] = edp[i] + alambda[i][k] * ed[k];
              }
           }
    /* --- end forces due to distributed loads --- */
       if (fabs(udl[2*n]) > 0 || fabs(udl[2*n+1]) > 0) {
	   ed[0] = 0;
	   ed[3] = 0;
	   ed[6] = 0;
	   ed[9] = 0;
	   ed[1] = -udl[2*n] * el / 2;
	   ed[7] = ed[1];
	   ed[5] = -udl[2*n] * el * el / 12;
	   ed[11] = -ed[5];
	   ed[2] = -udl[2*n+1] * el * el / 12;
	   ed[8] = -ed[2];
	   ed[4] = -udl[2*n+1] * el * el / 12;
	   ed[10] = -ed[4];
	   }
	else {
	   for (k = 0; k < 12; k++) {
		  ed[k] = 0;
		  }
	   }
       for (i = 0; i < 12; i++) {
	   ef[i] = ed[i];
	   for (k = 0; k < 12; k++) {
              ef[i] = ef[i] + sep[i][k] * edp[k];
              }
           }
	fprintf(fptr, " member # %d\n", n+1);
        for (i = 0; i < 2; i++) {
         ii = 6*i;
	 fprintf(fptr, "%11.4e %11.4e %11.4e ",ef[ii],ef[ii+1],ef[ii+2]);
	 fprintf(fptr, "%11.4e %11.4e %11.4e\n",ef[ii+3],ef[ii+4],ef[ii+5]);
	 }
     }
   fclose (fptr);
   printf("\n output is in file %s \n",file2);
   return(0);
}
/* ----- element stiffness ----- */
elstif(istf,se,alambda,sep,el1,x,noc,mat,pm,arin,nch,n)
int *noc,*mat,istf,nch,n;
float *pm,*x,*arin,se[][12],sep[][12],alambda[][12],*el1;
{
  int i1,i2,i3,m,i,j,k,ik;
  float x21,y21,z21,eal,eiyl,eizl,gjl,el,dcos[3][3];
  float c1,c2,c3,eip1,eip2,eip3,cc;
     i1 = 3*(noc[3*n]-1);
     i2 = 3*(noc[3*n+1]-1);
     i3 = 3*(noc[3*n+2]-1);
     m = mat[n]-1;
     x21 = x[i2] - x[i1];
     y21 = x[i2+1] - x[i1+1];
     z21 = x[i2+2] - x[i1+2];
     el = sqrt((double) (x21 * x21 + y21 * y21 + z21 * z21));
     *el1 = el;
     eal = pm[2*m] * arin[nch*n] / el;
     eiyl = pm[2*m] * arin[nch*n+1] / el;
     eizl = pm[2*m] * arin[nch*n+2] / el;
     gjl = pm[2*m+1] * arin[nch*n+3] / el;
     for (i = 0; i < 12; i++) {
         for (j = 0; j < 12; j++) {
            sep[i][j] = 0;
            }
	}     
	 sep[0][0] = eal;
	 sep[0][6] = -eal;
	 sep[6][6] = eal;
	 sep[3][3] = gjl;
	 sep[3][9] = -gjl;
	 sep[9][9] = gjl;
	 sep[1][1] = 12 * eizl / el / el;
	 sep[1][5] = 6 * eizl / el;
	 sep[1][7] = -sep[1][1];
	 sep[1][11] = sep[1][5];
	 sep[2][2] = 12 * eiyl / el / el;
	 sep[2][4] = -6 * eiyl / el;
	 sep[2][8] = -sep[2][2];
	 sep[2][10] = sep[2][4];
	 sep[4][4] = 4 * eiyl;
	 sep[4][8] = 6 * eiyl / el;
	 sep[4][10] = 2 * eiyl;
	 sep[5][5] = 4 * eizl;
	 sep[5][7] = -6 * eizl / el;
	 sep[5][11] = 2 * eizl;
	 sep[7][7] = 12 * eizl / el / el ;
	 sep[7][11] = -6 * eizl / el;
	 sep[8][8] = 12 * eiyl / el / el;
	 sep[8][10] = 6 * eiyl / el;
	 sep[10][10] = 4 * eiyl;
	 sep[11][11] = 4 * eizl;
    for (i = 0; i < 12; i++) {
	for (j = i; j < 12; j++) {
           sep[j][i] = sep[i][j];
           }
	    }
   /* --- convert element stiffness matrix to global system --- */
     dcos[0][0] = x21 / el;
     dcos[0][1] = y21 / el;
     dcos[0][2] = z21 / el;
     eip1 = x[i3] - x[i1];
     eip2 = x[i3 + 1] - x[i1 + 1];
     eip3 = x[i3 + 2] - x[i1 + 2];
     c1 = dcos[0][1] * eip3 - dcos[0][2] * eip2;
     c2 = dcos[0][2] * eip1 - dcos[0][0] * eip3;
     c3 = dcos[0][0] * eip2 - dcos[0][1] * eip1;
     cc = sqrt((double) c1 * c1 + c2 * c2 + c3 * c3);
     dcos[2][0] = c1 / cc;
     dcos[2][1] = c2 / cc;
     dcos[2][2] = c3 / cc;
     dcos[1][0] = dcos[2][1] * dcos[0][2] - dcos[0][1] * dcos[2][2];
     dcos[1][1] = dcos[0][0] * dcos[2][2] - dcos[2][0] * dcos[0][2];
     dcos[1][2] = dcos[2][0] * dcos[0][1] - dcos[0][0] * dcos[2][1];     
     for (i = 0; i < 12; i++) {
        for (j = 0; j < 12; j++) {
           alambda[i][j] = 0;
           }
	 }
     for (k = 0; k < 4; k++) {
	ik = 3 * k;
	for (i = 0; i < 3; i++) {
	   for (j = 0; j < 3; j++) {
	      alambda[i + ik][j + ik] = dcos[i][j];
	      }
	       }
        }
     if (istf == 1)
        return(0);
     for (i = 0; i < 12; i++) {
        for (j = 0; j < 12; j++) {
           se[i][j] = 0;
           for (k = 0; k < 12; k++) {
              se[i][j] = se[i][j] + sep[i][k] * alambda[k][j];
              }
           }
	    }
     for (i = 0; i < 12; i++) {
	    for (j = 0; j < 12; j++) {
	       sep[i][j] = se[i][j];
	       }
	    }
     for (i = 0; i < 12; i++) {
        for (j = 0; j < 12; j++) {
           se[i][j] = 0;
           for (k = 0; k < 12; k++) {
              se[i][j] = se[i][j] + alambda[k][i] * sep[k][j];
              }
           }
	    }
     return(0);
}
/* ----- band solver ----- */
bansol(s,f,nq,nbw)
  int nq,nbw;
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