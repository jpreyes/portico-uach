/*==============================*/
/*          frame2d.c           */
/*  2-d  frame analysis by fem  */
/*   chandrupatla & belegundu   */
/* =============================*/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr;
   int n,nq,i,j,k,m,i1,i2,ii,nrt,it,nr,jj,nct,jt,nc;
   char dummy[81], title[81], file1[81], file2[81];
   int ne,nn,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc,nch;
   int *noc,*nu,*mat,*mpc;
   float *x,*arin,*pm,*u,*s,*f,*beta,*udl,se[6][6],ed[6];
   float edp[6],alambda[6][6],sep[6][6],ef[6];
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
   fscanf(fptr,"%d %d %d %d %d %d\n", &nn, &ne, &nm, &ndim, &nen, &ndn);
   fgets(dummy, 80, fptr);
   fscanf(fptr,"%d %d %d %d %d\n", &nd, &nl, &nmpc);
   npr = 1;  /* Material property E */
   nch = 2;  /* Area, Inertia of Element */
   /* ----- total dof is nq ----- */
   nq = ndn*nn;
/* --------------- memory allocation ---------------- */
   x = (float *) calloc(nn*ndim, sizeof(float));
   noc = (int *) calloc(ne*nen, sizeof(int));
   u = (float *) calloc(nd, sizeof(float));
   nu = (int *) calloc(nd, sizeof(int));
   mat = (int *) calloc(ne,sizeof(int));
   arin = (float *) calloc(nch*ne, sizeof(float));
   f = (float *) calloc(nn*ndn, sizeof(float));
   pm = (float *) calloc(nm*npr, sizeof(float));
   mpc = (int *) calloc(2*nmpc, sizeof(int));
   beta = (float *) calloc(3*nmpc, sizeof(float));
   udl = (float *) calloc(ne, sizeof(float));
/* ------------------ read data --------------------- */
/* ----- coordinates ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < nn; i++) {
      fscanf(fptr, "%d", &n);
      for (j = 0; j < ndim; j++) {
	  fscanf(fptr, "%f\n",&c);
          x[ndim*(n-1)+j] = c;
          }
      }
/* --- connectivity, material, area, mom_inertia, distr_load --- */
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
         arin[2*(n-1)] = c;
         fscanf(fptr,"%f\n",&c);
         arin[2*(n-1)+1] = c;
	     fscanf(fptr,"%f\n",&c);
         udl[n-1] = c;
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
/* ----- global stiffness matrix ----- */
   for (n = 0; n < ne; n++) {
      printf("forming stiffness matrix of element... %d\n", n+1);
    elstif(2,se,alambda,sep,&el,x,noc,mat,pm,arin,nch,n);
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
   /* ----- loads due to uniformly distributed load on element ----- */
     for (n = 0; n < ne; n++) {
	if (fabs(udl[n]) > 0) {
	   elstif(1,se,alambda,sep,&el,x,noc,mat,pm,arin,nch,n);
	   ed[0] = 0;
	   ed[3] = 0;
	   ed[1] = udl[n] * el / 2;
	   ed[4] = ed[1];
	   ed[2] = udl[n] * el * el / 12;
	   ed[5] = -ed[2];
	   for (i = 0; i < 6; i++) {
	      edp[i] = 0;
	      for (k = 0; k < 6; k++) {
		 edp[i] = edp[i] + alambda[k][i] * ed[k];
		 }
	       }
	   i1 = 3*(noc[2*n]-1);
	   i2 = 3*(noc[2*n+1]-1);
	    for (i = 0; i < 3; i++) {
	      f[i1 + i] = f[i1 + i] + edp[i];
	      f[i2 + i] = f[i2 + i] + edp[i + 3];
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
   printf("node#    x-displ.     y-displ.     rotation\n");
   fprintf(fptr, "node#    x-displ.     y-displ.     rotation\n");
   for (i = 0; i < nn; i++) {
    i1 = 3*i;
    printf(" %3d  %11.4e  %11.4e  %11.4e\n",i+1,f[i1],f[i1+1],f[i1+2]);
    fprintf(fptr," %3d  %11.4e  %11.4e  %11.4e\n",i+1,f[i1],f[i1+1],f[i1+2]);
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
       i1 = 3*(noc[2*n]-1);
	   i2 = 3*(noc[2*n+1]-1);
       for (i = 0; i < 3; i++) {
	   ed[i] = f[i1 + i];
	       ed[i + 3] = f[i2 + i];
           }
       for (i = 0; i < 6; i++) {
           edp[i] = 0;
           for (k = 0; k < 6; k++) {
              edp[i] = edp[i] + alambda[i][k] * ed[k];
              }
           }
    /* --- end forces due to distributed loads --- */
       if (fabs(udl[n]) > 0) {
           ed[0] = 0;
	   ed[3] = 0;
	   ed[1] = -udl[n] * el / 2;
	   ed[4] = ed[1];
	   ed[2] = -udl[n] * el * el / 12;
	   ed[5] = -ed[2];
	   }
	else {
	   for (k = 0; k < 6; k++) {
		  ed[k] = 0;
		  }
	   }
       for (i = 0; i < 6; i++) {
	   ef[i] = ed[i];
	   for (k = 0; k < 6; k++) {
              ef[i] = ef[i] + sep[i][k] * edp[k];
              }
           }
	fprintf(fptr, " member # %d\n", n+1);
        for (i = 0; i < 2; i++) {
         ii = 3*i;
	 fprintf(fptr, "%11.4e %11.4e %11.4e\n",ef[ii],ef[ii+1],ef[ii+2]);
	 }
     }
   fclose (fptr);
   printf("\n output is in file %s \n",file2);
   return(0);
}
/* ----- element stiffness ----- */
elstif(istf,se,alambda,sep,el1,x,noc,mat,pm,arin,nch,n)
int *noc,*mat,istf,nch,n;
float *pm,*x,*arin,se[][6],sep[][6],alambda[][6],*el1;
{
  int i1,i2,m,i,j,k,ik;
  float x21,y21,eal,eizl,el,dcos[3][3];
     i1 = 2*(noc[2*n]-1);
     i2 = 2*(noc[2*n+1]-1);
     m = mat[n]-1;
     x21 = x[i2] - x[i1];
     y21 = x[i2+1] - x[i1+1];
     el = sqrt((double) (x21 * x21 + y21 * y21));
     *el1 = el;
     eal = pm[m] * arin[nch*n] / el;
     eizl = pm[m] * arin[nch*n+1] / el;
     for (i = 0; i < 6; i++) {
         for (j = 0; j < 6; j++) {
            sep[i][j] = 0;
            }
	}
     sep[0][0] = eal;
     sep[0][3] = -eal;
     sep[3][3] = eal;
     sep[1][1] = 12 * eizl / el / el;
     sep[1][2] = 6 * eizl / el;
     sep[1][4] = -sep[1][1];
     sep[1][5] = sep[1][2];
     sep[2][2] = 4 * eizl;
     sep[2][4] = -6 * eizl / el;
     sep[2][5] = 2 * eizl;
     sep[4][4] = 12 * eizl / el / el;
     sep[4][5] = -6 * eizl / el;
     sep[5][5] = 4 * eizl;
     for (i = 0; i < 6; i++) {
	for (j = i; j < 6; j++) {
           sep[j][i] = sep[i][j];
           }
	    }
   /* --- convert element stiffness matrix to global system --- */
     dcos[0][0] = x21 / el;
     dcos[0][1] = y21 / el;
     dcos[0][2] = 0;
     dcos[1][0] = -dcos[0][1];
     dcos[1][1] = dcos[0][0];
     dcos[1][2] = 0;
     dcos[2][0] = 0;
     dcos[2][1] = 0;
     dcos[2][2] = 1;
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
           alambda[i][j] = 0;
           }
	 }
     for (k = 0; k < 2; k++) {
	ik = 3 * k;
	for (i = 0; i < 3; i++) {
	   for (j = 0; j < 3; j++) {
	      alambda[i + ik][j + ik] = dcos[i][j];
	      }
	       }
        }
     if (istf == 1)
        return(0);
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
           se[i][j] = 0;
           for (k = 0; k < 6; k++) {
              se[i][j] = se[i][j] + sep[i][k] * alambda[k][j];
              }
           }
	    }
     for (i = 0; i < 6; i++) {
	    for (j = 0; j < 6; j++) {
	       sep[i][j] = se[i][j];
	       }
	    }
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
           se[i][j] = 0;
           for (k = 0; k < 6; k++) {
              se[i][j] = se[i][j] + alambda[k][i] * sep[k][j];
              }
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