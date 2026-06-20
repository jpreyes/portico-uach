/********************************************/
/*             program  tetra3d             */
/*     3-d  stress analysis using the       */
/*           tetrahedral element            *
/*   t.r.chandrupatla and a.d.belegundu     */
/********************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1;
   int n,i,j,k,m,ii,jj,m1,nmin,nmax,nrt,nct,it,jt,i1,i2;
   int nr,nc;
   char dummy[81], title[81], file1[81], file2[81], file3[81];
   int ne,nn,nq,nm,nd,nl,nen,ndn,ndim,npr,nbw,nmpc;
   int *noc, *nu, *mat, *mpc;
   float *x, *pm, *u, *tempr, *s, *f, *beta;
   float c,dj,al,pnu,tld,cnst,reaction,s1,s2,s3,pi;
   float b[6][12], d[6][6], db[6][12], se[12][12], str[6], tl[6];
   float ai1,ai21,ai22,ai2,ai31,ai32,ai3,c1,c2,c3,th,th2,p1,p2,p3;
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
   npr = 3;    /* material properties E, Nu, Alpha */
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
        for (j = 1; j < nen;j++) {
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
     for (n = 0; n < ne; n++) {
	printf("forming stiffness matrix of element %d\n", n+1);
	dbmat(n,b,d,db,mat,pm,npr,x,noc,&dj,&al,&pnu);
     /* --- element stiffness --- */
        for (i = 0; i < 12; i++) {
            for (j = 0; j < 12; j++) {
                c = 0.;
                for (k = 0; k < 6; k++) {
                    c = c + fabs(dj) * b[k][i] * db[k][j] / 6;
                    }
		se[i][j] = c;
		}
	    }

     /* --- temperature load vector --- */
	c = al * tempr[n];
	for (i = 0; i < 12; i++) {
	    tl[i] = c * fabs(dj) * (db[0][i] + db[1][i] + db[2][i]) / 6;
	    }
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
     pi = 3.141593;
     for (n = 0; n < ne; n++) {
	dbmat(n,b,d,db,mat,pm,npr,x,noc,&dj,&al,&pnu);
	stress(f,al,noc,tempr,n,d,db,str);
/* ----- principal stress calculations ----- */
	
    ai1 = str[0] + str[1] + str[2];
    ai21 = str[0] * str[1] + str[1] * str[2] + str[2] * str[0];
    ai22 = str[3] * str[3] + str[4] * str[4] + str[5] * str[5];
    ai2 = ai21 - ai22;
    ai31 = str[0] * str[1] * str[2] + 2 * str[3] * str[4] * str[5];
    ai32 = str[0]*str[3]*str[3]+str[1]*str[4]*str[4]+str[2]*str[5]*str[5];
    ai3 = ai31 - ai32;
    c1 = ai2 - ai1 * ai1 / 3;
    c2 = -2 * (ai1*ai1*ai1) / 27 + ai1 * ai2 / 3 - ai3;
    c3 = 2 * sqrt((double) -c1 / 3);
    th = -3 * c2 / (c1 * c3);
    th2 = sqrt((double) fabs(1 - th * th));
    if (th == 0)
       th = pi / 2;
    if (th > 0)
       th = atan((double) th2 / th);
    if (th < 0)
       th = pi - atan((double) th2 / th);
    th = th / 3;
   /* --- principal stresses --- */
    p1 = ai1 / 3 + c3 * cos(th);
    p2 = ai1 / 3 + c3 * cos(th + 2 * pi / 3);
    p3 = ai1 / 3 + c3 * cos(th + 4 * pi / 3);
    fprintf( fptr1, "stresses in element no. %4d\n", n+1);
    fprintf( fptr1, "  normal stresses sx,sy,sz\n");
    fprintf( fptr1, "  %11.4e  %11.4e  %11.4e\n", str[0],str[1],str[2]);
    fprintf( fptr1, "  shear stresses tyz,txz,txy\n");
    fprintf( fptr1, "  %11.4e  %11.4e  %11.4e\n", str[3],str[4],str[5]);
    fprintf( fptr1, "  principal stresses\n");
    fprintf( fptr1, "  %11.4e  %11.4e  %11.4e\n", p1,p2,p3);
     }
     fclose(fptr1);
     printf( "complete results are in file %s\n", file2);
     return(0);
}
/* ----- db[] matrix ----- */
dbmat(n,b,d,db,mat,pm,npr,x,noc,djj,al1,pnu1)
    int n,npr;
    float *djj,*pnu1,*al1;
    int *mat,*noc;
    float *x,*pm;
    float b[][12],db[][12],d[][6];
{
    int i,j,k,m,i1,i2,i3,i4;
    float e,c,c1,c2,c3,c4,dj,pnu,al,dj1,dj2,dj3;
    float x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    float x14,x24,x34,y14,y24,y34,z14,z24,z34;
    float a11,a21,a31,a12,a22,a32,a13,a23,a33;
/* ----- d(), b() and db() matrices ----- */
  /* --- first the d-matrix --- */
     m = mat[n]-1;
     e = pm[npr*m];
     pnu = pm[npr*m+1];
     *pnu1 = pnu;
     al = pm[npr*m+2];
     *al1 = al;
  /* --- d() matrix --- */
     c4 = e / ((1 + pnu) * (1 - 2 * pnu));
     c1 = c4 * (1 - pnu);
     c2 = c4 * pnu;
     c3 = .5 * e / (1 + pnu);
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
           d[i][j] = 0;
           }
        }
     d[0][0] = c1;
     d[0][1] = c2;
     d[0][2] = c2;
     d[1][0] = c2;
     d[1][1] = c1;
     d[1][2] = c2;
     d[2][0] = c2;
     d[2][1] = c2;
     d[2][2] = c1;
     d[3][3] = c3;
     d[4][4] = c3;
     d[5][5] = c3;
/* --- strain-displacement matrix b() --- */
     i1 = noc[4*n]-1;
     i2 = noc[4*n+1]-1;
     i3 = noc[4*n+2]-1;
     i4 = noc[4*n+3]-1;
     x14 = x[3*i1] - x[3*i4];
     x24 = x[3*i2] - x[3*i4];
     x34 = x[3*i3] - x[3*i4];
     y14 = x[3*i1+1] - x[3*i4+1];
     y24 = x[3*i2+1] - x[3*i4+1];
     y34 = x[3*i3+1] - x[3*i4+1];
     z14 = x[3*i1+2] - x[3*i4+2];
     z24 = x[3*i2+2] - x[3*i4+2];
     z34 = x[3*i3+2] - x[3*i4+2];
     dj1 = x14 * (y24 * z34 - z24 * y34);
     dj2 = y14 * (z24 * x34 - x24 * z34);
     dj3 = z14 * (x24 * y34 - y24 * x34);
     dj = dj1 + dj2 + dj3;          /* dj is determinant of jacobian */
     *djj = dj;
     a11 = (y24 * z34 - z24 * y34) / dj;
     a21 = (z24 * x34 - x24 * z34) / dj;
     a31 = (x24 * y34 - y24 * x34) / dj;
     a12 = (y34 * z14 - z34 * y14) / dj;
     a22 = (z34 * x14 - x34 * z14) / dj;
     a32 = (x34 * y14 - y34 * x14) / dj;
     a13 = (y14 * z24 - z14 * y24) / dj;
     a23 = (z14 * x24 - x14 * z24) / dj;
     a33 = (x14 * y24 - y14 * x24) / dj;     
/* --- definition of b() matrix --- */
     for (i = 0; i < 6; i++) {
        for (j = 0; j < 12; j++) {
           b[i][j] = 0;
           }
        }
     b[0][0] = a11;
	 b[0][3] = a12;
	 b[0][6] = a13;
	 b[0][9] = -a11 - a12 - a13;
     b[1][1] = a21;
	 b[1][4] = a22;
	 b[1][7] = a23;
	 b[1][10] = -a21 - a22 - a23;
     b[2][2] = a31;
	 b[2][5] = a32;
	 b[2][8] = a33;
	 b[2][11] = -a31 - a32 - a33;
     b[3][1] = a31;
	 b[3][2] = a21;
	 b[3][4] = a32;
	 b[3][5] = a22;
	 b[3][7] = a33;
     b[3][8] = a23;
	 b[3][10] = b[2][11];
	 b[3][11] = b[1][10];
     b[4][0] = a31;
	 b[4][2] = a11;
	 b[4][3] = a32;
	 b[4][5] = a12;
	 b[4][6] = a33;
     b[4][8] = a13;
	 b[4][9] = b[2][11];
	 b[4][11] = b[0][9];
     b[5][0] = a21;
	 b[5][1] = a11;
	 b[5][3] = a22;
	 b[5][4] = a12;
	 b[5][6] = a23;
     b[5][7] = a13;
	 b[5][9] = b[1][10];
	 b[5][10] = b[0][9];
/* --- db matrix db = d*b ---*/
     for (i = 0; i < 6; i++) {
         for (j = 0; j < 12; j++) {
             c = 0.;
             for (k = 0; k < 6; k++) {
		 c = c + d[i][k] * b[k][j];
		 }
	     db[i][j] = c;
	     }
	 }
    return(0);
}
/* ----- stress evaluation ----- */
stress(f,al,noc,tempr,n,d,db,str)
  int n, *noc;
  float *f,*tempr,d[][6],db[][12],str[];
  float al;
{  int i,j,in,ii,k;
   float c,c1,q[12];	
 /* --- stress evaluation (element nodal displacements stored in q() --- */
     for (i = 0; i < 4; i++) {
	in = 3 * (noc[4*n + i] - 1);
	    ii = 3 * i;
	for (j = 0; j < 3; j++) {
	   q[ii + j] = f[in + j];
	       }
	    }
     c1 = al * tempr[n-1];
     for (i = 0; i < 6; i++) {
	    c = 0.;
	    for (k = 0; k < 12; k++) {
	        c = c + db[i][k] * q[k];
            }
	     str[i] = c - c1 * (d[i][0] + d[i][1] + d[i][2]);
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