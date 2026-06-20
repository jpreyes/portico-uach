/*****       program jacobi        *****/
/*    generalized jacobi's method      */
/*       for symmetric matrices        */
/*  t.r.chandrupatla and a.d.belegundu */
/***************************************/
#include <stdio.h>
#include <math.h>
main()
{
   FILE *fptr1, *fptr2;
   float *s,*gm,*evl,*evc,c,c1,c2,tols,tolm,aa,bb,cc;
   float cab,bet,alp,sqc,si,sj,emi,emj,sii,sij,sjj,eii,eij,ejj;
   float evi,evj,gm2,omega,freq;
   int *nord,nq,nbw,i,j,k,n,jn,nswmax,k1,i1,nsw,ii,jj,j1,ij,ifl;
   float pi = 3.14159, tol = 0.000001;  
   char dummy[81], file1[81], file2[81];
   printf("\n");
   puts("Input file name < dr:fn.ext >: ");
   gets(file1);
   puts("Output file name < dr:fn.ext >: ");
   gets(file2);
   fptr1 = fopen(file1, "r");
   fptr2 = fopen(file2, "w");
   fgets(dummy,80,fptr1);
   fgets(dummy,80,fptr1);
/* --- read in number of equations and bandwidth --- */
   fscanf(fptr1,"%d %d\n", &nq, &nbw);
   printf("default tolerance is  1e-6\n");      
/* -----         memory allocation         ----- */
   s = (float *) calloc(nq*nq, sizeof(float));
   gm = (float *) calloc(nq*nq, sizeof(float));
   evl = (float *) calloc(nq, sizeof(float));
   evc = (float *) calloc(nq*nq, sizeof(float));
   nord = (int *) calloc(nq, sizeof(int));
/* --------------------------------------------- */

   for (i = 0; i < nq; i++) {
           nord[i] = i;
          }
/* ----- banded stiffness matrix into s(nq,nq) ----- */
   fgets(dummy,80,fptr1);
   for (i = 1; i <= nq; i++) {
         for (jn = 1; jn <= nbw; jn++) {
         fscanf(fptr1, "%f\n", &c);
            j = i + jn - 1;
        if (j <= nq) {
           s[nq*(i-1)+j-1] = c;
               s[nq*(j-1)+i-1] = c;
           }
         }
           }
/* ----- banded mass matrix into gm(nq,nq) ----- */
   fgets(dummy,80,fptr1);
   for (i = 1; i <= nq; i++) {
         for (jn = 1; jn <= nbw; jn++) {
         fscanf(fptr1, "%f\n", &c);
            j = i + jn - 1;
        if (j <= nq) {
           gm[nq*(i-1)+j-1] = c;
               gm[nq*(j-1)+i-1] = c;
           }
         }
           }
     fclose(fptr1);     
  /* ----- initialize eigenvector matrix ----- */
     for (i = 0; i < nq; i++) {
             for (j = 0; j < nq; j++) {
                evc[nq*i+j] = 0;
           }
             evc[nq*i+i] = 1;
            }
         c1 = s[0];
         c2 = gm[0];
     for (i = 0; i < nq; i++) {
        if (c1 > s[nq*i+i])
            c1 = s[nq*i+i];
        if (c2 < gm[nq*i+i])
            c2 = gm[nq*i+i];
        }
     tols = tol * c1;
         tolm = tol * c2;
         nswmax = 50;
     printf("maximum number of sweeps nswmax = 50  ");
     k1 = 1;
         i1 = 1;
         nsw = 0;
  do {
     nsw = nsw + 1;
     if (nsw > nswmax) {
        printf ("no convergence after %d sweeps.", nswmax);
        fprintf (fptr2, "no convergence after %d sweeps.", nswmax);
        fclose(fptr2);
        exit(1);
        }
     printf( "------  sweep number   %d\n", nsw);
     for (k = k1; k < nq; k++) {
         for (i = i1; i <= k; i++) {
         ii = i - 1;
         j = nq - k + i;
         jj = nq - k + i - 1;
       if (fabs(s[nq*ii+jj]) > tols || fabs(gm[nq*ii+jj]) > tolm) {
         aa = s[nq*ii+ii] * gm[nq*ii+jj] - gm[nq*ii+ii] * s[nq*ii+jj];
         bb = s[nq*jj+jj] * gm[nq*ii+jj] - gm[nq*jj+jj] * s[nq*ii+jj];
         cc = s[nq*ii+ii] * gm[nq*jj+jj] - gm[nq*ii+ii] * s[nq*jj+jj];
         cab = .25 * cc * cc + aa * bb;
     if (cab < 0) {
         printf("square root of negative term -- check matrices");
         fprintf(fptr2, "square root of negative term -- check matrices");
         fclose(fptr2);
         exit(2);
        }
     if (aa == 0) {
        bet = 0;
            alp = -s[nq*ii+jj] / s[nq*ii+ii];
        } else if (bb == 0) {
            alp = 0;
                bet = -s[nq*ii+jj] / s[nq*jj+jj];
            } else {
               sqc = sqrt((double) cab);
                   if (cc < 0)
                      sqc = -sqc;
               alp = (-.5 * cc + sqc) / aa;
               bet = -aa * alp / bb;
               }
   /* ----- only upper triangular part is used in diagonalization ----- */
     if (i > 1) {
        for (n = 1; n < i; n++) {
            si = s[nq*(n-1)+ii];
            sj = s[nq*(n-1)+jj];
            emi = gm[nq*(n-1)+ii];
            emj = gm[nq*(n-1)+jj];
            s[nq*(n-1)+ii] = si + bet * sj;
            s[nq*(n-1)+jj] = sj + alp * si;
            gm[nq*(n-1)+ii] = emi + bet * emj;
            gm[nq*(n-1)+jj] = emj + alp * emi;
           }
        }
     if (j < nq) {
        for (n = j + 1; n <= nq; n++) {
           si = s[nq*ii+n-1];
           sj = s[nq*jj+n-1];
           emi = gm[nq*ii+n-1];
           emj = gm[nq*jj+n-1];
           s[nq*ii+n-1] = si + bet * sj;
           s[nq*jj+n-1] = sj + alp * si;
           gm[nq*ii+n-1] = emi + bet * emj;
           gm[nq*jj+n-1] = emj + alp * emi;
           }
        }
     if (i < j) {
        for (n = i + 1; n < j; n++) {
           si = s[nq*ii+n-1];
           sj = s[nq*(n-1)+jj];
           emi = gm[nq*ii+n-1];
           emj = gm[nq*(n-1)+jj];
           s[nq*ii+n-1] = si + bet * sj;
           s[nq*(n-1)+jj] = sj + alp * si;
           gm[nq*ii+n-1] = emi + bet * emj;
           gm[nq*(n-1)+jj] = emj + alp * emi;
           }
         }
     sii = s[nq*ii+ii];
     sij = s[nq*ii+jj];
     sjj = s[nq*jj+jj];
     s[nq*ii+jj] = 0;
     s[nq*ii+ii] = sii + 2 * bet * sij + bet * bet * sjj;
     s[nq*jj+jj] = sjj + 2 * alp * sij + alp * alp * sii;
     eii = gm[nq*ii+ii];
     eij = gm[nq*ii+jj];
     ejj = gm[nq*jj+jj];
     gm[nq*ii+jj] = 0;
     gm[nq*ii+ii] = eii + 2 * bet * eij + bet * bet * ejj;
     gm[nq*jj+jj] = ejj + 2 * alp * eij + alp * alp * eii;
 /* ----- eigenvectors ----- */
     for (n = 0; n < nq; n++) {
         evi = evc[nq*n+ii];
         evj = evc[nq*n+jj];
         evc[nq*n+ii] = evi + bet * evj;
         evc[nq*n+jj] = evj + alp * evi;
        }
             }
          }
         }
     for (k = 1; k < nq; k++) {
         for (i = 1; i <= k; i++) {
             ii = i - 1;
             j = nq - k + i - 1;
         ifl = 0;
        if (fabs(s[nq*ii+j]) > tols || fabs(gm[nq*ii+j]) > tolm) {
           k1 = k;
           i1 = i;
           ifl = 1;
           }
           if (ifl == 1)
             break;
             }
           if (ifl == 1)
             break;
         }
     } while (ifl == 1);
  /* -----  calculation of eigenvalues ----- */
     for (i = 0; i < nq; i++) {
        if (fabs(gm[nq*i+i]) < tolm)
           gm[nq*i+i] = tolm;
        evl[i] = s[nq*i+i] / gm[nq*i+i];
        }
  /* ----- scaling of eigenvectors ----- */
     for (i = 0; i < nq; i++) {
        gm2 = sqrt((double) fabs(gm[nq*i+i]));
        for (j = 0; j < nq; j++) {
           evc[nq*j+i] = evc[nq*j+i] / gm2;
           }
        }
   /* -----   results   ----- */
   /* --- ascending order of eigenvalues --- */
     for (i = 0; i < nq-1; i++) {
        ii = nord[i];
        i1 = ii;
        c1 = evl[ii];
            j1 = i;
        for (j = i+1; j < nq; j++) {
           ij = nord[j];
           if (c1 > evl[ij]) {
              c1 = evl[ij];
                  i1 = ij;
                  j1 = j;
             }
           }
           nord[i] = i1;
           nord[j1] = ii;
        }
 printf("eigenvalues and eigenvectors for data in file %s\n)", file1);
 fprintf(fptr2, "eigenvalues and eigenvectors for data in file %s\n", file1);
     for (i = 0; i < nq; i++) {
        ii = nord[i];
        printf( "eigenvalue number  %d\n ", i+1);
        fprintf(fptr2, "eigenvalue number  %d\n ", i+1);
        c = evl[ii];
            omega = sqrt((double) c);
            freq = .5 * omega / pi;
        printf("eigenvalue = %11.4e",  c);
        fprintf(fptr2, "eigenvalue = %11.4e",  c);
        printf("   omega = %10.3e", omega);
        fprintf(fptr2, "   omega = %10.3e", omega);
        printf("   freq = %10.3e hz\n", freq);
        fprintf(fptr2, "   freq = %10.3e hz\n", freq);
        printf( "eigenvector \n");
        fprintf(fptr2, "eigenvector \n");
        ifl = 0;
        for (j = 0; j < nq; j++) {
           printf( "%10.3e ", evc[nq*j+ii]);
           fprintf(fptr2, "%10.3e ", evc[nq*j+ii]);
           ifl = ifl + 1;
           if (ifl > 5)
              ifl = 0;
           if (ifl == 0) {
             printf( "\n");
                 fprintf(fptr2, "\n");
                 }
           }
        printf("\n");
            fprintf(fptr2, "\n");
        }
        fclose(fptr2);
     return(0);
}

