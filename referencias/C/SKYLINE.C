  /********   program skyline  ***********/
  /*    skyline method for equations     */
  /* t.r.chandrupatla and a.d.belegundu  */
  /***************************************/
  #include <stdio.h>
  #include <math.h>
  main()
    {
     FILE *fptr;
     int i,j,k,n,ns,nt,ni,nj,nsum,ii,ij,ki,kj,kk,i1,j1,k1;
     int *id;
     char dummy[81],title[81],file1[81],file2[81];
     float *a,*b,c;
     printf("\n");
     puts("Input file name <fn.ext>:");
     gets(file1);
     puts("Output file name <fn.ext>:");
     gets(file2);
     fptr = fopen(file1, "r");
     fgets(title,80,fptr);
     fgets(dummy,80,fptr);
     fscanf(fptr,"%d\n",&n);
     id = (int *) calloc(n, sizeof(int));
     fgets(dummy,80,fptr);
/* --- read column heights then convert to diagonal pointers ----- */
     fscanf(fptr,"%d\n", &k);
     id[0] = k;
     for (i = 1; i < n; i++) {
	fscanf(fptr,"%d\n", &k);
	id[i] = k;
	id[i] = id[i] + id[i - 1];
	}
     nsum = id[n-1];
     a = (float *) calloc(nsum, sizeof(float));
     b = (float *) calloc(n, sizeof(float));
     nt = 0;
     ni = 1;
     fgets(dummy,80,fptr);
     for (i = 0; i < n; i++) {
	     if (i > 0) {
	        nt = id[i - 1];
	        ni = id[i] - id[i - 1];
	        }
	      for (j = 0; j < ni; j++) {
		  fscanf(fptr, "%f\n",&c);
	          a[nt + j] = c;
	        }
	      }
     fgets(dummy,80,fptr);
     for (i = 0; i < n; i++) {
		fscanf(fptr, "%f\n",&c);
	      b[i] = c;
	     }
     fclose(fptr);
     for (j = 2; j <= n; j++) {
	     nj = id[j-1] - id[j - 2];
	     if (nj > 1) {
	        k1 = 0;
	        nj = j - nj + 1;
	        for (k = nj; k < j; k++) {
	            k1 = k1 + 1;
	            kj = id[j - 2] + k1;
	            kk = id[k-1];
	            c = a[kj-1] / a[kk-1];
	            for (i = k + 1; i <= j; i++) {
		           ni = id[i-1] - id[i - 2];
		           if ((i - k + 1) <= ni) {
		              ij = id[j-1] - j + i;
		              ki = id[i-1] - i + k;
		              a[ij-1] = a[ij-1] - c * a[ki-1];
		              }
		           }
	            }
	        }
	    }
     for (k = 1; k < n; k++) {
	     kk = id[k-1];
	     c = b[k-1] / a[kk-1];
	     for (i = k + 1; i <= n; i++) {
	         ni = id[i-1] - id[i - 2];
	         if ((i - k + 1) <= ni) {
	            ki = id[i-1] - i + k;
	            b[i-1] = b[i-1] - c * a[ki-1];
	            }
	         }
	     }
    /* ----- backsubstitution ----- */
     ns = id[n-1];
     b[n-1] = b[n-1] / a[ns-1];
     for (i1 = 1; i1 < n; i1++) {
	     i = n - i1;
	     ii= id[i-1];
	     c = 1 / a[ii-1];
	     b[i-1] = c * b[i-1];
	     for (j = i + 1; j <= n; j++) {
	         j1 = j - i + 1;
	         nj = id[j-1] - id[j - 2];
	         if (j1 <= nj) {
	            ij = id[j-1] - j + i;
	            b[i-1] = b[i-1] - c * a[ij-1] * b[j-1];
	      }
	   }
	}
     fptr = fopen(file2, "w");
     fprintf(fptr, "\n   solution to eqns in file %s\n", file1);
     printf("Solution Vector\n");
     for (i = 0; i < n; i++) {
	fprintf(fptr, " %10.4e ", b[i]);
        printf("%10.4e ",b[i]);
	}
	fprintf(fptr,"\n");
	printf("\n");
	fclose(fptr);
     return(0);
  }