    /*---------   program gauss  ---------*/
    /*      gauss elimination method      */
    /*          general matrix            */
    /* t.r.chandrupatla and a.d.belegundu */
    /*------------------------------------*/
    #include <stdio.h>
    #include <math.h>
    main()
    {
     FILE *fptr;
    /* ----- Gauss Elimination Method to Solve ax = b ----- */
    int n,i,j,k;
    float c, *a, *b;
    char dummy[81], title[81], file1[81], file2[81];
    puts("Input file name < dr:fn.ext >: ");
    gets(file1);
    puts("Output file name < dr:fn.ext >: ");
    gets(file2);
    fptr = fopen(file1, "r");
    fgets(title,80,fptr);
    fgets(dummy,80,fptr);
    fscanf(fptr,"%d\n", &n);
 /* ----- memory allocation ----- */
    a = (float *) calloc(n*n, sizeof(float));
    b = (float *) calloc(n, sizeof(float));
/* ----- Read a[][] matrix ----- */
   fgets(dummy,80,fptr);
   for (i = 0; i < n; i++) {
     for (j = 0; j < n-1; j++) {
      fscanf(fptr, "%f",&c);
      a[n*i+j] = c;
     }
     fscanf(fptr, "%f\n",&c);
     a[n*i+n-1]=c;
   }
   fgets(dummy,80,fptr);
   for (j = 0; j < n-1; j++) {
     fscanf(fptr, "%f",&c);
     b[j] = c;
   }
    fscanf(fptr, "%f\n",&c);
    b[n-1]=c;
    fclose(fptr);
/* ----- forward elimination ----- */
     for (k = 0; k < n - 1; k++) {
	 for (i = k + 1; i < n; i++) {
	     c = a[n*i+k] / a[n*k+k];
	     for (j = k + 1; j < n; j++) {
		 a[n*i+j] = a[n*i+j] - c * a[n*k+j];
		 }
	     b[i] = b[i] - c * b[k];
	     }
	  }
/* ----- back-substitution ----- */
     b[n-1] = b[n-1] / a[n*(n-1)+n-1];
     for (i = n - 2; i >= 0; i--) {
	 c = 1 / a[n*i+i];
	 b[i] = c * b[i];
	 for (k = i + 1; k < n; k++) {
	     b[i] = b[i] - c * a[n*i+k] * b[k];
	     }
	 }
   fptr = fopen(file2, "w");
   printf("\n%s\n", title);
   fprintf(fptr, "\n%s\n", title);

     printf("\nx() =\n");
     fprintf(fptr,"\nx() =\n");
     for (i = 0; i < n; i++) {
	 printf("%7.4f ", b[i]);
	 fprintf(fptr, "%7.4f ", b[i]);
	 }
     printf("\n");
     fclose(fptr);
     return(0);
  }
