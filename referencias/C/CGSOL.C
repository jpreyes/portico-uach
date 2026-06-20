   /*-------  program cgsolve  ----------*/
   /*      conjugate gradient solver     */
   /* t.r.chandrupatla and a.d.belegundu */
   /*------------------------------------*/
   #include <stdio.h>
   #include <math.h>
   main()
    {
     FILE *fptr;
    /* ----- Conjugate Gradient Method to Solve ax = b ----- */
    int n,i,j,k,iter;
    float c, *a, *b, *x;
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
    x = (float *) calloc(n, sizeof(float));
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

    cgsolve(a, b, x, n, &iter);
    fptr = fopen(file2, "w");
    printf("\n%s\n", title);
    fprintf(fptr, "\n%s\n", title);
    printf("iterations = %d\n", iter);
    fprintf(fptr,"iterations = %d\n", iter);
    printf("\nx() =\n");
    fprintf(fptr,"\nx() =\n");
    for (i = 0; i < n; i++) {
	 printf("%7.4f ", x[i]);
	 fprintf(fptr, "%7.4f ", x[i]);
	 }
    printf("\n");
    fclose(fptr);
    return(0);
  }
    cgsolve (a, b, x, n, iter1)
      float *a, *b, *x;
      int n, *iter1;
    {
      float *g, *d, *ad, gg1, gg2;
      int iter, i, j;
      float c, al, bt, dad;
      g = (float *)calloc(n, sizeof(float));
      d = (float *)calloc(n, sizeof(float));
      ad = (float *)calloc(n, sizeof(float));      
      for (i = 0; i < n; i++) {
           x[i] = 0.;
           g[i] = -b[i];
           d[i] = b[i];
          }
          gg1 = 0;
        for (i = 0; i < n; i++) {
            gg1 = gg1 + g[i] * g[i];
            }
        iter = 0;
while(gg1 > .000001) {
        iter = iter + 1;
        dad = 0.;
        for (i = 0; i < n; i++) {
	   c = 0;
	   for (j = 0; j < n; j++) {
              c = c + a[n*i + j] * d[j];
              }
           ad[i] = c;
           dad = dad + c * d[i];
           }
        al = gg1 / dad;
        gg2 = 0.;
        for (i = 0; i < n; i++) {
           x[i] = x[i] + al * d[i];
           g[i] = g[i] + al * ad[i];
           gg2 = gg2 + g[i] * g[i];
           }
        bt = gg2 / gg1;
        for (i = 0; i < n; i++) {
           d[i] = -g[i] + bt * d[i];
           }
	gg1 = gg2;
      }
      *iter1 = iter;
      return(0);
  }

