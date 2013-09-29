#include "ed.h"

/********************************************************************************/

void print_dvec(double *vec, int size, char *str) {

   FILE *pf = fopen(str, "a+"); 
   int i; 

   for (i=0; i<size; i++) {
     fprintf(pf, "%f\n", vec[i]); 
   }
   fprintf(pf, "\n\n"); 
   fflush(pf); 

fclose(pf); 
return; 
}

/********************************************************************************/

void print_ivec(int *vec, int size, char *str) {

   FILE *pf = fopen(str, "a+");
   int i;

   for (i=0; i<size; i++) {
     fprintf(pf, "%d\n", vec[i]);
   }
   fprintf(pf, "\n\n");
   fflush(pf);

fclose(pf);
return;
}

/********************************************************************************/

void print_dmat(double *mat, int size1, int size2, char *str) {

    FILE *pf = fopen(str, "a+"); 
    int i, j; 

    for (i=0; i<size1; i++) {
      for (j=0; j<size2; j++) {
        fprintf(pf, "%f\t", mat[i*size2+j]); 
      }
      fprintf(pf, "\n"); 
    }
    fprintf(pf, "\n\n"); 
    fflush(pf); 

fclose(pf); 
return; 
}

/*********************************************************************************/
  
  
