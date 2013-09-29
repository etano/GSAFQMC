#include "afqmc.h"

/******************************************************************************************/

double compute_overlap(double *wf_1, double *wf_2, double *overlap_matrix_inverse, int size) {

   /*Computes the Overlap As Well as the Inverse of the Overlap Matrix*/
   int size_sq = size * size;
   double overlap;  
   double *overlap_marix; 

   overlap_matrix=(double*)calloc(size_sq,sizeof(double)); 

   /*First Obtain Overlap Matrix*/
   transpose_mat_mat(wf_1,wf_2,overlap_matrix_inverse,ist.n_sites,size,size);  
    
 



free(overlap_matrix); 
return(overlap); 
}

/******************************************************************************************/
