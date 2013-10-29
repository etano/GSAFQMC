#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_kinetic(double *kinetic_full,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,double *weights,int_st ist) {

   /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

   int i; 
   int n_sites_n_up = ist.n_sites * ist.n_up; 
   int n_sites_n_down = ist.n_sites * ist.n_down; 
   double previous_overlap_up, previous_overlap_down; 
   double *stored_product_up, *stored_product_down;    
   double *stored_product_up_2, *stored_product_down_2; 

   stored_product_up=(double*)calloc(ist.n_sites*ist.n_up,sizeof(double)); 
   stored_product_down=(double*)calloc(ist.n_sites*ist.n_down,sizeof(double));

   stored_product_up_2=(double*)calloc(ist.n_up*ist.n_up,sizeof(double));
   stored_product_down_2=(double*)calloc(ist.n_down*ist.n_down,sizeof(double)); 

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( weights[i] != 0 ) { 
  
      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i]; 
      previous_overlap_down = overlap_down[i]; 

      /*First Get Up Matrices*/
      mat_mat(kinetic_full,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
      copy_mat(&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites*ist.n_up); 

      /*Get Product of Trial and Actual*/
      transpose_mat_mat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up); 
      overlap_up[i]=inverse_det(stored_product_up_2,&overlap_inverse_up[i*n_sites_n_up],ist.n_up); 

      /*Get Down Matrices*/
      mat_mat(kinetic_full,&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
      copy_mat(&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites*ist.n_down);
   
      /*Get Product of Trial and Down Actual*/
      transpose_mat_mat(trial_wf_down,&wf_down[i*n_sites_n_down],stored_product_down_2,ist.n_sites,ist.n_down,ist.n_down);
      overlap_down[i]=inverse_det(stored_product_down_2,&overlap_inverse_down[i*n_sites_n_down],ist.n_down); 

      /*Get Overall Weight*/
      weights[i] *= overlap_up[i] * overlap_down[i] / ( previous_overlap_up * previous_overlap_down );  
 
    } /*weights walkers*/

   } /*i*/

free(stored_product_up); 
free(stored_product_down); 
free(stored_product_up_2); 
free(stored_product_down_2); 
return; 
}

/******************************************************************/

void propagate_half_backwards_kinetic(double *kinetic_backwards_half,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_up,double *overlap_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int_st ist){ 

    /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

    int i;
    int n_sites_n_up = ist.n_sites * ist.n_up; 
    int n_sites_n_down = ist.n_sites * ist.n_down; 
    double previous_overlap_up, previous_overlap_down; 
    double *stored_product_up, *stored_product_down;
    double *stored_product_up_2, *stored_product_down_2; 

    stored_product_up=(double*)calloc(ist.n_sites*ist.n_up,sizeof(double));
    stored_product_down=(double*)calloc(ist.n_sites*ist.n_down,sizeof(double));

    stored_product_up_2=(double*)calloc(ist.n_up*ist.n_up,sizeof(double));
    stored_product_down_2=(double*)calloc(ist.n_down*ist.n_down,sizeof(double));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( weights[i] != 0 ) { 

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i]; 
       previous_overlap_down = overlap_down[i]; 

       /*First Get Up Matrices*/
       mat_mat(kinetic_backwards_half,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       copy_mat(&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites*ist.n_up);

       /*Get Product of Trial and Actual*/
       transpose_mat_mat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up);
       overlap_up[i]=inverse_det(stored_product_up_2,&overlap_inverse_up[i*n_sites_n_up],ist.n_up);
 
       /*Get Down Matrices*/
       mat_mat(kinetic_backwards_half,&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       copy_mat(&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites*ist.n_down);

       /*Get Product of Trial and Actual*/
       transpose_mat_mat(trial_wf_down,&wf_down[i*n_sites_n_down],stored_product_down_2,ist.n_sites,ist.n_down,ist.n_down);
       overlap_down[i]=inverse_det(stored_product_down_2,&overlap_inverse_down[i*n_sites_n_down],ist.n_down);

       /*Get Overall Weights*/
       weights[i] *= overlap_up[i] * overlap_down[i] / ( previous_overlap_up * previous_overlap_down );   
  
      } /*Non Zero Weights*/ 

     }

free(stored_product_up);
free(stored_product_down);
free(stored_product_up_2); 
free(stored_product_down_2); 
return;   
}

/***********************************************************************/

void propagate_half_forwards_kinetic(double *kinetic_forwards_half,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_up,double *overlap_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Half Kinetic Time Slice*/

    int i;
    int n_sites_n_up = ist.n_sites * ist.n_up;
    int n_sites_n_down = ist.n_sites * ist.n_down;
    double previous_overlap_up, previous_overlap_down;
    double *stored_product_up, *stored_product_down;
    double *stored_product_up_2, *stored_product_down_2;

    stored_product_up=(double*)calloc(ist.n_sites*ist.n_up,sizeof(double));
    stored_product_down=(double*)calloc(ist.n_sites*ist.n_down,sizeof(double));

    stored_product_up_2=(double*)calloc(ist.n_up*ist.n_up,sizeof(double));
    stored_product_down_2=(double*)calloc(ist.n_down*ist.n_down,sizeof(double));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( weights[i] != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       mat_mat(kinetic_forwards_half,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       copy_mat(&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites*ist.n_up);

       /*Get Product of Trial and Actual*/
       transpose_mat_mat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up);
       overlap_up[i]=inverse_det(stored_product_up_2,&overlap_inverse_up[i*n_sites_n_up],ist.n_up);

       /*Get Down Matrices*/
       mat_mat(kinetic_forwards_half,&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       copy_mat(&wf_down[i*n_sites_n_down],stored_product_down,ist.n_sites*ist.n_down);

       /*Get Product of Trial and Actual*/
       transpose_mat_mat(trial_wf_down,&wf_down[i*n_sites_n_down],stored_product_down_2,ist.n_sites,ist.n_down,ist.n_down);
       overlap_down[i]=inverse_det(stored_product_down_2,&overlap_inverse_down[i*n_sites_n_down],ist.n_down);

       /*Get Overall Weights*/
       weights[i] *= overlap_up[i] * overlap_down[i] / ( previous_overlap_up * previous_overlap_down );
  
      }

     }

free(stored_product_up);
free(stored_product_down);
free(stored_product_up_2);
free(stored_product_down_2);
return;
} 
