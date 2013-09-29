#include "afqmc.h"

/*****************************************************************************/

void orthogonalize(double *wf_up,double *wf_down,double *overlap_up,double *overlap_down,double *weights,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker; 
    int n_sites_n_up = ist.n_sites * ist.n_up; 
    int n_sites_n_down = ist.n_sites * ist.n_down; 
    double *R_up, *R_down; 
    double determinant; 
    double old_overlap;

    R_up=(double *)calloc(ist.n_up*ist.n_up,sizeof(double)); 
    R_down=(double *)calloc(ist.n_down*ist.n_down,sizeof(double)); 

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( weights[walker] != 0 ) {

         modified_gram_schmidt(&wf_up[walker*n_sites_n_up],R_up,ist.n_sites,ist.n_up); 
         det(R_up,ist.n_up,&determinant);
         old_overlap = overlap_up[walker];  
         overlap_up[walker] *= determinant; 
         weights[walker] *= overlap_up[walker] / old_overlap; 

         modified_gram_schmidt(&wf_down[walker*n_sites_n_down],R_down,ist.n_sites,ist.n_down);
         det(R_down,ist.n_down,&determinant); 
         old_overlap = overlap_down[walker]; 
         overlap_down[walker] *= determinant;
         weights[walker] *= overlap_down[walker] / old_overlap;  

      } 

    }  

free(R_up); 
free(R_down); 
return; 
}

/*********************************************************************************/ 

void orthogonalize_without_weights(double *wf_up,double *wf_down,int_st ist) {

    /*Perform Gram-Schmidt Orthogonalization, But Without Weights*/
    double *R_up, *R_down;  
     
    R_up=(double *)calloc(ist.n_up*ist.n_up,sizeof(double));
    R_down=(double *)calloc(ist.n_down*ist.n_down,sizeof(double));

    modified_gram_schmidt(wf_up,R_up,ist.n_sites,ist.n_up);
    modified_gram_schmidt(wf_down,R_down,ist.n_sites,ist.n_down);

free(R_up); 
free(R_down); 
return; 
}

/********************************************************************************/ 
