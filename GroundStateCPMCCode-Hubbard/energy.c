#include "afqmc.h"

/***********************************************************************************************/

void compute_energy(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int *neighbors,int *number_neighbors,int_st ist,cns_st cns,char *str,int step) {

   /*Computes the Energy of the System*/ 
   int walker, i, j; 
   int n_sites_n_up = ist.n_sites * ist.n_up; 
   int n_sites_n_down = ist.n_sites * ist.n_down; 
   double *greens_function_up, *greens_function_down;
   double local_ke, local_pe; 
   double energy = 0.0, kinetic_energy = 0.0, potential_energy = 0.0;   
   double weight_denominator = 0.0;  
   FILE *pf = fopen(str, "a+"); 
   FILE *pf2 = fopen("errors.dat", "a+"); 

   greens_function_up=(double *)calloc(ist.n_sites_sq,sizeof(double)); 
   greens_function_down=(double *)calloc(ist.n_sites_sq,sizeof(double)); 
  
   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/    
     if ( weights[walker] != 0 ) { 

        /*Accumulate Weight in Denominator*/
        weight_denominator += weights[walker]; 

        /*First Compute the Green's Function*/
        compute_greens_function(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*n_sites_n_up],&overlap_inverse_down[walker*n_sites_n_down],greens_function_up,greens_function_down,ist);  

         print_dmat(&greens_function_up[walker*n_sites_n_up],ist.n_sites,ist.n_sites,"errors.dat"); 
         print_dmat(&greens_function_down[walker*n_sites_n_down],ist.n_sites,ist.n_sites,"errors.dat");

        /*Obtain Kinetic and Potential Energies*/ 
        local_ke = local_pe = 0.0; 
        for ( i=0; i<ist.n_sites; i++) {
         local_pe += (1 - greens_function_up[i*ist.n_sites+i] ) * (1 - greens_function_down[i*ist.n_sites+i]);  
         for (j=0; j<number_neighbors[i]; j++) {
            local_ke += greens_function_up[i*ist.n_sites+neighbors[4*i+j]] + greens_function_down[i*ist.n_sites+neighbors[4*i+j]]; 
         }
        }

        if ( weights[walker]< 0 ) {
          fprintf(pf2, "weights %f\n", weights[walker]); fflush(pf2); 
        } 
 
        /*Half Energy If Double Counting Occurs*/
        if ( ist.n_sites_one == 2 && ist.n_sites_two == 2 ) {
           local_ke *= .5; 
        } 

        local_ke *= cns.t * .5 * weights[walker]; 
        local_pe *= cns.U * weights[walker];  

        /*Add to Total*/
        kinetic_energy += local_ke; 
        potential_energy += local_pe;          
        energy += local_ke + local_pe;  


      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy /= weight_denominator; 
   potential_energy /= weight_denominator; 
   energy /= weight_denominator;  

   /*Print Desired Information*******************************************/
   fprintf(pf, "%d\t\t  %f\t\t   %f\t\t   %f\t\t   %f\t\n", step, step*cns.dtau, kinetic_energy, potential_energy, energy); fflush(pf);     

fclose(pf); 
fclose(pf2); 
free(greens_function_up); 
free(greens_function_down); 
return; 
}

/***********************************************************************************************/

void compute_greens_function(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *greens_function_up,double *greens_function_down,int_st ist) {

   /*Computes the Spin Up and Down Greens Functions*/
   int i, j; 
   double *stored_product1_up, *stored_product1_down; 
   double *stored_product2_up, *stored_product2_down; 
 
   stored_product1_up=(double *)calloc(ist.n_sites*ist.n_up,sizeof(double)); 
   stored_product1_down=(double *)calloc(ist.n_sites*ist.n_down,sizeof(double));  
     
   stored_product2_up=(double *)calloc(ist.n_sites_sq,sizeof(double)); 
   stored_product2_down=(double *)calloc(ist.n_sites_sq,sizeof(double));    

   /*Get Up Green's Function First*/
   mat_transpose_mat(overlap_inverse_up,trial_wf_up,stored_product1_up,ist.n_up,ist.n_sites,ist.n_up);
   mat_mat(wf_up,stored_product1_up,stored_product2_up,ist.n_sites,ist.n_up,ist.n_sites); 

   /*Get Down Green's Function*/
   mat_transpose_mat(overlap_inverse_down,trial_wf_down,stored_product1_down,ist.n_down,ist.n_sites,ist.n_down); 
   mat_mat(wf_down,stored_product1_down,stored_product2_down,ist.n_sites,ist.n_down,ist.n_sites);

    /*Now Get Greens Function By Subtracting From Unity*/  
   for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_sites; j++) {
      greens_function_up[i*ist.n_sites+j] = -1 * stored_product2_up[i*ist.n_sites+j];
      greens_function_down[i*ist.n_sites+j] = -1 * stored_product2_down[i*ist.n_sites+j];
    }
    greens_function_up[i*ist.n_sites+i] += 1.0;
    greens_function_down[i*ist.n_sites+i] += 1.0; 
   }


free(stored_product1_up); 
free(stored_product1_down);
free(stored_product2_up); 
free(stored_product2_down);  
return; 
}

/***********************************************************************************************/
