#include "afqmc.h"

/************************************************************************************/

double find_gamma(double U,double dtau) { 

   /*Finds Gamma Given U*********************************/
   int i;
   double gamma;  
   double gamma_interval = 0.00001; 
   double cosh_error = .0001; 
   double exponential_value = 2.0 * exp(dtau * U * .5);

   for (i=0; i<100000000; i++) {
      gamma = i * gamma_interval;
      if ( fabs(exp(gamma) + exp(-1.0 * gamma) - exponential_value) < cosh_error ) {
         return(gamma); 
      }
   }    

return(0); 
} 

/************************************************************************************/

void propagate_forwards_potential(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,double *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;   
   int n_sites_n_up = ist.n_sites * ist.n_up; 
   int n_sites_n_down = ist.n_sites * ist.n_down; 
   int n_up_sq = ist.n_up * ist.n_up; 
   int n_down_sq = ist.n_down * ist.n_down; 
   int flag_neg_overlap_up, flag_neg_overlap_down; 
   double probability_field_up, probability_field_down; 
   double spin_up = 1.0, spin_down = -1.0; 
   double field_up = 1.0, field_down = -1.0; 
   double field_selected; 
   double random_number, total_probability; 
   long tidum; 
   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {
    
        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) { 
 
          /*Check WEights*/
          if ( weights[walker] != 0 ) { 
 
            /*Ensure That Overlap Flags Are Positive*/
            flag_neg_overlap_up = flag_neg_overlap_down = 0;  
 
            /*Determine Probabilities for Each Scenario*/
            probability_field_up = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_up,sites);                 
            probability_field_up *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_up,sites);
            probability_field_up *= .5;

            probability_field_down = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_down,sites);
            probability_field_down *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_down,sites);
            probability_field_down *= .5;

            /*If Constrained Path Approximation Applied - With Mirror Correction To Check If Next Timeslice Will Cause Overlap and to Reduce Weights If Next Step Does*/
            /*if ( ist.flag_cp == 1 ) {
              if ( probability_field_up < 0 ) {
                 flag_neg_overlap_up = 1; 
                 probability_field_up = 0;
                 weights[walker] /= (1-probability_field_up*2.0);  
              }
              if ( probability_field_down < 0 ) {
                 flag_neg_overlap_down = 1;
                 probability_field_down = 0;  
                 weights[walker] /= (1-probability_field_down*2.0); 
              }
            }*/

            total_probability = probability_field_down + probability_field_up; 
 
            /*Now Propagate If Weight Is Not Negative************************************************/
             if ( total_probability <= 0 && ist.flag_cp == 1 ) {
                 weights[walker]=0; 
             }
             else {  
 
                /*Perform Monte Carlo Decision for Which Field to Pick*/
                random_number = ran_nrc(&tidum);  
                if ( random_number < probability_field_up/total_probability ) { /*If Spin Up Selected*/
                  weights[walker] *= total_probability; 
  
                  propagate_wave_functions(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],ist,cns,field_up,sites); 
                  update_overlaps(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*n_up_sq],&overlap_inverse_down[walker*n_down_sq],&overlap_up[walker],&overlap_down[walker],ist); 
  
                  field_selected = field_up; 
                }
                else { /*If Spin Down Selected*/  
                  weights[walker] *= total_probability; 
  
                  propagate_wave_functions(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],ist,cns,field_down,sites);
                  update_overlaps(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*n_up_sq],&overlap_inverse_down[walker*n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);
 
                  field_selected = field_down; 

               } 
            } /*If Not Traversing Boundary*/
            /*Done Propagating********************************************************************/

            /*Apply Mirror Correction Again Afterwards to Reduce Weight If Next Time Slice Will Result in Crossing*/
            /*if ( ist.flag_cp == 1 ) {

               if ( flag_neg_overlap_up == 0 && flag_neg_overlap_down == 0 ) {
                  if ( field_selected == field_up ) {
                      probability_field_up = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_up,sites);    
                      probability_field_up *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_up,sites);
                      if ( probability_field_up < 0 ) {
                         weights[walker] /= (1-probability_field_up);  
                      }
                  }
                  else {
                      probability_field_down = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_down,sites);
                      probability_field_down *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_down,sites);
                      if ( probability_field_down < 0 ) {
                         weights[walker] /= (1-probability_field_down); 
                      } 
                  } 
               }
                 
            }*/ /*Check Mirror Correction********************/  


         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return; 
}

/***************************************************************************************/
 
double determine_overlap_ratio(double *wf,double *trial_wf,double *overlap_inverse,int_st ist,cns_st cns,double spin_up,double field_up,int site) {

      /*Determines Overlap Ratio*/
      double overlap_ratio;  
      double *stored_product1, *stored_product2; 

      /*If An Up-Spin*/
      if (spin_up == 1.0 ) {
        stored_product1=(double *)calloc(ist.n_sites*ist.n_up,sizeof(double)); 
        stored_product2=(double *)calloc(ist.n_sites_sq,sizeof(double));  

        /*Get Matrix For Green's Function*/
        mat_transpose_mat(overlap_inverse,trial_wf,stored_product1,ist.n_up,ist.n_sites,ist.n_up);
        mat_mat(wf,stored_product1,stored_product2,ist.n_sites,ist.n_up,ist.n_sites);

        if ( field_up == 1.0 ) {
          overlap_ratio = 1 + (cns.factor_spin_up_field_up - 1.0 ) *  stored_product2[site * ist.n_sites + site]; 
        }
        else { 
          overlap_ratio = 1 + (cns.factor_spin_up_field_down - 1.0 ) *  stored_product2[site * ist.n_sites + site];
        }  

      }
      else { /*If Spin-Down*/
        stored_product1=(double *)calloc(ist.n_sites*ist.n_down,sizeof(double));
        stored_product2=(double *)calloc(ist.n_sites_sq,sizeof(double)); 

        /*Get Matrix For Green's Function*/
        mat_transpose_mat(overlap_inverse,trial_wf,stored_product1,ist.n_down,ist.n_sites,ist.n_down);
        mat_mat(wf,stored_product1,stored_product2,ist.n_sites,ist.n_down,ist.n_sites);

        if ( field_up == 1.0 ) {
          overlap_ratio = 1 + (cns.factor_spin_up_field_down - 1.0 ) *  stored_product2[site * ist.n_sites + site];
        }
        else {
           overlap_ratio = 1 + (cns.factor_spin_up_field_up - 1.0 ) *  stored_product2[site * ist.n_sites + site];
        }  

     }

free(stored_product1); 
free(stored_product2); 
return(overlap_ratio);  
}

/****************************************************************************************/  

void propagate_wave_functions(double *wf_up,double *wf_down,int_st ist,cns_st cns,double field_up,int sites) {

     /*Propagates the Wave Function*/
     int i; 

     /*For Up Wave Function*/ 
     for (i=0; i<ist.n_up; i++) {
        if ( field_up == 1.0 ) {
          wf_up[sites*ist.n_up+i] *= cns.factor_spin_up_field_up; 
        }
        else {
          wf_down[sites*ist.n_down+i] *= cns.factor_spin_up_field_down;   
        }
     }

     /*For Down Wave Function*/
     for (i=0; i<ist.n_down; i++) {
        if ( field_up == 1.0 ) {
          wf_down[sites*ist.n_down+i] *= cns.factor_spin_up_field_down; 
        }
        else {
          wf_down[sites*ist.n_down+i] *= cns.factor_spin_up_field_up; 
        }
     }

return; 
}

/****************************************************************************************/

void update_overlaps(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/ 
    double *stored_product_up; 
    double *stored_product_down; 

    stored_product_up=(double *)calloc(ist.n_up*ist.n_up,sizeof(double)); 
    stored_product_down=(double *)calloc(ist.n_down*ist.n_down,sizeof(double)); 

    /*Get Product of Trial and Actual*/
    transpose_mat_mat(trial_wf_up,wf_up,stored_product_up,ist.n_sites,ist.n_up,ist.n_up);
    (*overlap_up)=inverse_det(stored_product_up,overlap_inverse_up,ist.n_up);

    transpose_mat_mat(trial_wf_down,wf_down,stored_product_down,ist.n_sites,ist.n_down,ist.n_down);
    (*overlap_down)=inverse_det(stored_product_down,overlap_inverse_down,ist.n_down);
     
free(stored_product_up); 
free(stored_product_down); 
return; 
}

/*****************************************************************************************/
