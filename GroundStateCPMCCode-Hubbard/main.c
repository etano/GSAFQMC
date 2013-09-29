/*This Is an Example  Ground-State AFQMC Code Meant To Be Used on the Hubbard Model. Many Coding Aspects of It Can Be Improved; It is Meant to Exemplify the Physics. Started on September 18, 2013*/  

#include "afqmc.h"

int main() {

  FILE *pout, *pe;   
  int i,j; 
  int steps; 
  int *neighbors, *number_neighbors; 
  double *weights, *overlap_up, *overlap_down; 
  double *overlap_inverse_up, *overlap_inverse_down; 
  double *wf_up, *wf_down; 
  double *kinetic_full, *kinetic_backwards_half, *kinetic_forwards_half; 
  double *kinetic_eigs, *kinetic_eigvecs; 
  double *trial_wf_up, *trial_wf_down; 
  long idum; 
  cns_st cns; int_st ist;

  /*******************************************************************************************/

  pout = fopen("errors.dat", "a+"); 
  pe = fopen("energy.dat", "w+");  

  /********************************************************************************************/

  /*Obtain Parameters*/
  init(&ist,&cns); 
  Randomize(); idum = -random();

  /*******************************************************************************************/

  neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
  number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

  weights=(double*)calloc(ist.n_walkers,sizeof(double));

  overlap_up=(double*)calloc(ist.n_walkers,sizeof(double));
  overlap_down=(double*)calloc(ist.n_walkers,sizeof(double));

  overlap_inverse_up=(double*)calloc(ist.n_walkers*ist.n_sites*ist.n_sites,sizeof(double));
  overlap_inverse_down=(double*)calloc(ist.n_walkers*ist.n_sites*ist.n_sites,sizeof(double));

  wf_up=(double*)calloc(ist.n_walkers*ist.n_sites*ist.n_up,sizeof(double));
  wf_down=(double*)calloc(ist.n_walkers*ist.n_sites*ist.n_down,sizeof(double));

  trial_wf_up=(double*)calloc(ist.n_sites*ist.n_up,sizeof(double));
  trial_wf_down=(double*)calloc(ist.n_sites*ist.n_down,sizeof(double));

  kinetic_full=(double*)calloc(ist.n_sites_sq,sizeof(double));
  kinetic_backwards_half=(double*)calloc(ist.n_sites_sq,sizeof(double));
  kinetic_forwards_half=(double *)calloc(ist.n_sites_sq,sizeof(double)); 

  kinetic_eigs=(double*)calloc(ist.n_sites,sizeof(double));
  kinetic_eigvecs=(double*)calloc(ist.n_sites_sq,sizeof(double)); 

  /******************************************************************************************/
 
  /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
  init_neighbors(neighbors,number_neighbors,ist); 
  init_kinetic(kinetic_full,kinetic_backwards_half,kinetic_forwards_half,kinetic_eigs,kinetic_eigvecs,neighbors,number_neighbors,ist,cns);
  init_wf(trial_wf_up,trial_wf_down,kinetic_eigs,kinetic_eigvecs,ist,cns); 
  init_walkers(wf_up,wf_down,trial_wf_up,trial_wf_down,weights,overlap_up,overlap_down,ist); 

  /*******************************************************************************************/

  print_dmat(trial_wf_up,ist.n_sites,ist.n_up,"errors.dat"); 

  /*First Regress the Walkers By Half a Kinetic Propagator*/
  propagate_half_backwards_kinetic(kinetic_backwards_half,wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);  

  /*Equilibration Phase***********************************************************************/
  for (steps = 0; steps < ist.n_steps_equilibration; steps++) {

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic(kinetic_full,wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);      

  }

  /*Production Phase*************************************************************************/
  for (steps = 0; steps < ist.n_steps_production; steps++) {

    /*First Propagate a Full Step Forwards*/
    propagate_forwards_kinetic(kinetic_full,wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

    /*Then Propagate Forwards a Full U Step*/
    if ( cns.U != 0 ) {
     propagate_forwards_potential(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist,cns,&idum);
    }

    /*Now Collect Energies and Other Observables*/
    if ( steps != 0 && steps%ist.n_steps_energy==0 ) { 
   
      /*First Propagate One Half Step Forward*/
       propagate_half_forwards_kinetic(kinetic_forwards_half,wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

       /*Compute Energ*/
       compute_energy(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,weights,neighbors,number_neighbors,ist,cns,"energy.dat",steps); 

       /*Propagate Back Backwards*/
       propagate_half_backwards_kinetic(kinetic_backwards_half,wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist); 

     }   

     /*Orthogonalize Q and R if Necessary*/
     if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
        orthogonalize(wf_up,wf_down,overlap_up,overlap_down,weights,ist);   
     }

  }

  /******************************************************************************************/ 

  fclose(pout); 
  fclose(pe); 

  /*******************************************************************************************/

  free(weights); 
  free(overlap_up); free(overlap_down); 
  free(overlap_inverse_up); free(overlap_inverse_down); 
  free(wf_up); free(wf_down); 
  free(trial_wf_up); free(trial_wf_down); 
  free(kinetic_full); free(kinetic_backwards_half); free(kinetic_forwards_half);  
  free(kinetic_eigs); free(kinetic_eigvecs);  

return 0; 
}
