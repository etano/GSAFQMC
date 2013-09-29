#include "afqmc.h"

void init(int_st *ist,cns_st *cns) { 

  FILE *pf = fopen("afqmc.par", "r"); 
  FILE *po = fopen("afqmc-parameters.dat", "w+"); 

  /*Scan In Relevant Parameters*/

  fscanf (pf,"%d",&ist->n_sites_one); 
  fscanf (pf,"%d",&ist->n_sites_two); 
  fscanf (pf,"%d",&ist->n_up); 
  fscanf (pf,"%d",&ist->n_down); 
  fscanf (pf,"%d",&ist->n_walkers); 
  fscanf (pf,"%d",&ist->n_steps_energy); 
  fscanf (pf,"%d",&ist->n_steps_orthogonalize); 
  fscanf (pf,"%d",&ist->n_steps_free); 
  fscanf (pf,"%d",&ist->flag_trial); 
  fscanf (pf,"%d",&ist->flag_pbc);   
  fscanf (pf,"%d",&ist->flag_cp); 

  fscanf (pf,"%lg",&cns->beta);
  fscanf (pf,"%lg",&cns->beta_equilibration);  
  fscanf (pf,"%lg",&cns->dtau); 
  fscanf (pf,"%lg",&cns->U); 
  fscanf (pf,"%lg",&cns->t); 

  /*Steps*/
  ist->n_sites = ist->n_sites_one * ist->n_sites_two; 
  ist->n_steps = (int)(cns->beta/cns->dtau); 
  ist->n_steps_equilibration = (int)(cns->beta_equilibration/cns->dtau); 
  ist->n_steps_production = ist->n_steps - ist->n_steps_equilibration; 
  ist->n_sites_sq = ist->n_sites * ist->n_sites;  

  /*Print Out Parameters Scanned In As a Check*/ 

  fprintf (po,"Sites One = %d\n", ist->n_sites_one); 
  fprintf (po,"Sites Two = %d\n", ist->n_sites_two); 
  fprintf (po,"Sites = %d\n", ist->n_sites); 
  fprintf (po,"Number Electrons Up = %d\n", ist->n_up); 
  fprintf (po,"Number Electrons Down = %d\n", ist->n_down); 
  fprintf (po,"Number of Walkers = %d\n", ist->n_walkers); 
  fprintf (po,"Nsteps = %d\n", ist->n_steps);
  fprintf (po,"Nstepts Equilibration = %d\n", ist->n_steps_equilibration);
  fprintf (po,"Nsteps Energy = %d\n", ist->n_steps_energy); 
  fprintf (po,"Nsteps Orthogonalize = %d\n", ist->n_steps_orthogonalize); 
  fprintf (po,"Nsteps Free = %d\n", ist->n_steps_free);  
  fprintf (po,"Flag CP = %d\n", ist->flag_cp);  
  fprintf (po,"Trial Wavefunction Type = %d\n", ist->flag_trial); 
  fprintf (po,"Periodic Boundary Conditions = %d\n", ist->flag_pbc); 

  /*Find Gamma Given U and Delta Tau*/
  cns->gamma = find_gamma(cns->U,cns->dtau); 

  /*Find Exponentials for Spin Up and Spin Down and Related Fields*/
  cns->factor_spin_up_field_down = exp(cns->gamma); 
  cns->factor_spin_up_field_up = exp(-1*cns->gamma); 

  fprintf (po,"Beta = %g\n", cns->beta);
  fprintf (po,"Beta Equilib = %g\n", cns->beta_equilibration);  
  fprintf (po,"Dtau = %g\n", cns->dtau); 
  fprintf (po,"U = %g\n", cns->U); 
  fprintf (po,"t = %g\n", cns->t);
  fprintf (po,"gamma = %g\n", cns->gamma);
  fprintf (po,"factor_spin_up_field_up = %g\n", cns->factor_spin_up_field_up); 
  fprintf (po,"factor_spin_up_field_down = %g\n", cns->factor_spin_up_field_down);  
  fflush(po);    

  fclose(pf); 
  fclose(po);  

return; 
}

/*************************************************************/

void init_neighbors(int *neighbors,int *number_neighbors,int_st ist){

  int i;  

  /*Finds the Neighbors for Each Sites*/
  if (ist.flag_pbc == 1) {
   for (i=0; i<ist.n_sites; i++) {
      neighborsperiodicboundary(i,neighbors,number_neighbors,ist);
   }
  }
  else {
   for (i=0; i<ist.n_sites; i++) {
      neighborsopenboundary(i,neighbors,number_neighbors,ist);
   }
  }
 
return; 
}

/*************************************************************/

void init_walkers(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *weights,double *overlap_up,double *overlap_down,int_st ist) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2; 

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<ist.n_walkers; iw++) {
     weights[iw] = 1.0; 
     overlap_up[iw] = 1.0; 
     overlap_down[iw] = 1.0;

     iwp1 = iw * ist.n_sites * ist.n_up;      
     iwp2 = iw * ist.n_sites * ist.n_down;  

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {
              
       isp1 = is * ist.n_up;
       isp2 = is * ist.n_down;   

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) { 
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip]; 
       }
       for ( ip = 0 ; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp2 + isp2 + ip] = trial_wf_down[ isp2 + ip]; 
       }  

     } /*is*/

  } /*iw*/ 

return; 
}

/**************************************************************/

void init_wf(double *trial_wf_up,double *trial_wf_down,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_trial == 0 ) {
    /*Just an Identity Matrix*/

    trial_identity(trial_wf_up,trial_wf_down,ist); 

  }
  else if ( ist.flag_trial == 1 ) {
    /*The Free Case*/

    trial_free(trial_wf_up,trial_wf_down,kinetic_eigs,kinetic_eigvecs,ist,cns);  

  }
  else if ( ist.flag_trial == 2 ) {
     /*The Restricted Hartree Fock Case*/

  }
  else if ( ist.flag_trial == 3 ) { 
      /*Just a Random WF*/

     trial_random(trial_wf_up,trial_wf_down,ist); 

  }
  else { 

  }

return; 
}

/************************************************************/

void init_kinetic(double *kinetic_full,double *kinetic_backwards,double *kinetic_forwards_half,double *kinetic_eigs,double *kinetic_eigvecs,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

  /*Creates Kinetic Matrices********/

  int i, j, k;

  /*If Non-Zero T***************/
  if (cns.t!=0) {

    /*Creating Original Kinetic Matrix and Then Diagonalizing It*/
    for (i=0; i<ist.n_sites;  i++) {
      for (j=0; j<number_neighbors[i]; j++) {
        kinetic_full[i*ist.n_sites+neighbors[i*4+j]]=-cns.t;
        kinetic_full[neighbors[i*4+j]*ist.n_sites+i]=kinetic_full[i*ist.n_sites+neighbors[i*4+j]];
      }   
     }

    /*Now Compute Eigenvalues and Eigenvectors Using NRC Routines*/
    jacobireal(kinetic_full, kinetic_eigs, kinetic_eigvecs,ist);
    dzero_vec(kinetic_full,ist.n_sites_sq); 
    dzero_vec(kinetic_backwards,ist.n_sites_sq); 
    dzero_vec(kinetic_forwards_half,ist.n_sites_sq); 

    /*Find Various Kinetic Terms*/
    for (k=0; k<ist.n_sites; k++) {
     for (i=0; i<ist.n_sites; i++) {
      for (j=0; j<ist.n_sites; j++) {
        kinetic_full[i*ist.n_sites+j]+=exp(-1*kinetic_eigs[k]*cns.dtau)*kinetic_eigvecs[i*ist.n_sites+k]*kinetic_eigvecs[j*ist.n_sites+k];
        kinetic_backwards[i*ist.n_sites+j]+=exp(kinetic_eigs[k]*cns.dtau/2.0)*kinetic_eigvecs[i*ist.n_sites+k]*kinetic_eigvecs[j*ist.n_sites+k];
        kinetic_forwards_half[i*ist.n_sites+j]+=exp(-1*kinetic_eigs[k]*cns.dtau/2.0)*kinetic_eigvecs[i*ist.n_sites+k]*kinetic_eigvecs[j*ist.n_sites+k];
      }
     }
    }

  }
  else { /*If Zero T********************************************/

    /*Make Fermions Matrices Diagonal*/
    for (i=0; i<ist.n_sites; i++) {
     for (j=0; j<ist.n_sites; j++) {
       kinetic_full[i*ist.n_sites+j]=kinetic_backwards[i*ist.n_sites+j]=kinetic_forwards_half[i*ist.n_sites+j]=0; 
       kinetic_eigvecs[i*ist.n_sites+j]=0;
     }
     kinetic_eigs[i]=0;
     kinetic_full[i*ist.n_sites+i]=kinetic_backwards[i*ist.n_sites+i]=kinetic_forwards_half[i*ist.n_sites+i]=1; 
     kinetic_eigvecs[i*ist.n_sites+i]=1;   
   }

  } /**********************************************************/

return; 
}

/***********************************************************/

void init_kinetic_many(double *kinetic_many_full,double *kinetic_eigs,double *kinetic_eigvecs,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

    /*Create A Kinetic Matrix with Several Large Time Steps To Create Free Trial Wavefunction*/
   
    int i, j, k; 
    double full_time_step = cns.dtau * ist.n_steps_orthogonalize; 
    FILE *pf = fopen("errors.dat", "a+"); 
 
    /*If Non-Zero T*/
    if (cns.t != 0 ) {  
  
       fprintf(pf, "before dzer\n");  fflush(pf); 
       dzero_vec(kinetic_many_full,ist.n_sites_sq);
 
       /*Find Various Kinetic Terms*/
       for (k=0; k<ist.n_sites; k++) {
         for (i=0; i<ist.n_sites; i++) {
           for (j=0; j<ist.n_sites; j++) {
             kinetic_many_full[i*ist.n_sites+j]+=exp(-1*kinetic_eigs[k]*full_time_step)*kinetic_eigvecs[i*ist.n_sites+k]*kinetic_eigvecs[j*ist.n_sites+k];
           }
         }
       }  

      fprintf(pf,"in kinetic init %f\n", full_time_step); fflush(pf); 
      print_dmat(kinetic_many_full,ist.n_sites,ist.n_sites,"errors.dat"); 

    }
    else {  /*If T=0***************************************************************/
   
       /*Make Fermions Matrices Diagonal*/
       for (i=0; i<ist.n_sites; i++) {
        for (j=0; j<ist.n_sites; j++) {
         kinetic_many_full[i*ist.n_sites+j]=0.0; 
        }
       kinetic_many_full[i*ist.n_sites+i]=1.0; 
     } 

   }

fclose(pf); 
return; 
}

/**********************************************************/ 
