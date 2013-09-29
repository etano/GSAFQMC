#include "ed.h"

/*******************************************************************************/

void init(int_st *ist,cns_st *cns) {

  FILE *pf = fopen("ed.par", "r");
  FILE *po = fopen("ed-parameters.dat", "a+");

  /*Scan In Relevant Parameters*/

  fscanf (pf,"%d",&ist->dimension); 
  fscanf (pf,"%d",&ist->n_sites_one);
  fscanf (pf,"%d",&ist->n_sites_two);
  fscanf (pf,"%d",&ist->n_sites_three); 
  fscanf (pf,"%d",&ist->n_up);
  fscanf (pf,"%d",&ist->n_down);
  fscanf (pf,"%d",&ist->flag_pbc);
  fscanf (pf,"%d",&ist->flag_spin_polarized); 
  fscanf (pf,"%d",&ist->flag_ground_state); 

  fscanf (pf,"%lg",&cns->U);
  fscanf (pf,"%lg",&cns->t);
  fscanf (pf,"%lg",&cns->beta); 

  /*Steps*/
  ist->n_sites = ist->n_sites_one * ist->n_sites_two * ist->n_sites_three;
  ist->n_sites_sq = ist->n_sites * ist->n_sites;
  
  /*Find States Per Site*/
  if ( ist->flag_spin_polarized == 0 ) {
    ist->states_per_site = 4; 
  }
  else {
    ist->states_per_site = 2; 
    ist->n_down = 0; 
  }  
  cns->total_number_possible_states = pow(ist->states_per_site, ist->n_sites); 

  /*Print Out Parameters Scanned In As a Check*/

  fprintf (po,"Dimension = %d\n", ist->dimension); 
  fprintf (po,"Sites One = %d\n", ist->n_sites_one);
  fprintf (po,"Sites Two = %d\n", ist->n_sites_two);
  fprintf (po,"Sites Three = %d\n", ist->n_sites_three);  
  fprintf (po,"Sites = %d\n", ist->n_sites);
  fprintf (po,"States Per Site = %d\n", ist->states_per_site); 
  fprintf (po,"Number Electrons Up = %d\n", ist->n_up);
  fprintf (po,"Number Electrons Down = %d\n", ist->n_down);
  fprintf (po,"Periodic Boundary Conditions = %d\n", ist->flag_pbc);
  fprintf (po,"Spin Polarized = %d\n", ist->flag_spin_polarized); 
  fprintf (po,"Ground State = %d\n", ist->flag_ground_state);  

  fprintf (po,"U = %g\n", cns->U);
  fprintf (po,"t = %g\n", cns->t);
  fprintf (po,"Total Possible Number of States = %d\n", cns->total_number_possible_states);
  fprintf (po,"Beta = %g\n", cns->beta);  
  fflush(po);

  fclose(pf);
  fclose(po);

return;
}

/******************************************************************************/

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
