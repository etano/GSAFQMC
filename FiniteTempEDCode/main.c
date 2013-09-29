/*This is a Canonical Ensemble Code to Compute the Ground State Energies and Green's Functions for Electronic Systems to Check Against AFQMC*/

#include "ed.h"

/***************************************************************************************/

int main() {

  FILE *pout; 
  int i, j; 
  int total_overall_number_states = 0; 
  int *neighbors, *number_neighbors;
  int *electron_state, *temporary_electron_state; 
  int *possible_states, *inverse_possible_states; 
  double *full_matrix; 
  double *eigenvalues, *eigenvectors; 
  double *overall_eigenvalues; 
  cns_st cns; int_st ist;

  /**********************************************************************/

  pout = fopen("errors.dat", "a+");

  /**********************************************************************/

  /*Obtain Parameters*/
  init(&ist,&cns);

  /**********************************************************************/

  /*Initialize Relevant Vectors and Matrices*/
  neighbors=(int*)calloc(ist.dimension*2*ist.n_sites,sizeof(int)); 
  number_neighbors=(int*)calloc(ist.n_sites,sizeof(int)); 

  electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 
  temporary_electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 

  inverse_possible_states=(int*)calloc(cns.total_number_possible_states,sizeof(int)); 

  /*********************************************************************/ 

  /*Get the Lattice Site Neighbors*/
  init_neighbors(neighbors,number_neighbors,ist);

  /*Run GS Calculat If GS*/
  if ( ist.flag_ground_state == 1 ) { 

    /*Determine Total Number of States*/  
    total_number_states(ist,&cns);

    possible_states=(int*)calloc(cns.total_number_states,sizeof(int));
    full_matrix=(double*)calloc(cns.total_number_states_sq,sizeof(double));

    eigenvalues=(double*)calloc(cns.total_number_states,sizeof(double)); 
    eigenvectors=(double*)calloc(cns.total_number_states_sq,sizeof(double));

    /*Determine Possible States To Form Hubbard Matrix*/
    init_states(possible_states,inverse_possible_states,ist,cns);  

    /*Form the Hubbard Hamiltonian*/
    form_hubbard(full_matrix,possible_states,inverse_possible_states,neighbors,number_neighbors,ist,cns);  

    /*Diagonalize the Hamiltonian Formed*/
    jacobireal(full_matrix,eigenvalues,eigenvectors,ist,cns); 

    /*Determine the Energy and Print It*/
    determine_energy_ground_state(eigenvalues,ist,cns,"energy.dat"); 

    free(possible_states); free(inverse_possible_states);
    free(eigenvalues); free(eigenvectors);

  }
  else { /*If a Finite Temperature Calculation*/

    overall_eigenvalues=(double*)calloc(cns.total_number_possible_states,sizeof(double));

    /*Run Through All Possible Combinations of Fillings*/
    for ( ist.n_up = 0; ist.n_up <= ist.n_sites; ist.n_up++) {
     for ( ist.n_down = 0; ist.n_down <= ist.n_sites; ist.n_down++) { 

      /*Determine Total Number of States*/
      total_number_states(ist,&cns);

      possible_states=(int*)calloc(cns.total_number_states,sizeof(int));
      full_matrix=(double*)calloc(cns.total_number_states_sq,sizeof(double));

      eigenvalues=(double*)calloc(cns.total_number_states,sizeof(double));
      eigenvectors=(double*)calloc(cns.total_number_states_sq,sizeof(double));

      /*Determine Possible States To Form Hubbard Matrix*/
      init_states(possible_states,inverse_possible_states,ist,cns);

      /*Form the Hubbard Hamiltonian*/
      form_hubbard(full_matrix,possible_states,inverse_possible_states,neighbors,number_neighbors,ist,cns);

      /*Diagonalize the Hamiltonian Formed*/
      jacobireal(full_matrix,eigenvalues,eigenvectors,ist,cns);

      /*Copy Eigenvalues*/
      copy_eigenvalues(overall_eigenvalues,eigenvalues,total_overall_number_states,cns); 
      total_overall_number_states += cns.total_number_states; 

      free(possible_states); free(full_matrix);
      free(eigenvalues); free(eigenvectors);

     } /*down electrons*/
   } /*up electrons*/

   /*Determine Overall Energy*/
   determine_energy_finite_temperature(overall_eigenvalues,total_overall_number_states,cns,"energy.dat");  

  }

  /********************************************************************/  

  fclose(pout); 

  /*********************************************************************/

  free(neighbors); free(number_neighbors);
  free(electron_state); free(temporary_electron_state); 
  free(inverse_possible_states); 
  free(overall_eigenvalues);

  /********************************************************************/

return 0; 
} 

/***************************************************************************************/
