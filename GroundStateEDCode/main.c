/*This is a Canonical Ensemble Code to Compute the Ground State Energies and Green's Functions for Electronic Systems to Check Against AFQMC*/

#include "ed.h"

/***************************************************************************************/

int main() {

  FILE *pout; 
  int i, j; 
  int *neighbors, *number_neighbors;
  int *electron_state, *temporary_electron_state; 
  int *possible_states, *inverse_possible_states; 
  double *full_matrix; 
  double *eigenvalues, *eigenvectors; 
  cns_st cns; int_st ist;

  /**********************************************************************/

  pout = fopen("errors.dat", "a+");

  /**********************************************************************/

  /*Obtain Parameters*/
  init(&ist,&cns);
  total_number_states(ist,&cns);
 
  /**********************************************************************/

  /*Initialize Relevant Vectors and Matrices*/
  neighbors=(int*)calloc(ist.dimension*2*ist.n_sites,sizeof(int)); 
  number_neighbors=(int*)calloc(ist.n_sites,sizeof(int)); 

  electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 
  temporary_electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 

  full_matrix=(double*)calloc(cns.total_number_states_sq,sizeof(double)); 

  eigenvalues=(double*)calloc(cns.total_number_states,sizeof(double)); 
  eigenvectors=(double*)calloc(cns.total_number_states_sq,sizeof(double)); 

  inverse_possible_states=(int*)calloc(cns.total_number_possible_states,sizeof(int)); 
  possible_states=(int*)calloc(cns.total_number_states,sizeof(int)); 

  /*********************************************************************/ 

  /*Get the Lattice Site Neighbors*/
  init_neighbors(neighbors,number_neighbors,ist); 

  fprintf(pout, "before init states\n"); fflush(pout); 

  /*Determine Possible States To Form Hubbard Matrix*/
  init_states(possible_states,inverse_possible_states,ist,cns);  

   fprintf(pout, "after init states\n"); fflush(pout); 

  /*Form the Hubbard Hamiltonian*/
  form_hubbard(full_matrix,possible_states,inverse_possible_states,neighbors,number_neighbors,ist,cns);  

  print_dmat(full_matrix,cns.total_number_states,cns.total_number_states,"errors.dat"); 

  fprintf(pout, "form hubbard after \n"); fflush(pout); 

  /*Diagonalize the Hamiltonian Formed*/
  jacobireal(full_matrix,eigenvalues,eigenvectors,ist,cns); 

  print_dvec(eigenvalues,cns.total_number_states,"errors.dat"); 

  /*Determine the Energy and Print It*/
  determine_energy(eigenvalues,ist,cns,"energy.dat"); 

  /********************************************************************/  

  fclose(pout); 

  /*********************************************************************/

  free(neighbors); free(number_neighbors);
  free(electron_state); free(temporary_electron_state); 
  free(possible_states); free(inverse_possible_states);  
  free(eigenvalues); free(eigenvectors); 

  /********************************************************************/

return 0; 
} 

/***************************************************************************************/
