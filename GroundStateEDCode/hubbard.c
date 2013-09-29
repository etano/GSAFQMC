#include "ed.h"

/**************************************************************************/

void form_hubbard(double *full_matrix, int *possible_states, int *inverse_possible_states, int *neighbors, int *number_neighbors, int_st ist, cns_st cns){
  /*Makes the Full Hubbard Matrix*/ 

  /*Zero Both Matrices*/
  dzero_vec(full_matrix,cns.total_number_states_sq); 

  /*First Form Kinetic Part*/
  if ( cns.total_number_states > 1 ) {
    form_kinetic(full_matrix,possible_states,inverse_possible_states,neighbors,number_neighbors,ist,cns); 
  }

  if ( cns.U != 0 ) { 
    /*Then Form Potential Part*/
    form_potential(full_matrix,possible_states,ist,cns);  
  }

return; 
}

/*************************************************************************/

void form_kinetic(double *full_matrix, int *possible_states, int *inverse_possible_states, int *neighbors, int *number_neighbors, int_st ist, cns_st cns){

  /*Forms Matrix for Kinetic Portion*/
  int i, j, k;
  int row_number = 2 * ist.dimension;   
  int current_state, proposed_state;
  int *temporary_electron_state; 

  temporary_electron_state = (int*)calloc(ist.n_sites,sizeof(int)); 

  for (i=0; i<cns.total_number_states; i++) {

    /*Converts the State to a Lattice Representation*/
    convert_electron_state_to_lattice(possible_states[i],temporary_electron_state,ist); 
    current_state = convert_lattice_to_electron_state(temporary_electron_state,ist); 

    /*Run Through All Sites, Trying to Move Spins to Other Sites*/
    for (j=0; j<ist.n_sites; j++) {
      if (temporary_electron_state[j]!=0) {

        /*Not the Most Efficient Way, But Go Through Possible Hops*/
        if (temporary_electron_state[j]==1) {
          /*Run Through Nearest Neighbors*/
          for (k=0; k<number_neighbors[j]; k++) {  /*No Double Counting Here for SMALL Lattice of Length Less than or equal to 2*/
            if (temporary_electron_state[neighbors[j*row_number+k]]==2) {

              /*Convert to Proposed State and Allow Hop*/
              current_state=convert_lattice_to_electron_state(temporary_electron_state,ist); 

              temporary_electron_state[neighbors[j*row_number+k]]=3;
              temporary_electron_state[j]=0;

              proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist); 
              if (possible_state(proposed_state,ist)==1) {
               full_matrix[i*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 1,temporary_electron_state,ist)*cns.t; 
              }

             convert_electron_state_to_lattice(current_state,temporary_electron_state,ist); 
             }
             else if (temporary_electron_state[neighbors[j*row_number+k]]==0) {
              current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
              temporary_electron_state[neighbors[j*row_number+k]]=1;
              temporary_electron_state[j]=0;

              proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
              if (possible_state(proposed_state,ist)==1) {
                full_matrix[i*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 1,temporary_electron_state,ist)*cns.t;
              }

              convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
            }
            else if (temporary_electron_state[neighbors[j*row_number+k]]==3) {
               current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               temporary_electron_state[neighbors[j*row_number+k]]=1;
               temporary_electron_state[j]=3;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)!=0 ) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 2,temporary_electron_state,ist)*cns.t; 
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
            }
           }
         }
         else if (temporary_electron_state[j]==2) {
           for (k=0; k<number_neighbors[j]; k++) {

             if (temporary_electron_state[neighbors[j*row_number+k]]==0) {
               current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);

               temporary_electron_state[neighbors[j*row_number+k]]=2;
               temporary_electron_state[j]=0;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)==1) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 2,temporary_electron_state,ist)*cns.t;
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
             }
             else if (temporary_electron_state[neighbors[j*row_number+k]]==1) {
                current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
                temporary_electron_state[neighbors[j*row_number+k]]=3;
                temporary_electron_state[j]=0;

                proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
                if (possible_state(proposed_state,ist)==1) {
                 full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 2,temporary_electron_state,ist)*cns.t; 
                }

                convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
             }
             else if (temporary_electron_state[neighbors[j*row_number+k]]==3) {
                current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
                temporary_electron_state[neighbors[j*row_number+k]]=2;
                temporary_electron_state[j]=3;

                proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
                if (possible_state(proposed_state,ist)==1) {
                  full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 1,temporary_electron_state,ist)*cns.t;
                }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
             }
           }
         }
         else if (temporary_electron_state[j]==3) {
           for (k=0; k<number_neighbors[j]; k++) {
             if (temporary_electron_state[neighbors[j*row_number+k]]==0) {
               /*Two Choices of Up and Down Here*/
               current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               temporary_electron_state[neighbors[j*row_number+k]]=1;
               temporary_electron_state[j]=2;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)==1) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 1,temporary_electron_state,ist)*cns.t;
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);

               temporary_electron_state[neighbors[j*row_number+k]]=2;
               temporary_electron_state[j]=1;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)==1) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 2,temporary_electron_state,ist)*cns.t;
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
             }
            else if (temporary_electron_state[neighbors[j*row_number+k]]==1) {
               current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               temporary_electron_state[neighbors[j*row_number+k]]=3;
               temporary_electron_state[j]=1;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)==1) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 2,temporary_electron_state,ist)*cns.t;
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
            }
            else if (temporary_electron_state[neighbors[j*row_number+k]]==2) {
               current_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               temporary_electron_state[neighbors[j*row_number+k]]=3;
               temporary_electron_state[j]=2;

               proposed_state=convert_lattice_to_electron_state(temporary_electron_state,ist);
               if (possible_state(proposed_state,ist)==1) {
                full_matrix[inverse_possible_states[current_state]*cns.total_number_states+inverse_possible_states[proposed_state]]=sign(current_state, neighbors[j*row_number+k], j, 1,temporary_electron_state,ist)*cns.t;
               }

               convert_electron_state_to_lattice(current_state,temporary_electron_state,ist);
            }
           }
         }

      } /*site not in zero state*/
    }  /*j loop*/

  } /*i loop through states*/

free(temporary_electron_state); 
return; 
}

/************************************************************************/

void form_potential(double *full_matrix,int *possible_state,int_st ist,cns_st cns){

   /*Forms the Potential Portion of the Hubbard Matrix*/
   int site, state;
   int *temporary_electron_state; 

   temporary_electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 

   /*If Spin-Polarized, There Is a Hubbard U*/
   if ( ist.flag_spin_polarized == 0 ) {

     /*Add U Terms for all of the Different States*/
     for (state=0; state<cns.total_number_states; state++) {

      /*First Convert All Sites to Real Space*/
      convert_electron_state_to_lattice(possible_state[state],temporary_electron_state,ist);  

      for (site=0; site<ist.n_sites; site++) {
        if (temporary_electron_state[site]==3) {
          full_matrix[state*cns.total_number_states+state]+=cns.U; 
        }
      }

     } /*Running Through States*/

   }
     
free(temporary_electron_state); 
return; 
}

/************************************************************************/

double sign(double current_state, int first, int second, int updown,int *electron_state,int_st ist) {
  /*This Function Determines What Sign the Term Has in the Hopping Depending Upon the Current and Proposed States*/
  /*If Electrons Hop Around Periodic Boundary Conditions in a Ring, Then You Change Sign*/

  int i;
  double ultimatesign;
  int numberspinsgreaterthanfirst=0, numberspinsgreaterthansecond=0;

  if ( ist.n_up>1 || ist.n_down>1) {

     /*First Obtain Current States*/
     convert_electron_state_to_lattice(current_state,electron_state,ist); 

     /*Determine Number of Spins Greater Than First or Second Creation/Destruction Operator Site*/
     for (i=0; i<ist.n_sites; i++) {
       if (i>first && (electron_state[i]==3 || electron_state[i]==updown) ) {
           numberspinsgreaterthanfirst++;
       }
       if (i>second && (electron_state[i]==3 || electron_state[i]==updown) ) {
           numberspinsgreaterthansecond++;
        }
     }

     if (first<second) {
        ultimatesign=pow(-1.0, numberspinsgreaterthanfirst+numberspinsgreaterthansecond+1);
        return ultimatesign;
     }
     else if (first > second) {
        ultimatesign=pow(-1.0, numberspinsgreaterthanfirst+numberspinsgreaterthansecond);
        return ultimatesign;
     }

 } /*If Current Number Electrons*/

return(-1.0);
}

/************************************************************************/ 
