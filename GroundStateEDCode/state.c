#include "ed.h"

/********************************************************************************/

void convert_electron_state_to_lattice(double state, int *electron_state, int_st ist) {
    /*Converts the Number of the State to a Real Space Lattice Equivalent*/

    int remainder, divided, i;

    /*Convert Electron States*/
    for (i=ist.n_sites-1; i>-1; i--) {
       divided=(int)(state/((double) pow(ist.states_per_site, i)));
       remainder=(int)(state-(int) divided*pow(ist.states_per_site, i));
       state=remainder;
       electron_state[i]=divided;
     }

return; 
}

/********************************************************************************/

int convert_lattice_to_electron_state(int *electron_state, int_st ist) {
  /*Take an Electron State in Real Space and Converts It to the Electron State Number*/

  double state=0;
  int i;

  /*Go Through All Sites and Add Up State Number in Base of Possible Spin Configurations */
  for (i=0; i<ist.n_sites; i++) {
    switch(electron_state[i]) {
         case 0:
          break;
         case 1: {
           state += pow(ist.states_per_site, i);
           break;
         }
         case 2: {
           state += 2*pow(ist.states_per_site, i);
           break;
         }
         case 3: {
           state += 3*pow(ist.states_per_site, i);
           break;
         }
      }
    }

   state = (int) state; 

return(state);
}

/********************************************************************************/

int possible_state(double state, int_st ist) {

    /*Determines If Given State Is a Possible State Given Fixed Number of Electrons*/
 
    int i; 
    int number_electrons_up = 0, number_electrons_down = 0;   
    int *temporary_electron_state; 

    temporary_electron_state=(int*)calloc(ist.n_sites,sizeof(int)); 

    /*Convert To State*/
    convert_electron_state_to_lattice(state,temporary_electron_state,ist); 

    /*Figure Out How Many Electrons on Lattice in State*/
    if ( ist.flag_spin_polarized == 0 ) {

       for (i=0; i<ist.n_sites; i++) {
         switch(temporary_electron_state[i]) {
            case 0:
             break;
            case 1:
             number_electrons_up++; 
             break;
            case 2:
             number_electrons_down++; 
             break;
            case 3:
             number_electrons_up++; 
             number_electrons_down++; 
             break;
         }
         if (number_electrons_up > ist.n_up || number_electrons_down > ist.n_down) {
          return 0;
         }
       } /*Loop Through Sites*/

     }
     else {  /*If Spin Polarized******************************/

      for (i=0; i<ist.n_sites; i++) {
         if (temporary_electron_state[i]>1) {
          return 0;
         }

         switch(temporary_electron_state[i]) {
            case 0:
             break;
            case 1:
             number_electrons_up+=1;
             break;
          }
          if ( number_electrons_up > ist.n_up ) {
            return 0; 
          }

       } /*Loop Through Sites*/
       

    } /*If Spin Polarized******************************/

    /*Accept If Total Spin Is What It Should Be*/
    if ( number_electrons_up == ist.n_up && number_electrons_down == ist.n_down ) {
      print_ivec(temporary_electron_state,ist.n_sites,"errors.dat"); 
      return 1; 
    }

free(temporary_electron_state); 
return(0); 
}

/********************************************************************************/

void total_number_states(int_st ist, cns_st *cns) {
   /*Determines the Total Number of States for Fermions*/

   int i;
   int total_number_states = 0; 
   FILE *pf = fopen("ed-parameters.dat", "a+"); 

   for (i=0; i<cns->total_number_possible_states; i++) {
     if ( possible_state(i,ist) ) {    /*If of 0, means counts all possible spin states*/
       total_number_states+=1;
     }
   }
   cns->total_number_states = total_number_states; 
   cns->total_number_states_sq = cns->total_number_states * cns->total_number_states; 

   fprintf(pf, "Total Number States = %d\n", cns->total_number_states); fflush(pf); 

fclose(pf); 
return; 
}

/******************************************************************************/

void init_states(int *possible_states, int *inverse_possible_states, int_st ist, cns_st cns) {

   /*Maps Every Actual State Onto a Number and Vice Versa*/

   int i, j; 

   i=0;
   for (j=0; j<cns.total_number_possible_states; j++) {
     inverse_possible_states[j] = -1; 
     if (possible_state(j,ist)==1 ) {
       possible_states[i]=j;
       inverse_possible_states[j]=i;
       i++;
     }
   }

return; 
}

/*******************************************************************************/
