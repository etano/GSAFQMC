#include "ed.h"

/**************************************************************************/

void determine_energy_ground_state(double *eigenvalues,int_st ist,cns_st cns,char *str1) {
    /*Determines the Ground State Energy of a Fermionic System*/

    int i;
    double mineigenvalue=3*pow(10, 12);
    FILE *pf = fopen(str1, "a+"); 

    /*Run Through Eigenvalues to Find Lowest*/
    for (i=0; i<cns.total_number_states; i++) {
       if (eigenvalues[i]<mineigenvalue) {
         mineigenvalue=eigenvalues[i];
       }
    }

    /*Print Lowest Energy Eigenvalue*/
    fprintf(pf, "%f\t %f\t %d\t %d\t %f\n", cns.U, cns.t, ist.n_up, ist.n_down, mineigenvalue); fflush(pf);

fclose(pf); 
}

/**************************************************************************/

void determine_energy_finite_temperature(double *overall_eigenvalues,int total_overall_number_states,cns_st cns,char *str1) {
    /*Determines the Finite Temperature Energy*/

    int i; 
    double numerator = 0.0, denominator = 0.0; 
    FILE *pf = fopen(str1, "a+"); 

    for (i=0; i <total_overall_number_states; i++) {
      numerator += overall_eigenvalues[i] * exp(-cns.beta*overall_eigenvalues[i]); 
      denominator += exp(-cns.beta*overall_eigenvalues[i]); 
      fprintf(pf, "eigs %f\n", overall_eigenvalues[i]); fflush(pf); 
    }  
    fprintf(pf, "\n\n %f\t %f\t %f\t %f\n", cns.U, cns.t, cns.beta, numerator/denominator); fflush(pf);  

fclose(pf); 
return; 
}

/**************************************************************************/

void copy_eigenvalues(double *overall_eigenvalues,double *eigenvalues,int number_eigenvalues, cns_st cns) {
   /*Copy the Eigenvalues Over*/

   int i; 

   for (i=0; i<cns.total_number_states; i++) { 
     overall_eigenvalues[i+number_eigenvalues] = eigenvalues[i]; 
   }       
  
return; 
} 

/**************************************************************************/
