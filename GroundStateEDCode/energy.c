#include "ed.h"

/**************************************************************************/

void determine_energy(double *eigenvalues,int_st ist,cns_st cns,char *str1) {
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
    fprintf(pf, "%f\t %f\t %f\n", cns.U, cns.t, mineigenvalue); fflush(pf);

fclose(pf); 
}

/**************************************************************************/
