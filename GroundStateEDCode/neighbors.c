#include "ed.h"

/*Determines the Neighbors of a Given Lattice Site*/

/***************************************************************************************/

void neighborsopenboundary(int site, int *neighbors, int *number_neighbors,int_st ist) {
  /*Obtains the 2*dimensions nearest-neighbors directly surrounding a site*/
  /*Should be Done More Generally, But This Works for the Moment!!!*/

  int nsites_squared_one_two=ist.n_sites_one*ist.n_sites_two; 
  int nsites_squared_one_one=ist.n_sites_one*ist.n_sites_one; 
  int positionthirdside=(int)(site/nsites_squared_one_two);
  int positionsecondside=(int)(site-positionthirdside*nsites_squared_one_two)/ist.n_sites_one;
  int positionfirstside=(int)(site-positionthirdside*nsites_squared_one_two-positionsecondside*ist.n_sites_one);
  int row_number = 2 * ist.dimension; 

   /*Ordering Such That Double Counting Can Be Avoided In Certain Cases*/
   if (ist.dimension == 1) { /*********************************************************************/
     if (site!=ist.n_sites_one-1 && site!=0) {
       neighbors[site*row_number]=site-1;
       neighbors[site*row_number+1]=site+1;
       number_neighbors[site]=2;
     }
     else if (site==ist.n_sites_one-1){
       neighbors[site*row_number]=site-1;
       number_neighbors[site]=1;
     }
     else if (site==0) {
       neighbors[site*row_number]=site+1;
       number_neighbors[site]=1;
     }
   }
   else if (ist.dimension == 2 ) { /**************************************************************/ 
     if (ist.n_sites_one==ist.n_sites_two) {
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one; 
          number_neighbors[site]=3;
        }
        else if (site==0) {
           neighbors[row_number*site]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
        else if (site==ist.n_sites_one-1) {
           neighbors[row_number*site]=site-1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
       }
       else if (site<nsites_squared_one_one && site>=nsites_squared_one_one-ist.n_sites_one) {
         if (site!=nsites_squared_one_one-ist.n_sites_one && site!=nsites_squared_one_one-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=3;
         }
         else if (site==nsites_squared_one_one-ist.n_sites_one) {
          neighbors[row_number*site]=site+1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
        else if (site==nsites_squared_one_one-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
      }
      else if (site%ist.n_sites_one==0) {
         if (site!=0 && site!=nsites_squared_one_one-ist.n_sites_one) {
           neighbors[row_number*site+2]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site]=site-ist.n_sites_one;
           number_neighbors[site]=3;
         }
      }
      else if (site%ist.n_sites_one==ist.n_sites_one-1) {
         if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_one-1) {
            neighbors[row_number*site]=site-1;
            neighbors[row_number*site+1]=site+ist.n_sites_one;
            neighbors[row_number*site+2]=site-ist.n_sites_one;
            number_neighbors[site]=3;
       }
      }
      else {
           neighbors[row_number*site]=site-1;
           neighbors[row_number*site+2]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
      }
     } /*If LengthSideOne==LengthSideTwo*/
     else {  /*I Assume Length Side One Is Horizontal and LengthSide Two Is Vertical*/
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          number_neighbors[site]=3;
        }
        else if (site==0) {
           neighbors[row_number*site]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
        else if (site==ist.n_sites_one-1) {
           neighbors[site*ist.n_sites_one]=site-1;
           neighbors[site*ist.n_sites_one+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
       }
       else if (site<nsites_squared_one_two && site>=nsites_squared_one_two-ist.n_sites_one) {
        if (site!=nsites_squared_one_two-ist.n_sites_one && site!=nsites_squared_one_two-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=3;
        }
        else if (site==nsites_squared_one_two-ist.n_sites_one) {
          neighbors[row_number*site]=site+1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
        else if (site==nsites_squared_one_two-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
       }
       else if (site%ist.n_sites_one==0) {
         if (site!=0 && site!=nsites_squared_one_two-ist.n_sites_one) {
           neighbors[row_number*site+2]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site]=site-ist.n_sites_one;
           number_neighbors[site]=3;
         }
       }
       else if (site%ist.n_sites_one==ist.n_sites_one-1) {
         if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_two-1) {
            neighbors[row_number*site]=site-1;
            neighbors[row_number*site+1]=site+ist.n_sites_one;
            neighbors[row_number*site+2]=site-ist.n_sites_one;
            number_neighbors[site]=3;
        }
       }
       else {
           neighbors[row_number*site]=site-1;
           neighbors[row_number*site+2]=site+1;
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
       }
    } /*If Sides Are Unequal LengthSideOne!=LengthSideTwo*/

  }/*else if two dimensions*/
  else { /*If Three Dimensions*********************************************/

    /*If In Center of Cube*/
    if (positionfirstside!=0 && positionfirstside!=(ist.n_sites_one-1) && positionsecondside!=0 && positionsecondside!=(ist.n_sites_two-1) && positionthirdside!=0 && positionthirdside!=(ist.n_sites_three-1)) {
      neighbors[site*row_number]=site-1;
      neighbors[site*row_number+1]=site+1;
      neighbors[site*row_number+2]=site-ist.n_sites_one;
      neighbors[site*row_number+3]=site+ist.n_sites_one;
      neighbors[site*row_number+4]=site-nsites_squared_one_two; 
      neighbors[site*row_number+5]=site+nsites_squared_one_two; 
      number_neighbors[site]=6;
    }
    else if (positionfirstside==0) { /*If on Left Edge, But In Center*/
      if (positionsecondside!=0 && positionsecondside!=(ist.n_sites_two-1) && positionthirdside!=0 && positionthirdside!=(ist.n_sites_three-1)){ /*If Not On Any Other Edges*/
        neighbors[site*row_number+1]=site+1;
        neighbors[site*row_number+2]=site-ist.n_sites_one;
        neighbors[site*row_number+3]=site+ist.n_sites_one;
        neighbors[site*row_number+4]=site-nsites_squared_one_two;
        neighbors[site*row_number]=site+nsites_squared_one_two;
        number_neighbors[site]=5;
      }
      else if (positionsecondside==0) {
        if (positionthirdside==0) { /*Front Corner*/
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number]=site+ist.n_sites_one;
          neighbors[site*row_number+2]=site+nsites_squared_one_two;
          number_neighbors[site]=3;
        }
        else if (positionthirdside==ist.n_sites_three-1) {  /*Back Corner*/
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number]=site+ist.n_sites_one;
          neighbors[site*row_number+2]=site-nsites_squared_one_two;
          number_neighbors[site]=3;
        }
        else {
          neighbors[site*row_number]=site+1;
          neighbors[site*row_number+1]=site+ist.n_sites_one;
          neighbors[site*row_number+2]=site-nsites_squared_one_two; 
          neighbors[site*row_number+3]=site+nsites_squared_one_two;
          number_neighbors[site]=4;
        }
      }
      else if (positionsecondside==ist.n_sites_two-1) { /*So Top Left Edge*/
          if (positionthirdside==0) { /*Front Top Corner*/
            neighbors[site*row_number]=site+1;
            neighbors[site*row_number+1]=site-ist.n_sites_one;
            neighbors[site*row_number+2]=site+nsites_squared_one_two;
            number_neighbors[site]=3;
          }
          else if (positionthirdside==ist.n_sites_three-1) { /*Back Top Corner*/
            neighbors[site*row_number+1]=site+1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number]=site-nsites_squared_one_two;
            number_neighbors[site]=3;
          }
          else {
            neighbors[site*row_number+1]=site+1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number]=site-nsites_squared_one_two;
            neighbors[site*row_number+3]=site+nsites_squared_one_two;
            number_neighbors[site]=4;
          }
       }
       else if (positionthirdside==0) { /*Upfront Edge*/
           neighbors[site*row_number]=site+1;
           neighbors[site*row_number+1]=site-ist.n_sites_one;
           neighbors[site*row_number+2]=site+ist.n_sites_one;
           neighbors[site*row_number+3]=site+nsites_squared_one_two;
           number_neighbors[site]=4;
       }
       else { /*Back Side Edge*/
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number]=site-nsites_squared_one_two;
           number_neighbors[site]=4;
       }
     } /*If Position 1 Is Zero*/
     else if (positionfirstside==ist.n_sites_one-1) { /*If Right Side Now*/
       if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) {  /*So Not At Other Edges*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+1]=site+nsites_squared_one_two;
           number_neighbors[site]=5;
       }
       else if (positionsecondside==0) { /*Bottom Right Edge*/
         if (positionthirdside==0) {   /*Bottom, Right Corner*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+ist.n_sites_one;
           neighbors[site*row_number+2]=site+nsites_squared_one_two;
           number_neighbors[site]=3;
         }
         else if (positionthirdside==ist.n_sites_three-1) { /*Bottom, Back, Right Corner*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+ist.n_sites_one;
           neighbors[site*row_number+2]=site-nsites_squared_one_two;
           number_neighbors[site]=3;
         }
         else { /*If Back Edge*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+1]=site-nsites_squared_one_two;
           neighbors[site*row_number+2]=site+nsites_squared_one_two;
           number_neighbors[site]=4;
         }
        }
        else if (positionsecondside==ist.n_sites_two-1) { /*Top Right Edge*/
          if (positionthirdside==0) { /*Top Right, Front Corner*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+1]=site+nsites_squared_one_two;
            number_neighbors[site]=3;
          }
          else if (positionthirdside==ist.n_sites_three-1) { /*Top, Right, Back Corner*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+2]=site-(ist.n_sites_one);
            neighbors[site*row_number+1]=site-nsites_squared_one_two;
            number_neighbors[site]=3;
          }
          else {   /*If Top, Right Edge*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+1]=site-nsites_squared_one_two;
            neighbors[site*row_number+3]=site+nsites_squared_one_two;
            number_neighbors[site]=4;
          }
        }
        else { /*On Front or Back Edges*/
          if (positionthirdside==0) { /*If  Front, Right Edge*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+1]=site+ist.n_sites_one;
            neighbors[site*row_number+3]=site+nsites_squared_one_two;
            number_neighbors[site]=4;
          }
          else if (positionthirdside==ist.n_sites_three-1) {
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site+ist.n_sites_one;
            neighbors[site*row_number+1]=site-nsites_squared_one_two;
            number_neighbors[site]=4;
          }
         }
      } /*If On Right Side of Cube*/
      else if (positionsecondside==0) { /*If On Bottom of Cube*/
        if (positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two; 
           neighbors[site*row_number+2]=site+nsites_squared_one_two;
           number_neighbors[site]=5;
        }
        else if (positionthirdside==0) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+3]=site+ist.n_sites_one; 
           neighbors[site*row_number+2]=site+nsites_squared_one_two;
           number_neighbors[site]=4;
        }
        else if (positionthirdside==ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+2]=site-nsites_squared_one_two;
           number_neighbors[site]=4;
        }
      }
      else if (positionsecondside==ist.n_sites_two-1) { /*If On Top of Cube*/
        if (positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) { /*Top of Cube But Not On Edges*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two; 
           neighbors[site*row_number+3]=site+nsites_squared_one_two; 
           number_neighbors[site]=5;
        }
        else if (positionthirdside==0) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+nsites_squared_one_two;
           number_neighbors[site]=4;
        }
        else if (positionthirdside==ist.n_sites_three-1) {
          neighbors[site*row_number]=site-1;
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number+2]=site-ist.n_sites_one;
          neighbors[site*row_number+3]=site-nsites_squared_one_two;
          number_neighbors[site]=4;
        }
      }
      else if (positionthirdside==0) { /*If Front Face*/
        if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionfirstside!=0 && positionfirstside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site+nsites_squared_one_two;
           number_neighbors[site]=5;
        }
      }
      else if (positionthirdside==ist.n_sites_three-1) { /*If Back Face*/
         if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionfirstside!=0 && positionfirstside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           number_neighbors[site]=5;
         }
       }
   } /*If Three Dimensions*/

return; 
}

/*********************************************************************************************************/

void neighborsperiodicboundary(int site, int *neighbors, int *number_neighbors,int_st ist) {
  /*Obtains the 2*dimensions nearest-neighbors directly surrounding a site*/
  /*Should be Done More Generally, But This Works for the Moment!!!*/

  int nsites_squared_one_one = ist.n_sites_one * ist.n_sites_one; 
  int nsites_squared_one_two = ist.n_sites_one * ist.n_sites_two; 
  int nsites_cubed = ist.n_sites_one*ist.n_sites_two*(ist.n_sites_three-1);
  int row_number = ist.dimension * 2;  
  int positionthirdside=(int)(site/nsites_squared_one_two);
  int positionsecondside=(int)(site-positionthirdside*nsites_squared_one_two)/ist.n_sites_one;
  int positionfirstside=(int)(site-positionthirdside*nsites_squared_one_two-positionsecondside*ist.n_sites_one);

   /*Ordering Such That Double Counting Can Be Avoided In Certain Cases*/
  if ( ist.dimension == 1 ) {
    if (site!=ist.n_sites_one-1 && site!=0) {
      neighbors[site*row_number]=site-1;
      neighbors[site*row_number+1]=site+1;
      number_neighbors[site]=2;
    }
    else if (site==ist.n_sites_one-1){
      neighbors[site*row_number]=site-1;
      neighbors[site*row_number+1]=0;
      number_neighbors[site]=2;
    }
    else if (site==0) {
      neighbors[site*row_number]=ist.n_sites_one-1;
      neighbors[site*row_number+1]=site+1;
      number_neighbors[site]=2;
    }
  }
  else if ( ist.dimension == 2 )  { /*If 2 Dimensions*/ 

    if (ist.n_sites_one==ist.n_sites_two) {
      if (site<ist.n_sites_one) {
       if (site!=0 && site!=ist.n_sites_one-1){
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site+ist.n_sites_one;
         neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
         number_neighbors[site]=4;
       }
       else if (site==0) {
          neighbors[row_number*site]=ist.n_sites_one-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
          number_neighbors[site]=4;
       }
       else if (site==ist.n_sites_one-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=0;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
          number_neighbors[site]=4;
       }
      }
      else if (site<nsites_squared_one_one && site>=nsites_squared_one_one-ist.n_sites_one) {
        if (site!=nsites_squared_one_one-ist.n_sites_one && site!=nsites_squared_one_one-1) {
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
        }
        else if (site==nsites_squared_one_one-ist.n_sites_one) {
         neighbors[row_number*site]=nsites_squared_one_one-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_one-1) {
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=nsites_squared_one_one-ist.n_sites_one;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
     }
     else if (site%ist.n_sites_one==0) {
        if (site!=0 && site!=nsites_squared_one_one-ist.n_sites_one) {
          neighbors[row_number*site]=site+(ist.n_sites_one-1);
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
        }
     }
     else if (site%ist.n_sites_one==ist.n_sites_one-1) {
        if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_one-1) {
           neighbors[row_number*site]=site-1;
           neighbors[row_number*site+2]=site-(ist.n_sites_one-1);
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
      }
     }
     else {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
     }
    } /*If LengthSideOne==LengthSideTwo*/
    else {  /*I Assume Length Side One Is Horizontal and LengthSide Two Is Vertical*/
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site+ist.n_sites_one;
         neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
         number_neighbors[site]=4;
       }
       else if (site==0) {
          neighbors[row_number*site]=ist.n_sites_one-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
          number_neighbors[site]=4;
       }
       else if (site==ist.n_sites_one-1) {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=0;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
          number_neighbors[site]=4;
       }
      }
      else if (site<nsites_squared_one_two && site>=nsites_squared_one_two-ist.n_sites_one) {
       if (site!=nsites_squared_one_two-ist.n_sites_one && site!=nsites_squared_one_two-1) {
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_two-ist.n_sites_one) {
         neighbors[row_number*site]=nsites_squared_one_two-1;
         neighbors[row_number*site+2]=site+1;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one; 
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_two-1) {
         neighbors[row_number*site]=site-1;
         neighbors[row_number*site+2]=nsites_squared_one_two-ist.n_sites_one;
         neighbors[row_number*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[row_number*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
      }
      else if (site%ist.n_sites_one==0) {
        if (site!=0 && site!=nsites_squared_one_two-ist.n_sites_one) {
          neighbors[row_number*site]=site+(ist.n_sites_one-1);
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
        }
      }
      else if (site%ist.n_sites_one==ist.n_sites_one-1) {
        if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_two-1) {
           neighbors[row_number*site]=site-1;
           neighbors[row_number*site+2]=site-(ist.n_sites_one-1);
           neighbors[row_number*site+1]=site+ist.n_sites_one;
           neighbors[row_number*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
       }
      }
      else {
          neighbors[row_number*site]=site-1;
          neighbors[row_number*site+2]=site+1;
          neighbors[row_number*site+1]=site+ist.n_sites_one;
          neighbors[row_number*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
      }
    } /*If Sides Are Unequal LengthSideOne!=LengthSideTwo*/
 
  } /*If 2 Dimensions*/ 
  else { /*If 3 Dimensions***************************************/

    /*If In Center of Cube*/
    if (positionfirstside!=0 && positionfirstside!=(ist.n_sites_one-1) && positionsecondside!=0 && positionsecondside!=(ist.n_sites_two-1) && positionthirdside!=0 && positionthirdside!=(ist.n_sites_three-1)) {
      neighbors[site*row_number]=site-1;
      neighbors[site*row_number+1]=site+1;
      neighbors[site*row_number+2]=site-ist.n_sites_one;
      neighbors[site*row_number+3]=site+ist.n_sites_one;
      neighbors[site*row_number+4]=site-nsites_squared_one_two;
      neighbors[site*row_number+5]=site+nsites_squared_one_two;
      number_neighbors[site]=6;
    }
    else if (positionfirstside==0) { /*If on Left Edge, But In Center*/
      if (positionsecondside!=0 && positionsecondside!=(ist.n_sites_two-1) && positionthirdside!=0 && positionthirdside!=(ist.n_sites_three-1)){ /*If Not On Any Other Edges*/
        neighbors[site*row_number]=site+(ist.n_sites_one-1);
        neighbors[site*row_number+1]=site+1;
        neighbors[site*row_number+2]=site-ist.n_sites_one;
        neighbors[site*row_number+3]=site+ist.n_sites_one;
        neighbors[site*row_number+4]=site-nsites_squared_one_two;
        neighbors[site*row_number+5]=site+nsites_squared_one_two;
        number_neighbors[site]=6;
      }
      else if (positionsecondside==0) {
        if (positionthirdside==0) { /*Front Corner*/
          neighbors[site*row_number]=site+(ist.n_sites_one-1);
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
          neighbors[site*row_number+3]=site+ist.n_sites_one;
          neighbors[site*row_number+4]=site+nsites_cubed;
          neighbors[site*row_number+5]=site+nsites_squared_one_two;
          number_neighbors[site]=6;
        }
        else if (positionthirdside==ist.n_sites_three-1) {  /*Back Corner*/
          neighbors[site*row_number]=site+(ist.n_sites_one-1);
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
          neighbors[site*row_number+3]=site+ist.n_sites_one;
          neighbors[site*row_number+4]=site-nsites_squared_one_two;
          neighbors[site*row_number+5]=site-nsites_cubed;
          number_neighbors[site]=6;
        }
        else {
          neighbors[site*row_number]=site+(ist.n_sites_one-1);
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
          neighbors[site*row_number+3]=site+ist.n_sites_one;
          neighbors[site*row_number+4]=site-nsites_squared_one_two;
          neighbors[site*row_number+5]=site+nsites_squared_one_two;
          number_neighbors[site]=6;
        }
      }
      else if (positionsecondside==ist.n_sites_two-1) { /*So Top Left Edge*/
          if (positionthirdside==0) { /*Front Top Corner*/
            neighbors[site*row_number]=site+(ist.n_sites_one-1);
            neighbors[site*row_number+1]=site+1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site-(ist.n_sites_two-1)*ist.n_sites_one;
            neighbors[site*row_number+4]=site+nsites_cubed; 
            neighbors[site*row_number+5]=site+nsites_squared_one_two;
            number_neighbors[site]=6;
          }
          else if (positionthirdside==ist.n_sites_three-1) { /*Back Top Corner*/
            neighbors[site*row_number]=site+(ist.n_sites_one-1);
            neighbors[site*row_number+1]=site+1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site-(ist.n_sites_two-1)*ist.n_sites_one;
            neighbors[site*row_number+4]=site-nsites_squared_one_two;
            neighbors[site*row_number+5]=site-nsites_cubed;
            number_neighbors[site]=6;
          }
          else {
            neighbors[site*row_number]=site+(ist.n_sites_one-1);
            neighbors[site*row_number+1]=site+1;
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site-(ist.n_sites_two-1)*ist.n_sites_one;
            neighbors[site*row_number+4]=site-nsites_squared_one_two;
            neighbors[site*row_number+5]=site+nsites_squared_one_two;
            number_neighbors[site]=6;
          }
       }
       else if (positionthirdside==0) { /*Upfront Edge*/
           neighbors[site*row_number]=site+(ist.n_sites_one-1);
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site+nsites_cubed;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
       }
       else { /*Back Side Edge*/
           neighbors[site*row_number]=site+(ist.n_sites_one-1);
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site-nsites_cubed;
           number_neighbors[site]=6;
       }
     } /*If Position 1 Is Zero*/
     else if (positionfirstside==ist.n_sites_one-1) { /*If Right Side Now*/
       if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) {  /*So Not At Other Edges*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
       }
       else if (positionsecondside==0) { /*Bottom Right Edge*/
         if (positionthirdside==0) {   /*Bottom, Right Corner*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
           neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site+nsites_cubed;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
         }
         else if (positionthirdside==ist.n_sites_three-1) { /*Bottom, Back, Right Corner*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
           neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site-nsites_cubed;
           number_neighbors[site]=6;
         }
         else { /*If Back Edge*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
           neighbors[site*row_number+2]=site+(ist.n_sites_two-1)*ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
         }
        }
        else if (positionsecondside==ist.n_sites_two-1) { /*Top Right Edge*/
          if (positionthirdside==0) { /*Top Right, Front Corner*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
            neighbors[site*row_number+4]=site+nsites_cubed;
            neighbors[site*row_number+5]=site+nsites_squared_one_two;
            number_neighbors[site]=6;
          }
          else if (positionthirdside==ist.n_sites_three-1) { /*Top, Right, Back Corner*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
            neighbors[site*row_number+2]=site-(ist.n_sites_one);
            neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
            neighbors[site*row_number+4]=site-nsites_squared_one_two;
            neighbors[site*row_number+5]=site-nsites_cubed;
            number_neighbors[site]=6;
          }
          else {   /*If Top, Right Edge*/
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
            neighbors[site*row_number+4]=site-nsites_squared_one_two;
            neighbors[site*row_number+5]=site+nsites_squared_one_two;
            number_neighbors[site]=6;
          }
        }
        else { /*On Front or Back Edges*/
          if (positionthirdside==0) { /*If  Front, Right Edge*/
            neighbors[site*row_number+0]=site-1;
            neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site+ist.n_sites_one;
            neighbors[site*row_number+4]=site+nsites_cubed;
            neighbors[site*row_number+5]=site+nsites_squared_one_two;
            number_neighbors[site]=6;
          }
          else if (positionthirdside==ist.n_sites_three-1) {
            neighbors[site*row_number]=site-1;
            neighbors[site*row_number+1]=site-(ist.n_sites_one-1);
            neighbors[site*row_number+2]=site-ist.n_sites_one;
            neighbors[site*row_number+3]=site+ist.n_sites_one;
            neighbors[site*row_number+4]=site-nsites_squared_one_two;
            neighbors[site*row_number+5]=site-nsites_cubed;
            number_neighbors[site]=6;
          }
         }
      } /*If On Right Side of Cube*/
      else if (positionsecondside==0) { /*If On Bottom of Cube*/
        if (positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site+ist.n_sites_one*(ist.n_sites_two-1);
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
        }
        else if (positionthirdside==0) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site+ist.n_sites_one*(ist.n_sites_two-1);
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site+nsites_cubed;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
        }
        else if (positionthirdside==ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site+ist.n_sites_one*(ist.n_sites_two-1);
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site-nsites_cubed;
           number_neighbors[site]=6;
        }
      }
      else if (positionsecondside==ist.n_sites_two-1) { /*If On Top of Cube*/
        if (positionthirdside!=0 && positionthirdside!=ist.n_sites_three-1) { /*Top of Cube But Not On Edges*/
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
        }
        else if (positionthirdside==0) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
           neighbors[site*row_number+4]=site+nsites_cubed;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
        }
        else if (positionthirdside==ist.n_sites_three-1) {
          neighbors[site*row_number]=site-1;
          neighbors[site*row_number+1]=site+1;
          neighbors[site*row_number+2]=site-ist.n_sites_one;
          neighbors[site*row_number+3]=site-ist.n_sites_one*(ist.n_sites_two-1);
          neighbors[site*row_number+4]=site-nsites_squared_one_two;
          neighbors[site*row_number+5]=site-nsites_cubed;
          number_neighbors[site]=6;
        }
      }
      else if (positionthirdside==0) { /*If Front Face*/
        if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionfirstside!=0 && positionfirstside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site+nsites_cubed;
           neighbors[site*row_number+5]=site+nsites_squared_one_two;
           number_neighbors[site]=6;
        }
      }
      else if (positionthirdside==ist.n_sites_three-1) { /*If Back Face*/
         if (positionsecondside!=0 && positionsecondside!=ist.n_sites_two-1 && positionfirstside!=0 && positionfirstside!=ist.n_sites_three-1) {
           neighbors[site*row_number]=site-1;
           neighbors[site*row_number+1]=site+1;
           neighbors[site*row_number+2]=site-ist.n_sites_one;
           neighbors[site*row_number+3]=site+ist.n_sites_one;
           neighbors[site*row_number+4]=site-nsites_squared_one_two;
           neighbors[site*row_number+5]=site-nsites_cubed;
           number_neighbors[site]=6;
         }
       }
   } /*If Three Dimensions*/

return; 
}
