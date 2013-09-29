#include "ed.h"

/********************************************************************************************************/

void jacobireal(double *matrix, double *eigenvalues, double *eigenvectors, int_st ist, cns_st cns) {
     /*Find Eigenvalues of Input Matrix*/
 
     int j,iq,ip,i;
     double tresh,theta,tau,t,sm,s,h,g,c; 
     double *b,*z;
  
     b=(double *)calloc(cns.total_number_states,sizeof(double)); 
     z=(double *)calloc(cns.total_number_states,sizeof(double)); 

     dzero_vec(eigenvalues,cns.total_number_states); 
     dzero_vec(eigenvectors,cns.total_number_states_sq); 

     /*Numerical Recipes Jacobi Routine*/
     for (ip=1;ip<=cns.total_number_states;ip++) {
                for (iq=1;iq<=cns.total_number_states;iq++) eigenvectors[(ip-1)*cns.total_number_states+(iq-1)]=0.0;
                eigenvectors[(ip-1)*cns.total_number_states+(ip-1)]=1.0;
        }
        for (ip=1;ip<=cns.total_number_states;ip++) {
                b[ip-1]=eigenvalues[ip-1]=matrix[(ip-1)*cns.total_number_states+(ip-1)];
                z[ip-1]=0.0;
        }
        for (i=1;i<200;i++) {
                sm=0.0;
                for (ip=1;ip<=cns.total_number_states-1;ip++) {
                        for (iq=ip+1;iq<=cns.total_number_states;iq++)
                                sm += fabs(matrix[(ip-1)*cns.total_number_states+(iq-1)]);
                }
                if (sm == 0.0) {
                        free(b); 
                        free(z);
                        return;
                }
                if (i < 4)
                        tresh=0.2*sm/((float) cns.total_number_states_sq);
                else
                        tresh=0.0;
                for (ip=1;ip<=cns.total_number_states-1;ip++) {
                        for (iq=ip+1;iq<=cns.total_number_states;iq++) {
                                g=100.0*fabs(matrix[(ip-1)*cns.total_number_states+iq-1]);
                                if (i > 4 && (float)(fabs(eigenvalues[ip-1])+g) == (float)fabs(eigenvalues[ip-1])
                                        && (float)(fabs(eigenvalues[iq-1])+g) == (float)fabs(eigenvalues[iq-1]))
                                        matrix[(ip-1)*cns.total_number_states+iq-1]=0.0;
                                else if (fabs(matrix[(ip-1)*cns.total_number_states+iq-1]) > tresh) {
                                        h=eigenvalues[iq-1]-eigenvalues[ip-1];
                                        if ((float)(fabs(h)+g) == (float)fabs(h))
                                                t=(matrix[(ip-1)*cns.total_number_states+(iq-1)])/h;
                                        else {
                                                theta=0.5*h/(matrix[(ip-1)*cns.total_number_states+iq-1]);
                                                t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                                                if (theta < 0.0) t = -t;
}
                                        c=1.0/sqrt(1+t*t);
                                        s=t*c;
                                        tau=s/(1.0+c);
                                        h=t*matrix[(ip-1)*cns.total_number_states+iq-1];
                                        z[ip-1] -= h;
                                        z[iq-1] += h;
                                        eigenvalues[ip-1] -= h;
                                        eigenvalues[iq-1] += h;
                                        matrix[(ip-1)*cns.total_number_states+iq-1]=0.0;
                                        for (j=1;j<=ip-1;j++) {
                                           g=matrix[(j-1)*cns.total_number_states+ip-1];
                                           h=matrix[(j-1)*cns.total_number_states+iq-1];
                                           matrix[(j-1)*cns.total_number_states+ip-1]=g-s*(h+g*tau);
                                           matrix[(j-1)*cns.total_number_states+iq-1]=h+s*(g-h*tau);
                                        }
                                        for (j=ip+1;j<=iq-1;j++) {
                                                g=matrix[(ip-1)*cns.total_number_states+j-1];
                                                h=matrix[(j-1)*cns.total_number_states+iq-1];
                                                matrix[(ip-1)*cns.total_number_states+j-1]=g-s*(h+g*tau);
                                                matrix[(j-1)*cns.total_number_states+iq-1]=h+s*(g-h*tau);
                                        }
                                        for (j=iq+1;j<=cns.total_number_states;j++) {
                                                g=matrix[(ip-1)*cns.total_number_states+j-1];
                                                h=matrix[(iq-1)*cns.total_number_states+j-1];
                                                matrix[(ip-1)*cns.total_number_states+j-1]=g-s*(h+g*tau);
                                                matrix[(iq-1)*cns.total_number_states+j-1]=h+s*(g-h*tau);
                                        }
                                        for (j=1;j<=cns.total_number_states;j++) {
                                                g=eigenvectors[(j-1)*cns.total_number_states+ip-1];
                                                h=eigenvectors[(j-1)*cns.total_number_states+iq-1];
                                                eigenvectors[(j-1)*cns.total_number_states+ip-1]=g-s*(h+g*tau);
                                                eigenvectors[(j-1)*cns.total_number_states+iq-1]=h+s*(g-h*tau);
                                        }
                                }
                        }
                }
                for (ip=1;ip<=cns.total_number_states;ip++) {
                        b[ip-1] += z[ip-1];
                        eigenvalues[ip-1]=b[ip-1];
                        z[ip-1]=0.0;
                }
        }


free(b); 
free(z); 
}

/***************************************************************************************************************/

void mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/

    int i, j, k; 
    dzero_vec(product,size_1*size_3);   

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0; 
        
        for (k=0; k<size_2; k++) {
          product[i*size_3+j] += mat_1[i*size_2+k] * mat_2[k*size_3+j]; 
        }

      }
    }

return; 
}

/***************************************************************************************************************/

void mat_transpose_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/
    /*Note That I am Not Taking the Size of the Second Matrix As the Second Size and Third Size*/ 

    int i, j, k;
    dzero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_2; j++) {
        product[i*size_2+j]=0.0;

        for (k=0; k<size_3; k++) {
          product[i*size_2+j] += mat_1[i*size_1+k] * mat_2[j*size_3+k];
        }

      }
    }

return;
}

/***************************************************************************************************************/

void transpose_mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

     /*Multiplies the Transpose of the First Matrix with the Second, As Needed to Compute Overlaps*/

     int i, j, k; 

     for (i=0; i<size_2; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0;

        for (k=0; k<size_1; k++) {
          product[i*size_3+j] += mat_1[k*size_2+i] * mat_2[k*size_3+j];
        }
      }
    }

return; 
}

/***************************************************************************************************************/

void copy_mat(double *mat_1, double *mat_2, int size) {

   int i; 

   for (i=0; i<size; i++) {
      mat_1[i] = mat_2[i]; 
   }

return; 
}

/***************************************************************************************************************/

void det(double *matrix, int size, double *determinant) {

    int i; 
    int *indx;   
   
    indx=(int *)calloc(size,sizeof(int)); 

    /*First Perform LU Decomposition*/
    ludmp(matrix,size,indx,determinant); 
    
    /*Now Multiply By Diagaonl*/
    for (i=0; i<size; i++) {
      (*determinant) *= matrix[i*size+i]; 
    }

free(indx); 
return; 
}

/***************************************************************************************************************/

void inverse(double *matrix, double *inverse, int size) {

     int i,j; 
     int *indx; 
     double *col; 
     double d;    

     indx=(int *)calloc(size,sizeof(int));
     col=(double *)calloc(size,sizeof(double));

     ludmp(matrix,size,indx,&d); 
     for (j=0; j<size; j++) {
        for (i=0; i<size; i++) {
           col[i]=0.0; 
        }
        col[j]=1.0; 

        lubksb(matrix,size,indx,col); 
        for(i=0; i<size; i++) {
          inverse[i*size+j]=col[i]; 
        }
    } 

free(col); 
free(indx); 
return; 
}

/***************************************************************************************************************/

double inverse_det(double *matrix, double *inverse, int size) {

  int i,j;
  int *indx;
  double *col;
  double determinant;

  /*Gets Inverse and Determinant*/ 

  indx=(int *)calloc(size,sizeof(int));
  col=(double *)calloc(size,sizeof(double));

  ludmp(matrix,size,indx,&determinant);
  for (j=0; j<size; j++) {
    determinant *= matrix[j*size+j];  

    for (i=0; i<size; i++) {
        col[i]=0.0;
     }
     col[j]=1.0;
     lubksb(matrix,size,indx,col);

     for(i=0; i<size; i++) {
       inverse[i*size+j]=col[i];
      }
    }

free(col); 
free(indx); 
return(determinant);  
}

/***************************************************************************************************************/

void ludmp(double *matrix, int n, int *indx, double *d) { 

   /*Performs an LU Decomposition Based on Numerical Recipes*/

   int i,imax,j,k; 
   double big,dum,sum,temp; 
   double *vv; 
 
   vv=(double*)calloc(n,sizeof(double)); 

   *d=1.0; 
   for (i=1; i<=n; i++){
     big=0.0;
     for (j=1; j<=n; j++) {
        if ((temp==fabs(matrix[(i-1)*n+j-1])) > big) big=temp; 
     }
     vv[i-1]=1.0/big; 
    }

    for (j=1; j<=n; j++) {
      for (i=1; i<j; i++) {
         sum=matrix[(i-1)*n+j-1]; 
         for (k=1; k<i; k++) { 
           sum -= matrix[(i-1)*n+k-1] * matrix[(k-1)*n+j-1]; 
         }
         matrix[(i-1)*n+j-1]=sum; 
       }
    

    big=0.0; 
    for (i=j; i<=n; i++) {
       sum=matrix[(i-1)*n+j-1]; 
       for (k=1; k<j; k++) {
         sum -= matrix[(i-1)*n+k-1]*matrix[(k-1)*n+j-1]; 
       }
       matrix[(i-1)*n+j-1]=sum; 
      
       if ( (dum=vv[i-1]*fabs(sum)) >= big) {
          big = dum; 
          imax = i; 
       }
     }

     if ( j!= imax ) {
        for (k=1; k<=n; k++) {
           dum = matrix[(imax-1)*n+k-1]; 
           matrix[(imax-1)*n+k-1] = matrix[(j-1)*n+k-1]; 
           matrix[(j-1)*n+k-1]=dum; 
         }
         (*d) = -(*d); 
         vv[imax-1]=vv[j-1]; 
     }
     indx[j-1]=imax; 

     if (matrix[(j-1)*n+j-1] == 0.0 ) { 
         matrix[(j-1)*n+j-1] = TINY; 
     }

     if (j != n ) {
        dum=1.0/matrix[(j-1)*n+j-1]; 
        for (i=j+1; i<=n; i++) {
          matrix[(i-1)*n+j-1] *= dum; 
        }
      }

   } 

free(vv); 
return; 
}

/***************************************************************************************************************/

void lubksb(double *matrix, int n, int *indx, double *b) {

    /*Peforms LU Forward and Backward Substitution From NRC*/
    
    int i,ii=0,ip,j; 
    double sum; 

    for (i=1; i<=n; i++) {
       ip=indx[i-1]; 
       sum=b[ip-1]; 
       b[ip-1]=b[i-1]; 
      
        if (ii) 
              for (j=ii; j<=i-1; j++) sum -= matrix[(i-1)*n+j-1] * b[j-1]; 
        else if (sum) ii=i; 
        b[i-1]=sum; 
     }

     for (i=n; i>=1; i--) {
        sum=b[i-1]; 
        for (j=i+1; j<=n; j++) sum -= matrix[(i-1)*n+j-1]*b[j-1]; 
        b[i-1]=sum/matrix[(i-1)*n+i-1]; 
     }

return; 
}

/**************************************************************************************************************/

void modified_gram_schmidt(double *q, double *r, int size1, int size2) {

   /*The Modified Gram-Schmidt Routine Used to Stabilize Matrix Products*/
   int i, j, k;   
   double temporary, anorm; 
   double *d; 

   d=(double *)calloc(size2,sizeof(double)); 
   dzero_vec(r,size2*size2); 

   for (i=1; i<=size2; i++) {

     temporary = 0.0;
     for (j=1; j<=size1; j++) {
       temporary += q[(j-1)*size2+i-1]*q[(j-1)*size2+i-1]; 
     }
     d[i-1]=sqrt(temporary); 
     anorm = 1.0/d[i-1];   

     for (j=1; j<=size1; j++) {
       q[(j-1)*size2+i-1] *= anorm; 
     }

     for ( j=i+1; j<=size2; j++) {
       temporary = 0.0; 
       for (k=0; k<=size1; k++) {
         temporary += q[(k-1)*size2+i-1] * q[(k-1)*size2+j-1]; 
       }
 
       for (k=1; k<=size1; k++) {
         q[(k-1)*size2+j-1]-= temporary * q[(k-1)*size2+i-1]; 
       }
       r[(i-1)*size2+j-1]=temporary*anorm;  
      }

   } /*i loop*/ 

   /*Now Make V Into R Matrix*/
   for (i=0; i<size2; i++) {
    r[i*size2+i]=d[i]; 
   } 
   
   free(d); 

return; 
}

/*********************************************************************************************************************/
