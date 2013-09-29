#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <malloc.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <assert.h>

/*NRC Definitions*/
#define TINY 1.0e-20
#define IM1          2147483563
#define IM2          2147483399
#define AM           (1.0 / (double)IM1)
#define IMM1         (IM1 - 1)
#define IA1          40014
#define IA2          40692
#define IQ1          53668
#define IQ2          52774
#define IR1          12211
#define IR2          3791
#define NTAB         32
#define NDIV         (1 + IMM1 / NTAB)
#define REPS         1.2e-7
#define RNMX         (1.0 - REPS)

/**************************************************************/

typedef struct st1 {
  double beta, beta_equilibration, dtau; 
  double U, t;
  double gamma;   
  double factor_spin_up_field_up, factor_spin_up_field_down; 
  double trial_energy; 
} cns_st;

typedef struct st2 {
  int n_sites_one, n_sites_two; 
  int n_sites, n_sites_sq; 
  int n_up, n_down;  
  int n_walkers; 
  int flag_trial;              /*Type of Trial Wavefunction - 0 If Free, 1 If UHF*/ 
  int flag_pbc;                /*Periodic Boundary Conditions - 0 If Open, 1 If Periodic*/ 
  int flag_cp; 
  int n_steps, n_steps_equilibration, n_steps_production;  
  int n_steps_energy;          /*How Often You Collect Energy Measurements*/
  int n_steps_orthogonalize;   /*How Often You Orthogonalize*/
  int n_steps_free;            /*How Many Steps to Propagate the Free Matrix To Get It To Ground State*/ 
} int_st;

/***************************************************************/

void Randomize();
double ran_nrc(long *);
double ran();

/**************************************************************/

void dzero_vec(double *vec, int size); 

/***************************************************************/

void copy_mat(double *mat_1, double *mat_2, int size); 
void mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void inverse(double *matrix, double *inverse, int size); 
double inverse_det(double *matrix, double *inverse, int size); 
void transpose_mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void jacobireal(double *matrix, double *eigenvalues, double *eigenvectors, int_st ist); 
void ludmp(double *matrix, int n, int *indx, double *d); 
void lubksb(double *matrix, int n, int *indx, double *b); 
void det(double *matrix, int size, double *determinant); 
void mat_transpose_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void modified_gram_schmidt(double *q, double *r, int size1, int size2); 

/**************************************************************/

void init(int_st *ist,cns_st *cns);
void init_kinetic(double *kinetic_full,double *kinetic_backwards,double *kinetic_forwards_half,double *kinetic_eigs,double *kinetic_eigvecs,int *neighbors,int *number_neighbors,int_st ist,cns_st cns); 
void init_kinetic_many(double *kinetic_many_full,double *kinetic_eigs,double *kinetic_eigvecs,int *neighbors,int *number_neighbors,int_st ist,cns_st cns); 
void init_wf(double *trial_wf_up,double *trial_wf_down,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns);
void init_walkers(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *weights,double *overlap_up,double *overlap_down,int_st ist); 
void init_neighbors(int *neighbors,int *number_neighbors,int_st ist); 

/**************************************************************/

void trial_identity(double *trial_wf_up,double *trial_wf_down,int_st ist);  
void trial_random(double *trial_wf_up,double *trial_wf_down,int_st ist); 
void trial_free(double *trial_wf_up,double *trial_wf_down,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns); 

/**************************************************************/

void neighborsperiodicboundary(int site, int *neighbors, int *number_neighbors, int_st ist); 
void neighborsopenboundary(int site, int *neighbors, int *numberneighbors, int_st ist);

/**************************************************************/

void print_dmat(double *mat, int size1, int size2, char *str); 
void print_dvec(double *vec, int size, char *str); 

/*************************************************************/

void compute_greens_function(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *greens_function_up,double *greens_function_down,int_st ist); 
void compute_energy(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int *neighbors,int *number_neighbors,int_st ist,cns_st cns,char *str,int step); 

/************************************************************/

void orthogonalize(double *wf_up,double *wf_down,double *overlap_up,double *overlap_down,double *weights,int_st ist); 
void orthogonalize_without_weights(double *wf_up,double *wf_down,int_st ist); 

/************************************************************/

double find_gamma(double U,double dtau); 
double determine_overlap_ratio(double *wf,double *trial_wf,double *overlap_inverse,int_st ist,cns_st cns,double spin_up,double field_up,int site); 
void propagate_forwards_potential(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,double *weights,int_st ist,cns_st cns,long *idum); 
void update_overlaps(double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,int_st ist); 
void propagate_wave_functions(double *wf_up,double *wf_down,int_st ist,cns_st cns,double field_up,int sites); 

/************************************************************/

void propagate_forwards_kinetic(double *kinetic_full,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_inverse_up,double *overlap_inverse_down,double *overlap_up,double *overlap_down,double *weights,int_st ist); 
void propagate_half_backwards_kinetic(double *kinetic_backwards_half,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_up,double *overlap_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int_st ist);
void propagate_half_forwards_kinetic(double *kinetic_forwards_half,double *wf_up,double *wf_down,double *trial_wf_up,double *trial_wf_down,double *overlap_up,double *overlap_down,double *overlap_inverse_up,double *overlap_inverse_down,double *weights,int_st ist);  
