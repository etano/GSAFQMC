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

/*********************************************************************************************/

typedef struct st1 {
  double U, t;
  double beta; 
  int total_number_states, total_number_possible_states;  
  int total_number_states_sq; 
} cns_st;

typedef struct st2 {
  int n_sites_one, n_sites_two, n_sites_three; 
  int dimension; 
  int n_sites, n_sites_sq;
  int n_up, n_down;
  int states_per_site; 
  int flag_pbc;                /*Periodic Boundary Conditions - 0 If Open, 1 If Periodic*/
  int flag_spin_polarized; 
  int flag_ground_state; 
} int_st;

/*******************************************************************************************/

void init(int_st *ist,cns_st *cns); 
void init_neighbors(int *neighbors,int *number_neighbors,int_st ist); 

/******************************************************************************************/

void form_kinetic(double *full_matrix, int *possible_states, int *inverse_possible_states, int *neighbors, int *number_neighbors, int_st ist, cns_st cns); 
void form_hubbard(double *full_matrix, int *possible_states, int *inverse_possible_states, int *neighbors, int *number_neighbors, int_st ist, cns_st cns); 
void form_potential(double *full_matrix,int *possible_state,int_st ist,cns_st cns); 
double sign(double current_state, int first, int second, int updown,int *electron_state,int_st ist); 

/******************************************************************************************/

void neighborsperiodicboundary(int site, int *neighbors, int *number_neighbors,int_st ist); 
void neighborsopenboundary(int site, int *neighbors, int *number_neighbors,int_st ist); 

/******************************************************************************************/

void dzero_vec(double *vec, int size); 

/******************************************************************************************/

void print_dvec(double *vec, int size, char *str); 
void print_ivec(int *vec, int size, char *str); 
void print_dmat(double *mat, int size1, int size2, char *str); 

/******************************************************************************************/

void determine_energy_finite_temperature(double *overall_eigenvalues,int total_overall_number_states,cns_st cns,char *str1); 
void determine_energy_ground_state(double *eigenvalues,int_st ist,cns_st cns,char *str1); 
void copy_eigenvalues(double *overall_eigenvalues,double *eigenvalues,int number_eigenvalues, cns_st cns); 

/******************************************************************************************/

void convert_electron_state_to_lattice(double state, int *electron_state, int_st ist); 
int convert_lattice_to_electron_state(int *electron_state, int_st ist); 

/******************************************************************************************/

void jacobireal(double *matrix, double *eigenvalues, double *eigenvectors, int_st ist, cns_st cns); 
void mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void mat_transpose_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void transpose_mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3); 
void copy_mat(double *mat_1, double *mat_2, int size); 
void det(double *matrix, int size, double *determinant); 
void inverse(double *matrix, double *inverse, int size); 
double inverse_det(double *matrix, double *inverse, int size); 
void ludmp(double *matrix, int n, int *indx, double *d); 
void lubksb(double *matrix, int n, int *indx, double *b); 
void modified_gram_schmidt(double *q, double *r, int size1, int size2); 

/*****************************************************************************************/

int possible_state(double state, int_st ist); 
void total_number_states(int_st ist, cns_st *cns); 
void init_states(int *possible_states, int *inverse_possible_states, int_st ist, cns_st cst); 

/*****************************************************************************************/ 
