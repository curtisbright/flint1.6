#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "flint.h"
#include "mpfr.h"
#include "F_mpz.h"
#include "mpz_mat.h"
#include "gmp.h"
#include "long_extras.h"
#include "mpq_mat.h"
#include "mpfr_mat.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL.h"
#include "memory-manager.h"
#include "test-support.h"

void _mpfr_vec_clean_scalar_product2(mpfr_t sp, __mpfr_struct * vec1, __mpfr_struct * vec2, int n, mp_prec_t prec)
{
   if (n <= 1){
      mpfr_mul(sp, vec1, vec2, GMP_RNDN);
      return;
   }
   mpfr_t tmp;
   mpfr_init2(tmp, prec);
  
   mpfr_mul(sp, vec1, vec2, GMP_RNDN);

   for (long i = 1; i < n; i++)
   {
      mpfr_mul(tmp, vec1 + i, vec2 + i, GMP_RNDN);
      mpfr_add(sp, sp, tmp, GMP_RNDN);
   }

   mpfr_clear(tmp);

   return;
}

void F_mpz_mat_RQ_factor_mgso(F_mpz_mat_t B, __mpfr_struct ** R, __mpfr_struct ** Q, long r, long c, mp_prec_t prec){

//Doing modified GSO will convert Q from B to Q row by row
   mpfr_t tmp;
   mpfr_init2(tmp, prec);

   long i, k, j;
   for (i = 0; i < r; i++)
	   _F_mpz_vec_to_mpfr_vec(Q[i], B->rows[i], c); 

   //iteration k should start with Q = q1, ..., q_(k-1), a_k', ... a_r' then convert a_k to q_k by a subtraction and the other a_j' should be modded by a_k  
   for (k = 0; k < r; k++){
      _mpfr_vec_norm2(R[k]+k, Q[k], c, prec);
      mpfr_sqrt(R[k]+k, R[k]+k, GMP_RNDN);
      for (i = 0; i < c; i++)
         mpfr_div(Q[k]+i, Q[k]+i, R[k]+k, GMP_RNDN);
      for (j = k+1; j < r; j++){
         _mpfr_vec_clean_scalar_product2(R[j]+k, Q[k], Q[j], c, prec);
         for (i = 0; i < c; i++){
            mpfr_mul(tmp, R[j]+k, Q[k]+i, GMP_RNDN);
            mpfr_sub(Q[j]+i, Q[j]+i, tmp, GMP_RNDN);
         }
      }
   }

   mpfr_clear(tmp);
   return;
}

int mpfr_mat_R_reduced(__mpfr_struct ** R, long d, double delta, double eta, mp_prec_t prec)
{

   if (d == 1)
      return 1;

   mpfr_t tmp1;
   mpfr_t tmp2;
   mpfr_init2(tmp1, prec);
   mpfr_init2(tmp2, prec);
   int reduced = 1;

   long i;
   for (i = 0; (i < d -1) && (reduced == 1); i++)
   {
      mpfr_pow_ui(tmp1, R[i+1] + i, 2L, GMP_RNDN);
      mpfr_pow_ui(tmp2, R[i+1] + i + 1, 2L, GMP_RNDN);
      mpfr_add(tmp1, tmp1, tmp2, GMP_RNDN);

      mpfr_pow_ui(tmp2, R[i] + i, 2L, GMP_RNDN);
      mpfr_mul_d(tmp2, tmp2, (double) delta, GMP_RNDN);

      mpfr_sub(tmp1, tmp1, tmp2, GMP_RNDN);
//      mpfr_add_d(tmp1, tmp1, .001, GMP_RNDN);
      if (mpfr_sgn(tmp1) < 0) 
      {
         reduced = 0;
         printf(" happened at index i = %ld\n", i);
         break;
      }
      long j;
      for (j = 0; (j <= i) && (reduced == 1); j++)
      {
         mpfr_mul_d(tmp2, R[i] + i, (double) eta, GMP_RNDN);
         fprintf(stderr, "i,j = %ld, %ld\n", i,j);
         if (mpfr_cmpabs(R[j] + i + 1, tmp2) > 0)
         {
            reduced = 0;
            printf(" size red problem at index i = %ld, j = %ld\n", i, j);
            break;
         }
      }
   }

   mpfr_clear(tmp1);
   mpfr_clear(tmp2);
   return reduced;
}

void F_mpz_mat_R_factor_house(F_mpz_mat_t B, __mpfr_struct ** squareR, __mpfr_struct ** R, long r, long c, mp_prec_t prec){

   mpfr_t tmp;
   mpfr_init2(tmp, prec);

   mpfr_t tmp2;
   mpfr_init2(tmp2, prec);

   __mpfr_struct ** V;

// V will store the r householder transforms with row j having c-j entries for j = 0... r-1
   V = mpfr_mat_init2(r, c, prec);

//for Householder QR we will only compute the R but from a rectangular R := B in the first iteration
   long i, k, j;
   for (i = 0; i < r; i++)
	   _F_mpz_vec_to_mpfr_vec(R[i], B->rows[i], c); 

   //iteration k should start with first k-1 rows already triangularized and diagonal = length of gs vector  
   for (i = 0; i < r; i++){
      //we have to update row i by applying the householder matrices from previous iterations, each matrix is stored as a vector v
      for (j = 0; j < i; j++){
         _mpfr_vec_clean_scalar_product2(tmp, R[i]+j, V[j], c - j, prec);
         for (k = 0; k < c - j; k++){
            mpfr_mul(tmp2, tmp, V[j]+k, GMP_RNDN);
            mpfr_sub(R[i] + j + k, R[i] + j + k, tmp2, GMP_RNDN);
         }
      }

      for (k = 0; k < c - i; k++){
         mpfr_set(V[i] + k, R[i] + i + k, GMP_RNDN);
      }

      _mpfr_vec_clean_scalar_product2(tmp, R[i] + i, R[i] + i, c - i, prec);
      mpfr_sqrt(tmp, tmp, GMP_RNDN);

      if (c - i > 1){
         int sgn = mpfr_sgn(R[i] + i);
         if (sgn < 0){
            mpfr_neg(tmp, tmp, GMP_RNDN);
         }
         //Now temp = +/-|| r || for the subvector we want to make (*,0,0,...,0)
         _mpfr_vec_clean_scalar_product2(tmp2, R[i] + i + 1, R[i] + i + 1, c - i - 1, prec);
         mpfr_neg(V[i], tmp2, GMP_RNDN);
         mpfr_add(tmp2, R[i] + i, tmp, GMP_RNDN);
         mpfr_div(V[i], V[i], tmp2, GMP_RNDN);
      //now first entry of V[i] is r[1] - sigma

         sgn = mpfr_sgn(V[i]);
         if (sgn != 0){
            mpfr_mul(tmp2, tmp, V[i], GMP_RNDN);
            mpfr_neg(tmp2, tmp2, GMP_RNDN);
            mpfr_sqrt(tmp2, tmp2, GMP_RNDN);
            for (k = 0; k < c - i; k++){
               mpfr_div(V[i] + k, V[i] + k, tmp2, GMP_RNDN);
            }
         }
//now just write out r = size, 0, 0, 0...
         mpfr_abs(R[i] + i, tmp, GMP_RNDN);
         for (k = 1; k < c - i; k++){
            mpfr_set_ui(R[i] + i + k, 0, GMP_RNDN);
         }
      }
      else{
         mpfr_abs(R[i] + i, tmp, GMP_RNDN);
      }
   }

   for (i = 0; i < r; i++)
      for (j = 0; j < r; j++)
   	   mpfr_set(squareR[i] + j, R[i] + j, GMP_RNDN); 

   mpfr_mat_clear(V, r, c);
   mpfr_clear(tmp);
   mpfr_clear(tmp2);
   return;
}


int mpfr_mat_R_best_reduction(__mpfr_struct ** R, long d, mp_prec_t prec)
{

   if (d == 1){
      printf("trivially reduced\n");
      return 1;
   }

   mpfr_t tmp1;
   mpfr_t tmp2;
   mpfr_init2(tmp1, prec);
   mpfr_init2(tmp2, prec);

   mpfr_t delta;
   mpfr_t eta;
   mpfr_t theta;
   mpfr_init2(delta, prec);
   mpfr_init2(eta, prec);
   mpfr_init2(theta, prec);

   int reduced = 1;

   mpfr_set_ui(eta, 0, GMP_RNDN);
   mpfr_set_ui(theta, 0, GMP_RNDN);

   long i,j;
   for (i = 0; (i < d -1) && (reduced == 1); i++)
   {
      mpfr_pow_ui(tmp1, R[i+1] + i, 2L, GMP_RNDN);
      mpfr_pow_ui(tmp2, R[i+1] + i + 1, 2L, GMP_RNDN);
      mpfr_add(tmp1, tmp1, tmp2, GMP_RNDN);
//here tmp1 is the r_i,i+1^2 + r_i+1,i+1^2 which should compare to r_i,i^2

      mpfr_pow_ui(tmp2, R[i] + i, 2L, GMP_RNDN);
//here tmp2 is r_i,i^2

      mpfr_div(tmp1, tmp1, tmp2, GMP_RNDN);
//now tmp1 is the comparable ratio if tmp1 < delta replace delta
      if (mpfr_number_p(delta) == 0){
         mpfr_set(delta, tmp1, GMP_RNDN);
      }
      else{
         int cmp = mpfr_cmp(tmp1, delta);
         if (cmp < 0)
            mpfr_set(delta, tmp1, GMP_RNDN);
      }

      if (mpfr_cmp_d(delta, .65) < 0){
         printf("delta dropped below .65 at index i = %ld\n", i);
         //reduced = 0;
      }

      for (j = 0; (j < i + 1) && (reduced == 1); j++)
      {
//         printf("check R[%ld]+%ld\n", i + 1, j);
         mpfr_sub(tmp1 , R[i + 1] + j, R[j] + j, GMP_RNDN);
         mpfr_div(tmp2 , tmp1, R[j] + j, GMP_RNDN);
         //tmp2 is now a theta which would work for eta == 1
         int cmp2 = mpfr_cmp(tmp2, theta);
         if (cmp2 > 0)
             mpfr_set(theta, tmp2, GMP_RNDN);

         mpfr_div(tmp1 , R[i + 1] + j, R[j] + j, GMP_RNDN);
         mpfr_abs(tmp1 , tmp1, GMP_RNDN);
         cmp2 = mpfr_cmp(tmp1, eta);
         if (cmp2 > 0)
             mpfr_set(eta, tmp1, GMP_RNDN);

      }
   }

   mpfr_printf("Worst-case Delta is %.12Rf\n", delta); 
   mpfr_printf("Worst-case theta is %.12Rf\n", theta); 
   mpfr_printf("Worst-case eta is %.12Rf\n", eta); 

   mpfr_clear(tmp1);
   mpfr_clear(tmp2);
   mpfr_clear(theta);
   mpfr_clear(eta);
   mpfr_clear(delta);
   return reduced;
}



int main(int argc, char * argv[]){

   F_mpz_mat_t M;
   F_mpz_mat_init(M, 0, 0);

   F_mpz_mat_fread_pretty(M, stdin);

//   double x = 25.250293560534;
//   long exp = 42;

//   printf("%.12f is x, %.12f is ldexp\n", x, ldexp(x, exp));

   long r, c, i, j;
   mp_prec_t prec;

   if (argc == 1){
      prec = 50;
   }
   else{
      prec = atoi(argv[1]);
   }

   r = M->r;
   c = M->c;
   
   __mpfr_struct ** Q, ** R;

   Q = mpfr_mat_init2(r, c, prec);
   R = mpfr_mat_init2(r, r, prec);

//   LLL(M);

//   F_mpz_mat_print_pretty(M);


   F_mpz_mat_R_factor_house(M, R, Q, r, c, prec); 

/*   for (j = 0; j < r; j++)
   { 
      mpfr_printf("%.12Rf was R[i][i] for i = %ld\n", R[j] + j, j); 
   }

   for (i = 0; i < r; i++)
   for (j = 0; j < r; j++)
   { 
      mpfr_printf("%.12Rf was R[i][j] for i = %ld, j = %ld\n", R[i] + j, i, j); 
   }

   printf("was the Householder R factor\n");
*/

   int red = mpfr_mat_R_best_reduction(R, r, prec);

   mpfr_mat_clear(Q, r, c);
   mpfr_mat_clear(R, r, r);

   F_mpz_mat_clear(M);

   flint_stack_cleanup();

   return 0;
}

/*      _mpfr_vec_norm2(R[k]+k, Q[k], c, prec);
      mpfr_sqrt(R[k]+k, R[k]+k, GMP_RNDN);
      for (i = 0; i < c; i++)
         mpfr_div(Q[k]+i, Q[k]+i, R[k]+k, GMP_RNDN);
*/
