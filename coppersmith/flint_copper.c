#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "F_mpz.h"
#include "F_mpz_poly.h"

void F_mpz_poly_power(F_mpz_poly_t output, const F_mpz_poly_t poly, const unsigned long exp){

   fmpz_poly_t foutput, fpoly;
   fmpz_poly_init(foutput);
   fmpz_poly_init(fpoly);

   F_mpz_poly_to_fmpz_poly(fpoly, poly);

   fmpz_poly_power(foutput, fpoly, exp);

   fmpz_poly_to_F_mpz_poly(output, foutput);

   fmpz_poly_clear(fpoly);
   fmpz_poly_clear(foutput);
}

int main(int argc, char * argv[]){

   long alpha, Xpow, Ypow, i, j;

   if (argc == 1){
      alpha = 3;
      Xpow = 20;
      Ypow = 220;
   }
   else{
      alpha = atoi(argv[1]);
      Xpow = atoi(argv[2]);
      Ypow = atoi(argv[3]);
   }
   
   F_mpz_poly_t f;
   F_mpz_poly_init(f);
   FILE * inst = fopen("pol_in","r");
   F_mpz_poly_fread(f, inst);
   fclose(inst);

   F_mpz_t Modulus;
   F_mpz_init(Modulus);

   F_mpz_read(Modulus);

//   F_mpz_print(Modulus); printf(" was modulus\n");
//   F_mpz_poly_print(f); printf(" was f\n");
   long d = f->length - 1;

   long ypow = d*(alpha) + 1;
   long maxlength = ypow*alpha + 1;
//   printf("ypow = %ld\n", ypow);

   F_mpz_poly_t temp;
   F_mpz_poly_init(temp);
   F_mpz_poly_t P[alpha + 1];
   F_mpz_mat_t M;
   F_mpz_mat_init(M, ((alpha+1)*(alpha+2))/2 , maxlength);
   long row = 0;

   F_mpz_poly_t m_arr; //an array of powers of the modulus, m_arr->coeffs + i should be Modulus^i
   F_mpz_poly_init(m_arr);
   F_mpz_poly_set_coeff_ui(m_arr, 0, 1);

   F_mpz_t temp_m;
   F_mpz_init(temp_m);

   F_mpz_set_ui(temp_m, 1);

   for (j = 1; j < alpha + 1; j++)
   {
      F_mpz_mul2(temp_m, temp_m, Modulus);
      F_mpz_poly_set_coeff(m_arr, j, temp_m);
   } 

   F_mpz_mat_t row_scale;
   F_mpz_mat_init(row_scale, M->r, M->r);

   for (j = 0; j < alpha + 1; j++)
   {
      if (j == 0){
         F_mpz_poly_init(P[0]);
         F_mpz_poly_set_coeff_ui(P[0], 0, 1);
      }
      else if (j == 1){
         F_mpz_poly_init(P[1]);
         F_mpz_poly_set(P[1], f);
//      printf("i is %ld and j is %ld\n", i, j);
         F_mpz_poly_set_coeff_si(P[1], ypow, -1);
      }
      else if (j > 1){
         F_mpz_poly_init(P[j]);
         F_mpz_poly_mul(P[j], P[j-1], P[1]);
      }

      for (i = 0; i < alpha + 1 - j; i++){
         F_mpz_poly_left_shift(temp, P[j], i);
//         F_mpz_poly_print(temp); printf(" was P with i = %ld and j = %ld\n", i, j);
         long k;
         for (k = 0; k < maxlength; k++){
            if (k < temp->length)
               F_mpz_set(M->rows[row]+k, temp->coeffs + k);
            else
               F_mpz_zero(M->rows[row]+k);
         }
         F_mpz_set(row_scale->rows[row] + row, m_arr->coeffs + alpha - j);          
         row++;
      }
   } 
   F_mpz_clear(temp_m);
   F_mpz_poly_clear(m_arr);

   F_mpz_t X;
   F_mpz_t Y;
   F_mpz_init(X);
   F_mpz_init(Y);

   F_mpz_set_ui(X, 2);
   F_mpz_pow_ui(X, X, Xpow);
   F_mpz_set_ui(Y, 2);
   F_mpz_pow_ui(Y, Y, Ypow);

   F_mpz_mat_t col_scale;
   F_mpz_mat_init(col_scale, maxlength, maxlength);

   F_mpz_set_ui(col_scale->rows[0], 1);
   for (i = 1; i < ypow; i++){
      F_mpz_mul2(col_scale->rows[i] + i, col_scale->rows[i-1] + i-1, X);
   }

   for (j = 1; j < alpha; j++){
      for (i = ypow*(j-1); i < j*ypow; i++){
         F_mpz_mul2(col_scale->rows[i + ypow] + i + ypow, col_scale->rows[i] + i, Y);
      }
   }
//using X bounds and Y bounds as powers of 2 but not taking advantage of 
//that via cheap multiplication (or to save in LLL costs) and we should


   F_mpz_mul2(col_scale->rows[maxlength-1] + maxlength -1 , col_scale->rows[maxlength - ypow - 1] + maxlength - ypow - 1, Y);

//   F_mpz_mat_print_pretty(col_scale);

   F_mpz_mat_mul_classical(M, M, col_scale);

   F_mpz_mat_mul_classical(M, row_scale, M);

   FILE * fpre = fopen("pre_LLL","w");
   F_mpz_mat_fprint_pretty(M, fpre);
   fclose(fpre);

   LLL(M);

   FILE * fpost = fopen("post_LLL","w");
   F_mpz_mat_fprint_pretty(M, fpost);
   fclose(fpost);

   F_mpz_poly_clear(temp);
   F_mpz_poly_init(temp);

   for (i = 0; i < maxlength; i++){
      F_mpz_poly_set_coeff(temp, i, M->rows[0]+i);
      F_mpz_divexact(temp->coeffs + i, temp->coeffs + i, col_scale->rows[i] + i); //undo col_scale
   }

   FILE * fpoly1 = fopen("poly1","w");
   F_mpz_poly_fprint(temp, fpoly1);
   fclose(fpoly1);
//   F_mpz_poly_print(temp); printf(" was the first row after LLL as a polynomial\n");

   F_mpz_poly_clear(temp);
   F_mpz_poly_init(temp);

   for (i = 0; i < maxlength; i++){
      F_mpz_poly_set_coeff(temp, i, M->rows[1]+i);
      F_mpz_divexact(temp->coeffs + i, temp->coeffs + i, col_scale->rows[i] + i);
   }
//   F_mpz_poly_print(temp); printf(" was the second row after LLL as a polynomial\n"); //undo col_scale
   FILE * fpoly2 = fopen("poly2","w");
   F_mpz_poly_fprint(temp, fpoly2);
   fclose(fpoly2);

   F_mpz_poly_clear(temp);
   F_mpz_poly_init(temp);

   for (i = 0; i < maxlength; i++){
      F_mpz_poly_set_coeff(temp, i, M->rows[2]+i);
      F_mpz_divexact(temp->coeffs + i, temp->coeffs + i, col_scale->rows[i] + i); //undo col_scale
   }

   FILE * fpoly3 = fopen("poly3","w");
   F_mpz_poly_fprint(temp, fpoly3);
   fclose(fpoly3);
//   F_mpz_poly_print(temp); printf(" was the first row after LLL as a polynomial\n");

   F_mpz_poly_clear(temp);
   F_mpz_poly_init(temp);

   for (i = 0; i < maxlength; i++){
      F_mpz_poly_set_coeff(temp, i, M->rows[3]+i);
      F_mpz_divexact(temp->coeffs + i, temp->coeffs + i, col_scale->rows[i] + i);
   }
//   F_mpz_poly_print(temp); printf(" was the second row after LLL as a polynomial\n"); //undo col_scale
   FILE * fpoly4 = fopen("poly4","w");
   F_mpz_poly_fprint(temp, fpoly4);
   fclose(fpoly4);

   F_mpz_clear(X);
   F_mpz_clear(Y);
   F_mpz_clear(Modulus);

   for (j = 0; j < alpha + 1; j++){
      F_mpz_poly_clear(P[j]);
   }
   F_mpz_poly_clear(temp);
   F_mpz_mat_clear(col_scale);
   F_mpz_mat_clear(row_scale);
   F_mpz_mat_clear(M);
   F_mpz_poly_clear(f);
   flint_stack_cleanup();
   return 0;
}
