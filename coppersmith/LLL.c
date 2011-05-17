#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "flint.h"
#include "F_mpz_LLL.h"

int main(int argc, char * argv[]){

   F_mpz_mat_t M;
   F_mpz_mat_init(M, 0, 0);

   F_mpz_mat_fread_pretty(M, stdin);

   clock_t t1,t2;

   t1=clock();
   LLL(M);
   t2 = clock();

   printf("took %ld ticks = %f seconds\n", t2 - t1, ( (double) t2 - (double) t1 ) / ( (double) CLOCKS_PER_SEC ) );

   F_mpz_mat_clear(M);

   flint_stack_cleanup();

   return 0;
}

