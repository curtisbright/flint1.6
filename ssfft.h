#ifndef SSFFT_H
#define SSFFT_H

#include <gmp.h>

void ssfft_fft(mp_limb_t** x, unsigned long t,
               unsigned long m, unsigned long z, unsigned long g,
               unsigned long ru, unsigned long rU,
               unsigned long n, mp_limb_t** scratch);

void ssfft_fft_threaded(
               mp_limb_t** x, unsigned long t,
               unsigned long m, unsigned long z, unsigned long g,
               unsigned long ru, unsigned long rU,
               unsigned long n);

void ssfft_ifft(mp_limb_t** x, unsigned long t,
               unsigned long m, unsigned long z, unsigned long g, int e,
               unsigned long ru, unsigned long rU, unsigned long n,
               mp_limb_t** scratch);

// exported for testing only:

void ssfft_signed_add_1(mp_limb_t* x, unsigned long count,
                        mp_limb_signed_t limb);

void basic_fast_reduce(mp_limb_t* x, unsigned long n);

void basic_unrotate_bits(mp_limb_t* b, mp_limb_t* a,
                         unsigned long s, unsigned long n);

void basic_rotate_limbs(mp_limb_t* b, mp_limb_t* a,
                        unsigned long s, unsigned long n);

void basic_unrotate_limbs(mp_limb_t* b, mp_limb_t* a,
                          unsigned long s, unsigned long n);

void basic_sub_rotate_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                            unsigned long s, unsigned long n);

void basic_rotate_sub_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                            unsigned long s, unsigned long n);

void basic_unrotate_sub_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                              unsigned long s, unsigned long n);

void basic_unrotate_add_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                              unsigned long s, unsigned long n);

void basic_copy(mp_limb_t* b, mp_limb_t* a, unsigned long n);

void basic_negate(mp_limb_t* b, mp_limb_t* a, unsigned long n);

void coeff_forward_simple_butterfly(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch, unsigned long n);

void coeff_forward_butterfly_limbs(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n);

void coeff_inverse_butterfly_limbs(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n);

void coeff_sub_rotate_bits(mp_limb_t** c, mp_limb_t** a, mp_limb_t** b,
                           unsigned long s, unsigned long n);

void coeff_forward_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n);

void coeff_inverse_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n);

void coeff_rotate_bits(mp_limb_t** b, mp_limb_t** a,
                       unsigned long s, unsigned long n);

void coeff_cross_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n);

void coeff_cross_butterfly2(mp_limb_t** a, mp_limb_t** b, unsigned long n);

void coeff_sqrt2_helper(mp_limb_t** b, mp_limb_t** a, unsigned long n);

void coeff_rotate_arbitrary(mp_limb_t** b, mp_limb_t** a, mp_limb_t** scratch,
                            unsigned long s, unsigned long n);


void ssfft_fft_size4_bits(
       mp_limb_t** x, unsigned long t, unsigned long z, unsigned long g,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch);

void ssfft_fft_size4_z4g4_limbs(
       mp_limb_t** x, unsigned long t,
       unsigned long ru_limbs, unsigned long rU_limbs, unsigned long n,
       mp_limb_t** scratch);

void ssfft_fft_size8_bits(
       mp_limb_t** x, unsigned long t, unsigned long z, unsigned long g,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch);

void ssfft_fft_size8_z8g8_limbs(
       mp_limb_t** x, unsigned long t,
       unsigned long ru_limbs, unsigned long rU_limbs, unsigned long n,
       mp_limb_t** scratch);


void ssfft_ifft_size4_bits(
       mp_limb_t** x, unsigned long t,
       unsigned long z, unsigned long g, int e,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch);

void ssfft_ifft_size4_z4g4_limbs(
       mp_limb_t** x, unsigned long t,
       unsigned long ru_limbs, unsigned long rU_limbs, unsigned long n,
       mp_limb_t** scratch);

void ssfft_ifft_size8_bits(
       mp_limb_t** x, unsigned long t,
       unsigned long z, unsigned long g, int e,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch);

void ssfft_ifft_size8_z8g8_limbs(
       mp_limb_t** x, unsigned long t,
       unsigned long ru_limbs, unsigned long rU_limbs, unsigned long n,
       mp_limb_t** scratch);
#endif
