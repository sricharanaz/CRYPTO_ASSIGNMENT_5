#ifndef __H_ELGAMAL_ECC_H
#define __H_ELGAMAL_ECC_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <stdlib.h>

/* format specifier for different cosole colors */
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KWHT  "\x1B[37m"

typedef struct point {
  mpz_t x, y;
} point;


typedef struct elliptic_curve {
  mpz_t a, b, p;
  point *base;
} elliptic_curve;


typedef struct elgam_ctx {
  mpz_t dom_par_q, dom_par_p, dom_par_g;
  mpz_t priv_x, pub_h, eph_k;
} elgam_ctx;


typedef struct ciphertext {
  mpz_t c1, c2;
} ciphertext;


typedef struct cipherec {
	point *c1, *c2;
} cipherec;


typedef struct elgam_ecc_ctx {
  mpz_t priv_key, eph_k;
  point *pub_key;
  elliptic_curve *ec;
} elgam_ecc_ctx;


void init_point(point **);
void destroy_point(point *);

point* ecc_scalar_mul(elliptic_curve *, mpz_t, point *);
point* ecc_addition(elliptic_curve *, point *, point *);
point* ecc_doubling(elliptic_curve *, point *);

point* ecc_scalar_mul2(elliptic_curve *, mpz_t, point *);
void ecc_addition2(elliptic_curve *, point *, point *, point *);
void ecc_doubling2(elliptic_curve *, point *, point *);
#endif
