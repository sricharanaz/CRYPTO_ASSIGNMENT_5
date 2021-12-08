#include "elg_ecc.h"

void init_point(point **p) 
{
	*p = (point*) malloc(sizeof(point));
	mpz_init((*p)->x);
	mpz_init((*p)->y);
}


void destroy_point(point *p)
{
	if (p) {
		mpz_clears(p->x, p->y, NULL);
		free(p);
		p = NULL;
	}
}

point* ecc_scalar_mul(elliptic_curve *ec, mpz_t m, point *p) {
	if (mpz_cmp_ui(m, 1) == 0) {
		return p;
	} else if (mpz_even_p(m) == 0) {
		mpz_sub_ui(m, m, 1);
		return ecc_addition(ec, p, ecc_scalar_mul(ec,m,p));
	} else {
		mpz_div_ui(m, m, 2);
		return ecc_scalar_mul(ec, m, ecc_doubling(ec, p));
	}
}


point* ecc_doubling(elliptic_curve *ec, point *p)
{
	point *r = malloc(sizeof(point));
	mpz_init(r->x);
	mpz_init(r->y);
	mpz_mod(p->x, p->x, ec->p);
	mpz_mod(p->y, p->y, ec->p);
	mpz_t temp, slope;
	mpz_init(temp);
	mpz_init_set_ui(slope, 0);

	// temp = 2*y1
	mpz_mul_ui(temp, p->y, 2);
	// temp = temp^-1 mod p
	mpz_invert(temp, temp, ec->p);
	// slope = x1*x1 = x1^2
	mpz_mul(slope, p->x, p->x);
	// slope = slope * 3
	mpz_mul_ui(slope, slope, 3);

	// slope = slope + a
	mpz_add(slope, slope, ec->a);

	// slope = slope * temp (numinator * denuminator)
	mpz_mul(slope, slope, temp);
	// slope = slope mod p
	mpz_mod(slope, slope, ec->p);

	// x3 = slope * slope
	mpz_mul(r->x, slope, slope);
	mpz_sub(r->x, r->x, p->x);
	mpz_sub(r->x, r->x, p->x);
	mpz_mod(r->x, r->x, ec->p);
	mpz_sub(temp, p->x, r->x);
	mpz_mul(r->y, slope, temp);
	mpz_sub(r->y, r->y, p->y);
	mpz_mod(r->y, r->y, ec->p);

	mpz_clears(temp, slope, NULL);
	return r;
}

point* ecc_addition(elliptic_curve *ec, point *p, point *q)
{
	point *r = malloc(sizeof(point));
	mpz_init(r->x);
	mpz_init(r->y);
	mpz_mod(p->x, p->x, ec->p);
	mpz_mod(p->y, p->y, ec->p);
	mpz_mod(q->x, q->x, ec->p);
	mpz_mod(q->y, q->y, ec->p);
	mpz_t temp,slope;
	mpz_init(temp);
	mpz_init_set_ui(slope, 0);

	// temp = x1-x2
	mpz_sub(temp, p->x, q->x);
	// temp = temp mod p
	mpz_mod(temp, temp, ec->p);
	// temp^-1 mod p
	mpz_invert(temp, temp, ec->p);
	// slope = y1-y2
	mpz_sub(slope, p->y, q->y);
	// slope = slope * temp
	mpz_mul(slope, slope, temp);
	// slope = slope mod p
	mpz_mod(slope, slope, ec->p);

	// x3 = slope * slope = alpha^2
	mpz_mul(r->x, slope, slope);

	// x3 = x3 - x1
	mpz_sub(r->x, r->x, p->x);
	// x3 = x3 - x2
	mpz_sub(r->x, r->x, q->x);
	// x3 = x3 mod p
	mpz_mod(r->x, r->x, ec->p);

	// temp = x1 - x3
	mpz_sub(temp, p->x, r->x);
	// y3 = slope * temp
	mpz_mul(r->y, slope, temp);
	// y3 = y3 - y1
	mpz_sub(r->y, r->y, p->y);
	// y3 = y3 mod p
	mpz_mod(r->y, r->y, ec->p);
	mpz_clears(temp, slope, NULL);
	return r;
}


/*
  Sets r to a random GMP integer with the specified number
  of bits.
*/
void get_random_n_bits(mpz_t r, size_t bits)
{
	size_t size = (size_t) ceilf(bits/8);
	char *buffer = (char*) malloc(sizeof(char)*size);
	int prg = open("/dev/random", O_RDONLY);
	read(prg, buffer, size);
	close(prg);
	mpz_import (r, size, 1, sizeof(char), 0, 0, buffer);
	free(buffer);
}

/*
  Sets r to a random GMP *prime* integer, smaller than max.
*/
void get_random_n_prime(mpz_t r, mpz_t max) 
{
	do {
		get_random_n_bits(r, mpz_sizeinbase(max, 2));
		mpz_nextprime(r, r);
	} while (mpz_cmp(r, max) >= 0);
}


/*
  Sets r to a random GMP integer smaller than max.
*/
void get_random_n(mpz_t r, mpz_t max) 
{
	do {
		get_random_n_bits(r, mpz_sizeinbase(max, 2));
	} while (mpz_cmp(r, max) >= 0);
}

void init_elgam_ecc(elgam_ecc_ctx **eec_ctx)
{
    *eec_ctx = (elgam_ecc_ctx*) malloc(sizeof(elgam_ecc_ctx));
    elliptic_curve *ecc = malloc(sizeof(elliptic_curve));
    (*eec_ctx)->ec = ecc;

	//Set elliptic curve.
    mpz_set_str(ecc->a, "340E7BE2A280EB74E2BE61BADA745D97E8F7C300", 16); 
    mpz_set_str(ecc->b, "1E589A8595423412134FAA2DBDEC95C8D8675E58", 16); 
    mpz_set_str(ecc->p, "E95E4A5F737059DC60DFC7AD95B3D8139515620F", 16); 

	mpz_init((*eec_ctx)->priv_key);
	init_point(&(ecc->base));
	init_point(&((*eec_ctx)->pub_key));

	// Set the base point.
	mpz_set_str(ecc->base->x, "BED5AF16EA3F6A4F62938C4631EB5AF7BDBCDBC3", 16); 
	mpz_set_str(ecc->base->y, "1667CB477A1A8EC338F94741669C976316DA6321", 16); 
	gmp_printf("\np = %Zd\n", ecc->p);

	// Choose a random private key.
	get_random_n((*eec_ctx)->priv_key, ecc->p);
	gmp_printf("x = %Zd\n", (*eec_ctx)->priv_key);

	mpz_t tmp;
	mpz_init_set(tmp, (*eec_ctx)->priv_key);
	(*eec_ctx)->pub_key = ecc_scalar_mul((*eec_ctx)->ec, tmp, ecc->base);
	mpz_clears(tmp, NULL);
	gmp_printf("Base point P = (%Zd,%Zd)\n", ecc->base->x, ecc->base->y);
	gmp_printf("Public key xP =  (%Zd,%Zd)\n\n", ((*eec_ctx)->pub_key)->x, ((*eec_ctx)->pub_key)->y);
}

void destroy_elgam_ecc(elgam_ecc_ctx *eec_ctx) 
{
	if (eec_ctx) {
 		mpz_clears(eec_ctx->priv_key, eec_ctx->eph_k, NULL);
		mpz_clears(eec_ctx->ec->a, eec_ctx->ec->b, eec_ctx->ec->p, NULL);
		destroy_point(eec_ctx->ec->base);
		destroy_point(eec_ctx->pub_key);
		if (eec_ctx->ec) {
			free(eec_ctx->ec);
			eec_ctx->ec = NULL;
		}
		free(eec_ctx);
		eec_ctx = NULL;
	}
}

cipherec* encrypt_ecc(elgam_ecc_ctx *eec, point *pm)
{
  printf
    ("\n----------------------------------------------------------------------------------------\n\n");
  printf ("%sElgamal Elliptic Curve Cryptography Encryption:%s\n", KYEL, KWHT);
	gmp_printf("plaintext: (%Zd,%Zd)\n", pm->x, pm->y);  

	mpz_init(eec->eph_k);
	get_random_n(eec->eph_k, eec->ec->p);
	gmp_printf("\nEphemeral key = %Zd\n", eec->eph_k);

	cipherec *cipher = malloc(sizeof(cipherec));
	init_point(&cipher->c1);
	init_point(&cipher->c2);
	mpz_t tmp;
	mpz_init_set(tmp, eec->eph_k);
	cipher->c1 = ecc_scalar_mul(eec->ec, tmp, eec->ec->base);
	mpz_clears(tmp, NULL);

	mpz_init_set(tmp, eec->eph_k);
	cipher->c2 = ecc_scalar_mul(eec->ec, tmp, eec->pub_key);
	mpz_clears(tmp, NULL);
	gmp_printf("Cipher c1: (%Zd,%Zd)\n", cipher->c1->x, cipher->c1->y);
	gmp_printf("Cipher c2 without msg: (%Zd,%Zd)\n", cipher->c2->x, cipher->c2->y);
	cipher->c2 = ecc_addition(eec->ec, cipher->c2, pm);
	gmp_printf("Cipher c2 with msg: (%Zd,%Zd)\n", cipher->c2->x, cipher->c2->y);
	mpz_clears(eec->eph_k, NULL);
	return cipher;
}

point* decrypt_ecc(elgam_ecc_ctx *eec, cipherec *c)
{
	printf
	("\n\n----------------------------------------------------------------------------------------\n\n");
	printf ("%sElgamal Elliptic Curve Cryptography Deccryption:%s\n", KYEL, KWHT);
  	point *d1, *d2;
  	init_point(&d1);
  	init_point(&d2);
	mpz_t tmp;
  	mpz_init_set(tmp, eec->priv_key);
  	d1 = ecc_scalar_mul(eec->ec, tmp, c->c1);

  	mpz_clears(tmp, NULL);
  	gmp_printf("\nD1=(%Zd,%Zd)\n", d1->x, d1->y);
	gmp_printf("Before neg: (%Zd,%Zd)\n", d1->x, d1->y);
	mpz_neg(d1->y, d1->y);
  	gmp_printf("After neg: (%Zd,%Zd)\n", d1->x, d1->y);
  	d2 = ecc_addition(eec->ec, c->c2, d1);
  	gmp_printf("Decrypted: (%Zd,%Zd)\n", d2->x, d2->y);
	destroy_point(d1);
	return d2;
}


void destroy_cipherecc(cipherec *c)
{
	if(!c)
		return;
	destroy_point(c->c1);
	destroy_point(c->c2);
	free(c);
	c = NULL;
}

int main() 
{
	// ElGamal-EC
	elgam_ecc_ctx *eec;
	init_elgam_ecc(&eec);
	printf("\nElGamal Elliptic Curve Cryptography\n");
	printf("Elliptic Curve General Form \t y2 = x3 + ax + b, y^2 mod p = (x^3  + A*x + B) mod p \n");
	point *p;
	cipherec *c;
	init_point(&p);
	gmp_printf("\nEnter Plain text in the form of point P(x,y):");
	gmp_scanf("%Zd,%Zd", p->x,p->y);
	c = encrypt_ecc(eec, p);
	destroy_point(p);
	
	init_point(&p);
	p = decrypt_ecc(eec, c);
	destroy_point(p);
	destroy_cipherecc(c);

	destroy_elgam_ecc(eec);
}
