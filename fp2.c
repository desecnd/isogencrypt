#include "fp2.h"

static fp_t fp_char;

void fp_init(fp_t res) {
    mpz_init(res);
}

void fp_clear(fp_t res) {
    mpz_clear(res);
}

// set the field characteristic to p
void fp_set_char(fp_t p) {
    mpz_set(fp_char, p);
}

// set: res = a 
void fp_set(fp_t res, const fp_t a) {
    mpz_set(res, a);
}

void fp_set_uint(fp_t res, unsigned long int a) {
    mpz_set_ui(res, a);
}

// add: result = a + b (mod p)
void fp_add(fp_t res, const fp_t a, const fp_t b) {
    mpz_add(res, a, b);
    mpz_mod(res, res, fp_char);
} 

// add unsigned int: result = a + (unsigned int) b (mod p)
void fp_add_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_add_ui(res, a, b);
    mpz_mod(res, res, fp_char);
}

// sub: res = a - b (mod p)
void fp_sub(fp_t res, const fp_t a, const fp_t b) {
    mpz_sub(res, a, b);
    mpz_mod(res, res, fp_char);
}

// sub uint: res = a - (unsigned int) b (mod p)
void fp_sub_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_sub_ui(res, a, b);
    mpz_mod(res, res, fp_char);
}

// mul: res = a * b (mod p)
void fp_mul(fp_t res, const fp_t a, const fp_t b) {
    mpz_mul(res, a, b);
    mpz_mod(res, res, fp_char);
}

// mul int: res = a * (int) b (mod p)
void fp_mul_int(fp_t res, const fp_t a, long int b) {
    mpz_mul_si(res, a, b);
    mpz_mod(res, res, fp_char);
}

// modular inverse: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a) {
    mpz_invert(res, a, fp_char);
}

// div: a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b) {
    mpz_invert(res, b, fp_char);
    mpz_mul(res, res, a);
    mpz_mod(res, res, fp_char);
}

// neg: a = -a (mod p)
void fp_neg(fp_t res, const fp_t a) {
    mpz_sub(res, fp_char, a);
    mpz_mod(res, res, fp_char);
}

// sqrt: 
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a) {
    mpz_t exp; 
    mpz_init(exp);

    // calculate a^((p + 1)/4)

    // step 1: e = p + 1
    mpz_add_ui(exp, fp_char, 1);

    // step 2: e = (p + 1)/4
    // can use divexact only if we know the divisor: d = 4
    mpz_divexact_ui(exp, exp, 4);

    // step 3: res = a ^ e
    mpz_powm(res, a, exp, fp_char);

    // step 4: clear memory for e variable
    mpz_clear(exp);
}
