#include "fp2.h"
#include <stdlib.h>

static fp_t fp_char;

// init: allocate memory and set value = 0
void fp_init(fp_t res) {
    mpz_init(res);
}

// clear: deallocate memory
void fp_clear(fp_t res) {
    mpz_clear(res);
}

// set the field characteristic to p
void fp_set_global_char(fp_t p) {
    mpz_set(fp_char, p);
}

// set: result <- a
void fp_set(fp_t res, const fp_t a) {
    mpz_set(res, a);
}

// set uint: result <- (uint) a
void fp_set_uint(fp_t res, unsigned long int a) {
    mpz_set_ui(res, a);
}

// add: result = a + b (mod p)
void fp_add(fp_t res, const fp_t a, const fp_t b) {
    mpz_add(res, a, b);
    mpz_mod(res, res, fp_char);
} 

// add uint: result = a + (unsigned int) b (mod p)
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

// =========== 
// FP2 methods
// ===========

// init: allocate memory and set value to 0 + 0i
void fp2_init(fp2_t *res) {
    // allocate memory for the fp2 structure
    *res = (fp2_t) malloc(sizeof(fp2));
    // initialize both variables
    fp_init((*res)->a);
    fp_init((*res)->b);
}

// clear: free memory and set pointer to NULL 
void fp2_clear(fp2_t *res) {
    // clear both variables
    fp_clear((*res)->a);
    fp_clear((*res)->b);
    // free memory for res pointer
    free(*res);
    *res = NULL;
}

// fill: set individual fields as FP elements
void fp2_fill(fp2_t res, fp_t a, fp_t b) {
    fp_set(res->a, a);
    fp_set(res->b, b);
}

// set: result <- (unsigned int) rhs 
void fp2_set_uint(fp2_t res, unsigned long int rhs) {
    fp_set_uint(res->a, rhs);
    fp_set_uint(res->b, 0);
}

// set: result = a
void fp2_set(fp2_t res, const fp2_t arg) {
    fp_set(res->a, arg->a);
    fp_set(res->b, arg->b);
}

// add: result <- lhs[a + bi] + rhs[a + bi]
void fp2_add(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    fp_add(res->a, lhs->a, rhs->a);
    fp_add(res->b, lhs->b, rhs->b);
}

// add: result <- lhs[a + bi] + (unsigned int) rhs
void fp2_add_uint(fp2_t res, const fp2_t lhs, unsigned long int rhs) {
    fp_add_uint(res->a, lhs->a, rhs);
    fp_set(res->b, lhs->b);
}

// sub: result <- lhs[a + bi] - rhs[a + bi]
void fp2_sub(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    fp_sub(res->a, lhs->a, rhs->a);
    fp_sub(res->b, lhs->b, rhs->b);
}

// sub: result <- lhs[a + bi] - (unsigned int) rhs 
void fp2_sub_uint(fp2_t res, const fp2_t lhs, unsigned long int rhs) {
    fp_sub_uint(res->a, lhs->a, rhs);
    fp_set(res->b, lhs->b);
}

// mul: result <- lhs[a + bi] * rhs[a + bi]
void fp2_mul(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    // add temprary field elements for storing values
    fp_t tmp;
    fp_init(tmp); // allocate and set to 0

    // (a + bi) * (c + di) = [ac - bd] + [ad + bc]i

    // calculate the first variable (real part); ac - bd 
    // res[0] <- a * c
    fp_mul(res->a, lhs->a, rhs->a);
    // temp <- b * d
    fp_mul(tmp, lhs->b, rhs->b);
    // res[0] <- [res0](a * c) - [tmp](b * d)
    fp_sub(res->a, res->a, tmp);

    // calculate the second variable (real part); ad - bc
    // res[1] <- a * d
    fp_mul(res->b, lhs->a, rhs->b);
    // TODO: make sure we can set the register without UB
    // i.e-> there is no garbage left from other operations
    // so no need to set value to 0?

    // tmp <- b * c
    fp_mul(tmp, lhs->b, rhs->a);
    fp_add(res->b, res->b, tmp);

    fp_clear(tmp);
}

// sq: result <- arg^2
void fp2_sq(fp2_t res, const fp2_t arg) {
    fp2_mul(res, arg, arg);
}