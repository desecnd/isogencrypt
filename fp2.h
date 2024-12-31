#ifndef FP2_H
#define FP2_H

#include <gmp.h>

typedef mpz_t fp_t;


typedef union {
    struct {
        fp_t a;
        fp_t b;
    };
    fp_t elem[2];
} fp2_t;


// ------ FP Methods

// init memory for variable
void fp_init(fp_t res); 

// clear memory for variable 
void fp_clear(fp_t res); 

// set the field characteristic to p
void fp_set_global_char(fp_t p); 

// set: res = a 
void fp_set(fp_t res, const fp_t a); 

// set uint: res = (unsigned int) a
void fp_set_uint(fp_t res, unsigned long int a);

// add: result = a + b (mod p)
void fp_add(fp_t res, const fp_t a, const fp_t b);

// add unsigned int: result = a + (unsigned int) b (mod p)
void fp_add_uint(fp_t res, const fp_t a, unsigned long int b); 

// sub: res = a - b (mod p)
void fp_sub(fp_t res, const fp_t a, const fp_t b); 

// sub uint: res = a - (unsigned int) b (mod p)
void fp_sub_uint(fp_t res, const fp_t a, unsigned long int b); 

// mul: res = a * b (mod p)
void fp_mul(fp_t res, const fp_t a, const fp_t b); 

// mul int: res = a * (int) b (mod p)
void fp_mul_int(fp_t res, const fp_t a, long int b); 

// modular inverse: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a); 

// div: res = a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b); 

// neg: res = -a (mod p)
void fp_neg(fp_t res, const fp_t a); 

// sqrt: res = a^(1/2) (mod p)
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a); 

// --- FP2 Methods

void fp2_init(fp2_t res);

void fp2_clear(fp2_t res);

// fill: set individual fields as FP elements
void fp2_fill(fp2_t res, fp_t a, fp_t b);

// set: result <- a
void fp2_set(fp2_t res, const fp2_t arg);

// add: result <- lhs[a + bi] + rhs[a + bi]
void fp2_add(fp2_t res, const fp2_t lhs, const fp2_t rhs);

// add: result <- lhs[a + bi] + (unsigned int) rhs
void fp2_add_uint(fp2_t res, const fp2_t lhs, unsigned long int rhs);

// sub: result <- lhs[a + bi] - rhs[a + bi] (mod p)
void fp2_sub(fp2_t res, const fp2_t lhs, const fp2_t rhs);

#endif