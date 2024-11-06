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
void fp_set_char(fp_t p); 

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

// div: a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b); 

// neg: a = -a (mod p)
void fp_neg(fp_t res, const fp_t a); 

// sqrt: 
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a); 

// --- FP2 Methods



#endif