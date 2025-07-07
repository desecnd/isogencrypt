#pragma once

#include <assert.h>
#include <gmp.h>

typedef mpz_t fp_t;

int fpchar_setup(fp_t p); 
int fpchar_setup_uint(unsigned int p); 
int fpchar_clear();
int fpchar_clear_if_set();
int fpchar_check();

// init memory for variable
void fp_init(fp_t res); 

// clear memory for variable 
void fp_clear(fp_t res); 

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

// inv: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a); 

// div: res = a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b); 

// neg: res = -a (mod p)
void fp_neg(fp_t res, const fp_t a); 

// sqrt: res = a^(1/2) (mod p)
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a); 

// Return 1 if fp is zero, 0 otherwise
int fp_is_zero(const fp_t a);

/*
 * @brief Return 1 if FiniteField element is equal to ulong number
 */
int fp_equal_uint(fp_t a, unsigned long int b);

/*
 * @brief Return 1 if two FiniteField elements are equal, 0 otherwise
 */
int fp_equal(fp_t a, fp_t b);


/*
 * @brief Compare fp_t element with fp_t element represented as a str type. 
 * Return 1 if equal, 0 otherwise.
 */
int fp_equal_str(fp_t a, const char *b_str);

/* 
 * @brief Output fp element to the stdout in format: "name: a"
 */
void fp_print(fp_t a, const char* name);
