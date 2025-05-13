#ifndef FP2_H
#define FP2_H

#include <assert.h>
#include <gmp.h>

typedef mpz_t fp_t;


// TODO: Maybe use a typedef to 1-dim array?
// https://stackoverflow.com/a/69731101

typedef union {
     struct {
         fp_t a;
         fp_t b;
     };
     fp_t raw[2];
} fp2_elem, *fp2_t;

// ------ FP Methods

// set the global static field characteristic to p
// TODO: add check for valid modulus x^2 + 1 == 0
// prime should be select such that x^2 + 1 is irreducible over Fp^2
int global_fpchar_setup(fp_t p); 
int global_fpchar_setup_uint(unsigned int p); 
int global_fpchar_clear();

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

// modular inverse: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a); 

// div: res = a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b); 

// neg: res = -a (mod p)
void fp_neg(fp_t res, const fp_t a); 

// sqrt: res = a^(1/2) (mod p)
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a); 

// Return 1 if fp is zero, 0 otherwise
static int fp_is_zero(const fp_t a) {
    return (int) (mpz_sgn(a) == 0);
}

/*
 * @brief Return 1 if FiniteField element is equal to ulong number
 */
int fp_equal_uint(fp_t a, unsigned long int b);

/*
 * @brief Return 1 if two FiniteField elements are equal, 0 otherwise
 */
int fp_equal(fp_t a, fp_t b);

// --- FP2 Methods

void fp2_init(fp2_t *res);

void fp2_clear(fp2_t *res);


// set: result <- a + bi
void fp2_set(fp2_t res, const fp2_t arg);

// set: result <- (unsigned int) rhs 
void fp2_set_uint(fp2_t res, unsigned long int rhs);

// set: result <- x[a + bi]
int fp2_set_str(fp2_t res, const char *x);

// fill: set individual fields as FP elements
void fp2_fill(fp2_t res, fp_t a, fp_t b);

// fill: res = (a, b)
void fp2_fill_uint(fp2_t res, unsigned long int a, unsigned long int b);

// fill: res = (a, b)
// if string are given with prefix (0x, 0, 0b) 
// corresponding base will be detected and applied to conversion
void fp2_fill_str(fp2_t res, const char *a, const char *b);

int fp2_equal(fp2_t x, fp2_t y);

/*
 * @brief Compare fp2_t element with fp2 element represented as str type. 
 * Returns 1 if equal, 0 otherwise. 
 * Allocates additional fp2_t element for comparison.
 */
int fp2_equal_str(fp2_t res, const char * x_str);

int fp2_equal_uint(fp2_t x, unsigned long int a);


// add: result <- lhs[a + bi] + rhs[a + bi]
void fp2_add(fp2_t res, const fp2_t lhs, const fp2_t rhs);

// add: result <- lhs[a + bi] + (unsigned int) rhs
void fp2_add_uint(fp2_t res, const fp2_t lhs, unsigned long int rhs);

// sub: result <- lhs[a + bi] - rhs[a + bi] (mod p)
void fp2_sub(fp2_t res, const fp2_t lhs, const fp2_t rhs);

// sub: result <- lhs[a + bi] - (unsigned int) rhs 
void fp2_sub_uint(fp2_t res, const fp2_t lhs, unsigned long int rhs);

/*
 * @brief Multiply fp2_t by a fp-rational long int number
 * @param[out]  r   result of the operation x * y
 * @param[in]   x   fp2 element to multiply, can safely point to `r`
 * @param[in]   y   long int element to multiply
 */
void fp2_mul_int(fp2_t r, const fp2_t x, long int y);

// mul: result <- lhs[a + bi] * rhs[a + bi] (mod p)
void fp2_mul_unsafe(fp2_t res, const fp2_t lhs, const fp2_t rhs);
void fp2_mul_safe(fp2_t x, const fp2_t y);

void fp2_inv_unsafe(fp2_t res, const fp2_t arg);
void fp2_inv_safe(fp2_t x);

// sq: result <- arg^2
// This function cannot be used when res == arg
static inline void fp2_sq_unsafe(fp2_t res, const fp2_t arg) {
    assert(res != arg && "fp2_sq cannot be called with res = arg");
    fp2_mul_unsafe(res, arg, arg);
}

// sq: result <- arg^2
// This function can be used safely with res == arg
static inline void fp2_sq_safe(fp2_t x) {
    fp2_mul_safe(x, x);
}

// Return true if fp2 is zero
static inline int fp2_is_zero(const fp2_t arg) {
    return (int) (fp_is_zero(arg->a) && fp_is_zero(arg->b));
}

// Calculate res: lhs / rhs
static inline void fp2_div_unsafe(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    assert(res != rhs && res != lhs && "Fp^2 division is not arg-safe, cannot be called with r = rhs or r = lhs");
    fp2_inv_unsafe(res, rhs);
    // Uses safe in order to not allocate additional memory
    fp2_mul_safe(res, lhs);
}

/* 
 * @brief Output fp2 element to the stdout in format: "name: a*i + b"
 */
void fp2_print(fp2_t x, const char* name);

#endif