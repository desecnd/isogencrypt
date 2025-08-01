#pragma once

#include <assert.h>
#include <gmp.h>

#include "fp.h"

// TODO: Maybe use a typedef to 1-dim array?
// https://stackoverflow.com/a/69731101
typedef struct {
    fp_t a;
    fp_t b;
} fp2_elem, *fp2_t;

void fp2_init(fp2_t *x);

void fp2_clear(fp2_t *x);

// set: result = x[a + bi]
void fp2_set(fp2_t res, const fp2_t x);

// set: result = (ul) a + 0i
void fp2_set_uint(fp2_t res, unsigned long int a);

/*
 * @brief Set fp2_t into value represented by string x.
 * Possible formats are: a + b*i, b*i + a, a, b*i.
 * Whitespaces cannot be inserted only between "*" and "*" - in all other cases
 * they are permited.
 */
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
 * @brief Compare fp2_t element with fp2_t element represented as str type.
 * Returns 1 if equal, 0 otherwise.
 * Allocates additional fp2_t element for comparison.
 */
int fp2_equal_str(fp2_t res, const char *x_str);

int fp2_equal_uint(fp2_t x, unsigned long int a);

// add: result <- x[a + bi] + y[a + bi]
void fp2_add(fp2_t res, const fp2_t x, const fp2_t y);

// add: result <- x[a + bi] + (unsigned int) y
void fp2_add_uint(fp2_t res, const fp2_t x, unsigned long int y);

// sub: result <- x[a + bi] - y[a + bi] (mod p)
void fp2_sub(fp2_t res, const fp2_t x, const fp2_t y);

// sub: result <- x[a + bi] - (unsigned int) y
void fp2_sub_uint(fp2_t res, const fp2_t x, unsigned long int y);

/*
 * @brief Multiply fp2_t by a fp-rational long int number
 * @param[out]  r   result of the operation x * y
 * @param[in]   x   fp2 element to multiply, can safely point to `r`
 * @param[in]   y   long int element to multiply
 */
void fp2_mul_int(fp2_t r, const fp2_t x, long int y);

/*
 * @brief Calculate mul: result = x[a + bi] * y[a + bi].
 * This function is not argument-safe, i.e. calling fp2_mul(r, x, y) for r = x
 * or r = y will provide incorrect results.
 */
void fp2_mul_unsafe(fp2_t res, const fp2_t x, const fp2_t y);

/*
 * @brief Calculate mul: x[a + bi] *= y[a + bi].
 * This function is argument-safe, i.e. calling fp2_mul(r, x, y) for r = x or r
 * = y is allowed.
 */
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
static inline void fp2_sq_safe(fp2_t x) { fp2_mul_safe(x, x); }

// Return true if fp2 is zero
static inline int fp2_is_zero(const fp2_t arg) {
    return (int)(fp_is_zero(arg->a) && fp_is_zero(arg->b));
}

// Calculate res: x / y
static inline void fp2_div_unsafe(fp2_t res, const fp2_t x, const fp2_t y) {
    assert(
        res != y && res != x &&
        "Fp^2 division is not arg-safe, cannot be called with r = y or r = x");
    fp2_inv_unsafe(res, y);
    // Uses safe in order to not allocate additional memory
    fp2_mul_safe(res, x);
}

/*
 * @brief Output fp2 element to the stdout in format: "name: a*i + b"
 */
void fp2_print(fp2_t x, const char *name);

/*
 * @brief Calculate size required to store fp2 in character buffer, used by
 * fp2_write
 */
size_t fp2_write_size(fp2_t x);

/*
 * @brief Write fp2 element x into character buffer adding NULL-termination
 * '\0'. X is stored in decimal format.
 */
void fp2_write(fp2_t x, char *buffer);
