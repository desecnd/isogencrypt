#include <assert.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp2.h"

// init: allocate memory and set value to 0 + 0i
void fp2_init(fp2_t *res) {
    // allocate memory for the fp2 structure
    *res = (fp2_t)malloc(sizeof(fp2_elem));
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

// set: result = x
void fp2_set(fp2_t r, const fp2_t x) {
    fp_set(r->a, x->a);
    fp_set(r->b, x->b);
}

// set: result = (uint) x
void fp2_set_uint(fp2_t res, unsigned long int x) {
    fp_set_uint(res->a, x);
    fp_set_uint(res->b, 0);
}

int fp2_set_str(fp2_t res, const char *x) {
    int two_parts = (int)(strstr(x, "+") != NULL);
    int imag_part = (int)(strstr(x, "*i") != NULL);

    int only_real = !two_parts && !imag_part;
    int only_imag = !two_parts && imag_part;

    // Only a part was passed as argument
    if (only_real) {
        mpz_set_ui(res->b, 0);
        // "This function returns 0 if the entire string is a valid number in
        // base base. Otherwise it returns −1."
        return mpz_set_str(res->a, x, 0);
    }

    // Allocate string on the heap to allow for strsep and truncating
    unsigned long n_chars = strlen(x);
    char *buffer = malloc(sizeof(char) * (n_chars + 1));
    memcpy(buffer, x, n_chars);
    buffer[n_chars] = '\0';

    // Only b*i part was passed as argument
    if (only_imag) {
        char *temp = buffer;
        // "Cut" the string before *i, imag will point to the start of the
        // buffer, "*" will be replaced with "\0"
        char *imag = strsep(&temp, "*i");
        mpz_set_ui(res->a, 0);
        int ret_imag = mpz_set_str(res->b, imag, 0);
        free(buffer);
        return ret_imag;
    }

    char *imag = buffer;
    char *real = strsep(&imag, "+");
    char *temp = imag;

    // swap imag and real
    if (strstr(imag, "*i") == NULL) {
        imag = real;
        real = temp;
        temp = imag;
    } else if (strstr(real, "*i") != NULL) {
        // Found i in both real and imag
        free(buffer);
        return 1;
    }

    // Strip the "*i" part
    imag = strsep(&temp, "*i");

    // TODO: change inconsistent API
    int ret_real = mpz_set_str(res->a, real, 0);
    int ret_imag = mpz_set_str(res->b, imag, 0);

    free(buffer);

    // Return -1 if any errors occured
    return (ret_real == 0 && ret_imag == 0) ? 0 : -1;
}

// fill: set individual fields as FP elements
void fp2_fill(fp2_t res, fp_t a, fp_t b) {
    fp_set(res->a, a);
    fp_set(res->b, b);
}

// fill_str: result = a + bi
void fp2_fill_str(fp2_t res, const char *a, const char *b) {
    // TODO: change inconsistent API
    mpz_set_str(res->a, a, 0);
    mpz_set_str(res->b, b, 0);
}

// fill_uint: result = a + bi
void fp2_fill_uint(fp2_t res, unsigned long int a, unsigned long int b) {
    fp_set_uint(res->a, a);
    fp_set_uint(res->b, b);
}

// add: result <- x[a + bi] + y[a + bi]
void fp2_add(fp2_t res, const fp2_t x, const fp2_t y) {
    fp_add(res->a, x->a, y->a);
    fp_add(res->b, x->b, y->b);
}

// add: result <- x[a + bi] + (unsigned int) y
void fp2_add_uint(fp2_t res, const fp2_t x, unsigned long int y) {
    fp_add_uint(res->a, x->a, y);
    fp_set(res->b, x->b);
}

// sub: result <- x[a + bi] - y[a + bi]
void fp2_sub(fp2_t res, const fp2_t x, const fp2_t y) {
    fp_sub(res->a, x->a, y->a);
    fp_sub(res->b, x->b, y->b);
}

// sub: result <- x[a + bi] - (unsigned int) y
void fp2_sub_uint(fp2_t res, const fp2_t x, unsigned long int y) {
    fp_sub_uint(res->a, x->a, y);
    fp_set(res->b, x->b);
}

void fp2_mul_unsafe(fp2_t res, const fp2_t x, const fp2_t y) {
    assert(res != x && res != y &&
           "fp2_mul cannot be called with res = x or res = y");
    // add temprary field elements for storing values
    fp_t t;
    fp_init(t); // allocate and set to 0

    // (a + bi) * (c + di) = [ac - bd] + [ad + bc]i

    // calculate the first variable (real part); ac - bd
    // res[0] <- a * c
    fp_mul(res->a, x->a, y->a);
    // temp <- b * d
    fp_mul(t, x->b, y->b);
    // res[0] <- [res0](a * c) - [tmp](b * d)
    fp_sub(res->a, res->a, t);

    // calculate the second variable (real part); ad - bc
    // res[1] <- a * d
    fp_mul(res->b, x->a, y->b);
    // TODO: make sure we can set the register without UB
    // i.e-> there is no garbage left from other operations
    // so no need to set value to 0?

    // tmp <- b * c
    fp_mul(t, x->b, y->a);
    fp_add(res->b, res->b, t);

    fp_clear(t);
}

void fp2_mul_int(fp2_t r, const fp2_t x, long int y) {
    // Multiply each of the "coordinates"
    // xy = (a + bi)y = ay + byi
    fp_mul_int(r->a, x->a, y);
    fp_mul_int(r->b, x->b, y);
}

void fp2_mul_safe(fp2_t x, const fp2_t y) {
    // x = a + bi, y = c + di
    // r = (ac - bd) + (ad + bc)i = e + fi
    // cost: 4m + 2a

    fp_t t0, t1, tc;
    fp_init(t0);
    fp_init(t1);
    fp_init(tc);

    // tc is used for storing value c in case
    // it will be overwritten in step: e = ac
    // if fp2_mul is called with r = x = y
    fp_set(tc, y->a);       // tc = c
    fp_mul(t0, x->b, y->b); // t0 = bd
    fp_mul(t1, x->a, y->b); // t1 = ad

    fp_mul(x->a, x->a, y->a); // e = ac
    fp_sub(x->a, x->a, t0);   // e = ac - bd

    fp_mul(x->b, x->b, tc); // f = bc
    fp_add(x->b, x->b, t1); // f = bc + ad

    fp_clear(t0);
    fp_clear(t1);
    fp_clear(tc);
}

// Calculate inverse of element in Fp^2
// Cannot be used with r = arg
void fp2_inv_unsafe(fp2_t r, const fp2_t x) {
    assert(
        r != x &&
        "Fp^2 inversion is not argument safe, cannot be called with r = arg");
    assert(!fp2_is_zero(x) && "Fp^2 inversion does not accept x = 0");

    fp_mul(r->a, x->a, x->a); // r.a = a^2
    fp_mul(r->b, x->b, x->b); // r.b = b^2

    fp_add(r->a, r->a, r->b); // r.a = a^2 + b^2
    fp_inv(r->a, r->a);       // r.a = (a^2 + b^2)^-1

    fp_mul(r->b, x->b, r->a); // r.b = b * (a^2 + b^2)^-1
    fp_neg(r->b, r->b);       // r.b = -b * (a^2 + b^2)^-1

    fp_mul(r->a, r->a, x->a); // r.a = a * (a^2 + b^2)^-1

    // r = a(a^2 + b^2)^-1 + b(a^2 + b^2)^-1 * i
}

void fp2_inv_safe(fp2_t x) {
    // Calculates the inverse of fp^2 element in a argument-safe manner (i.e.
    // valid when r = x) input: x = a + bi output: r = [a * (a^2 + b^2)^-1] + [b
    // * (a^2 + b^2)^-1]i cost:

    assert(!fp2_is_zero(x) && "Fp^2 inversion does not accept x = 0");

    fp_t t0, t1;
    fp_init(t0);
    fp_init(t1);

    fp_mul(t0, x->a, x->a); // t0 = a^2
    fp_mul(t1, x->b, x->b); // t1 = b^2

    fp_add(t0, t0, t1); // t0 = a^2 + b^2
    fp_inv(t0, t0);     // t0 = (a^2 + b^2)^-1

    fp_mul(x->b, x->b, t0); // x.b = b * (a^2 + b^2)^-1
    fp_neg(x->b, x->b);     // x.b = -b * (a^2 + b^2)^-1

    fp_mul(x->a, x->a, t0); // x.a = a * (a^2 + b^2)^-1

    // Clear the allocations
    fp_clear(t0);
    fp_clear(t1);
}

/*
 * @brief Output fp2 element to the stdout in format: "name: a*i + b"
 */
void fp2_print(fp2_t x, const char *name) {
    if (fp_is_zero(x->b)) {
        gmp_printf("%s: %Zd\n", name, x->a);
    } else if (fp_is_zero(x->a)) {
        gmp_printf("%s: %Zd*i\n", name, x->b);
    } else {
        gmp_printf("%s: %Zd*i + %Zd\n", name, x->b, x->a);
    }
}

int fp2_equal_uint(fp2_t x, unsigned long int a) {
    return fp_equal_uint(x->a, a) && fp_is_zero(x->b);
}

int fp2_equal(fp2_t x, fp2_t y) {
    return !!(fp_equal(x->a, y->a) && fp_equal(x->b, y->b));
}

int fp2_equal_str(fp2_t res, const char *x_str) {
    fp2_t x;
    fp2_init(&x);

    fp2_set_str(x, x_str);

    int are_equal = fp2_equal(res, x);
    fp2_clear(&x);

    return are_equal;
}

size_t fp2_write_size(fp2_t x) {
    // Calculate buffer size required before calling "fp2_write"

    // TODO: Maybe store as hex in the future to save space?
    size_t size_a = mpz_sizeinbase(x->a, 10);
    size_t size_b = mpz_sizeinbase(x->b, 10);
    size_t size_mid = 5; // "*i + "

    // Add +1 for null termination
    size_t size_sum = size_a + size_mid + size_b + 1;
    return size_sum;
}

void fp2_write(fp2_t x, char *buffer) {
    // Use GMP sprintf to store the fp2 value in format respected by fp2_set_str
    gmp_sprintf(buffer, "%Zd*i + %Zd", x->b, x->a);
}
