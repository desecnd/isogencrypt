#include "fp2.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

static fp_t fpchar;
static int is_fpchar_set = 0;

int global_fpchar_setup_uint(unsigned int p) {
    fp_t pp;
    fp_init(pp);
    fp_set_uint(pp, p);
    int retval = global_fpchar_setup(pp);
    fp_clear(pp);
    return retval;
}

// Initialize and set the field characteristic of Fp and Fp^2 arithmetic.
// For now, only p = 3 (mod 4) is allowed
int global_fpchar_setup(fp_t p) {
    if (is_fpchar_set) {
        printf("[!] trying to initialize the characteristic second time!\n");
        return -1;
    }

    int invalid_mod4 = 0;
    
    // Check for modulo 4 operans
    {
        fp_t t; fp_init(t);
        // TODO: change inconsistent API
        mpz_mod_ui(t, p, 4);
        // invalid := p mod 4 is different than 3
        invalid_mod4 = mpz_cmp_ui(t, 3);
        fp_clear(t);
    }

    if (invalid_mod4) {
        printf("[!] Trying to set prime field characteristic with prime != 3 (mod 4)\n");
        return -2;
    }

    fp_init(fpchar);
    fp_set(fpchar, p);
    is_fpchar_set = 1;
    return 0;
}

int global_fpchar_clear() {
    if (!is_fpchar_set) {
        printf("[!] trying to clear not-initialized characteristic!\n");
        return -1;
    }
    fp_clear(fpchar);
    is_fpchar_set = 0;
    return 0;
}

// init: allocate memory and set value = 0
void fp_init(fp_t res) {
    mpz_init(res);
}

// clear: deallocate memory
void fp_clear(fp_t res) {
    mpz_clear(res);
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
    mpz_mod(res, res, fpchar);
} 

// add uint: result = a + (unsigned int) b (mod p)
void fp_add_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_add_ui(res, a, b);
    mpz_mod(res, res, fpchar);
}

// sub: res = a - b (mod p)
void fp_sub(fp_t res, const fp_t a, const fp_t b) {
    mpz_sub(res, a, b);
    mpz_mod(res, res, fpchar);
}

// sub uint: res = a - (unsigned int) b (mod p)
void fp_sub_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_sub_ui(res, a, b);
    mpz_mod(res, res, fpchar);
}

// mul: res = a * b (mod p)
// This function is argument-safe and can be called 
// with: fp_mul(n, n, n), where n is the same variable
void fp_mul(fp_t res, const fp_t a, const fp_t b) {
    mpz_mul(res, a, b);
    mpz_mod(res, res, fpchar);
}

// mul int: res = a * (int) b (mod p)
void fp_mul_int(fp_t res, const fp_t a, long int b) {
    mpz_mul_si(res, a, b);
    mpz_mod(res, res, fpchar);
}

// modular inverse: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a) {
    mpz_invert(res, a, fpchar);
}

// div: a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b) {
    mpz_invert(res, b, fpchar);
    mpz_mul(res, res, a);
    mpz_mod(res, res, fpchar);
}


// neg: a = -a (mod p)
void fp_neg(fp_t res, const fp_t a) {
    mpz_sub(res, fpchar, a);
    mpz_mod(res, res, fpchar);
}

// sqrt: 
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a) {
    mpz_t exp; 
    mpz_init(exp);

    // calculate a^((p + 1)/4)

    // step 1: e = p + 1
    mpz_add_ui(exp, fpchar, 1);

    // step 2: e = (p + 1)/4
    // can use divexact only if we know the divisor: d = 4
    mpz_divexact_ui(exp, exp, 4);

    // step 3: res = a ^ e
    mpz_powm(res, a, exp, fpchar);

    // step 4: clear memory for e variable
    mpz_clear(exp);
}

// =========== 
// FP2 methods
// ===========

// init: allocate memory and set value to 0 + 0i
void fp2_init(fp2_t *res) {
    // allocate memory for the fp2 structure
    *res = (fp2_t) malloc(sizeof(fp2_elem));
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
// This function is **not** argument-safe, i.e.
// calling fp_mul(a, a, a) will provide incorrect results.
// Therefore, res != lhs and res != rhs, but calling (a, n, n) is ok.
void fp2_mul_unsafe(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    assert(res != lhs && res != rhs && "fp2_mul cannot be called with res = lhs or res = rhs");
    // add temprary field elements for storing values
    fp_t t;
    fp_init(t); // allocate and set to 0

    // (a + bi) * (c + di) = [ac - bd] + [ad + bc]i

    // calculate the first variable (real part); ac - bd 
    // res[0] <- a * c
    fp_mul(res->a, lhs->a, rhs->a);
    // temp <- b * d
    fp_mul(t, lhs->b, rhs->b);
    // res[0] <- [res0](a * c) - [tmp](b * d)
    fp_sub(res->a, res->a, t);

    // calculate the second variable (real part); ad - bc
    // res[1] <- a * d
    fp_mul(res->b, lhs->a, rhs->b);
    // TODO: make sure we can set the register without UB
    // i.e-> there is no garbage left from other operations
    // so no need to set value to 0?

    // tmp <- b * c
    fp_mul(t, lhs->b, rhs->a);
    fp_add(res->b, res->b, t);

    fp_clear(t);
}

void fp2_mul_safe(fp2_t res, const fp2_t lhs, const fp2_t rhs) {
    // x = a + bi, y = c + di
    // r = (ac - bd) + (ad + bc)i = e + fi
    // cost: 4m + 2a

    fp_t t0, t1, tc;
    fp_init(t0); fp_init(t1); fp_init(tc);

    // tc is used for storing value c in case 
    // it will be overwritten in step: e = ac 
    // if fp2_mul is called with r = x = y
    fp_set(tc, rhs->raw[0]);             // tc = c
    fp_mul(t0, lhs->raw[1], rhs->raw[1]);     // t0 = bd
    fp_mul(t1, lhs->raw[0], rhs->raw[1]);     // t1 = ad

    fp_mul(res->raw[0], lhs->raw[0], rhs->raw[0]); // e = ac
    fp_sub(res->raw[0], res->raw[0], t0);     // e = ac - bd

    fp_mul(res->raw[1], lhs->raw[1], tc);     // f = bc 
    fp_add(res->raw[1], res->raw[1], t1);     // f = bc + ad

    fp_clear(t0); fp_clear(t1); fp_clear(tc);
}

// sq: result <- arg^2
// This function cannot be used when res == arg
void fp2_sq_unsafe(fp2_t res, const fp2_t arg) {
    assert(res != arg && "fp2_sq cannot be called with res = arg");
    fp2_mul_unsafe(res, arg, arg);
}

// sq: result <- arg^2
// This function can be used safely with res == arg
void fp2_sq_safe(fp2_t res, const fp2_t arg) {
    fp2_mul_safe(res, arg, arg);
}

void fp2_print_uint(fp2_t arg, const char* name) {
    printf("%s = %ld + %ld * i\n", name, mpz_get_ui(arg->a), mpz_get_ui(arg->b));
}