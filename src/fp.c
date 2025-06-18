#include <gmp.h>
#include <stdio.h>
#include <assert.h>

#include "fp.h"

static fp_t g_fpchar;
static int g_is_fpchar_set = 0;

int fpchar_clear_if_set() {
    if (g_is_fpchar_set) { 
        fpchar_clear();
        return 1;
    } else {
        return 0;
    }
}

int fpchar_check() {
    return g_is_fpchar_set;
}

int fpchar_setup_uint(unsigned int p) {
    fp_t pp;
    fp_init(pp);
    fp_set_uint(pp, p);
    int retval = fpchar_setup(pp);
    fp_clear(pp);
    return retval;
}

// Initialize and set the field characteristic of Fp and Fp^2 arithmetic.
// For now, only p = 3 (mod 4) is allowed
int fpchar_setup(fp_t p) {
    if (g_is_fpchar_set) {
        fprintf(stderr, "[!] trying to initialize the characteristic second time!\n");
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
        fprintf(stderr, "[!] Trying to set prime field characteristic with prime != 3 (mod 4)\n");
        assert(0);
    }

    fp_init(g_fpchar);
    fp_set(g_fpchar, p);
    g_is_fpchar_set = 1;
    return 0;
}

int fpchar_clear() {
    if (!g_is_fpchar_set) {
        fprintf(stderr, "[!] trying to clear not-initialized characteristic!\n");
        return -1;
    }
    fp_clear(g_fpchar);
    g_is_fpchar_set = 0;
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
    mpz_mod(res, res, g_fpchar);
} 

// add uint: result = a + (unsigned int) b (mod p)
void fp_add_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_add_ui(res, a, b);
    mpz_mod(res, res, g_fpchar);
}

// sub: res = a - b (mod p)
void fp_sub(fp_t res, const fp_t a, const fp_t b) {
    mpz_sub(res, a, b);
    mpz_mod(res, res, g_fpchar);
}

// sub uint: res = a - (unsigned int) b (mod p)
void fp_sub_uint(fp_t res, const fp_t a, unsigned long int b) {
    mpz_sub_ui(res, a, b);
    mpz_mod(res, res, g_fpchar);
}

// mul: res = a * b (mod p)
// This function is argument-safe and can be called 
// with: fp_mul(n, n, n), where n is the same variable
void fp_mul(fp_t res, const fp_t a, const fp_t b) {
    mpz_mul(res, a, b);
    mpz_mod(res, res, g_fpchar);
}

// mul int: res = a * (int) b (mod p)
void fp_mul_int(fp_t res, const fp_t a, long int b) {
    mpz_mul_si(res, a, b);
    mpz_mod(res, res, g_fpchar);
}

// modular inverse: res = a^-1 (mod p)
void fp_inv(fp_t res, const fp_t a) {
    mpz_invert(res, a, g_fpchar);
}

// div: a / b (mod p) = a * b^-1 (mod p)
void fp_div(fp_t res, const fp_t a, const fp_t b) {
    mpz_invert(res, b, g_fpchar);
    mpz_mul(res, res, a);
    mpz_mod(res, res, g_fpchar);
}


// neg: a = -a (mod p)
void fp_neg(fp_t res, const fp_t a) {
    mpz_sub(res, g_fpchar, a);
    mpz_mod(res, res, g_fpchar);
}

// sqrt: 
// assumes that prime is in form: p = 3 (mod 4)
void fp_sqrt(fp_t res, const fp_t a) {
    mpz_t exp; 
    mpz_init(exp);

    // calculate a^((p + 1)/4)

    // step 1: e = p + 1
    mpz_add_ui(exp, g_fpchar, 1);

    // step 2: e = (p + 1)/4
    // can use divexact only if we know the divisor: d = 4
    mpz_divexact_ui(exp, exp, 4);

    // step 3: res = a ^ e
    mpz_powm(res, a, exp, g_fpchar);

    // step 4: clear memory for e variable
    mpz_clear(exp);
}

int fp_equal_uint(fp_t a, unsigned long int b) {
    return !mpz_cmp_ui(a, b);
}

int fp_equal(fp_t a, fp_t b) {
    return !mpz_cmp(a, b);
}

int fp_equal_str(fp_t a, const char *b_str) {
    fp_t b; 
    mpz_init_set_str(b, b_str, 0);
    int equal = fp_equal(a, b);
    fp_clear(b);
    return equal;
}

void fp_print(fp_t a, const char* name) {
    gmp_printf("%s: %Zd\n", name, a);
}
