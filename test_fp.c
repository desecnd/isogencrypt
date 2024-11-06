#include <stdio.h>
#include <assert.h>

#include "fp2.h"

int main() {
    // init characteristic
    fp_t p, r;
    fp_init(p);
    fp_init(r);

    fp_set_uint(p, 431);
    fp_set_char(p);

    fp_set_uint(r, 16); // r == 16
    assert(!mpz_cmp_ui(r, 16));

    fp_sqrt(r, r); // r == 4
    assert(!mpz_cmp_ui(r, 4));

    fp_add_uint(r, r, 10); // r == 14
    assert(!mpz_cmp_ui(r, 14));

    fp_sub_uint(r, r, 15); // r == -1 == 430
    assert(!mpz_cmp_ui(r, 430));

    fp_inv(r, r); // r == r ^-1 == 430
    assert(!mpz_cmp_ui(r, 430));

    fp_clear(p);
    fp_clear(r);

    return 0;
}