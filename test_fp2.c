#include <stdio.h>
#include <assert.h>

#include "fp2.h"

int main() {
    // init characteristic
    fp_t p;
    fp_init(p);
    fp_set_uint(p, 431);
    fp_set_char(p);

    fp2_t x;
    fp2_init(x);

    printf("x: %ld, %ldi\n", mpz_get_ui(x.a), mpz_get_ui(x.b));
    printf("x: %ld, %ldi\n", mpz_get_ui(x.elem[0]), mpz_get_ui(x.elem[1]));


    return 0;
}