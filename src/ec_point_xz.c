#include <stdlib.h>

#include "fp2.h"
#include "ec_point_xz.h"

// Elliptic Curve coefficient related conversion methods
void A24p_from_A(fp2_t A24p, fp2_t C24, const fp2_t A, const fp2_t C) {
    // Set A24p := A + 2C and C24 := 4C
    fp2_set(A24p, A);               // A24p = A
    fp2_add(C24, C, C);                 // C24 = 2C
    fp2_add(A24p, A24p, C24);   // A24p = A + 2C
    fp2_add(C24, C24, C24);             // C24 = 4C
}

void A_from_A24p(fp2_t A, fp2_t C, const fp2_t A24p, const fp2_t C24) {
    // Set A := 4*A24p - 2*C24 and C := C24
    fp2_set(C, C24);

    // A = 4A'
    fp2_mul_int(A, A24p, 4);

    // A = 4A' - 2C'
    fp2_sub(A, A, C);
    fp2_sub(A, A, C);
}

// Elliptic Curve Point related methods

void point_init(point_t *P) {
    *P = (point_t) malloc(sizeof(struct point_xz));
    // Recursively allocate memory for the member pointers
    fp2_init(&(*P)->X);
    fp2_init(&(*P)->Z);
}

void point_clear(point_t *P) {
    if (P) {
        fp2_clear(&(*P)->X);
        fp2_clear(&(*P)->Z);
        free(*P);
        P = NULL;
    }
}

// copy coordinates of point: R <- P
void point_set(point_t R, const point_t P) {
    // TODO: Not a deepcopy -> only copies the pointers
    // memcpy(R, P, sizeof(struct point_xz));
    fp2_set(R->X, P->X);
    fp2_set(R->Z, P->Z);
}

// Set: x(P) = x, and z(P) = 1
void point_set_str_x(point_t P, const char *x) {
    fp2_set_str(P->X, x);
    fp2_set_uint(P->Z, 1);
}


// Set: x(P) = x, and z(P) = 1
void point_set_fp2_x(point_t P, fp2_t x) {
    fp2_set(P->X, x);
    fp2_set_uint(P->Z, 1);
}

void point_printx(point_t P, const char* name) {
    // Make sure that the coordinate is normalized, otherwise we get false results when Z != 1: x = (X : Z)
    assert(point_is_normalized(P));
    fp2_print(P->X, name);
}


void point_printx_normalized(point_t P, const char* name) {
    // Make sure that the coordinate is normalized, otherwise we get false results when Z != 1: x = (X : Z)
    point_normalize_coords(P);
    fp2_print(P->X, name);
}

int point_equal_str_x(point_t P, const char* str) {
    point_t P_str;
    point_init(&P_str);

    // Make sure that the Point is normalized, i.e: P->Z = 1
    point_normalize_coords(P);

    point_set_str_x(P_str, str);
    int equal = fp2_equal(P->X, P_str->X) && fp2_equal(P->Z, P_str->Z);

    point_clear(&P_str);
    return equal;
}

int point_is_normalized(point_t P) {
    return fp2_equal_uint(P->Z, 1);
}

// Normalize: P = (X : Z) -> (X' : 1) with X' = X/Z
void point_normalize_coords(point_t P) {
    assert(!fp2_is_zero(P->Z) && "Normalized Point cannot have Z = 0");
    // Registers: 1

    fp2_t t; fp2_init(&t);

    fp2_div_unsafe(t, P->X, P->Z);
    fp2_set(P->X, t);
    fp2_set_uint(P->Z, 1);

    fp2_clear(&t);
}

