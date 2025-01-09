#include <stdio.h>
#include "ec_mont.h"
#include <assert.h>

int main() {

    // Set field characteristic
    assert(!global_fpchar_setup_uint(431));

    // Initialize the Elliptic Curve E: y^2 = x^3 + 6x^2 + x
    // (A : C) are projective coordinates of (a : 1) = (6 : 1)
    // We additionally define A24+, C24 and A24- as:
    // A24+: A + 2C
    // C24:  4C
    // A24-: A - 2C
    // This constants are used to speed up the calculations inside 
    // certain functions e.g. xDBL algorithm

    // (A, C) = (6, 1)
    fp2_t A, C;
    fp2_init(&A);
    fp2_init(&C);
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    fp2_print_uint(A, "A");
    fp2_print_uint(C, "C");

    // (A24+ : C24) ~ (A + 2C : 4C)
    // This projective pair is the representant of the variable (a + 2)/4
    fp2_t A24_plus, C24;
    fp2_init(&A24_plus);
    fp2_init(&C24);

    // Set A24p := A + 2C and C24 := 4C
    fp2_set(A24_plus, A); // A24p = A
    fp2_add(C24, C, C);   // C24 = 2C
    fp2_add(A24_plus, A24_plus, C24); // A24p = A + 2C
    fp2_add(C24, C24, C24); // C24 = 4C

    fp2_print_uint(A24_plus, "A24p");
    fp2_print_uint(C24, "C24");

    // x-coordinate for P: x(P)

    // fp2_t xP;
    // fp2_init(&xP);

    // P = (15*i + 292 : 235*i + 281 : 1)
    point_t P;
    point_init(&P);

    // Px = 292 + 15i
    fp_set_uint(P->X->a, 292);
    fp_set_uint(P->X->b, 15);
    // Pz = 1
    fp_set_uint(P->Z->a, 1);

    fp2_print_uint(P->X, "xP");
    fp2_print_uint(P->Z, "zP");

    // Q = [2]P
    { 
        point_t Q;
        point_init(&Q);

        xDBL(Q, P, A24_plus, C24);
        
        fp2_print_uint(Q->X, "xQ");
        fp2_print_uint(Q->Z, "zQ");

        // Answer in affine coordinates: x(Q) = 61 + 184 * i
        // Which is equal to projective coordinates (this implementation)
        // X(Q) = 157 + 180 * i
        // Z(Q) = 65 + 358 * i
        // X(Q)/Z(Q) = 61 + 184 * i = x(Q)

        assert(!mpz_cmp_ui(Q->X->a, 157) && !mpz_cmp_ui(Q->X->b, 180));
        assert(!mpz_cmp_ui(Q->Z->a, 65) && !mpz_cmp_ui(Q->Z->b, 358));

        point_clear(&Q);
    }

    fp2_clear(&A);
    fp2_clear(&C);
    fp2_clear(&A24_plus);
    fp2_clear(&C24);
    point_clear(&P);

    assert(!global_fpchar_clear());

    return 0;
}