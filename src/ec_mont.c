#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "ec_mont.h"
#include "fp2.h"


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

void xDBL(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24) {
    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

    // P = (X : Z)
    fp2_sub(t0, P->X, P->Z);    // t0 = X - Z
    fp2_add(t1, P->X, P->Z);    // t1 = X + Z

    // Warning! Must use "safe" function for square-ing
    fp2_sq_safe(t0);             // t0 = (X - Z)^2 
    fp2_sq_safe(t1);             // t1 = (X + Z)^2

    fp2_mul_unsafe(R->Z, t0, C24);     // Z' = (X - Z)^2 * C24
    fp2_mul_unsafe(R->X, R->Z, t1);    // X' = (X + Z)^2 * (X - Z)^2 * C24

    fp2_sub(t1, t1, t0);        // t1 = (X + Z)^2 - (X - Z)^2 
    fp2_mul_unsafe(t0, A24p, t1);  // t0 = A24p * [(X + Z)^2 - (X - Z)^2]

    fp2_add(R->Z, R->Z, t0);    // Z' = A24p * [(X + Z)^2 - (X - Z)^2] + (X - Z)^2 * C24
    fp2_mul_safe(R->Z, t1);    // Z' = { A24p * [(X + Z)^2 - (X - Z)^2] + C24 (X - Z)^2] } [(X + Z)^2 - (X - Z)^2]

    fp2_clear(&t0);
    fp2_clear(&t1);
}

void xDBLe(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24, const int e) {
    point_set(R, P);
    // Repeat the step of doubling multiple times
    for (int i = 0; i < e; i++) {
        xDBL(R, R, A24p, C24);
    }
}

void xADD(point_t PQsum, const point_t P, const point_t Q, const point_t PQdiff) {
    // Function is argument safe for calling PQSum = Q or P or PQdiff
    // Given P, Q and PQdiff = P-Q output P+Q

    // xADD formula in projective coordinates: P + Q = (X' : Z')
    // X' = ZR * (XP * XQ - ZP * ZQ)^2 
    // Z' = XR * (XP * ZQ - XQ * ZP)^2
    // 
    // Optimised version:
    // X' = ZR * [(XP - ZP)(XQ + ZQ) + (XP + ZP)(XQ - ZQ)]^2 
    // Z' = XR * [(XP - ZP)(XQ + ZQ) - (XP + ZP)(XQ - ZQ)]^2 

    fp2_t t0, t1, t2, t3;
    fp2_init(&t0); fp2_init(&t1); fp2_init(&t2); fp2_init(&t3); 

    // Defining additional variables:
    // a: XP + ZP
    // b: XP - ZP
    // c: XQ + ZQ
    // d: XQ - ZQ
    //
    // X' = ZR * (b * c + a * d)^2 
    // Z' = XR * (b * c - a * d)^2 
    fp2_add(t0, P->X, P->Z);        // t0 = a: XP + ZP
    fp2_sub(t1, P->X, P->Z);        // t1 = b: XP - ZP
    fp2_add(t2, Q->X, Q->Z);        // t2 = c: xQ + zQ
    fp2_sub(t3, Q->X, Q->Z);        // t3 = d: xQ - zQ

    fp2_mul_safe(t0, t3);           // t0 = t0 * t3: a * d
    fp2_mul_safe(t1, t2);           // t1 = t1 * t2: b * c

    fp2_add(t2, t0, t1);            // t2 = t0 + t1: ad + bc
    fp2_sub(t3, t0, t1);            // t3 = t0 - t1: ad - bc

    fp2_sq_unsafe(t0, t2);          // t0 = t2^2: (ad + bc)^2
    fp2_sq_unsafe(t1, t3);          // t1 = t3^2: (ad - bc)^2

    // We must use t2 to not override the used later PQdiff->X in scenario when PQsum = PQdiff 
    // t2 = t0 * Z(P-Q): (ad + bc)^2 * Z(P-Q)
    fp2_mul_unsafe(t2, t0, PQdiff->Z);          
    // Z' = t1 * X(P-Q): (ad - bc)^2 * Z(P-Q)
    fp2_mul_unsafe(PQsum->Z, t1, PQdiff->X);   
    // X' = t2
    fp2_set(PQsum->X, t2);

    fp2_clear(&t0); fp2_clear(&t1); fp2_clear(&t2); fp2_clear(&t3);
}

// Out: P = 2P, Q = P + Q
// TODO: check argument-safeness
// TODO: think about expressing function 
// as set of arithmetic instructions 
// instead of combinatin of 2 function calls
void xDBLADD(point_t P, point_t Q, const point_t PQdiff, const fp2_t A24p, const fp2_t C24) {
    xADD(Q, P, Q, PQdiff);
    xDBL(P, P, A24p, C24);
}

void criss_cross(fp2_t lsum, fp2_t rdiff, const fp2_t x, const fp2_t y, const fp2_t z, const fp2_t w) {
    // Calculate (xw + yz, xw - yz) given (x, y, z, w)
    // Argument-safe: Yes
    // Registers: 2
    // Cost: 2M + 2a

    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1);

    // t0 = x * w
    fp2_mul_unsafe(t0, x, w);
    // t1 = y * z
    fp2_mul_unsafe(t1, y, z);
    
    // rdiff = t0 - t1: xw - yz
    fp2_sub(rdiff, t0, t1);
    // lsum = t0 + t1: xw + yz
    fp2_add(lsum, t0, t1);
    
    fp2_clear(&t0); fp2_clear(&t1);
}

void xLADDER(point_t R0, const point_t P, const mpz_t m, const fp2_t A24p, const fp2_t C24) {
    assert(mpz_sgn(m) > 0 && "Given scalar m must be nonnegative");
    
    mpz_t r;
    mpz_init(r);

    point_t R1;
    point_init(&R1);

    // R0 = P, R1 = [2]R
    point_set(R0, P);
    xDBL(R1, P, A24p, C24);

    // Get number of "active" bits
    int n_bits = mpz_sizeinbase(m, 2);

    // Iterate over bits downwards (leading bit not included).
    // Invariant of the algorithm R1 - R0 = P
    for (int bit = n_bits - 2; bit >= 0; bit--) {
        // Bit is equal to 1
        if (mpz_tstbit(m, bit)) {
            // R1 = [2]R1; R0 = R0 + R1
            xDBLADD(R1, R0, P, A24p, C24);
        } else {
            // R0 = [2]R0; R1 = R0 + R1
            xDBLADD(R0, R1, P, A24p, C24);
        }
    }

    mpz_clear(r);
    point_clear(&R1);
}

void xLADDER_int(point_t R0, const point_t P, long int m, const fp2_t A24p, const fp2_t C24) {
    assert(m > 0 && "Given scalar m must be nonnegative");

    point_t R1;
    point_init(&R1);

    // R0 = P, R1 = [2]R
    point_set(R0, P);
    xDBL(R1, P, A24p, C24);

    // Get number of "active" bits (count until leading bit is found)
    int bits = 0;
    for (long int x = m; x > 0; x /= 2) bits++;

    // Iterate over bits downwards (leading bit not included).
    // Invariant of the algorithm R1 - R0 = P
    for (int bit = bits - 2; bit >= 0; bit--) {

        // bit = 1
        if (m & (1 << bit)) {
            // R1 = [2]R1; R0 = R0 + R1
            xDBLADD(R1, R0, P, A24p, C24);
        } else {
            // R0 = [2]R0; R1 = R0 + R1
            xDBLADD(R0, R1, P, A24p, C24);
        }
    }

    point_clear(&R1);
}


// calculate P = P + [m]Q
void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m, const fp2_t A24p, const fp2_t C24) {
    assert(m >= 0 && "Given scalar m must be nonnegative");

    while (m > 0) {
        if (m & 1) 
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else       
            xDBLADD(Q, PQdiff, P, A24p, C24);
        m /= 2;
    }
}

// calculate P = P + [m]Q
void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, const mpz_t m, const fp2_t A24p, const fp2_t C24) {
    assert(mpz_sgn(m) >= 0 && "Given scalar m must be nonnegative");

    mpz_t r, n;
    mpz_init(r);
    mpz_init_set(n, m);

    while (mpz_sgn(n) > 0) {
        // m = m//2; r = m % 2
        mpz_fdiv_qr_ui(n, r, n, 2);

        if (mpz_sgn(r) > 0)
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else 
            xDBLADD(Q, PQdiff, P, A24p, C24);
    }

    mpz_clear(r);
    mpz_clear(n);
}

void KPS(point_t * kpts, size_t n, const point_t K, const fp2_t A24p, const fp2_t C24) {
    // "deepcopy" generator point as the first [1]K point
    point_set(kpts[0], K);

    // Calculate the second as simple [2]K
    if (n >= 2) {
        xDBL(kpts[1], K, A24p, C24);
    }

    // Calculate next using: [i]K = [i - 1]K + K
    // To get difference: [i - 1]K - K = [i - 2]K 
    for (size_t i = 2; i < n; i++) {
        xADD(kpts[i], kpts[i - 1], kpts[0], kpts[i - 2]);
    }
}

void prepare_kernel_points(point_t *kpoints, size_t n) {
    fp2_t temp; fp2_init(&temp);
    for (size_t i = 0; i < n; i++) {
        // (X : Z) -> (X + Z, X - Z)
        fp2_add(temp, kpoints[i]->X, kpoints[i]->Z);
        fp2_sub(kpoints[i]->Z, kpoints[i]->X, kpoints[i]->Z);
        fp2_set(kpoints[i]->X, temp);
    }
    fp2_clear(&temp);
}

void xISOG_odd(point_t Q, const point_t *prep_kpts, size_t n, const point_t P) {
    assert(n > 0 && prep_kpts != NULL && "List of kernel points cannot be empty");
    // Registers: 4

    fp2_t t0, t1, t2, t3;
    fp2_init(&t0); fp2_init(&t1);
    fp2_init(&t2); fp2_init(&t3);

    // Prepare point P = (X : Z) => P^ = (X + Z : X - Z)
    // Unfortunately we cannot modify P - so it requires additional 2 registers
    // t2 = XP + ZP: XP^ (hat)
    fp2_add(t2, P->X, P->Z);
    // t3 = XP - ZP: ZP^ (hat)
    fp2_sub(t3, P->X, P->Z);

    // By XP^ we represent the "prepared" variant
    // X' = [(XK + ZK)(XP - ZP) + (XK - ZK)(XP + ZP)] = [XK^ * ZP^ + ZK^ * XP^]
    // Z' = [(XK + ZK)(XP - ZP) - (XK - ZK)(XP + ZP)] = [XK^ * ZP^ - ZK^ * XP^]
    criss_cross(Q->X, Q->Z, prep_kpts[0]->X, prep_kpts[0]->Z, t2, t3);

    // Multiply X' and Z' by same formula for Ki
    for (size_t i = 1; i < n; i++) {
        criss_cross(t0, t1, prep_kpts[i]->X, prep_kpts[i]->Z, t2, t3);
        fp2_mul_safe(Q->X, t0);
        fp2_mul_safe(Q->Z, t1);
    }

    // t0 = XQ^2: prod(i: XKi_ * ZP_ + ZKi_ * XP_)^2
    fp2_sq_unsafe(t0, Q->X);
    // t1 = ZQ^2: prod(i: XKi_ * ZP_ - ZKi_ * XP_)^2
    fp2_sq_unsafe(t1, Q->Z);
    
    // XQ = XP * t0: XP * prod(i: XKi_ * ZP_ - ZKi_ * XP_)^2
    fp2_mul_unsafe(Q->X, t0, P->X);
    // ZQ = ZP * t1: ZP * prod(i: XKi_ * ZP_ - ZKi_ * XP_)^2
    fp2_mul_unsafe(Q->Z, t1, P->Z);

    fp2_clear(&t0);
    fp2_clear(&t1);
    fp2_clear(&t2);
    fp2_clear(&t3);
}

void aISOG_curve_KPS(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24, const point_t * kpts, size_t n) {

    fp2_t sigma, sigma_inv, pi;
    fp2_init(&sigma); fp2_init(&sigma_inv); fp2_init(&pi);

    // pi is equal to product of points x-coordinates, therefore it must be initialized with 1
    fp2_set_uint(pi, 1);

    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1);

    for (size_t i = 0; i < n; i++) {
        // x(P) = X/Z
        fp2_div_unsafe(t0, kpts[i]->X, kpts[i]->Z);
        // x(P)^-1 = (X/Z)^-1 = (Z/X)
        fp2_div_unsafe(t1, kpts[i]->Z, kpts[i]->X);
        
        // sigma += x([i]K)
        fp2_add(sigma, sigma, t0);
        // pi *= x([i]K)
        fp2_mul_safe(pi, t0);

        // sigma_inv += x([i]K)^-1
        fp2_add(sigma_inv, sigma_inv, t1);
    }

    // Obtain original coordinates (A:C) from (A24:C24)
    // use (t0 : t1) as registers
    A_from_A24p(t0, t1, A24p, C24);

    // Use "C_" as register to hold "a" = (a:1) = (A:C) value
    // C_ = t0 / t1: A / C = a
    fp2_div_unsafe(C_, t0, t1);

    // t0 = sigma_inv - sigma
    fp2_sub(t0, sigma_inv, sigma);

    // t0 = t0 * 6: A_ = 6(sigma_inv - sigma) = 6sigma_inv - 6sigma
    fp2_mul_int(t0, t0, 6);

    // t0:: t0 + C_: 6sigma_inv - 6sigma + A/C
    // use value stored in C_ register := A/C = a
    fp2_add(t0, t0, C_);

    // t1 = pi^2
    fp2_sq_unsafe(t1, pi);

    // A_ = t0 * t1: (6sigma_inv - 6sigma + A/C) * pi^2
    fp2_mul_unsafe(A_, t0, t1);
    // C_ = 1
    fp2_set_uint(C_, 1);

    fp2_clear(&sigma); fp2_clear(&sigma_inv); fp2_clear(&pi);
    fp2_clear(&t0); fp2_clear(&t1);
}

void aISOG_curve(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24, const point_t K, int degree) {

    size_t n = KPS_DEG2SIZE(degree);
    point_t * kpts = calloc(sizeof(point_t), n);

    // Initialize kernel points
    for (size_t i = 0; i < n; i++) {
        point_init(&kpts[i]);
    }

    // Calculate [1]K, [2]K, [3]K, ...
    KPS(kpts, n, K, A24p, C24);

    aISOG_curve_KPS(A_, C_, A24p, C24, kpts, n);

    for (size_t i = 0; i < n; i++) {
        point_clear(&kpts[i]);
    }
    free(kpts);
}

void xISOG2_unsafe(point_t Q, const point_t K, const point_t P) {
    // Formula works only for K = (x, y=0) where x != 0
    assert(!fp2_is_zero(K->X));

    fp2_t t0, t1, t2, t3;
    fp2_init(&t0); fp2_init(&t1); fp2_init(&t2); fp2_init(&t3);
    fp2_sub(t0, P->X, P->Z); // t0: XP - ZP
    fp2_add(t1, P->X, P->Z); // t1: XP + ZP
    fp2_sub(t2, K->Z, K->X); // t2: ZK - XK
    fp2_add(t3, K->Z, K->X); // t3: ZK + XK

    // Z = (XP - ZP)(ZK + XK) + (XP + ZP)(ZK - XK)
    // X = (XP - ZP)(ZK + XK) - (XP + ZP)(ZK - XK)
    criss_cross(Q->Z, Q->X, t0, t1, t2, t3);

    fp2_mul_safe(Q->X, P->X);
    fp2_mul_safe(Q->Z, P->Z);

    fp2_clear(&t0); fp2_clear(&t1); fp2_clear(&t2); fp2_clear(&t3);
}

void prepare_isog2_kernel(point_t K) {
    // Formula works only for K = (x, y=0) where x != 0
    assert(!fp2_is_zero(K->X));

    fp2_t t;
    fp2_init(&t);

    // t = ZK + XK
    fp2_add(t, K->Z, K->X);
    // ZK = ZK - XK
    fp2_sub(K->Z, K->Z, K->X);
    // XK = t: XK + ZK
    fp2_set(K->X, t);
    // K = (ZK + XK : ZK - XK)

    fp2_clear(&t);
}

void xISOG2_prep(point_t Q, const point_t prep_K, const point_t P) {
    // Assertion for K.x != 0 is present in `prepare_isog2_kernel`.
    // At this point we cannot tell whether K.x = 0.

    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1); 
    fp2_sub(t0, P->X, P->Z); // t0: XP - ZP
    fp2_add(t1, P->X, P->Z); // t1: XP + ZP

    // Z = (XP - ZP)(ZK + XK) + (XP + ZP)(ZK - XK)
    // X = (XP - ZP)(ZK + XK) - (XP + ZP)(ZK - XK)
    criss_cross(Q->Z, Q->X, t0, t1, prep_K->Z, prep_K->X);

    fp2_mul_safe(Q->X, P->X);
    fp2_mul_safe(Q->Z, P->Z);

    fp2_clear(&t0); fp2_clear(&t1);
}

void aISOG2_24p(fp2_t A24p_, fp2_t C24_, const point_t K) {
    // Formula works only for K = (x, 0) where x != 0
    assert(!fp2_is_zero(K->X));
    // (A + 2: 4) = (XK^2 - ZK^2 : ZK^2)

    // A = XK^2
    fp2_sq_unsafe(A24p_, K->X);
    // C = ZK^2
    fp2_sq_unsafe(C24_, K->Z);
    // A = C - A: ZK^2 - XK^2 
    fp2_sub(A24p_, C24_, A24p_);
}

void aISOG2(fp2_t A_, fp2_t C_, const point_t K) {
    // Formula works only for K = (x, 0) where x != 0
    assert(!fp2_is_zero(K->X));

    // A = XK^2
    fp2_sq_unsafe(A_, K->X);
    // A = 2*XK^2
    fp2_add(A_, A_, A_);

    // C = ZK^2
    fp2_sq_unsafe(C_, K->Z);
    
    // A = C - A: ZK^2 - 2XK^2
    fp2_sub(A_, C_, A_);

    // A = 2A: 2(ZK^2 - 2XK^2)
    fp2_add(A_, A_, A_);
}

void ISOG2e(fp2_t A24p, fp2_t C24, const fp2_t A24p_init, const fp2_t C24_init, const point_t K, uint32_t e, point_t* push_points) {

    point_t K0, T, R;
    point_init(&K0); point_init(&T); point_init(&R);

    // K0 = K
    point_set(K0, K);

    // Copy initial curve parameters
    fp2_set(A24p, A24p_init);
    fp2_set(C24, C24_init);

    for (uint32_t i = 0; i < e; i++) {
        // Calculate "local" kernel - T.order() == 2
        if (i + 1 < e) {
            xDBLe(T, K0, A24p, C24, e - i - 1);
        } else {
            point_set(T, K0);
        }

        assert(!fp2_is_zero(T->X) && "Kernel point cannot lay above (0, 0)");

        // Push every point on the list through the partial isogeny
        for (point_t *pptr = push_points; *pptr != NULL; pptr++) {
            xISOG2_unsafe(R, T, *pptr);
            point_set(*pptr, R);
        }

        // Calculate codomain
        aISOG2_24p(A24p, C24, T);
        
        // Push the kernel of [2^e] degree isogeny, only if not last iteration.
        // During the last iteration K0 will go to E(0) 
        if (i + 1 < e) {
            xISOG2_unsafe(R, T, K0);
            point_set(K0, R);
        }
    }
    point_clear(&K0); point_clear(&T); point_clear(&R);
}

void ISOG_chain(fp2_t A24p, fp2_t C24, const fp2_t A24p_init, const fp2_t C24_init, const point_t K, pprod_t isog_degree, point_t *push_points) { 

    fp2_t A24p_next, C24_next;
    fp2_init(&A24p_next); fp2_init(&C24_next);

    fp2_set(A24p, A24p_init);
    fp2_set(C24, C24_init);

    point_t K0, Q, T, R;
    point_init(&K0); point_init(&Q); point_init(&T); point_init(&R);

    point_set(K0, K);

    // Push K0 into the push_points list, as first occured "NULL"
    point_t *pp = push_points;
    while (*pp != NULL) pp++;
    // Make sure that push_points end with 2x NULL
    assert(*pp == NULL && *(pp+1) == NULL);
    // Replace the first NULL with K0
    *pp = K0;

    // TODO: optimize: calculate MAX out of degree->div and allocate space
    // maybe store it inside pprod structure? 
    unsigned int max_div = 0;
    for (unsigned int i = 0; i < isog_degree->n_primes; i++) {
        max_div = isog_degree->primes[i] > max_div ? isog_degree->primes[i] : max_div;
    }

    size_t max_n = KPS_DEG2SIZE(max_div);
    point_t * kpts = calloc(sizeof(point_t), max_n);

    // Initialize kernel points
    for (size_t i = 0; i < max_n; i++) {
        point_init(&kpts[i]);
    }

    // Iterate over all distinct degrees that produce final isogeny
    for (unsigned int i = 0; i < isog_degree->n_primes; i++) {

        // divisor
        unsigned int div = isog_degree->primes[i];
        
        // Calculate the kernel of ith-degree isogeny
        // T = [deg/div]K0 is a point of order "div"
        point_set(T, K0);
        for (unsigned int j = i + 1; j < isog_degree->n_primes; j++) {
            // Ki = [m]Ki;  Ki *= m
            xLADDER_int(Q, T, isog_degree->primes[j], A24p, C24);
            point_set(T, Q);
        }

        // TODO: For now we assume that every component can be 'even'
        // In the future we can store "power_of_two" inside the number
        // We could "prepare" the kernel points if it is required more than once?
        //
        // TODO: Move it as "power2" isogeny function
        if (div % 2 == 0) {
            assert(i == 0 && "Only first number can be a power of 2");
            int log2 = 0;
            while (div > 1) {
                log2++;
                div /= 2;
            }

            // K0 was already appended into the list of push_points
            ISOG2e(A24p_next, C24_next, A24p, C24, T, log2, push_points);
            fp2_set(A24p, A24p_next);
            fp2_set(C24, C24_next);
            continue;
        }

        // Calculate [1]T, [2]T, [3]T ... [div//2]T
        size_t n = KPS_DEG2SIZE(div);
        KPS(kpts, n, T, A24p, C24);

        // Calculate coefficients of the next curve in the isogeny chain
        aISOG_curve_KPS(A24p_next, C24_next, A24p, C24, kpts, n);
        A24p_from_A(A24p, C24, A24p_next, C24_next); 

        // Step required for the multiple calculations of the points
        prepare_kernel_points(kpts, n);

        // K0 is included in the push_points
        for (point_t *pp = push_points; *pp != NULL; pp++) {
            xISOG_odd(Q, kpts, n, *pp);
            point_set(*pp, Q);
        }
    }

    assert(fp2_is_zero(K0->Z) && "Kernel of the isogeny should end-up as E(0) - E0->Z = 0");

    for (size_t i = 0; i < max_n; i++) {
        point_clear(&kpts[i]);
    }
    free(kpts);

    fp2_clear(&A24p_next); fp2_clear(&C24_next);
    point_clear(&K0); point_clear(&Q); point_clear(&T); point_clear(&R);
}

void j_invariant(fp2_t j_inv, const fp2_t A, const fp2_t C) {
    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1);

    // TODO: compare this approach with trying to stay in projective coords until last division

    // Use jinv as register until last line:
    // jinv = A/C: a
    fp2_div_unsafe(j_inv, A, C); 

    // t0 = jinv^2: a^2 
    fp2_sq_unsafe(t0, j_inv);

    // jinv = t0 - 3: a^2 - 3 
    fp2_sub_uint(j_inv, t0, 3);
    
    // t0 = jinv^3: (a^2 - 3)^3
    fp2_sq_unsafe(t0, j_inv);
    fp2_mul_safe(t0, j_inv);
    
    // t1 = jinv - 1: a^2 - 4
    fp2_sub_uint(t1, j_inv, 1);
    assert(!fp2_is_zero(t1) && "A was equal 2 or -2");
    
    // jinv = t0 / t1: (a^2 - 3)^3/(a^2 - 4)
    fp2_div_unsafe(j_inv, t0, t1);

    // jinv = 256 * jinv: 256(a^2 - 3)^3/(a^2 - 4)
    // fp2_mul_int(j_inv, j_inv, 256);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);

    fp2_clear(&t0); fp2_clear(&t1);
}