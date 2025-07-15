#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "isog_mont.h"
#include "ec_mont.h"

void criss_cross(fp2_t lsum, fp2_t rdiff, const fp2_t x, const fp2_t y,
                 const fp2_t z, const fp2_t w) {
    // Calculate (xw + yz, xw - yz) given (x, y, z, w)
    // Argument-safe: Yes
    // Registers: 2
    // Cost: 2M + 2a

    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

    // t0 = x * w
    fp2_mul_unsafe(t0, x, w);
    // t1 = y * z
    fp2_mul_unsafe(t1, y, z);

    // rdiff = t0 - t1: xw - yz
    fp2_sub(rdiff, t0, t1);
    // lsum = t0 + t1: xw + yz
    fp2_add(lsum, t0, t1);

    fp2_clear(&t0);
    fp2_clear(&t1);
}

void KPS(point_t *kpts, size_t n, const point_t K, const fp2_t A24p,
         const fp2_t C24) {
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
    fp2_t temp;
    fp2_init(&temp);
    for (size_t i = 0; i < n; i++) {
        // (X : Z) -> (X + Z, X - Z)
        fp2_add(temp, kpoints[i]->X, kpoints[i]->Z);
        fp2_sub(kpoints[i]->Z, kpoints[i]->X, kpoints[i]->Z);
        fp2_set(kpoints[i]->X, temp);
    }
    fp2_clear(&temp);
}

void xISOG_odd(point_t Q, const point_t *prep_kpts, size_t n, const point_t P) {
    assert(n > 0 && prep_kpts != NULL &&
           "List of kernel points cannot be empty");
    // Registers: 4

    fp2_t t0, t1, t2, t3;
    fp2_init(&t0);
    fp2_init(&t1);
    fp2_init(&t2);
    fp2_init(&t3);

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

void aISOG_curve_KPS(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24,
                     const point_t *kpts, size_t n) {

    fp2_t sigma, sigma_inv, pi;
    fp2_init(&sigma);
    fp2_init(&sigma_inv);
    fp2_init(&pi);

    // pi is equal to product of points x-coordinates, therefore it must be
    // initialized with 1
    fp2_set_uint(pi, 1);

    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

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

    fp2_clear(&sigma);
    fp2_clear(&sigma_inv);
    fp2_clear(&pi);
    fp2_clear(&t0);
    fp2_clear(&t1);
}

void aISOG_curve(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24,
                 const point_t K, int degree) {

    size_t n = KPS_DEG2SIZE(degree);
    point_t *kpts = calloc(n, sizeof(point_t));

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
    fp2_init(&t0);
    fp2_init(&t1);
    fp2_init(&t2);
    fp2_init(&t3);
    fp2_sub(t0, P->X, P->Z); // t0: XP - ZP
    fp2_add(t1, P->X, P->Z); // t1: XP + ZP
    fp2_sub(t2, K->Z, K->X); // t2: ZK - XK
    fp2_add(t3, K->Z, K->X); // t3: ZK + XK

    // Z = (XP - ZP)(ZK + XK) + (XP + ZP)(ZK - XK)
    // X = (XP - ZP)(ZK + XK) - (XP + ZP)(ZK - XK)
    criss_cross(Q->Z, Q->X, t0, t1, t2, t3);

    fp2_mul_safe(Q->X, P->X);
    fp2_mul_safe(Q->Z, P->Z);

    fp2_clear(&t0);
    fp2_clear(&t1);
    fp2_clear(&t2);
    fp2_clear(&t3);
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
    fp2_init(&t0);
    fp2_init(&t1);
    fp2_sub(t0, P->X, P->Z); // t0: XP - ZP
    fp2_add(t1, P->X, P->Z); // t1: XP + ZP

    // Z = (XP - ZP)(ZK + XK) + (XP + ZP)(ZK - XK)
    // X = (XP - ZP)(ZK + XK) - (XP + ZP)(ZK - XK)
    criss_cross(Q->Z, Q->X, t0, t1, prep_K->Z, prep_K->X);

    fp2_mul_safe(Q->X, P->X);
    fp2_mul_safe(Q->Z, P->Z);

    fp2_clear(&t0);
    fp2_clear(&t1);
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

void ISOG2e(fp2_t A24p, fp2_t C24, const fp2_t A24p_init, const fp2_t C24_init,
            const point_t K, uint32_t e, point_t *push_points) {

    point_t K0, T, R;
    point_init(&K0);
    point_init(&T);
    point_init(&R);

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
    point_clear(&K0);
    point_clear(&T);
    point_clear(&R);
}

void ISOG_chain(fp2_t A24p, fp2_t C24, const fp2_t A24p_init,
                const fp2_t C24_init, const point_t K, pprod_t isog_degree,
                point_t *push_points) {

    fp2_t A24p_next, C24_next;
    fp2_init(&A24p_next);
    fp2_init(&C24_next);

    fp2_set(A24p, A24p_init);
    fp2_set(C24, C24_init);

    point_t K0, Q, T, R;
    point_init(&K0);
    point_init(&Q);
    point_init(&T);
    point_init(&R);

    point_set(K0, K);

    // Push K0 into the push_points list, as first occured "NULL"
    point_t *pp = push_points;
    while (*pp != NULL)
        pp++;
    // Make sure that push_points end with 2x NULL
    assert(*pp == NULL && *(pp + 1) == NULL);
    // Replace the first NULL with K0
    *pp = K0;

    // TODO: optimize: calculate MAX out of degree->div and allocate space
    // maybe store it inside pprod structure?
    unsigned int max_div = 0;
    for (unsigned int i = 0; i < isog_degree->n_primes; i++) {
        max_div =
            isog_degree->primes[i] > max_div ? isog_degree->primes[i] : max_div;
    }

    size_t max_n = KPS_DEG2SIZE(max_div);
    point_t *kpts = calloc(max_n, sizeof(point_t));

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
        // We could "prepare" the kernel points if it is required more than
        // once?
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

    assert(fp2_is_zero(K0->Z) &&
           "Kernel of the isogeny should end-up as E(0) - E0->Z = 0");

    for (size_t i = 0; i < max_n; i++) {
        point_clear(&kpts[i]);
    }
    free(kpts);

    fp2_clear(&A24p_next);
    fp2_clear(&C24_next);
    point_clear(&K0);
    point_clear(&Q);
    point_clear(&T);
    point_clear(&R);
}