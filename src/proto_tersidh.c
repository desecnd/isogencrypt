#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ec_mont.h"
#include "ec_point_xz.h"
#include "isog_mont.h"
#include "proto_tersidh.h"

// Catch! The first number is not a prime number => 2^2, this is a special case
// to obtain p == 3 (mod 4)
static unsigned int PRIMES_ALICE[TERSIDH_TMAX] = {
    4,    5,    11,   17,   23,   31,   41,   47,   59,   67,   73,   83,
    97,   103,  109,  127,  137,  149,  157,  167,  179,  191,  197,  211,
    227,  233,  241,  257,  269,  277,  283,  307,  313,  331,  347,  353,
    367,  379,  389,  401,  419,  431,  439,  449,  461,  467,  487,  499,
    509,  523,  547,  563,  571,  587,  599,  607,  617,  631,  643,  653,
    661,  677,  691,  709,  727,  739,  751,  761,  773,  797,  811,  823,
    829,  853,  859,  877,  883,  907,  919,  937,  947,  967,  977,  991,
    1009, 1019, 1031, 1039, 1051, 1063, 1087, 1093, 1103, 1117, 1129, 1153,
    1171, 1187, 1201, 1217, 1229, 1237, 1259, 1279, 1289, 1297, 1303, 1319,
    1327, 1367, 1381, 1409, 1427, 1433, 1447, 1453, 1471, 1483, 1489, 1499,
    1523, 1543, 1553, 1567, 1579, 1597, 1607, 1613, 1621, 1637, 1663, 1669,
    1697, 1709, 1723, 1741, 1753, 1777, 1787, 1801, 1823, 1847, 1867, 1873,
    1879, 1901, 1913, 1933, 1951, 1979, 1993, 1999, 2011, 2027, 2039, 2063,
    2081, 2087, 2099, 2113, 2131, 2141, 2153, 2179, 2207, 2221, 2239, 2251,
    2269, 2281, 2293, 2309, 2333, 2341, 2351, 2371, 2381, 2389, 2399, 2417,
    2437, 2447, 2467, 2477, 2521, 2539, 2549, 2557, 2591, 2609, 2621, 2647,
    2659, 2671, 2683, 2689, 2699, 2711, 2719, 2731, 
};

static unsigned int PRIMES_BOB[TERSIDH_TMAX] = {
    3,    7,    13,   19,   29,   37,   43,   53,   61,   71,   79,   89,
    101,  107,  113,  131,  139,  151,  163,  173,  181,  193,  199,  223,
    229,  239,  251,  263,  271,  281,  293,  311,  317,  337,  349,  359,
    373,  383,  397,  409,  421,  433,  443,  457,  463,  479,  491,  503,
    521,  541,  557,  569,  577,  593,  601,  613,  619,  641,  647,  659,
    673,  683,  701,  719,  733,  743,  757,  769,  787,  809,  821,  827,
    839,  857,  863,  881,  887,  911,  929,  941,  953,  971,  983,  997,
    1013, 1021, 1033, 1049, 1061, 1069, 1091, 1097, 1109, 1123, 1151, 1163,
    1181, 1193, 1213, 1223, 1231, 1249, 1277, 1283, 1291, 1301, 1307, 1321,
    1361, 1373, 1399, 1423, 1429, 1439, 1451, 1459, 1481, 1487, 1493, 1511,
    1531, 1549, 1559, 1571, 1583, 1601, 1609, 1619, 1627, 1657, 1667, 1693,
    1699, 1721, 1733, 1747, 1759, 1783, 1789, 1811, 1831, 1861, 1871, 1877,
    1889, 1907, 1931, 1949, 1973, 1987, 1997, 2003, 2017, 2029, 2053, 2069,
    2083, 2089, 2111, 2129, 2137, 2143, 2161, 2203, 2213, 2237, 2243, 2267,
    2273, 2287, 2297, 2311, 2339, 2347, 2357, 2377, 2383, 2393, 2411, 2423,
    2441, 2459, 2473, 2503, 2531, 2543, 2551, 2579, 2593, 2617, 2633, 2657,
    2663, 2677, 2687, 2693, 2707, 2713, 2729, 2741, 
};

/*
 * @brief Generate random secret for the tersidh and compute the kernel points. 
 * @reads: t, is_bob, secret (if skip_secret = 1), A24_start, C24_start, PQ_self
 * @modifies: secret (if skip_secret = 0), KP_deg, KQ_deg, KP, KQ
*/
void tersidh_generate_kernel_points(struct tersidh_state* tersidh, int skip_secret) {
    // Interpret secret as ternary number of length `t`.
    mpz_t r, n, cP, cQ;
    mpz_init(r); 
    mpz_init(n);
    mpz_init(cP);
    mpz_init(cQ);

    if (!skip_secret) {
        // Calculate upper bound for the secret: ternary string of length `t` 
        mpz_ui_pow_ui(n, 3, tersidh->t);  // n = 3^t
        // Sample random integer from the set: [0, 3^t)
        mpz_urandomm(tersidh->secret, tersidh->randstate, n);
    }

    mpz_set(n, tersidh->secret);
    mpz_set_ui(cP, 1);
    mpz_set_ui(cQ, 1);

    unsigned int *primes = tersidh->is_bob ? PRIMES_BOB : PRIMES_ALICE;

    unsigned int *kp_primes = malloc(sizeof(unsigned int) * tersidh->t);
    unsigned int *kq_primes = malloc(sizeof(unsigned int) * tersidh->t);
    int kp_size = 0, kq_size = 0;

    for (int i = 0; i < tersidh->t; i++) {
        // n = n//3; r = n % 3
        mpz_fdiv_qr_ui(n, r, n, 3);
        // r = {-1, 0, 1}
        mpz_sub_ui(r, r, 1);
        int sgn = mpz_sgn(r);

        // current prime number
        int p = primes[i];

        switch (sgn) {

            // Increase KP order, Decrease KQ order
            case -1: { 
                kp_primes[kp_size++] = p;   // ord(KP) *= p
                mpz_mul_ui(cQ, cQ, p);      // KQ = [p]KQ
            } break;

            // Increase KQ order, Decrease KP order
            case  0: {
                mpz_mul_ui(cP, cP, p);      // KP = [p]KP
                kq_primes[kq_size++] = p;   // ord(KQ) *= p
            } break;

            // Decrease both orders
            case  1: {
                mpz_mul_ui(cP, cP, p);      // KP = [p]KP
                mpz_mul_ui(cQ, cQ, p);      // KQ = [p]KQ
            } break;
            default: assert(0 && "Something went wrong during TerSIDH secret calculation.");
        }
    }

    pprod_set_array(tersidh->KP_deg, kp_primes, kp_size);
    pprod_set_array(tersidh->KQ_deg, kq_primes, kq_size);

    // KP = [cP]P
    xLADDER(tersidh->KP, tersidh->PQ_self.P, cP, tersidh->A24p_start, tersidh->C24_start);
    // KQ = [cQ]Q
    xLADDER(tersidh->KQ, tersidh->PQ_self.Q, cQ, tersidh->A24p_start, tersidh->C24_start);

    free(kp_primes);
    free(kq_primes);
    mpz_clear(r);
    mpz_clear(n);
    mpz_clear(cP);
    mpz_clear(cQ);
}


static inline int _apply_and_test_cofactor(mpz_t result, const mpz_t base,
                                           int f) {
    // test if result = base * f - 1 is a prime number
    mpz_mul_ui(result, base, f);
    mpz_sub_ui(result, result, 1);
    int is_prime = mpz_probab_prime_p(result, 100);
    // Value 2 means "confident prime", value 1 means "probably prime".
    // Return 0 if is prime or probably prime, -1 otherwise
    return (is_prime == 1 || is_prime == 2) ? 0 : -1;
}

// Find number f such that base * f - 1 is a prime number
// Return -1 if max number of tries exceeded
static int find_cofator(mpz_t result, const mpz_t base) {
    for (int f = 1; f < 1000; f++) {
        int errors = _apply_and_test_cofactor(result, base, f);
        if (!errors) {
            return f;
        }
    }

    return -1;
}

int tersidh_calc_pub_params(mpz_t p, pprod_t A, pprod_t B, int t, int f) {
    if (t > TERSIDH_TMAX || t < TERSIDH_TMIN || f < 0) {
        return -1;
    }

    // Generate composite numbers A, B
    pprod_set_array(A, PRIMES_ALICE, t);
    pprod_set_array(B, PRIMES_BOB, t);

    // Find cofactor f: p = fAB - 1
    mpz_t AB;
    mpz_init(AB);

    // temp = A * B
    mpz_mul(AB, A->value, B->value);

    // calculate p = ABf - 1, return 0 if prime, -1 otherwise
    int ret = _apply_and_test_cofactor(p, AB, f);
    mpz_clear(AB);

    return ret;
}

int tersidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, int t) {
    if (t > TERSIDH_TMAX || t < TERSIDH_TMIN) {
        return -1;
    }

    // Generate composite numbers A, B
    pprod_set_array(A, PRIMES_ALICE, t);
    pprod_set_array(B, PRIMES_BOB, t);

    // Find cofactor f: p = fAB - 1
    mpz_t AB;
    mpz_init(AB);

    // temp = A * B
    mpz_mul(AB, A->value, B->value);

    // f = -1 if cannot find such prime number p = fAB - 1
    int f = find_cofator(p, AB);

    mpz_clear(AB);

    return f;
}

void tersidh_state_reset(struct tersidh_state *tersidh) {
    assert(tersidh->status == TERSIDH_STATUS_PREPARED ||
           tersidh->status == TERSIDH_STATUS_INITIALIZED ||
           tersidh->status == TERSIDH_STATUS_EXCHANGED);

    // state_prepare" initializes global fpchar
    if (tersidh->status != TERSIDH_STATUS_INITIALIZED) {
        fpchar_clear_if_set();
    }

    // We dont deallocate the variables - simply change the status so "prepare"
    // can be called
    tersidh->status = TERSIDH_STATUS_INITIALIZED;
}

void tersidh_state_prepare(struct tersidh_state *tersidh,
                         const struct tersidh_data *params, int is_bob) {
    assert(tersidh->status == TERSIDH_STATUS_INITIALIZED);

    // `a = 2` is invalid in montgomery model
    assert(!fp2_equal_uint(params->a, 2) &&
           "Curve coefficient cannot be equal to 2");
    assert(params->t >= TERSIDH_TMIN && params->t <= TERSIDH_TMAX && "Invalid t-parameter size");
    // assert(params->f > 0 && "TERSIDH Params cofactor f must be a positive
    // integer");

    // Middle node - elliptic curve between both isogenies (codomain of KP's isogeny)
    fp2_t A24p_mid, C24_mid;
    fp2_init(&A24p_mid);
    fp2_init(&C24_mid);

    // Image of the KQ point under KP's isogeny
    point_t phi_KQ;
    point_init(&phi_KQ);

    // Convert xP, xQ and xR to torsion basis
    struct tors_basis PQ;
    tors_basis_init(&PQ);

    // Generate public protocol params: p, A, B given t;
    tersidh->is_bob = is_bob;
    tersidh->t = params->t;
    tersidh->f = params->f;

    // Do not seek the cofactor f, only calculate p = AB-f and check for
    // primality
    int ret =
        tersidh_calc_pub_params(tersidh->p, tersidh->A, tersidh->B, tersidh->t, tersidh->f);
    assert(ret == 0 && "TERSIDH cannot calculate public params");

    // Initialize global characteristic if its not set
    fpchar_clear_if_set();
    ret = fpchar_setup(tersidh->p);
    assert(ret == 0 &&
           "TERSIDH cannot work properly if global characteristic is invalid");

    // fill torsion basis P,Q = E[n] data based on given params
    point_set_fp2_x(PQ.P, params->xP);
    point_set_fp2_x(PQ.Q, params->xQ);
    point_set_fp2_x(PQ.PQd, params->xR);
    mpz_add_ui(PQ.n, tersidh->p, 1);

    // Initalize starting Ellitptic Curve: y^2 = x^3 + ax^2 + x
    fp2_set(tersidh->A24p_start, params->a);
    fp2_set_uint(tersidh->C24_start, 1);
    A24p_from_A(tersidh->A24p_start, tersidh->C24_start, tersidh->A24p_start,
                tersidh->C24_start);

    pprod_t *deg_self  = is_bob ? &tersidh->B : &tersidh->A;
    pprod_t *deg_other = is_bob ? &tersidh->A : &tersidh->B;

    // Generate my torsion basis PA, QA = E0[A]
    tors_basis_get_subgroup(&tersidh->PQ_self, (*deg_self)->value, &PQ,
                            tersidh->A24p_start, tersidh->C24_start);

    // Generate other torsion basis: PB, QB = E0[B]
    tors_basis_get_subgroup(&tersidh->PQ_pubkey, (*deg_other)->value, &PQ,
                            tersidh->A24p_start, tersidh->C24_start);


    // Draft random secret; Generate kernel points: KP, KQ
    // Sample new random secret value only if it's equal to 0 => otherwise it was set by the user (unit tests)
    int skip_secret = !!(mpz_sgn(tersidh->secret));
    tersidh_generate_kernel_points(tersidh, skip_secret);
    point_set(phi_KQ, tersidh->KQ);

    // -- Calculate both isogenies from KP and KQ
    // We push all Bob's torsion points + the KQ kernel point to calculate the next isogeny based on it
    point_t push_points[] = {tersidh->PQ_pubkey.P, tersidh->PQ_pubkey.Q, tersidh->PQ_pubkey.PQd, phi_KQ, NULL, NULL};

    // Calculate first isogeny phi_KP
    ISOG_chain(A24p_mid, C24_mid, tersidh->A24p_start, tersidh->C24_start, tersidh->KP, tersidh->KP_deg, push_points);

    // Remove the KQ kernel point from the push points list
    push_points[3] = NULL;

    // Calculate second isogeny phi_KQ
    ISOG_chain(tersidh->A24p_pubkey, tersidh->C24_pubkey, A24p_mid, C24_mid,phi_KQ, tersidh->KQ_deg, push_points);

    // Normalize for further access
    point_normalize_coords(tersidh->PQ_pubkey.P);
    point_normalize_coords(tersidh->PQ_pubkey.Q);
    point_normalize_coords(tersidh->PQ_pubkey.PQd);

    tors_basis_clear(&PQ);

    point_clear(&phi_KQ);

    fp2_clear(&A24p_mid);
    fp2_clear(&C24_mid);
    tersidh->status = TERSIDH_STATUS_PREPARED;
}

void tersidh_key_exchange(struct tersidh_state *tersidh,
                        const struct tersidh_data *pk_other) {
    assert(tersidh->status == TERSIDH_STATUS_PREPARED);

    fp2_t A24p_final, C24_final, A24p_mid, C24_mid;
    fp2_init(&A24p_final);
    fp2_init(&C24_final);
    fp2_init(&A24p_mid);
    fp2_init(&C24_mid);

    // TODO: Copy the point into my own torsion, to not discard the original one
    // TODO: add verification of the pairing
    assert(tersidh->t == pk_other->t);

    // Update my torsion basis (it will be destroyed during the key_exchange
    // process)
    point_set_fp2_x(tersidh->PQ_self.P, pk_other->xP);
    point_set_fp2_x(tersidh->PQ_self.Q, pk_other->xQ);
    point_set_fp2_x(tersidh->PQ_self.PQd, pk_other->xR);

    pprod_t *deg_self = tersidh->is_bob ? &tersidh->B : &tersidh->A;

    // Order should not change from previous iteration
    assert(0 == mpz_cmp(tersidh->PQ_self.n, (*deg_self)->value));

    // Reconstruct the Elliptic Curve given by pk_other, into new E0 for tersidh state
    // We need to overwrite E0 that is used in tersidh generate_kernel_points
    fp2_set(tersidh->A24p_start, pk_other->a);
    fp2_set_uint(tersidh->C24_start, 1);
    A24p_from_A(tersidh->A24p_start, tersidh->C24_start, tersidh->A24p_start, tersidh->C24_start);

    // Use already calculated secret, new basis <PA,QA> = EB[A] and E0 := EB to generate KP and KQ
    tersidh_generate_kernel_points(tersidh, 1); 

    point_t phi_KQ;
    point_init(&phi_KQ);
    point_set(phi_KQ, tersidh->KQ);

    // -- Calculate both isogenies from KP and KQ
    point_t push_points[] = { phi_KQ, NULL, NULL};

    // First isogeny; phi_KP
    ISOG_chain(A24p_mid, C24_mid, tersidh->A24p_start, tersidh->C24_start, tersidh->KP, tersidh->KP_deg, push_points);

    // Remove the KQ kernel from the list of points
    push_points[0] = NULL;

    // Second isogeny; phi_KQ
    ISOG_chain(A24p_final, C24_final, A24p_mid, C24_mid, phi_KQ, tersidh->KQ_deg, push_points);

    // Calculate j_invariant of the curve
    A_from_A24p(A24p_final, C24_final, A24p_final, C24_final);
    j_invariant(tersidh->j_inv, A24p_final, C24_final);

    point_clear(&phi_KQ);
    fp2_clear(&A24p_final);
    fp2_clear(&C24_final);
    fp2_clear(&A24p_mid);
    fp2_clear(&C24_mid);

    tersidh->status = TERSIDH_STATUS_EXCHANGED;
}

void tersidh_get_pubkey(const struct tersidh_state *tersidh,
                      struct tersidh_data *pk_self) {
    assert(tersidh->status == TERSIDH_STATUS_PREPARED);

    // Point have to be normalized, otherwise we will get false results
    assert(point_is_normalized(tersidh->PQ_pubkey.P));
    assert(point_is_normalized(tersidh->PQ_pubkey.Q));
    assert(point_is_normalized(tersidh->PQ_pubkey.PQd));

    pk_self->t = tersidh->t;
    pk_self->f = tersidh->f;

    fp2_set(pk_self->xP, tersidh->PQ_pubkey.P->X);
    fp2_set(pk_self->xQ, tersidh->PQ_pubkey.Q->X);
    fp2_set(pk_self->xR, tersidh->PQ_pubkey.PQd->X);

    // Calculate the a coefficient
    // TODO: use A:C inside TERSIDH state instead of converting
    fp2_t A, C;
    fp2_init(&A);
    fp2_init(&C);
    A_from_A24p(A, C, tersidh->A24p_pubkey, tersidh->C24_pubkey);
    assert(!fp2_is_zero(C));

    fp2_div_unsafe(pk_self->a, A, C);
    fp2_clear(&A);
    fp2_clear(&C);
}

void tersidh_data_init(struct tersidh_data *md) {
    fp2_init(&md->a);
    fp2_init(&md->xP);
    fp2_init(&md->xQ);
    fp2_init(&md->xR);
}

void tersidh_data_clear(struct tersidh_data *md) {
    fp2_clear(&md->a);
    fp2_clear(&md->xP);
    fp2_clear(&md->xQ);
    fp2_clear(&md->xR);
}

void tersidh_state_init(struct tersidh_state *tersidh) {
    gmp_randinit_mt(tersidh->randstate);

    mpz_init(tersidh->p);
    pprod_init(&tersidh->A);
    pprod_init(&tersidh->B);

    tors_basis_init(&tersidh->PQ_self);
    tors_basis_init(&tersidh->PQ_pubkey);

    point_init(&tersidh->KP);
    point_init(&tersidh->KQ);

    pprod_init(&tersidh->KP_deg);
    pprod_init(&tersidh->KQ_deg);

    fp2_init(&tersidh->A24p_start);
    fp2_init(&tersidh->C24_start);

    fp2_init(&tersidh->A24p_pubkey);
    fp2_init(&tersidh->C24_pubkey);

    fp2_init(&tersidh->j_inv);
    mpz_init(tersidh->secret);
    assert(mpz_sgn(tersidh->secret) == 0 && "TerSIDH state secret should be set to 0 after init.");

    tersidh->status = TERSIDH_STATUS_INITIALIZED;
}

void tersidh_state_clear(struct tersidh_state *tersidh) {
    gmp_randclear(tersidh->randstate);

    mpz_clear(tersidh->p);
    pprod_clear(&tersidh->A);
    pprod_clear(&tersidh->B);

    tors_basis_clear(&tersidh->PQ_self);
    tors_basis_clear(&tersidh->PQ_pubkey);

    point_clear(&tersidh->KP);
    point_clear(&tersidh->KQ);

    pprod_clear(&tersidh->KP_deg);
    pprod_clear(&tersidh->KQ_deg);

    fp2_clear(&tersidh->A24p_start);
    fp2_clear(&tersidh->C24_start);

    fp2_clear(&tersidh->A24p_pubkey);
    fp2_clear(&tersidh->C24_pubkey);

    fp2_clear(&tersidh->j_inv);
    mpz_clear(tersidh->secret);

    tersidh->status = TERSIDH_STATUS_UNINITIALIZED;
}
