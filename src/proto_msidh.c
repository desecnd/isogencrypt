#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>

#include "fp2.h"
#include "pprod.h"
#include "ec_mont.h"
#include "proto_msidh.h"


/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q) of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, pprod_t n, const struct tors_basis *PQ, const fp2_t A24p, const fp2_t C24) {
    // Temporarily set the order as [N/n]
    mpz_divexact(PQsub->n->value, PQ->n->value, n->value); 
    // Multiply all points in the subbasis by [N/n]
    xLADDER(PQsub->P, PQ->P, PQsub->n->value, A24p, C24);
    xLADDER(PQsub->Q, PQ->Q, PQsub->n->value, A24p, C24);
    xLADDER(PQsub->PQd, PQ->PQd, PQsub->n->value, A24p, C24);
    // Set the proper order of the subgroup torsion basis
    pprod_set(PQsub->n, n);
}

// Sample an element `x` from ``Z/mZ`` where ``x^2 = 1 (mod m)``.
// Return != 0 if something goes wrong
int sample_quadratic_root_of_unity(mpz_t result, const pprod_t modulus) {
    // TODO: CRT: check if it will overflow?
    // Use fp_t instead of mpz_t if so

    int ret = 0;
    mpz_set_ui(result, 0);

    mpz_t m, inv;
    mpz_init(m);
    mpz_init(inv);

    for (unsigned i = 0; i < modulus->n_primes; i++) {
        unsigned int p = modulus->primes[i];

        // Toss a coin and choose 1 or -1 (== p - 1)
        // TODO: change rand to something more serious in the future
        unsigned int root = 1;
        if (rand() % 2 == 0) {
            root = p - 1;
        }

        // Solve CRT Congruence: root^2 == 1 (mod p)

        // m = M / p 
        mpz_div_ui(m, modulus->value, p);

        // set temporary: inv = p
        // inv = (M/p)^-1 (mod p)
        // if returns 0: inverse does not exist
        mpz_set_ui(inv, p);
        int inv_ret = mpz_invert(inv, m, inv);
        if (inv_ret == 0) {
            ret = -1;
            break;
        }

        // m = root * (M / p)
        mpz_mul_ui(m, m, root);

        // result += root * (M/p) * (M/p)^-1
        mpz_addmul(result, m, inv);
        mpz_mod(result, result, modulus->value);
    }
    mpz_clear(m);
    mpz_clear(inv);

    return ret;
}

// TODO: From MSIDH paper: 128bit prime: p = 2^2 * l1 ... l571 * 10 - 1 
// So at least 571//2 primes per side is required
#define MSIDH_TMAX 600
#define MSIDH_TMAX_HALF MSIDH_TMAX/2

// TODO:
// Catch! The first number is not a prime number => 2^2, this is a special case
unsigned int PRIMES_ALICE[MSIDH_TMAX_HALF] = {
    4, 5, 11, 17, 23, 31, 41, 47, 59, 67, 73, 83, 97, 103, 109, 127, 137, 149, 157, 167, 179, 191, 197, 211, 227, 
    233, 241, 257, 269, 277, 283, 307, 313, 331, 347, 353, 367, 379, 389, 401, 419, 431, 439, 449, 461, 467, 487, 
    499, 509, 523, 547, 563, 571, 587, 599, 607, 617, 631, 643, 653, 661, 677, 691, 709, 727, 739, 751, 761, 773, 
    797, 811, 823, 829, 853, 859, 877, 883, 907, 919, 937, 947, 967, 977, 991, 1009, 1019, 1031, 1039, 1051, 1063, 
    1087, 1093, 1103, 1117, 1129, 1153, 1171, 1187, 1201, 1217, 1229, 1237, 1259, 1279, 1289, 1297, 1303, 1319, 1327, 
    1367, 1381, 1409, 1427, 1433, 1447, 1453, 1471, 1483, 1489, 1499, 1523, 1543, 1553, 1567, 1579, 1597, 1607, 1613,
    1621, 1637, 1663, 1669, 1697, 1709, 1723, 1741, 1753, 1777, 1787, 1801, 1823, 1847, 1867, 1873, 1879, 1901, 1913, 
    1933, 1951, 1979, 1993, 1999, 2011, 2027, 2039, 2063, 2081, 2087, 2099, 2113, 2131, 2141, 2153, 2179, 2207, 2221, 
    2239, 2251, 2269, 2281, 2293, 2309, 2333, 2341, 2351, 2371, 2381, 2389, 2399, 2417, 2437, 2447, 2467, 2477, 2521, 
    2539, 2549, 2557, 2591, 2609, 2621, 2647, 2659, 2671, 2683, 2689, 2699, 2711, 2719, 2731, 2749, 2767, 2789, 2797, 
    2803, 2833, 2843, 2857, 2879, 2897, 2909, 2927, 2953, 2963, 2971, 3001, 3019, 3037, 3049, 3067, 3083, 3109, 3121, 
    3163, 3169, 3187, 3203, 3217, 3229, 3253, 3259, 3299, 3307, 3319, 3329, 3343, 3359, 3371, 3389, 3407, 3433, 3457, 
    3463, 3469, 3499, 3517, 3529, 3539, 3547, 3559, 3581, 3593, 3613, 3623, 3637, 3659
};

unsigned int PRIMES_BOB[MSIDH_TMAX_HALF] = {
    3, 7, 13, 19, 29, 37, 43, 53, 61, 71, 79, 89, 101, 107, 113, 131, 139, 151, 163, 173, 181, 193, 199, 223, 229, 
    239, 251, 263, 271, 281, 293, 311, 317, 337, 349, 359, 373, 383, 397, 409, 421, 433, 443, 457, 463, 479, 491, 
    503, 521, 541, 557, 569, 577, 593, 601, 613, 619, 641, 647, 659, 673, 683, 701, 719, 733, 743, 757, 769, 787, 
    809, 821, 827, 839, 857, 863, 881, 887, 911, 929, 941, 953, 971, 983, 997, 1013, 1021, 1033, 1049, 1061, 1069, 
    1091, 1097, 1109, 1123, 1151, 1163, 1181, 1193, 1213, 1223, 1231, 1249, 1277, 1283, 1291, 1301, 1307, 1321, 1361,
    1373, 1399, 1423, 1429, 1439, 1451, 1459, 1481, 1487, 1493, 1511, 1531, 1549, 1559, 1571, 1583, 1601, 1609, 1619,
    1627, 1657, 1667, 1693, 1699, 1721, 1733, 1747, 1759, 1783, 1789, 1811, 1831, 1861, 1871, 1877, 1889, 1907, 1931, 
    1949, 1973, 1987, 1997, 2003, 2017, 2029, 2053, 2069, 2083, 2089, 2111, 2129, 2137, 2143, 2161, 2203, 2213, 2237, 
    2243, 2267, 2273, 2287, 2297, 2311, 2339, 2347, 2357, 2377, 2383, 2393, 2411, 2423, 2441, 2459, 2473, 2503, 2531, 
    2543, 2551, 2579, 2593, 2617, 2633, 2657, 2663, 2677, 2687, 2693, 2707, 2713, 2729, 2741, 2753, 2777, 2791, 2801, 
    2819, 2837, 2851, 2861, 2887, 2903, 2917, 2939, 2957, 2969, 2999, 3011, 3023, 3041, 3061, 3079, 3089, 3119, 3137, 
    3167, 3181, 3191, 3209, 3221, 3251, 3257, 3271, 3301, 3313, 3323, 3331, 3347, 3361, 3373, 3391, 3413, 3449, 3461, 
    3467, 3491, 3511, 3527, 3533, 3541, 3557, 3571, 3583, 3607, 3617, 3631, 3643, 3671
};


// Find number f such that base * f - 1 is a prime number
// Return -1 if max number of tries exceeded
int find_cofator(mpz_t result, const mpz_t base) {

    for (int f = 1; f < 1000; f++) {
        mpz_mul_ui(result, base, f);
        mpz_sub_ui(result, result, 1);
        int is_prime = mpz_probab_prime_p(result, 100);
        // 2 means "confident", 1 means "probably prime"
        if (is_prime == 1 || is_prime == 2) {
            return f;
        }
    }

    return -1;
}

int msidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, unsigned int t) {
    if (t > MSIDH_TMAX) {
        return -1;
    }

    // Generate composite numbers A, B
    pprod_set_array(A, PRIMES_ALICE, (t + 1)/2);
    pprod_set_array(B, PRIMES_BOB, t/2);

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

void msidh_state_reset(struct msidh_state *msidh) {
    assert(
        msidh->status == MSIDH_STATUS_PREPARED || 
        msidh->status == MSIDH_STATUS_INITIALIZED || 
        msidh->status == MSIDH_STATUS_EXCHANGED
    );

    // We dont deallocate the variables - simply change the status so "prepare" can be called
    msidh->status = MSIDH_STATUS_INITIALIZED;
}

void msidh_key_exchange(struct msidh_state *msidh, const fp2_t A24p_other, const fp2_t C24_other, const struct tors_basis* PQ_masked) {
    assert(msidh->status == MSIDH_STATUS_PREPARED);

    fp2_t A24p_final, C24_final;
    fp2_init(&A24p_final);
    fp2_init(&C24_final);

    // TODO: Copy the point into my own torsion, to not discard the original one
    // TODO: add verification of the pairing

    // Update my torsion basis (it will be destroyed during the key_exchange process)
    assert(!mpz_cmp(msidh->PQ_self.n->value, PQ_masked->n->value));
    point_set(msidh->PQ_self.P, PQ_masked->P);
    point_set(msidh->PQ_self.Q, PQ_masked->Q);
    point_set(msidh->PQ_self.PQd, PQ_masked->PQd);

    _msidh_key_exchange_alice(msidh->j_inv, A24p_final, C24_final, A24p_other, C24_other, &msidh->PQ_self, msidh->secret);

    fp2_clear(&A24p_final);
    fp2_clear(&C24_final);
    
    msidh->status = MSIDH_STATUS_EXCHANGED;
}

void msidh_state_prepare(struct msidh_state *msidh, unsigned int t, const fp2_t a, const struct tors_basis *PQ, int is_bob) {
    // TODO: what about the seed?
    assert(msidh->status == MSIDH_STATUS_INITIALIZED);
    // `a = 2` is invalid in montgomery model
    assert(!fp2_equal_uint(a, 2));

    // Generate public protocol params: p, A, B given t;
    int ret = msidh_gen_pub_params(msidh->p, msidh->A, msidh->B, t);
    assert(ret > 0);

    // TODO: Make sure that the global arithmetic is set to FP
    // Initialize global characteristic if its not set
    ret = global_fpchar_clear();
    ret = global_fpchar_setup(msidh->p);

    mpz_t mask;
    mpz_init(mask);

    // Initalize starting Ellitptic Curve: y^2 = x^3 + ax^2 + x
    fp2_set(msidh->A24p_start, a);
    fp2_set_uint(msidh->C24_start, 1); 
    A24p_from_A(msidh->A24p_start, msidh->C24_start, msidh->A24p_start, msidh->C24_start);

    if (is_bob) {
        tors_basis_get_subgroup(&msidh->PQ_self, msidh->B, PQ, msidh->A24p_start, msidh->C24_start); 
        tors_basis_get_subgroup(&msidh->PQ_pubkey, msidh->A, PQ, msidh->A24p_start, msidh->C24_start); 
        // Generate random secret s in range [0, msidh->B)
        mpz_urandomm(msidh->secret, msidh->randstate, msidh->B->value);
        sample_quadratic_root_of_unity(mask, msidh->A); 
    } else {
        tors_basis_get_subgroup(&msidh->PQ_self, msidh->A, PQ, msidh->A24p_start, msidh->C24_start); 
        tors_basis_get_subgroup(&msidh->PQ_pubkey, msidh->B, PQ, msidh->A24p_start, msidh->C24_start); 
        // Generate random secret s in range [0, msidh->B)
        mpz_urandomm(msidh->secret, msidh->randstate, msidh->A->value);
        sample_quadratic_root_of_unity(mask, msidh->B); 
    }

    _msidh_gen_pubkey_alice(msidh->A24p_pubkey, msidh->C24_pubkey, &msidh->PQ_self, &msidh->PQ_pubkey, msidh->A24p_start, msidh->C24_start, msidh->secret, mask);

    mpz_clear(mask);

    msidh->status = MSIDH_STATUS_PREPARED;
}

/* 
 * @brief Generate MSIDH public key from Alice perspective
*/
void _msidh_gen_pubkey_alice(fp2_t A24p_alice, fp2_t C24_alice, struct tors_basis* PQ_alice, struct tors_basis* PQ_bob, const fp2_t A24p_base, const fp2_t C24_base, const mpz_t secret, const mpz_t mask) {
    // P, Q is a torsion basis for deg

    // assert(global_fpchar_setup(p));

    // 1. Calculate the kernel of the Alice isogeny
    // PA = PA + [s]QA
    xLADDER3PT(PQ_alice->P, PQ_alice->Q, PQ_alice->PQd, secret, A24p_base, C24_base);

    point_t push_points[] = {PQ_bob->P, PQ_bob->Q, PQ_bob->PQd, NULL, NULL};

    ISOG_chain(A24p_alice, C24_alice, A24p_base, C24_base, PQ_alice->P, PQ_alice->n, push_points);

    // Mask the Bob torsion basis using quadratic root - alpha.
    // mpz_t alpha;
    // mpz_init(alpha);
    // int ret = sample_quadratic_root_of_unity(alpha, deg_B);
    // assert(!ret);

    // 3. Apply masking
    // We multiply all the points (PB, QB, PQBd) by `alpha` 
    // We use QA as temporary register for holding the point result
    xLADDER(PQ_alice->Q, PQ_bob->P, mask, A24p_alice, C24_alice);
    point_set(PQ_bob->P, PQ_alice->Q);
    xLADDER(PQ_alice->Q, PQ_bob->Q, mask, A24p_alice, C24_alice);
    point_set(PQ_bob->Q, PQ_alice->Q);
    xLADDER(PQ_alice->Q, PQ_bob->PQd, mask, A24p_alice, C24_alice);
    point_set(PQ_bob->PQd, PQ_alice->Q);
}


// TODO: Note that BPQA get destroyed
void _msidh_key_exchange_alice(fp2_t j_inv, fp2_t A24p_final, fp2_t C24_final, const fp2_t A24p_bob, const fp2_t C24_bob, struct tors_basis * BPQA, const mpz_t A_sec) {

    // 1. Calculate the kernel of the Alice isogeny
    // PA = PA + [s]QA
    xLADDER3PT(BPQA->P, BPQA->Q, BPQA->PQd, A_sec, A24p_bob, C24_bob);

    point_t push_points[] = {NULL, NULL};

    ISOG_chain(A24p_final, C24_final, A24p_bob, C24_bob, BPQA->P, BPQA->n, push_points);

    fp2_t A, C;
    fp2_init(&A); fp2_init(&C);

    A_from_A24p(A, C, A24p_final, C24_final);
    j_invariant(j_inv, A, C);

    fp2_clear(&A); fp2_clear(&C);
}

void msidh_state_init(struct msidh_state* msidh) {
    gmp_randinit_mt(msidh->randstate);

    mpz_init(msidh->p);
    pprod_init(&msidh->A);
    pprod_init(&msidh->B);

    point_init(&msidh->PQ_self.P);
    point_init(&msidh->PQ_self.Q);
    point_init(&msidh->PQ_self.PQd);
    pprod_init(&msidh->PQ_self.n);

    point_init(&msidh->PQ_pubkey.P);
    point_init(&msidh->PQ_pubkey.Q);
    point_init(&msidh->PQ_pubkey.PQd);
    pprod_init(&msidh->PQ_pubkey.n);

    fp2_init(&msidh->A24p_start);
    fp2_init(&msidh->C24_start);

    fp2_init(&msidh->A24p_pubkey);
    fp2_init(&msidh->C24_pubkey);

    fp2_init(&msidh->j_inv);
    mpz_init(msidh->secret);

    msidh->status = MSIDH_STATUS_INITIALIZED;
}

void msidh_state_clear(struct msidh_state* msidh) {
    gmp_randclear(msidh->randstate);

    mpz_clear(msidh->p);
    pprod_clear(&msidh->A);
    pprod_clear(&msidh->B);

    point_clear(&msidh->PQ_self.P);
    point_clear(&msidh->PQ_self.Q);
    point_clear(&msidh->PQ_self.PQd);
    pprod_clear(&msidh->PQ_self.n);

    point_clear(&msidh->PQ_pubkey.P);
    point_clear(&msidh->PQ_pubkey.Q);
    point_clear(&msidh->PQ_pubkey.PQd);
    pprod_clear(&msidh->PQ_pubkey.n);

    fp2_clear(&msidh->A24p_start);
    fp2_clear(&msidh->C24_start);

    fp2_clear(&msidh->A24p_pubkey);
    fp2_clear(&msidh->C24_pubkey);

    fp2_clear(&msidh->j_inv);
    mpz_clear(msidh->secret);

    msidh->status = MSIDH_STATUS_UNINITIALIZED;
}