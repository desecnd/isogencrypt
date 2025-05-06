#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "ec_mont.h"
#include "fp2.h"
#include "testing.h"
#include "proto_msidh.h"
#include "pprod.h"

fp2_t A24p, C24;
point_t P, Q, PQd, K;
point_t PA, QA, PQAd;
point_t PB, QB, PQBd;
pprod_t A_deg, B_deg;
mpz_t p, m;

void init_test_variables() {
    fp2_init(&A24p);
    fp2_init(&C24);
    point_init(&P);
    point_init(&Q);
    point_init(&PQd);
    point_init(&K);
    point_init(&PA);
    point_init(&QA);
    point_init(&PQAd);
    point_init(&PB);
    point_init(&QB);
    point_init(&PQBd);
    pprod_init(&A_deg);
    pprod_init(&B_deg);
    mpz_init(p);
    mpz_init(m);
}

void clear_test_variables() {
    fp2_clear(&A24p);
    fp2_clear(&C24);
    point_clear(&P);
    point_clear(&Q);
    point_clear(&PQd);
    point_clear(&K);
    point_clear(&PA);
    point_clear(&QA);
    point_clear(&PQAd);
    point_clear(&PB);
    point_clear(&QB);
    point_clear(&PQBd);
    pprod_clear(&A_deg);
    pprod_clear(&B_deg);
    mpz_clear(p);
    mpz_clear(m);
}

void test_pprod_init() {
    // product: 223'092'870 < 2^32
    unsigned int primes[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    pprod_t M;
    pprod_init(&M);

    pprod_set_array(M, primes, 9);
    CHECK_MSG(!mpz_cmp_ui(M->value, 223092870), "Incorrect product of primes value");

    pprod_clear(&M);
}

void test_random_unit_sampling_large() {
    srand(0xdeafbeef);

    // 100 prime numbers
    unsigned int primes_large[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 
        83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 
        199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 
        347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 
        479, 487, 491, 499, 503, 509, 521, 523, 541
    };

    pprod_t M;
    pprod_init(&M);
    pprod_set_array(M, primes_large, sizeof(primes_large)/sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret, "sample_quadratic_root_of_unity returned non-zero value");
        if (ret) break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);
        
        CHECK_MSG(!mpz_cmp_ui(result, 1), "Sampling result is not a square root of unity");
    }

    mpz_clear(result);
    pprod_clear(&M);
}

void test_random_unit_sampling_small() {
    // seed random number generators
    // TODO: make sure still using rand in implementation
    srand(0xdeafbeef);

    unsigned int primes_small[5] = {2, 3, 5, 7, 11};
    pprod_t M;
    pprod_init(&M);
    pprod_set_array(M, primes_small, sizeof(primes_small)/sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret, "sample_quadratic_root_of_unity returned non-zero value");
        if (ret) break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);
        
        CHECK_MSG(!mpz_cmp_ui(result, 1), "Sampling result is not a square root of unity");
    }

    mpz_clear(result);
    pprod_clear(&M);
}

void test_msidh_gen_pub_params() {
    int t, f;

    t = 4;
    printf("t: %d\n", t);

    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    CHECK_MSG(f>0, "msidh_gen_pub_params returned < 0");

    gmp_printf("p: %Zd\n", p);
    CHECK(!mpz_cmp_ui(p, 419));

    gmp_printf("A: %Zd\n", A_deg->value);
    CHECK(!mpz_cmp_ui(A_deg->value, 20));

    gmp_printf("B: %Zd\n", B_deg->value);
    CHECK(!mpz_cmp_ui(B_deg->value, 21));

    printf("f: %d\n", f);
    CHECK(f == 1);

    // char *p_str = mpz_get_str(NULL, 10, p);
    // CHECK_MSG(!strcmp(p_str, "419"), "Incorrect str comparison");

    t = 9;
    printf("t: %d\n", t);

    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    CHECK_MSG(f>0, "msidh_gen_pub_params returned < 0");

    gmp_printf("p: %Zd\n", p);
    CHECK(!mpz_cmp_ui(p, 892371479));

    gmp_printf("A: %Zd\n", A_deg->value);
    CHECK(!mpz_cmp_ui(A_deg->value, 86020));

    gmp_printf("B: %Zd\n", B_deg->value);
    CHECK(!mpz_cmp_ui(B_deg->value, 5187));

    printf("f: %d\n", f);
    CHECK(f == 2);
}


void setup_params_t4() {
    // p + 1 = 420 = 4 * 3 * 5 * 7
    int f = msidh_gen_pub_params(p, A_deg, B_deg, 4);
    CHECK(f==1);

    // Clear just to be sure
    global_fpchar_clear();
    global_fpchar_setup(p);

    // Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of size 419^2
    fp2_t A, C;
    fp2_init(&A); fp2_init(&C);

    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);
    
    // Convert to (A+2 : 4) form used in xDBL
    A24p_from_A(A24p, C24, A, C);

    fp2_clear(&A); fp2_clear(&C);
}

void test_msidh_secret_zero() {
    point_set_str_x(P, "295*i + 398");
    fp2_print(P->X, "xP");

    point_set_str_x(Q, "314*i + 149");
    fp2_print(Q->X, "xQ");

    point_set_str_x(PQd, "29*i + 395");
    fp2_print(PQd->X, "xPQd");

    mpz_set_ui(m, 0);
    xLADDER3PT(P,  Q, PQd, m, A24p, C24);

    // P should not change
    // xP: 295*i + 398
    point_printx(P, "xK");
    CHECK(!mpz_cmp_ui(P->X->a, 398) && !mpz_cmp_ui(P->X->b, 295));

    xLADDER3PT_int(P, Q, PQd, 0, A24p, C24);

    // P should not change
    // xP: 295*i + 398
    point_printx(P, "xK");
    CHECK(!mpz_cmp_ui(P->X->a, 398) && !mpz_cmp_ui(P->X->b, 295));
}


void test_msidh_internals() {
    point_set_str_x(P, "295*i + 398");
    point_printx(P, "xP");

    point_set_str_x(Q, "314*i + 149");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "29*i + 395");
    point_printx(PQd, "xPQd");

    struct tors_basis PQ = { .P=P, .Q=Q, .PQd=PQd  };
    pprod_init(&PQ.n);

    // n = p + 1
    // TODO: WARNING! not a valid pprod type! We just use it
    mpz_set(PQ.n->value, p);
    mpz_add_ui(PQ.n->value, PQ.n->value, 1);

    // Constuct Alice Basis (PA, QA) = [n//A](P, Q).
    struct tors_basis PQA = { .P=PA, .Q=QA, .PQd=PQAd };
    pprod_init(&PQA.n);
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24); 


    point_printx(PA, "xPA");
    point_printx(QA, "xQA");
    point_printx(PQAd, "xPQAd");

    // Constuct Bob Basis (PB, QB) = [n//B](P, Q)
    struct tors_basis PQB = { .P=PB, .Q=QB, .PQd=PQBd };
    pprod_init(&PQB.n);
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24); 

    point_printx(PB, "xPB");
    point_printx(QB, "xQB");
    point_printx(PQBd, "xPQBd");


    mpz_t a_sec, a_mask, b_sec, b_mask;
    mpz_init(a_sec);
    mpz_init(b_sec);
    mpz_init(a_mask);
    mpz_init(b_mask);
    
    fp2_t A24p_alice, C24_alice, aE_alice;
    fp2_init(&A24p_alice);
    fp2_init(&C24_alice);
    fp2_init(&aE_alice);

    fp2_t A24p_final, C24_final, aE_final;
    fp2_init(&A24p_final);
    fp2_init(&C24_final);
    fp2_init(&aE_final);

    fp2_t j_inv;
    fp2_init(&j_inv);

    // 1. --- Testcase sec = 2, mask = 1
    mpz_set_ui(a_sec, 2);
    mpz_set_ui(a_mask, 1);
    mpz_set_ui(b_sec, 5);
    mpz_set_ui(b_mask, 9);
    gmp_printf("a_sec: %Zd\n", a_sec);
    gmp_printf("a_mask: %Zd\n", a_mask);
    gmp_printf("b_sec: %Zd\n", b_sec);
    gmp_printf("b_mask: %Zd\n", b_mask);

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24,  a_sec, a_mask); 

    // aφ(E)(24p): 271*i + 111
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice); 
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 111) && !mpz_cmp_ui(aE_alice->b, 271));

    // xφ(PB): 363*i + 349
    point_printx(PB, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PQB.P->X->a, 349) && !mpz_cmp_ui(PB->X->b, 363));

    // xφ(QB): 303*i + 391
    point_printx(QB, "xφ(QB)");
    CHECK(!mpz_cmp_ui(PQB.Q->X->a, 391) && !mpz_cmp_ui(QB->X->b, 303));

    // xφ(PB-QB): 239*i + 90
    point_printx(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 90) && !mpz_cmp_ui(PQBd->X->b, 239));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice, C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p): 356*i + 219
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(!mpz_cmp_ui(aE_final->a, 219) && !mpz_cmp_ui(aE_final->b, 356));

    // jτ(φ(E)): 398*i + 166
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(!mpz_cmp_ui(j_inv->a, 166) && !mpz_cmp_ui(j_inv->b, 398));

    // 2. --- testcase sec = 13, mask = 8
    mpz_set_ui(a_sec, 13);
    mpz_set_ui(a_mask, 8);
    mpz_set_ui(b_sec, 7);
    mpz_set_ui(b_mask, 1);
    gmp_printf("a_sec: %Zd\n", a_sec);
    gmp_printf("a_mask: %Zd\n", a_mask);
    gmp_printf("b_sec: %Zd\n", b_sec);
    gmp_printf("b_mask: %Zd\n", b_mask);

    // Recalculate the PQA and PQB bases
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24); 
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24); 

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24,  a_sec, a_mask); 

    // aφ(E)(24p): 408*i + 332
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice); 
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 332) && !mpz_cmp_ui(aE_alice->b, 408));

    // xφ(PB): 192*i + 82
    point_printx(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 82) && !mpz_cmp_ui(PB->X->b, 192));

    // xφ(QB): 299*i + 211
    point_printx(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 211) && !mpz_cmp_ui(QB->X->b, 299));

    // xφ(PQBd): 127*i + 181
    point_printx(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 181) && !mpz_cmp_ui(PQBd->X->b, 127));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice, C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p): 204*i + 395
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(!mpz_cmp_ui(aE_final->a, 395) && !mpz_cmp_ui(aE_final->b, 204));

    // jτ(φ(E)): 407
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(!mpz_cmp_ui(j_inv->a, 407) && !mpz_cmp_ui(j_inv->b, 0));

    // 3. --- Testcase sec = 3, mask = 13
    mpz_set_ui(a_sec, 3);
    mpz_set_ui(a_mask, 13);
    mpz_set_ui(b_sec, 17);
    mpz_set_ui(b_mask, 19);
    gmp_printf("a_sec: %Zd\n", a_sec);
    gmp_printf("a_mask: %Zd\n", a_mask);
    gmp_printf("b_sec: %Zd\n", b_sec);
    gmp_printf("b_mask: %Zd\n", b_mask);

    // Recalculate the PQA and PQB bases
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24); 
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24); 

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24,  a_sec, a_mask); 

    // aφ(E)(24p): 109*i + 386
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice); 
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 386) && !mpz_cmp_ui(aE_alice->b, 109));

    // xφ(PB): 298*i + 413
    point_printx(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 413) && !mpz_cmp_ui(PB->X->b, 298));

    // xφ(QB): 323*i + 222
    point_printx(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 222) && !mpz_cmp_ui(QB->X->b, 323));

    // xφ(PQBd): 323*i + 66
    point_printx(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 66) && !mpz_cmp_ui(PQBd->X->b, 323));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice, C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p): 244*i + 279
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(!mpz_cmp_ui(aE_final->a, 279) && !mpz_cmp_ui(aE_final->b, 244));

    // jτ(φ(E)): 175*i + 351
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(!mpz_cmp_ui(j_inv->a, 351) && !mpz_cmp_ui(j_inv->b, 175));

    // 4. --- Testcase sec = 15, mask = 20
    mpz_set_ui(a_sec, 15);
    mpz_set_ui(a_mask, 20);
    mpz_set_ui(b_sec, 1);
    mpz_set_ui(b_mask, 11);
    gmp_printf("a_sec: %Zd\n", a_sec);
    gmp_printf("a_mask: %Zd\n", a_mask);
    gmp_printf("b_sec: %Zd\n", b_sec);
    gmp_printf("b_mask: %Zd\n", b_mask);

    // Recalculate the PQA and PQB bases
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24); 
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24); 

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24,  a_sec, a_mask); 

    // aφ(E)(24p): 16*i + 353
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice); 
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 353) && !mpz_cmp_ui(aE_alice->b, 16));

    // xφ(PB): 340
    point_printx(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 340) && !mpz_cmp_ui(PB->X->b, 0));

    // xφ(QB): 97*i + 356
    point_printx(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 356) && !mpz_cmp_ui(QB->X->b, 97));

    // xφ(PQBd): 330*i + 244
    point_printx(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 244) && !mpz_cmp_ui(PQBd->X->b, 330));


    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice, C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p): 302
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(!mpz_cmp_ui(aE_final->a, 302) && !mpz_cmp_ui(aE_final->b, 0));

    // jτ(φ(E)): 98
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(!mpz_cmp_ui(j_inv->a, 98) && !mpz_cmp_ui(j_inv->b, 0));

    fp2_clear(&A24p_alice);
    fp2_clear(&C24_alice);
    fp2_clear(&aE_alice);

    fp2_clear(&A24p_final);
    fp2_clear(&C24_final);
    fp2_clear(&aE_final);

    fp2_clear(&j_inv);

    mpz_clear(a_sec);
    mpz_clear(b_sec);
    mpz_clear(a_mask);
    mpz_clear(b_mask);

    pprod_clear(&PQ.n);
    pprod_clear(&PQA.n);
    pprod_clear(&PQB.n);
}


void test_msidh_non_deterministic() {

    struct msidh_state m1, m2, m3;
    msidh_state_init(&m1);
    msidh_state_init(&m2);
    msidh_state_init(&m3);

    struct tors_basis PQ;
    PQ.P = P;
    PQ.Q = Q;
    PQ.PQd = PQd;
    pprod_init(&PQ.n);
    mpz_add_ui(PQ.n->value, p, 1);

    point_set_str_x(P, "209*i + 332");
    point_printx(P, "xP");

    point_set_str_x(Q, "345*i + 223");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "98*i + 199");
    point_printx(PQd, "xPQd");

    msidh_state_prepare(&m1, A24p, C24, &PQ, A_deg, B_deg, 0);
    msidh_state_prepare(&m2, A24p, C24, &PQ, A_deg, B_deg, 0);
    msidh_state_prepare(&m3, A24p, C24, &PQ, A_deg, B_deg, 0);

    // At least one has to be different -> extremally low chance for false test
    CHECK(m1.secret != m2.secret || m2.secret == m3.secret || m3.secret == m1.secret);

    msidh_state_clear(&m1);
    msidh_state_clear(&m2);
    msidh_state_clear(&m3);
    pprod_clear(&PQ.n);
}

void test_msidh_monte_carlo() {

    struct msidh_state alice, bob;
    msidh_state_init(&alice);
    msidh_state_init(&bob);

    struct tors_basis PQ;
    PQ.P = P;
    PQ.Q = Q;
    PQ.PQd = PQd;
    pprod_init(&PQ.n);
    mpz_add_ui(PQ.n->value, p, 1);

    point_set_str_x(P, "209*i + 332");
    point_printx(P, "xP");

    point_set_str_x(Q, "345*i + 223");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "98*i + 199");
    point_printx(PQd, "xPQd");

    // Test is not deterministic - we try many different combinations of (secret/mask) for both parties
    for (int iter = 0; iter < 100; iter++) {
        msidh_state_prepare(&alice, A24p, C24, &PQ, A_deg, B_deg, 0);
        msidh_state_prepare(&bob, A24p, C24, &PQ, A_deg, B_deg, 1);

        msidh_key_exchange(&alice, bob.pk_A24p, bob.pk_C24, &bob.PQ_other);
        msidh_key_exchange(&bob, alice.pk_A24p, alice.pk_C24, &alice.PQ_other);

        // Make sure that the computed shared secret is the same for both parties
        CHECK(!mpz_cmp(alice.sk_jinv->a, bob.sk_jinv->a) && !mpz_cmp(alice.sk_jinv->b, bob.sk_jinv->b));

        msidh_state_reset(&alice);
        msidh_state_reset(&bob);
    }

    msidh_state_clear(&alice);
    msidh_state_clear(&bob);
}


int main() {
    init_test_variables();

    // General Functions
    TEST_RUN(test_pprod_init());
    TEST_RUN(test_random_unit_sampling_small());
    TEST_RUN(test_random_unit_sampling_large());
    TEST_RUN(test_msidh_gen_pub_params());

    // t = 4 for MSIDH
    setup_params_t4();

    TEST_RUN(test_msidh_internals());
    TEST_RUN(test_msidh_secret_zero());
    TEST_RUN(test_msidh_non_deterministic());
    TEST_RUN(test_msidh_monte_carlo());

    clear_test_variables();

    TEST_RUNS_END;
}