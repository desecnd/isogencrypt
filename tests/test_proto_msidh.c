#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "ec_mont.h"
#include "fp2.h"
#include "pprod.h"
#include "proto_msidh.h"
#include "testing.h"

fp2_t A24p, C24, a0;
point_t P, Q, PQd, K;
point_t PA, QA, PQAd;
point_t PB, QB, PQBd;
pprod_t A_deg, B_deg;
mpz_t p, m;
unsigned int g_t;

void init_test_variables() {
    fp2_init(&A24p);
    fp2_init(&C24);
    fp2_init(&a0);
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
    fp2_clear(&a0);
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
    CHECK_MSG(!mpz_cmp_ui(M->value, 223092870),
              "Incorrect product of primes value");

    pprod_clear(&M);
}

void test_random_unit_sampling_large() {
    srand(0xdeafbeef);

    // 100 prime numbers
    unsigned int primes_large[] = {
        2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,
        43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101,
        103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
        173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
        241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
        317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
        401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
        479, 487, 491, 499, 503, 509, 521, 523, 541};

    pprod_t M;
    pprod_init(&M);
    pprod_set_array(M, primes_large,
                    sizeof(primes_large) / sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret,
                  "sample_quadratic_root_of_unity returned non-zero value");
        if (ret)
            break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);

        CHECK_MSG(!mpz_cmp_ui(result, 1),
                  "Sampling result is not a square root of unity");
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
    pprod_set_array(M, primes_small,
                    sizeof(primes_small) / sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret,
                  "sample_quadratic_root_of_unity returned non-zero value");
        if (ret)
            break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);

        CHECK_MSG(!mpz_cmp_ui(result, 1),
                  "Sampling result is not a square root of unity");
    }

    mpz_clear(result);
    pprod_clear(&M);
}

void test_msidh_gen_pub_params() {
    int t, f;

    t = 4;
    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    assert(f > 0);
    printf("t: %d\n", t);
    fp_print(p, "p");
    fp_print(A_deg->value, "A");
    fp_print(B_deg->value, "B");
    printf("f: %d\n", f);

    // (1115881660253397921934830779, 20312793853220, 54934917782199, 1)
    CHECK(fp_equal_str(p, "419"));
    CHECK(fp_equal_str(A_deg->value, "20"));
    CHECK(fp_equal_str(B_deg->value, "21"));
    CHECK(f == 1);

    t = 20;
    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    assert(f > 0);
    printf("t: %d\n", t);
    fp_print(p, "p");
    fp_print(A_deg->value, "A");
    fp_print(B_deg->value, "B");
    printf("f: %d\n", f);

    // (1115881660253397921934830779, 20312793853220, 54934917782199, 1)
    CHECK(fp_equal_str(p, "1115881660253397921934830779"));
    CHECK(fp_equal_str(A_deg->value, "20312793853220"));
    CHECK(fp_equal_str(B_deg->value, "54934917782199"));
    CHECK(f == 1);

    t = 100;
    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    assert(f > 0);
    printf("t: %d\n", t);
    fp_print(p, "p");
    fp_print(A_deg->value, "A");
    fp_print(B_deg->value, "B");
    printf("f: %d\n", f);

    CHECK(fp_equal_str(
        p, "8575714055829256614755727859263673968077446087605609446743315408101"
           "8759112379622825936264472040913211898391694579544275450839721172773"
           "3337046938979806269649383591155750277998320880157855312245874492523"
           "686094863130587658379"));
    CHECK(
        fp_equal_str(A_deg->value,
                     "374801510478146654390916851616342588666386301694775232095"
                     "65100266538230879013060067601181882604603923651275420"));
    CHECK(
        fp_equal_str(B_deg->value,
                     "251436062458500732842543840100611612168366079689161565071"
                     "573990525611946403712566989284087642734230643996241279"));
    CHECK(f == 91);

    t = 170;
    f = msidh_gen_pub_params(p, A_deg, B_deg, t);
    assert(f > 0);
    printf("t: %d\n", t);
    fp_print(p, "p");
    fp_print(A_deg->value, "A");
    fp_print(B_deg->value, "B");
    printf("f: %d\n", f);

    CHECK(fp_equal_str(
        p, "0x76af20c40b3a9608503c0a977a5aa4b0166b507d7acae27c768f24e70445d8be4"
           "e8b28bc0d3e0cac0099671b40de4af52aeb44c481dfa10cd022053fced2df57c0ac"
           "ae54418987d9249919b738514f78cf1ea9fd8229b05e6084e53d5c603fdd373ab00"
           "bd558005996240cff9c20701a1b8a6da432ea20b2fc988d5889afb10d27bb55a19f"
           "39d71d3a7c45732d6bd265291103266d3b8909b4a1f1a0b1b17cab3c477aa428a5e"
           "e13fbf2d896e7b4e777"));
    CHECK(fp_equal_str(
        A_deg->value,
        "0x65cff008b4f1c2198f9588e205714bcf98a8d48a43a17950b50ef9348bb24543d426"
        "061ff7ac66d01171762b758a9db7cb8aa493cfe138ad8e8a6b74a2c1f38100b6a3ccaf"
        "116f6cb98817202d181d15521848734125df4"));
    CHECK(fp_equal_str(
        B_deg->value,
        "0x3a3a98e20a937f99bec60b3d3f9609b13c9c88cc6db1f40270e1e023752b58a25b4b"
        "b6722532c60ca32b3c66f62587003c5462a3dff7d6e0afd455e2e22c69d4167335e99b"
        "dbf9d8f84375d7a552fa4d8852951eb3688a63"));
    CHECK(f == 82);
}

void setup_params_t4() {
    // p + 1 = 420 = 4 * 3 * 5 * 7
    g_t = 4;
    int f = msidh_gen_pub_params(p, A_deg, B_deg, g_t);
    CHECK(f == 1);

    // Clear just to be sure
    fpchar_clear_if_set();
    fpchar_setup(p);

    // Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of
    // size 419^2
    fp2_set_uint(a0, 6);

    fp2_set(A24p, a0);
    fp2_set_uint(C24, 1);

    // Convert to (A+2 : 4) form used in xDBL
    // Function is safe to call with self-args
    A24p_from_A(A24p, C24, A24p, C24);
}

void test_msidh_secret_zero() {
    point_set_str_x(P, "295*i + 398");
    point_printx(P, "xP");

    point_set_str_x(Q, "314*i + 149");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "29*i + 395");
    point_printx(PQd, "xPQd");

    mpz_set_ui(m, 0);
    xLADDER3PT(P, Q, PQd, m, A24p, C24);

    // P should not change
    // xP: 295*i + 398
    point_printx_normalized(P, "xK");
    CHECK(!mpz_cmp_ui(P->X->a, 398) && !mpz_cmp_ui(P->X->b, 295));

    xLADDER3PT_int(P, Q, PQd, 0, A24p, C24);

    // P should not change
    // xP: 295*i + 398
    point_printx_normalized(P, "xK");
    CHECK(!mpz_cmp_ui(P->X->a, 398) && !mpz_cmp_ui(P->X->b, 295));
}

void test_msidh_internals() {
    point_set_str_x(P, "295*i + 398");
    point_printx(P, "xP");

    point_set_str_x(Q, "314*i + 149");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "29*i + 395");
    point_printx(PQd, "xPQd");

    struct tors_basis PQ = {.P = P, .Q = Q, .PQd = PQd};
    pprod_init(&PQ.n);

    // n = p + 1
    // TODO: WARNING! not a valid pprod type! We just use it
    mpz_set(PQ.n->value, p);
    mpz_add_ui(PQ.n->value, PQ.n->value, 1);

    // Constuct Alice Basis (PA, QA) = [n//A](P, Q).
    struct tors_basis PQA = {.P = PA, .Q = QA, .PQd = PQAd};
    pprod_init(&PQA.n);
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24);

    point_printx_normalized(PA, "xPA");
    point_printx_normalized(QA, "xQA");
    point_printx_normalized(PQAd, "xPQAd");

    // Constuct Bob Basis (PB, QB) = [n//B](P, Q)
    struct tors_basis PQB = {.P = PB, .Q = QB, .PQd = PQBd};
    pprod_init(&PQB.n);
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24);

    point_printx_normalized(PB, "xPB");
    point_printx_normalized(QB, "xQB");
    point_printx_normalized(PQBd, "xPQBd");

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

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p): 271*i + 111
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 111) && !mpz_cmp_ui(aE_alice->b, 271));

    // xφ(PB): 363*i + 349
    point_printx_normalized(PB, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PQB.P->X->a, 349) && !mpz_cmp_ui(PB->X->b, 363));

    // xφ(QB): 303*i + 391
    point_printx_normalized(QB, "xφ(QB)");
    CHECK(!mpz_cmp_ui(PQB.Q->X->a, 391) && !mpz_cmp_ui(QB->X->b, 303));

    // xφ(PB-QB): 239*i + 90
    point_printx_normalized(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 90) && !mpz_cmp_ui(PQBd->X->b, 239));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
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

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p): 408*i + 332
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 332) && !mpz_cmp_ui(aE_alice->b, 408));

    // xφ(PB): 192*i + 82
    point_printx_normalized(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 82) && !mpz_cmp_ui(PB->X->b, 192));

    // xφ(QB): 299*i + 211
    point_printx_normalized(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 211) && !mpz_cmp_ui(QB->X->b, 299));

    // xφ(PQBd): 127*i + 181
    point_printx_normalized(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 181) && !mpz_cmp_ui(PQBd->X->b, 127));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
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

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p): 109*i + 386
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 386) && !mpz_cmp_ui(aE_alice->b, 109));

    // xφ(PB): 298*i + 413
    point_printx_normalized(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 413) && !mpz_cmp_ui(PB->X->b, 298));

    // xφ(QB): 323*i + 222
    point_printx_normalized(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 222) && !mpz_cmp_ui(QB->X->b, 323));

    // xφ(PQBd): 323*i + 66
    point_printx_normalized(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 66) && !mpz_cmp_ui(PQBd->X->b, 323));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
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

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p): 16*i + 353
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(!mpz_cmp_ui(aE_alice->a, 353) && !mpz_cmp_ui(aE_alice->b, 16));

    // xφ(PB): 340
    point_printx_normalized(PQB.P, "xφ(PB)");
    CHECK(!mpz_cmp_ui(PB->X->a, 340) && !mpz_cmp_ui(PB->X->b, 0));

    // xφ(QB): 97*i + 356
    point_printx_normalized(PQB.Q, "xφ(QB)");
    CHECK(!mpz_cmp_ui(QB->X->a, 356) && !mpz_cmp_ui(QB->X->b, 97));

    // xφ(PQBd): 330*i + 244
    point_printx_normalized(PQBd, "xφ(PQBd)");
    CHECK(!mpz_cmp_ui(PQBd->X->a, 244) && !mpz_cmp_ui(PQBd->X->b, 330));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
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

    point_set_str_x(P, "209*i + 332");
    point_printx_normalized(P, "xP");

    point_set_str_x(Q, "345*i + 223");
    point_printx_normalized(Q, "xQ");

    point_set_str_x(PQd, "98*i + 199");
    point_printx_normalized(PQd, "xPQd");

    struct msidh_data md = {
        .t = g_t, .a = a0, .xP = P->X, .xQ = Q->X, .xR = PQd->X};

    msidh_state_prepare(&m1, &md, 0);
    msidh_state_prepare(&m2, &md, 0);
    msidh_state_prepare(&m3, &md, 0);

    // At least one has to be different -> extremally low chance for false test
    CHECK(m1.secret != m2.secret || m2.secret == m3.secret ||
          m3.secret == m1.secret);

    msidh_state_clear(&m1);
    msidh_state_clear(&m2);
    msidh_state_clear(&m3);
}

void test_msidh_monte_carlo() {

    struct msidh_state alice, bob;
    msidh_state_init(&alice);
    msidh_state_init(&bob);

    point_set_str_x(P, "209*i + 332");
    point_printx(P, "xP");

    point_set_str_x(Q, "345*i + 223");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "98*i + 199");
    point_printx(PQd, "xPQd");

    struct msidh_data md = {
        .t = g_t, .a = a0, .xP = P->X, .xQ = Q->X, .xR = PQd->X};

    struct msidh_data alice_pk, bob_pk;
    msidh_data_init(&alice_pk);
    msidh_data_init(&bob_pk);

    // Test is not deterministic - we try many different combinations of
    // (secret/mask) for both parties
    for (int iter = 0; iter < 100; iter++) {
        msidh_state_prepare(&alice, &md, 0);
        msidh_state_prepare(&bob, &md, 1);

        msidh_get_pubkey(&alice, &alice_pk);
        msidh_get_pubkey(&bob, &bob_pk);

        msidh_key_exchange(&alice, &bob_pk);
        msidh_key_exchange(&bob, &alice_pk);

        CHECK(fp2_equal(alice.j_inv, bob.j_inv));

        msidh_state_reset(&alice);
        msidh_state_reset(&bob);
    }

    msidh_data_clear(&alice_pk);
    msidh_data_clear(&bob_pk);

    msidh_state_clear(&alice);
    msidh_state_clear(&bob);
}

void setup_params_t30() {
    g_t = 30;
    int f = msidh_gen_pub_params(p, A_deg, B_deg, g_t);
    assert(f > 0);

    // Clear just to be sure
    fpchar_clear_if_set();
    fpchar_setup(p);

    // Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of
    // size 419^2
    fp2_set_uint(a0, 6);

    fp2_set(A24p, a0);
    fp2_set_uint(C24, 1);

    // Convert to (A+2 : 4) form used in xDBL
    // Function is safe to call with self-args
    A24p_from_A(A24p, C24, A24p, C24);
}

void test_msidh_internals_large() {

    point_set_str_x(P, "32381872305678490404833289490608450363172933700*i + "
                       "38566570518230924614310417068523536638018699310");
    point_set_str_x(Q, "29454235622145096109316297773070819902970029047*i + "
                       "17242937661247998401353436361850378272505949076");
    point_set_str_x(PQd, "29746668073241433805825980414131965658693152428*i + "
                         "17760541352524573821929254886145960108078787218");

    point_printx(P, "xP");
    point_printx(Q, "xQ");
    point_printx(PQd, "xPQd");

    // Setup Torsion basis
    struct tors_basis PQ = {.P = P, .Q = Q, .PQd = PQd};

    pprod_init(&PQ.n);
    mpz_add_ui(PQ.n->value, p, 1);

    // Constuct Alice Basis (PA, QA) = [n//A](P, Q).
    struct tors_basis PQA = {.P = PA, .Q = QA, .PQd = PQAd};
    pprod_init(&PQA.n);
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24);

    point_printx_normalized(PA, "xPA");
    point_printx_normalized(QA, "xQA");
    point_printx_normalized(PQAd, "xPQAd");

    CHECK(point_equal_str_x(
        PA, "29699839172659064070785630633025116007167774796*i + "
            "2504062641765751109972582500320053624879816710"));
    CHECK(point_equal_str_x(
        QA, "33271862144249447985137496787936923322336452945*i + "
            "61213944166497621905178667109855815871670127212"));
    CHECK(point_equal_str_x(
        PQAd, "16754712571370394119457689413046119537544122176*i + "
              "25057305442840430077505557836192648418816321985"));

    // Constuct Bob Basis (PB, QB) = [n//B](P, Q)
    struct tors_basis PQB = {.P = PB, .Q = QB, .PQd = PQBd};
    pprod_init(&PQB.n);
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24);

    point_printx_normalized(PB, "xPB");
    point_printx_normalized(QB, "xQB");
    point_printx_normalized(PQBd, "xPQBd");

    CHECK(point_equal_str_x(
        PB, "49114309482119019874119284972332462593113788145*i + "
            "4761433153046662093084589855451688595958230559"));
    CHECK(point_equal_str_x(
        QB, "30431121762650905008454289312347872723430751013*i + "
            "2343460809488608363031424186950047556179737763"));
    CHECK(point_equal_str_x(
        PQBd, "23391135694148403625145265581278108641686825034*i + "
              "26976742350009751812962081032240052720831329725"));

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

    // 1. --- Testcase 1
    mpz_set_str(a_sec, "11349330453264090324638", 10);
    mpz_set_str(a_mask, "131339333569636892083385", 10);
    mpz_set_str(b_sec, "72593563569111258510719", 10);
    mpz_set_str(b_mask, "17334875921447289106681", 10);
    fp_print(a_sec, "A_sec");
    fp_print(a_mask, "A_mask");
    fp_print(b_sec, "B_sec");
    fp_print(b_mask, "B_mask");

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p)
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(fp2_equal_str(aE_alice,
                        "58406864750297447375941148275890126964932964189*i + "
                        "60400200391283286845520983845656557221610499288"));

    // xφ(PB)
    point_printx_normalized(PQB.P, "xφ(PB)");
    CHECK(point_equal_str_x(
        PQB.P, "61460438692052277410377272364685473376640776051*i + "
               "55514416620914970543945665178326301897479630385"));

    // xφ(QB)
    point_printx_normalized(PQB.Q, "xφ(QB)");
    CHECK(point_equal_str_x(
        PQB.Q, "59420147577812126043911530847798020708549199424*i + "
               "45153371532976521708986975584064598899456381726"));

    // xφ(PQBd)
    point_printx_normalized(PQB.PQd, "xφ(PQBd)");
    CHECK(point_equal_str_x(
        PQB.PQd, "12338170382935113028682312433735367927625301703*i + "
                 "47299240837639061178684012620848485232866182686"));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p)
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(fp2_equal_str(aE_final,
                        "30563038667104601337368965988575062985170253283*i + "
                        "22273408037389398982670537391321135761023812334"));

    // jτ(φ(E))
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(fp2_equal_str(j_inv,
                        "9661335542795285403316671192384434122935976866*i + "
                        "5474694925468985143408009855884208752167844644"));

    // Testcase 2 -----
    // Reset the second testscase
    tors_basis_get_subgroup(&PQA, A_deg, &PQ, A24p, C24);
    tors_basis_get_subgroup(&PQB, B_deg, &PQ, A24p, C24);

    mpz_set_str(a_sec, "102959793904566459067664", 10);
    mpz_set_str(a_mask, "234614427428290727103469", 10);
    mpz_set_str(b_sec, "271374038831227351065368", 10);
    mpz_set_str(b_mask, "93531163257802324154959", 10);
    fp_print(a_sec, "A_sec");
    fp_print(a_mask, "A_mask");
    fp_print(b_sec, "B_sec");
    fp_print(b_mask, "B_mask");

    _msidh_gen_pubkey_alice(A24p_alice, C24_alice, &PQA, &PQB, A24p, C24, a_sec,
                            a_mask);

    // aφ(E)(24p)
    fp2_div_unsafe(aE_alice, A24p_alice, C24_alice);
    fp2_print(aE_alice, "aφ(E)(24p)");
    CHECK(fp2_equal_str(aE_alice,
                        "32721529036899972732279468408471975069661362292*i + "
                        "17269325977833293966138407821489567068853821114"));

    // xφ(PB)
    point_printx_normalized(PQB.P, "xφ(PB)");
    CHECK(point_equal_str_x(
        PQB.P, "16875856946578546258578008123468502154460718295*i + "
               "38847645774315915434007838296654577037686340296"));

    // xφ(QB)
    point_printx_normalized(PQB.Q, "xφ(QB)");
    CHECK(point_equal_str_x(
        PQB.Q, "16054449845823153742619853253831785845856263425*i + "
               "15381476167260289678857731510026166833683607342"));

    // xφ(PQBd)
    point_printx_normalized(PQB.PQd, "xφ(PQBd)");
    CHECK(point_equal_str_x(
        PQB.PQd, "51405460425642783842251929023025336502575438127*i + "
                 "11702887079976981107716609038554493137279265981"));

    _msidh_key_exchange_alice(j_inv, A24p_final, C24_final, A24p_alice,
                              C24_alice, &PQB, b_sec);
    fp2_div_unsafe(aE_final, A24p_final, C24_final);

    // aτ(φ(E))(24p)
    fp2_print(aE_final, "aτ(φ(E))(24p)");
    CHECK(fp2_equal_str(aE_final,
                        "995903660540748720123956841665787112221678284*i + "
                        "52717497346814873517492805505778428959241213199"));

    // jτ(φ(E))
    fp2_print(j_inv, "jτ(φ(E))");
    CHECK(fp2_equal_str(j_inv,
                        "41919959196404182164998159522259445819384604188*i + "
                        "46346974086814678100459816796141382648465154561"));

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

int main() {
    init_test_variables();

    // General Functions
    TEST_RUN_SILENT(test_pprod_init());
    TEST_RUN_SILENT(test_random_unit_sampling_small());
    TEST_RUN_SILENT(test_random_unit_sampling_large());
    TEST_RUN(test_msidh_gen_pub_params());

    // t = 4 for MSIDH
    setup_params_t4();

    TEST_RUN(test_msidh_internals());
    TEST_RUN(test_msidh_secret_zero());
    TEST_RUN(test_msidh_non_deterministic());

    // This test will "override" the characteristic
    TEST_RUN(test_msidh_monte_carlo());

    // t = 30 for MSIDH
    setup_params_t30();

    TEST_RUN(test_msidh_internals_large());

    clear_test_variables();

    TEST_RUNS_END;
}