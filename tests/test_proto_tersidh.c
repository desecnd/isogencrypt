#include <gmp.h>
#include <stdio.h>

#include "ec_point_xz.h"
#include "ec_tors_basis.h"
#include "fp.h"
#include "fp2.h"
#include "pprod.h"
#include "proto_tersidh.h"
#include "testing.h"

fp2_t A24p, C24, g_a;
point_t P, Q, PQd, K;
point_t PA, QA, PQAd;
point_t PB, QB, PQBd;
pprod_t A_deg, B_deg;
mpz_t p, g_m;
int g_t, g_f;

void init_test_variables() {
    fp2_init(&A24p);
    fp2_init(&C24);
    fp2_init(&g_a);
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
    mpz_init(g_m);
}

void clear_test_variables() {
    fp2_clear(&A24p);
    fp2_clear(&C24);
    fp2_clear(&g_a);
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
    mpz_clear(g_m);
    fpchar_clear_if_set();
}

void test_tersidh_gen_pub_params() {
    int t, f;

    t = 2;
    f = tersidh_gen_pub_params(p, A_deg, B_deg, t);
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

    t = 10;
    f = tersidh_gen_pub_params(p, A_deg, B_deg, t);
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

    t = 50;
    f = tersidh_gen_pub_params(p, A_deg, B_deg, t);
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

    t = 85;
    f = tersidh_gen_pub_params(p, A_deg, B_deg, t);
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

void setup_params_t2() {
    // p + 1 = 420 = 4 * 3 * 5 * 7
    g_t = 2;
    g_f = tersidh_gen_pub_params(p, A_deg, B_deg, g_t);
    CHECK(g_f == 1);

    // Clear just to be sure
    fpchar_clear_if_set();
    fpchar_setup(p);

    // Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of
    // size 419^2
    fp2_set_uint(g_a, 6);

    fp2_set(A24p, g_a);
    fp2_set_uint(C24, 1);

    // Convert to (A+2 : 4) form used in xDBL
    // Function is safe to call with self-args
    A24p_from_A(A24p, C24, A24p, C24);
}

void test_tersidh_generate_kernel_points() {
    // Testcases uses:  
    // * tersidh_state_init 
    // * tors_basis_got_subgroup

    // Initialize the points
    point_set_str_x(P, "124*i + 4");
    point_printx(P, "xP");

    point_set_str_x(Q, "234*i + 290");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "75*i + 345");
    point_printx(PQd, "xPQd");

    // Constuct Alice (PA, QA) = E0[A]
    struct tors_basis T = { .P = P, .Q = Q, .PQd = PQd};
    mpz_init(T.n);
    mpz_add_ui(T.n, p, 1);

    struct tersidh_state tersidh;
    tersidh_state_init(&tersidh);

    // Set params: t param and Alice perspective
    tersidh.t = g_t;
    tersidh.is_bob = 0;

    // Calculate the torsion basis PA, QA = E0[A]
    tors_basis_get_subgroup(&tersidh.PQ_self, A_deg->value, &T, A24p, C24);

    // Print the values
    point_printx_normalized(tersidh.PQ_self.P, "xPA");
    CHECK(fp2_equal_str(tersidh.PQ_self.P->X, "153*i + 357"));
    point_printx_normalized(tersidh.PQ_self.Q, "xQA");
    CHECK(fp2_equal_str(tersidh.PQ_self.Q->X, "115*i + 221"));

    // Set the starting curve as E0:
    fp2_set(tersidh.A24p_start, A24p);
    fp2_set(tersidh.C24_start,C24);

    // === Test Cases ===

    // 1. s: 0 = '00' -> only KQ gets reduced
    mpz_set_ui(tersidh.secret, 0);
    fp_print(tersidh.secret, "s");

    // KP = PA, KQ = E(0)
    tersidh_generate_kernel_points(&tersidh, 1);
    point_printx_normalized(tersidh.KP, "xKP");
    CHECK(point_equal_str_x(tersidh.KP, "153*i + 357"));
    fp2_print(tersidh.KQ->Z, "ZKQ");
    CHECK(fp2_equal_uint(tersidh.KQ->Z, 0));

    // 2. s: 4 = '11' -> only KP gets reduced
    mpz_set_ui(tersidh.secret, 4);
    fp_print(tersidh.secret, "s");

    // KP = E(0), KQ = QA
    tersidh_generate_kernel_points(&tersidh, 1);
    fp2_print(tersidh.KP->Z, "ZKP");
    CHECK(fp2_equal_uint(tersidh.KP->Z, 0));
    point_printx_normalized(tersidh.KQ, "xKQ");
    CHECK(point_equal_str_x(tersidh.KQ, "115*i + 221"));

    // 3. s: 8 = '22' -> both get reduced to neutral element
    mpz_set_ui(tersidh.secret, 8);
    fp_print(tersidh.secret, "s");

    // Both points get reduced to E(0)
    tersidh_generate_kernel_points(&tersidh, 1);
    fp2_print(tersidh.KP->Z, "ZKP");
    CHECK(fp2_equal_uint(tersidh.KP->Z, 0));
    fp2_print(tersidh.KQ->Z, "ZKQ");
    CHECK(fp2_equal_uint(tersidh.KQ->Z, 0));

    // 4. s: 3 = '10' -> KP = [4]PA, KQ = [5]QA 
    mpz_set_ui(tersidh.secret, 3);
    fp_print(tersidh.secret, "s");

    tersidh_generate_kernel_points(&tersidh, 1);
    point_printx_normalized(tersidh.KP, "xKP");
    CHECK(point_equal_str_x(tersidh.KP, "70*i + 249"));
    point_printx_normalized(tersidh.KQ, "xKQ");
    CHECK(point_equal_str_x(tersidh.KQ, "376*i + 379"));
    
    // 5. s: 5 = '12' 
    mpz_set_ui(tersidh.secret, 5);
    fp_print(tersidh.secret, "s");

    tersidh_generate_kernel_points(&tersidh, 1);
    fp2_print(tersidh.KP->Z, "ZKP");
    CHECK(fp2_equal_uint(tersidh.KP->Z, 0));
    point_printx_normalized(tersidh.KQ, "xKQ");
    CHECK(point_equal_str_x(tersidh.KQ, "376*i + 379"));

    tersidh_state_clear(&tersidh);
    mpz_clear(T.n);
}


void setup_params_t15() {
    g_t = 15;
    g_f = tersidh_gen_pub_params(p, A_deg, B_deg, g_t);
    CHECK(g_f == 1);

    // Clear just to be sure
    fpchar_clear_if_set();
    fpchar_setup(p);

    // Elliptic Curve defined by y^2 = x^3 + 6*x^2 + x over Finite Field in i of
    // size 419^2
    fp2_set_uint(g_a, 6);

    fp2_set(A24p, g_a);
    fp2_set_uint(C24, 1);

    // Convert to (A+2 : 4) form used in xDBL
    // Function is safe to call with self-args
    A24p_from_A(A24p, C24, A24p, C24);
}

void test_tersidh_state_prepare() {

    // Initialize the points
    point_set_str_x(P, "29584706441725156430045723882658927461711344409*i + 8723605936412145621728287379910510577218745101");
    point_printx(P, "xP");

    point_set_str_x(Q, "37774940619389439252585951117910616415029401650*i + 49431932972083700521783639676599147784719313148");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "27017898316266456978360989490543647986396295282*i + 24459173475706578813044034099488276350833204471");
    point_printx(PQd, "xPQd");

    struct tersidh_state tersidh;
    tersidh_state_init(&tersidh);

    // Set static contant secret for the test purpose
    mpz_set_ui(tersidh.secret, 12722590);
    gmp_printf("secret: %Zd\n", tersidh.secret);

    struct tersidh_data td = {
        .t = g_t, .f = g_f, .a = g_a, 
        .xP = P->X, .xQ = Q->X, .xR = PQd->X
    };

    // Run the fist handshake stage
    tersidh_state_prepare(&tersidh, &td, 0);

    gmp_printf("ord(KP): %Zd\n", tersidh.KP_deg->value);
    gmp_printf("ord(KQ): %Zd\n", tersidh.KQ_deg->value);
    CHECK(!mpz_cmp_ui(tersidh.KP_deg->value, 1362295));
    CHECK(!mpz_cmp_ui(tersidh.KQ_deg->value, 1535803748));

    fp2_t a_final;
    fp2_init(&a_final);

    A_from_A24p(tersidh.A24p_pubkey, tersidh.C24_pubkey, tersidh.A24p_pubkey, tersidh.C24_pubkey);
    fp2_div_unsafe(a_final, tersidh.A24p_pubkey, tersidh.C24_pubkey);

    fp2_print(a_final, "aφ(E)");
    CHECK(fp2_equal_str(a_final, "20823762012321546008743277796759767280614403163*i + 6374112813778050444004853308318914927205109921"));

    point_printx_normalized(tersidh.PQ_pubkey.P, "xφ(PB)");
    CHECK(point_equal_str_x(tersidh.PQ_pubkey.P, "20876846029638301559648453625045466417726680845*i + 6210091974667176375575623846556352273019744322"));
    point_printx_normalized(tersidh.PQ_pubkey.Q, "xφ(QB)");
    CHECK(point_equal_str_x(tersidh.PQ_pubkey.Q, "50225951711665673265080604151407102500355509060*i + 25666354126923489688218829391982119080587425206"));

    fp2_clear(&a_final);
    tersidh_state_clear(&tersidh);
}


void test_tersidh_key_exchange() {

    // Initialize the points
    point_set_str_x(P, "45255132863296035939428643087923170526055812335*i + 35207532789640029607392085315164843785886696913");
    point_printx(P, "xP");

    point_set_str_x(Q, "62188135383560125911606431706677411561756802948*i + 52478616152221885238374224345805957897724858098");
    point_printx(Q, "xQ");

    point_set_str_x(PQd, "31403620116220219651357966569215397638854000763*i + 20465179760444544011039140556083357241775723149");
    point_printx(PQd, "xPQd");

    struct tersidh_data td = {
        .t = g_t, .f = g_f, .a = g_a, 
        .xP = P->X, .xQ = Q->X, .xR = PQd->X
    };

    struct tersidh_state alice, bob;
    tersidh_state_init(&alice);
    tersidh_state_init(&bob);

    // Set static contant secret for the test purpose
    mpz_set_ui(alice.secret, 6631513);
    gmp_printf("a_secret: %Zd\n", alice.secret);

    mpz_set_ui(bob.secret, 4980130);
    gmp_printf("b_secret: %Zd\n", bob.secret);

    // Run the fist handshake stage
    tersidh_state_prepare(&alice, &td, 0);

    // Copy public key values
    struct tersidh_data a_pk, b_pk;
    tersidh_data_init(&a_pk); 
    tersidh_data_init(&b_pk); 

    tersidh_get_pubkey(&alice, &a_pk); 

    // Verify calculated Alice public key
    fp2_print(a_pk.a, "aφA(E0)");
    fp2_print(a_pk.xP, "xφA(PB)");
    fp2_print(a_pk.xQ, "xφA(QB)");
    CHECK(fp2_equal_str(a_pk.a, "54801376089903970925339307196694608662317200341*i + 3248665947417223701882116168791158651436364321"));
    CHECK(fp2_equal_str(a_pk.xP, "43868479477697879809566585639320978671149440254*i + 7059168512921389348457558705746610439407082700"));
    CHECK(fp2_equal_str(a_pk.xQ, "59025181202054368531235216709432152582925746866*i + 32671210671678364256011079742850149244716416643"));

    tersidh_state_prepare(&bob, &td, 1);
    tersidh_get_pubkey(&bob, &b_pk); 

    // Verify calculated Bob public key
    fp2_print(b_pk.a, "aφB(E0)");
    fp2_print(b_pk.xP, "xφB(PA)");
    fp2_print(b_pk.xQ, "xφB(QA)");
    CHECK(fp2_equal_str(b_pk.a, "3181532679381964117351769805003588383567570283*i + 8002906843868878916100601361840467966138396994"));
    CHECK(fp2_equal_str(b_pk.xP, "19693539544337260969862178248809154682696090848*i + 32057768515926752112257154360884331701455434315"));
    CHECK(fp2_equal_str(b_pk.xQ, "49661089078839342568423577829915630875462416634*i + 52781414101202328567169277344100576338432492581"));

    tersidh_key_exchange(&alice, &b_pk);
    tersidh_key_exchange(&bob, &a_pk);

    // Verify that calculated j-invariant is valid
    fp2_print(alice.j_inv, "j(EBA)");
    fp2_print(bob.j_inv, "j(EAB)");
    CHECK(fp2_equal(alice.j_inv, bob.j_inv));
    CHECK(fp2_equal_str(alice.j_inv, "16477822473601326380854754948643703255616912674*i + 12600493404726659034043219507047767318775160385"));

    tersidh_state_clear(&alice);
    tersidh_state_clear(&bob);

    tersidh_data_clear(&a_pk); 
    tersidh_data_clear(&b_pk); 
}

int main() {
    init_test_variables();

    // General Functions
    TEST_RUN(test_tersidh_gen_pub_params());

    // t = 2 for TerSIDH
    setup_params_t2();
    TEST_RUN(test_tersidh_generate_kernel_points());
    
    // t = 15 for TerSIDH
    setup_params_t15();
    TEST_RUN(test_tersidh_state_prepare());
    TEST_RUN(test_tersidh_key_exchange());
    
    clear_test_variables();

    TEST_RUNS_END;
}