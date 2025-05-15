
from lib.isogeny import sample_torsion_basis_smooth
from lib.protocols import MSIDH
from sage.all import set_random_seed, GF, EllipticCurve

template_code = """
#define N_BENCHMARKS {}
const struct bench_task BENCH_TASKS[N_BENCHMARKS] = {{
    {}
}};
"""

template_c = """    {{
        .t        = {0},
        .a        = {1},
        .xP_str   = "{2}",
        .xQ_str   = "{3}", 
        .xPQd_str = "{4}",
    }}, 
"""

template_sage = """    # MSIDH params for t = {t}
    p = {p}
    a = {a}
    F = GF(p**2, names=('i',), modulus=[1, 0, 1])
    (i,) = F._first_ngens(1)
    E = EllipticCurve(F, [0, a, 0, 1, 0])
    P = E({xP}, {yP})
    Q = E({xQ}, {yQ})
"""

if __name__ == "__main__":

    SEED = 0
    a = 6
    # for t in range(10, 601, 10):
    for t in range(10, 270, 10):

        p, A, B, _ = MSIDH.gen_pub_params(t)

        F = GF(p**2, names=('i',), modulus=[1, 0, 1])
        (i,) = F._first_ngens(1)

        # Montgomery Curve with 'a' as x^2 coefficient
        E = EllipticCurve(F, [0, a, 0, 1, 0])

        set_random_seed(SEED)
        P, Q = sample_torsion_basis_smooth(E, p + 1, True)

        # Add new bench_struct to the array
        print(template_c.format(t, a, P.x(), Q.x(), (P - Q).x()), flush=True, end='')
        # print(template_sage.format(t=t, p=p, a=a, xP=P.x(), yP=P.y(), xQ=Q.x(), yQ=Q.y()))

    # print(template_code.format(len(structs), ''.join(structs)))


        

