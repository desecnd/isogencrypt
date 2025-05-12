# from isogencrypt_sage import 

from lib.protocols import MSIDH
from lib.isogeny import sample_torsion_basis_smooth
from sage.all import GF, EllipticCurve, proof, set_random_seed
import timeit
import sys

# Speed up the calculations
proof.all("false")

# Benchmark designed "timeit" python module - no function arguments, global context
def run_msidh_exchange():
    global p, A, B, E, P, Q

    # Prepare the interfaces (step 0)
    A_mul = (p + 1)//A
    Alice = MSIDH(p, A, B, E0=E, P=P * A_mul, Q=Q * A_mul, is_bob=False)
    B_mul = (p + 1)//B
    Bob = MSIDH(p, A, B, E0=E, P=P * B_mul, Q=Q * B_mul, is_bob=True)

    # Calculate public key information (step 1)
    pubkey_alice = Alice.gen_pubkey(*Bob.basis)
    pubkey_bob = Bob.gen_pubkey(*Alice.basis)

    # Exchange - Calculate shared key information (step 2)
    jinv_alice = Alice.key_exchange(*pubkey_bob)
    jinv_bob = Bob.key_exchange(*pubkey_alice)

    # This will break the benchmark computations if result is incorrect
    assert jinv_alice == jinv_bob

def prep_environment(t: int, a: int, basis_P = None, basis_Q = None):
    global p, A, B, i, E, P, Q

    set_random_seed(0)
    p, A, B, _ = MSIDH.gen_pub_params(t)

    F = GF(p**2, names=('i',))
    (i,) = F._first_ngens(1)
    # Montgomery Curve with 'a' as x^2 coefficient
    E = EllipticCurve(F, [0, a, 0, 1, 0])

    P = basis_P
    Q = basis_Q

    # Generate torsion basis:
    # Uses sage's randomness, with set_random_seed called earlier 
    # in the global context should provide deterministic results
    if P is None or Q is None:
        # Create montgomery-friendly basis -> for integration with the .c code
        P, Q = sample_torsion_basis_smooth(E, (p+1), montgomery_basis=True)
        print("Generated Basis P, Q", file=sys.stderr)
        print(f"{P.xy() = }", file=sys.stderr)
        print(f"{Q.xy() = }", file=sys.stderr)
        print(f"{(P-Q).xy() = }", file=sys.stderr, flush=True)
    else:
        P, Q = E(P), E(Q)


def run_benchmark():
    global p, A, B, i, E, P, Q

    N_REPS = 10
    COL_SEP = ";"
    RUN_PARAMS = [
        {"t": 10, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 20, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 30, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 40, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 50, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 60, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 70, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 80, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 90, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 100, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 110, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 120, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 130, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 140, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 150, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 160, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 170, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 180, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 190, "a": 6, "basis_P": None, "basis_Q": None},
        {"t": 200, "a": 6, "basis_P": None, "basis_Q": None},
    ]

    print("# Running Benchmark for MSIDH protocol")
    headers = [ "n", "t_param", "p_bits", "time_avg", "time_median", "time_sum", "n_reps"]
    print(COL_SEP.join(headers))

    for n, params in enumerate(RUN_PARAMS):

        # Unpack the params to the function, run the context
        prep_environment(**params)

        t = params["t"]

        # Run The benchmark
        times = timeit.repeat(
            "run_msidh_exchange()", 
            setup="from __main__ import run_msidh_exchange", 
            # globals={"p": p, "A": A, "B": B, "E": E, "P": P, "Q": Q},
            number=1, 
            repeat=N_REPS
        )

        # Make sure that the number of results match
        assert len(times) == N_REPS

        # Prepare the statistics:
        times = sorted(times)
        p_bits = int(p).bit_length()
        time_sum = sum(times)
        time_avg = time_sum / len(times)
        time_median = times[N_REPS//2] if N_REPS % 2 else (times[N_REPS//2] + times[(N_REPS + 1)//2])/2

        # Flush will send the result immediately to the File / stream
        print(COL_SEP.join(str(x) for x in [n + 1, t, p_bits, f"{time_avg:.2f}", f"{time_median:.2f}", f"{time_sum:.2f}", N_REPS]), flush=True)


if __name__ == "__main__":
    run_benchmark()