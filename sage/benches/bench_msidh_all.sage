#!/usr/bin/sage -python

from sage.all import GF, EllipticCurve, proof, set_random_seed, order_from_multiple
from isogencrypt_sage.msidh import MSIDH, MSIDHBenchTask
# https://stackoverflow.com/a/25823885
from timeit import default_timer as timer
import sys
import math

# Speed up the calculations
proof.all(False)

from bench_tasks import BENCH_TASKS

def prep_environment(bt: MSIDHBenchTask, check: bool = False):
    global p, A, B, f, i, E, P, Q

    set_random_seed(0)
    p, A, B = MSIDH.calc_pub_params(bt.t, bt.f)
    f = bt.f

    F = GF(p**2, names=('i',), modulus=[1, 0, 1])
    (i,) = F._first_ngens(1)
    # Montgomery Curve with 'a' as x^2 coefficient
    a = F(bt.a)
    E = EllipticCurve(F, [0, a, 0, 1, 0])
    E.set_order((p+1)**2)
    if check:
        assert E.is_supersingular()

    P = E(F(bt.xP), F(bt.yP))
    Q = E(F(bt.xQ), F(bt.yQ))
    R = E(F(bt.xR), F(bt.yR))
    assert P - Q  == R

    if check:
        assert order_from_multiple(P.weil_pairing(Q, p + 1), p + 1, operation='*') == p + 1

def run_benchmark():
    global p, A, B, i, E, P, Q

    N_REPS = 5
    COL_SEP = "\t"

    print("# Sage Benchmark results for MSIDH protocol")
    headers = [ "n", "t", "p_bitsize", "avg", "stddev", "n_reps"]
    print(COL_SEP.join(headers))

    # Optional list of only-checked t values
    tvals_whitelist = [ ]

    try:
        for n, bt_dict in enumerate(BENCH_TASKS):
            # Unpack the params to the function, run the context
            bt = MSIDHBenchTask.from_dict(bt_dict)

            # Apply whitelisting filter
            if tvals_whitelist and bt.t not in tvals_whitelist:
                print(f"Skipping: n={n+1}, t={bt.t}. Not on whitelist")
                continue

            prep_environment(bt)

            times = []

            # Run The benchmark N_REPS times
            for j in range(N_REPS):
                A_mul = (p + 1)//A
                B_mul = (p + 1)//B
                B_basis = (P * B_mul, Q * B_mul)
                Alice = MSIDH(p, A, B, f, E0=E, P=P * A_mul, Q=Q * A_mul, is_bob=False)


                # Calculate public key information (step 1)
                pubkey_alice = Alice.gen_pubkey(*B_basis)

                # Measure single exchange by one of the parties
                start = timer()
                Bob = MSIDH(p, A, B, f, E0=E, P=P * B_mul, Q=Q * B_mul, is_bob=True)
                pubkey_bob = Bob.gen_pubkey(*Alice.basis)
                jinv_bob = Bob.key_exchange(*pubkey_alice)
                end = timer()

                # Exchange - Calculate shared key information (step 2)
                jinv_alice = Alice.key_exchange(*pubkey_bob)
                assert jinv_alice == jinv_bob

                time_s = end - start 
                print(f"[t={bt.t}][{j+1}/{N_REPS}]: Single M-SIDH exchange took {time_s:.3f} seconds to execute", file=sys.stderr)
                times.append(time_s)

            # Make sure that the number of results match
            assert len(times) == N_REPS

            # Prepare the statistics:
            times = sorted(times)
            time_sum = sum(times)
            time_avg = time_sum / len(times)
            time_median = times[N_REPS//2] if N_REPS % 2 else (times[N_REPS//2] + times[(N_REPS + 1)//2])/2
            time_stddev = math.sqrt(sum([ (x - time_avg)**2 for x in times ]) / len(times))

            # Flush will send the result immediately to the File / stream
            print(COL_SEP.join(str(x) for x in [n + 1, bt.t, p.bit_length(), f"{time_avg:.2f}", f"{time_stddev:.2f}", N_REPS]), flush=True)
    except KeyboardInterrupt:
        print("Interrupted. Ending.", file=sys.stderr)


if __name__ == "__main__":
    run_benchmark()
