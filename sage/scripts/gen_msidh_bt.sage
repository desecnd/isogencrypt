import argparse

from isogencrypt_sage.utils import print_error_and_exit, print_ok, print_info, print_run
from isogencrypt_sage.isogeny import sample_torsion_basis_smooth
from isogencrypt_sage.msidh import MSIDH, MSIDHBenchTask, load_msidh_bench_tasks, store_msidh_bench_tasks
from sage.all import set_random_seed, GF, EllipticCurve

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="gen_msidh_bench_tasks",
        description="Generate MSIDH benchmarks testcases in JSON format"
    )

    parser.add_argument("-r", "--restore", help="Continue generating from previous json file")
    parser.add_argument("-o", "--output", default="msidh_bench_tasks.json", help="JSON file to store the generated benchmark tasks")
    parser.add_argument("--seed", type=int, default=0, help="Setup random seed for the generated tests")
    parser.add_argument("-i", "--init", type=int, default=10, help="Initial value of t parameter")
    parser.add_argument("-e", "--end", type=int, default=301, help="Ending value of t (excluded)")
    parser.add_argument("-s", "--step", type=int, default=10, help="Step value of t in loop iteration")
    args = parser.parse_args()

    bench_tasks = []

    t_values = list(range(args.init, args.end, args.step))

    if args.restore:

        try:
            bench_tasks = load_msidh_bench_tasks(args.restore)
        except ValueError as e:
            print_error_and_exit(str(e))
        except FileNotFoundError as e:
            print_error_and_exit(str(e))

        # Remove the values we already calculated
        for bt in bench_tasks:
            if bt.t in t_values:
                t_values.remove(bt.t)

        print_ok(f"Correctly restored {len(bench_tasks)} bench tasks from: {args.restore}")

    print_run(f"Starting bench tasks generation for t in range({args.init}, {args.end}, {args.step})")
    try:
        for t in t_values:
            print_info(f"Running bench task t={t}")

            a = 6
            p, A, B, _ = MSIDH.gen_pub_params(t)

            F = GF(p**2, names=('i',), modulus=[1, 0, 1])
            (i,) = F._first_ngens(1)

            # Montgomery Curve with 'a' as x^2 coefficient
            E = EllipticCurve(F, [0, a, 0, 1, 0])
            assert E.is_supersingular()

            set_random_seed(args.seed)
            P, Q = sample_torsion_basis_smooth(E, p + 1, 'Q')
            R = P - Q

            bt = MSIDHBenchTask(
                t=t, a=str(a),
                xP=str(P.x()), yP=str(P.y()), 
                xQ=str(Q.x()), yQ=str(Q.y()),
                xR=str(R.x()), yR=str(R.y()),
            )

            bench_tasks.append(bt)
            print_ok(f"Generated MSIDH bench tasks for t={t}")
    except KeyboardInterrupt:
        print_info(f"Interrupted with t={t}")

    bench_tasks = sorted(bench_tasks, key=lambda bt: bt.t)
    store_msidh_bench_tasks(args.output, bench_tasks)
    print_ok(f"Successfully stored {len(bench_tasks)} bench tasks")
