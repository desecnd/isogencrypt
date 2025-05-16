from dataclasses import dataclass, asdict
import json
import argparse

from lib.isogeny import sample_torsion_basis_smooth
from lib.protocols import MSIDH
from sage.all import set_random_seed, GF, EllipticCurve

@dataclass
class MSIDHBenchTask:
    t: int
    a: str
    xP: str
    yP: str
    xQ: str
    yQ: str
    xR: str
    yR: str

TEMPLATE_C_PREFIX = """// Auto-generated code using gen_msidh_bench_tasks.sage script
struct bench_task {{
    int t;
    const char *a_str, *xP_str, *xQ_str, *xPQd_str;
}};
#define N_BENCHMARKS {}
const struct bench_task BENCH_TASKS[N_BENCHMARKS] = {{
"""

TEMPLATE_C_SUFFIX = """};
"""

TEMPLATE_C_ITEM = """    {{
        .t        = {t},
        .a_str    = "{a}",
        .xP_str   = "{xP}",
        .xQ_str   = "{xQ}", 
        .xPQd_str = "{xPQd}",
    }}, 
"""

SEED = 0
JSON_FIELD = "bench_tasks"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="gen_msidh_bench_tasks",
        description="Generate MSIDH benchmarks testcases in JSON format"
    )

    parser.add_argument("-r", "--restore", help="Continue generating from previous json file")
    parser.add_argument("-o", "--output", default="msidh_bench_tasks.json", help="JSON file to store the generated benchmark tasks")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Setup random seed for the generated tests")
    parser.add_argument("-c", "--code", action="store_true", default=False, help="Generate C code to include in the sources using template")
    args = parser.parse_args()

    SEED = args.seed
    bench_tasks = []

    t_values = list(range(10, 301, 10))

    if args.restore:

        with open(args.restore) as in_json:
            file_dict = json.loads(in_json.read())

        try:
            bench_tasks = [ MSIDHBenchTask(**obj) for obj in file_dict[JSON_FIELD]]
        except KeyError as _:
            print(f"Given restore file '{args.restore}' JSON object does not contain '{JSON_FIELD}' object")
            exit(1)
        except ValueError as e:
            print(f"Error occured while parsing the given restore file: {args.restore} - {e}")
            exit(1)

        # Remove the values we already calculated
        for bt in bench_tasks:
            if bt.t in t_values:
                t_values.remove(bt.t)

        print(f"[OK]: Correctly restored bench_tasks with {len(bench_tasks)} items.")
    try:
        for t in t_values:

            a = 6
            p, A, B, _ = MSIDH.gen_pub_params(t)

            F = GF(p**2, names=('i',), modulus=[1, 0, 1])
            (i,) = F._first_ngens(1)

            # Montgomery Curve with 'a' as x^2 coefficient
            E = EllipticCurve(F, [0, a, 0, 1, 0])

            set_random_seed(SEED)
            P, Q = sample_torsion_basis_smooth(E, p + 1, True)
            R = P - Q

            bt = MSIDHBenchTask(
                t=t, a=str(a),
                xP=str(P.x()), yP=str(P.y()), 
                xQ=str(Q.x()), yQ=str(Q.y()),
                xR=str(R.x()), yR=str(R.y()),
            )

            bench_tasks.append(bt)
            print(f"[t={t}]: Finished generating MSIDH bench task.")
    except KeyboardInterrupt:
        print(f"Interrupted when t={t}")

    if args.code:
        print(f"Dumping {len(bench_tasks)} bench tasks into C code: {args.output}.c")
        with open(args.output + ".c", "w") as bt_code:
            bt_code.write(TEMPLATE_C_PREFIX.format(len(bench_tasks)))
            for bt in bench_tasks:
                bt_code.write(TEMPLATE_C_ITEM.format(
                    t=bt.t, a=bt.a, xP=bt.xP, xQ=bt.xQ, xPQd=bt.xR
                ))
            bt_code.write(TEMPLATE_C_SUFFIX)

    with open(args.output, 'w') as bt_json:
        print(f"Dumping {len(bench_tasks)} bench tasks into JSON output file: {args.output}")
        bench_tasks = sorted(bench_tasks, key=lambda bt: bt.t)
        bench_tasks = [ asdict(bt) for bt in bench_tasks ]
        bt_json.write(json.dumps({JSON_FIELD: bench_tasks}, indent=4))

