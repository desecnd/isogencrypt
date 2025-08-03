#/usr/bin/sage

import argparse
import pathlib

from isogencrypt_sage.msidh import load_msidh_bench_tasks
from isogencrypt_sage.utils import print_error_and_exit, print_info, print_ok

TEMPLATE_C_PREFIX = """// Auto-generated code using conv_msidh_bt.sage script
struct bench_task {{
    int t, f;
    const char *a_str, *xP_str, *xQ_str, *xPQd_str;
}};
#define N_BENCHMARKS {}
const struct bench_task BENCH_TASKS[N_BENCHMARKS] = {{
"""

TEMPLATE_C_SUFFIX = """};
"""

TEMPLATE_C_ITEM = """    {{
        .t        = {t},
        .f        = {f},
        .a_str    = "{a}",
        .xP_str   = "{xP}",
        .xQ_str   = "{xQ}", 
        .xPQd_str = "{xPQd}",
    }}, 
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="conv_msidh_bt",
        description="Convert MSIDH benchmarks testcases into format valid for C language"
    )

    parser.add_argument("bench_tasks")
    parser.add_argument("-o", "--output", help="Output file or directory for generated C code")
    args = parser.parse_args()

    try:
        bench_tasks = load_msidh_bench_tasks(args.bench_tasks)
    except FileNotFoundError as e:
        print_error_and_exit(str(e))
    
    print_ok(f"Loaded {len(bench_tasks)} bench tasks from: {args.bench_tasks}")

    dest_file = pathlib.Path(args.bench_tasks).with_suffix(".c")

    if args.output is not None:
        out = pathlib.Path(args.output)
        if out.is_dir():
            dest_file = out / dest_file.name
        else:
            dest_file = out

    print_info(f"Dumping converted C code into: {dest_file}")

    with open(dest_file, "w") as bt_code:
        bt_code.write(TEMPLATE_C_PREFIX.format(len(bench_tasks)))
        for bt in bench_tasks:
            bt_code.write(TEMPLATE_C_ITEM.format(
                t=bt.t, f=bt.f, a=bt.a, xP=bt.xP, xQ=bt.xQ, xPQd=bt.xR
            ))
        bt_code.write(TEMPLATE_C_SUFFIX)

    print_ok("Success")

