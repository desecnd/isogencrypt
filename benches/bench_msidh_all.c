#include "bench_msidh.h"
#include <stdio.h>

int main() {

    printf("# C Benchmark results for MSIDH protocol\n");
    printf("n\tt\tp_bitsize\tavg\tstddev\tn_reps\n");

    struct benchmark_data bd;

    for (int j = 0; j < N_BENCHMARKS; j++) {
        run_benchmark(&BENCH_TASKS[j], &bd);
        printf("%d\t%d\t%d\t%0.2lf\t%0.2lf\t%d\n", j + 1, BENCH_TASKS[j].t,
               bd.p_bitsize, bd.average, bd.stddev, N_REPS);
        fflush(stdout);
    }
}