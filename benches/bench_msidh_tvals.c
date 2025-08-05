#include "bench_msidh.h"
#include <stdio.h>
#include "fp2.h"

int main() {
    fp2_layer_ctx_init();

    printf("# C Benchmark results for MSIDH protocol\n");
    printf("n\tt\tp_bitsize\tavg\tstddev\tn_reps\n");

    int t_values[] = {100, 200, 300, 400};
    const int N_RUNS = sizeof(t_values) / sizeof(int);

    struct benchmark_data bd;

    for (int i = 0; i < N_RUNS; i++) {
        int t_found = 0;
        for (int j = 0; !t_found && j < N_BENCHMARKS; j++) {
            if (BENCH_TASKS[j].t != t_values[i])
                continue;

            t_found = 1;
            run_benchmark(&BENCH_TASKS[j], &bd);
            printf("%d\t%d\t%d\t%0.2lf\t%0.2lf\t%d\n", i + 1, t_values[i],
                   bd.p_bitsize, bd.average, bd.stddev, N_REPS);
            fflush(stdout);
        }

        if (!t_found) {
            fprintf(stderr, "Cannot find BenchTask for MSIDH param t=%d\n",
                    t_values[i]);
        }
    }

    fp2_layer_ctx_init();
}