#ifndef TESTING_H
#define TESTING_H

#include <stdio.h>

int g_any_test_failed = 0;
int g_curr_test_error = 0;

#define __CHECK_COND(__cond, __message)                          \
    do {                                                        \
        if (!(__cond)) {                                      \
            fprintf(stderr, "\x1b[31m[!] Error\x1b[0m: %s:%d -> %s\n",                    \
                __FILE__, __LINE__, __message);   \
            g_curr_test_error = 1;                                  \
        }                                                       \
    } while(0)                                                 \

#define CHECK_MSG(__cond, __message) __CHECK_COND(__cond, __message)
#define CHECK(__cond) __CHECK_COND(__cond, "Condition failed.")

#define __TEST_RUN(__func, __test_name) \
    do { \
        g_curr_test_error = 0; \
        /* printf("\x1b[33m[*] Go\x1b[0m: %s\n", __test_name); */ \
        /* print to stdout instead of stderr in order to be present in teest vector */ \
        printf("---: %s\n", __test_name); \
        __func; \
        if (g_curr_test_error) { \
            g_any_test_failed |= 1; \
            fprintf(stderr, "\x1b[31m[!] Error\x1b[0m: %s\n", __test_name); \
        } else { \
            fprintf(stderr, "\x1b[32m[+] Ok\x1b[0m: %s\n", __test_name); \
        } \
    } while (0) \


#define TEST_RUN(__func) __TEST_RUN(__func, #__func)

#define TEST_RUNS_END \
    do { \
        if (g_any_test_failed) {   \
            fprintf(stderr, "\x1b[31m[!] Error\x1b[0m: Tests finished with errors.\n"); \
            return 1; \
        } else { \
            fprintf(stderr, "\x1b[32m[+] Ok\x1b[0m: Tests finished successfully\n"); \
            return 0; \
        } \
    } while (0) \

#endif
