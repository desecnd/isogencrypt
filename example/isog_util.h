#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <openssl/sha.h>

#define BUFFER_SIZE 1024
#define IV_SIZE 32

const char *PREFIX_INFO = "\x1b[34m[.]:\x1b[0m";
const char *PREFIX_RUN = "\x1b[33m[%]:\x1b[0m";

// Pretty-print for cmdline context
#define COLCTX(str) "\x1b[36m" str "\x1b[0m"

/*
 * @brief Derive the symmetric shared key for both sides using shared secret from the M-SIDH. This function frees the memory
 */
void derive_key(unsigned char shared_key[SHA256_DIGEST_LENGTH], char **shared_secret, size_t *shared_secret_len) {

    // Derive the shared key using SHA256
    SHA256((unsigned char *) *shared_secret, *shared_secret_len, shared_key);

    // Free allocated memory for the shared secret
    free(*shared_secret);
    *shared_secret = NULL;
    *shared_secret_len = 0;
    printf("%s Derived shared key using SHA256.\n", PREFIX_INFO);
}