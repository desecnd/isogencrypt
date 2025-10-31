#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <openssl/kdf.h>
#include <openssl/sha.h>
#include <openssl/evp.h>
#include <openssl/core_names.h>

#define BUFFER_SIZE 1024
#define IV_SIZE 32

const char *PREFIX_INFO = "\x1b[34m[.]:\x1b[0m";
const char *PREFIX_RUN = "\x1b[33m[%]:\x1b[0m";

// Pretty-print for cmdline context
#define COLCTX(str) "\x1b[36m" str "\x1b[0m"

/*
 * @brief Derive the symmetric shared key for both sides using shared secret from the M-SIDH. This function frees the memory
 */
int derive_key(unsigned char encryption_key[SHA256_DIGEST_LENGTH], char **shared_secret, size_t *shared_secret_len) {
    OSSL_PARAM params[4], *p = params;
    EVP_KDF *kdf = EVP_KDF_fetch(NULL, "HKDF", NULL);
    EVP_KDF_CTX *kctx = EVP_KDF_CTX_new(kdf);
    EVP_KDF_free(kdf);
    int ret = 0;

    *p++ = OSSL_PARAM_construct_utf8_string(OSSL_KDF_PARAM_DIGEST, 
                                            SN_sha256, strlen(SN_sha256));
    *p++ = OSSL_PARAM_construct_octet_string(OSSL_KDF_PARAM_KEY,
                                            *shared_secret, *shared_secret_len);
    *p++ = OSSL_PARAM_construct_octet_string(OSSL_KDF_PARAM_INFO,
                                            "encryption key", (size_t)14);
    *p = OSSL_PARAM_construct_end();

    if (EVP_KDF_derive(kctx, encryption_key, SHA256_DIGEST_LENGTH, params) <= 0) {
        perror("EVP_KDF_derive");
        ret = -1;
        goto err;
    }

    printf("%s Derived shared key using HKDF-SHA256.\n", PREFIX_INFO);

err:
    // Free HKDF context
    EVP_KDF_CTX_free(kctx);

    // Free allocated memory for the shared secret
    free(*shared_secret);
    *shared_secret = NULL;
    *shared_secret_len = 0;
    return ret;
}