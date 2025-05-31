/*
 * This is a simple 'client' demo presenting the capabilites of isogeny-based
 * key exchange. Code is only a demonstration and should not be used in
 * production systems.
 */

#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/sha.h>

#include "sock_msidh.h"

#define IV_SIZE 32
#define BUFFER_SIZE 1024

const char *PREFIX_INFO = "\x1b[34m[.]:\x1b[0m";
const char *PREFIX_RUN = "\x1b[33m[%]:\x1b[0m";

// Pretty-print for cmdline context
#define COLCTX(str) "\x1b[36m" str "\x1b[0m"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <SERVER_IP> <PORT>\n", argv[0]);
        exit(1);
    }

    const char *server_ip = argv[1];
    int port = atoi(argv[2]);

    int sock_fd;
    struct sockaddr_in server_addr;

    // Obtain the new socket file descriptor, handle the errors
    sock_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (sock_fd < 0) {
        perror("socket");
        exit(1);
    }

    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(port);
    inet_pton(AF_INET, server_ip, &server_addr.sin_addr);

    // Connect to the specified address on given socket
    if (connect(sock_fd, (struct sockaddr *)&server_addr, sizeof(server_addr)) <
        0) {
        perror("connect");
        close(sock_fd);
        exit(1);
    }

    printf("%s Connected to server at %s:%d\n", PREFIX_INFO, server_ip, port);

    // Generate the Initialization Vector for AES-CTR
    unsigned char iv[IV_SIZE];
    // Fill IV_SIZE cryptographically secure random bytes from the system.
    // Returns 1 on success, 0 otherwise.
    if (RAND_bytes(iv, IV_SIZE) != 1) {
        perror("RAND_bytes");
        exit(1);
    }

    // Send the IV data to the server
    if (IV_SIZE != write(sock_fd, iv, IV_SIZE)) {
        perror("write(iv)");
        exit(1);
    }

    printf("%s Sent IV to server.\n", PREFIX_INFO);
    printf("%s Isogeny Handshake...\n", PREFIX_RUN);

    // Run MSIDH handshake with specified level
    unsigned char shared_key[SHA256_DIGEST_LENGTH];
    int status = msidh_handshake(sock_fd, 1, shared_key, MSIDH_T150);
    if (status < 0) {
        fprintf(stderr, "MSIDH handshake returned with errors.\n");
        close(sock_fd);
        exit(1);
    }

    printf("%s Handshake Completed.\n", PREFIX_INFO);
    printf(COLCTX("--- Begin Encrypted Channel ---\n"));

    // Initialize the AES-CTR encryption context
    EVP_CIPHER_CTX *enc_ctx = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(enc_ctx, EVP_aes_256_ctr(), NULL, shared_key, iv);

    // User input
    char input[BUFFER_SIZE];
    // Encrypted input sent to the server
    char enc_buffer[BUFFER_SIZE];

    while (1) {
        printf(COLCTX("A> "));

        if (!fgets(input, BUFFER_SIZE, stdin)) {
            break;
        }
        size_t len = strlen(input);
        if (input[len - 1] == '\n')
            input[len - 1] = '\0';

        int n_encrypted;
        EVP_EncryptUpdate(enc_ctx, (unsigned char *)enc_buffer, &n_encrypted,
                          (unsigned char *)input, len);

        if (n_encrypted != write(sock_fd, enc_buffer, n_encrypted)) {
            break;
        }
    }

    printf(COLCTX("\n------------- End -------------\n"));

    EVP_CIPHER_CTX_free(enc_ctx);
    close(sock_fd);

    return 0;
}
