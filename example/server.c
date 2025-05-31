/*
 * This is a simple 'server' demo presenting the capabilites of isogeny-based
 * key exchange. Code is only a demonstration and should not be used in
 * production systems.
 */

#include <arpa/inet.h>
#include <netinet/in.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <openssl/evp.h>
#include <openssl/sha.h>

#include "sock_msidh.h"

#define BUFFER_SIZE 1024
#define IV_SIZE 32

const char *PREFIX_INFO = "\x1b[34m[.]:\x1b[0m";
const char *PREFIX_RUN = "\x1b[33m[%]:\x1b[0m";

// Pretty-print for cmdline context
#define COLCTX(str) "\x1b[36m" str "\x1b[0m"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <IP_ADDRESS> <PORT>\n", argv[0]);
        exit(1);
    }

    const char *ip_address = argv[1];
    int port = atoi(argv[2]);

    int server_fd, client_fd;
    struct sockaddr_in server_addr, client_addr;
    socklen_t client_len = sizeof(client_addr);

    // Obtain a new socket and handle errors
    server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd < 0) {
        perror("socket");
        exit(1);
    }

    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(port);
    inet_pton(AF_INET, ip_address, &server_addr.sin_addr);

    // Bind to specified address, handle error if something gone wrong
    if (bind(server_fd, (struct sockaddr *)&server_addr, sizeof(server_addr)) <
        0) {
        perror("bind");
        close(server_fd);
        exit(1);
    }

    // Listen for incoming clients on server's sock
    if (listen(server_fd, 1) < 0) {
        perror("listen");
        close(server_fd);
        exit(1);
    }
    printf("%s Server listening on %s:%d\n", PREFIX_RUN, ip_address, port);

    // Accept connection from incomming client
    client_fd = accept(server_fd, (struct sockaddr *)&client_addr, &client_len);
    // Client connected, handle errors
    if (client_fd < 0) {
        perror("accept");
        close(server_fd);
        exit(1);
    }

    printf("%s Client connected\n", PREFIX_INFO);

    // Read the Initialization Vector obtained from client to use in AES-CTR
    unsigned char iv[IV_SIZE];
    if (IV_SIZE != read(client_fd, iv, IV_SIZE)) {
        fprintf(stderr, "Cannot read IV from client\n");
        close(server_fd);
        close(client_fd);
        exit(1);
    }

    printf("%s Received IV data\n", PREFIX_INFO);
    printf("%s Isogeny Handshake...\n", PREFIX_RUN);

    // Run MSIDH handshake with specified level
    unsigned char shared_key[SHA256_DIGEST_LENGTH];
    int status = msidh_handshake(client_fd, 0, shared_key, MSIDH_T150);
    if (status < 0) {
        fprintf(stderr, "MSIDH handshake returned with errors.\n");
        close(server_fd);
        close(client_fd);
        exit(1);
    }

    printf("%s Handshake Completed.\n", PREFIX_INFO);
    printf(COLCTX("--- Begin Encrypted Channel ---\n"));

    // Initialize the AES-CTR decryption context
    EVP_CIPHER_CTX *dec_ctx = EVP_CIPHER_CTX_new();
    EVP_DecryptInit_ex(dec_ctx, EVP_aes_256_ctr(), NULL, shared_key, iv);

    // Encrypted input comming from client
    char enc_input[BUFFER_SIZE];
    // Decrypted client plaintext
    char buffer[BUFFER_SIZE];

    while (1) {
        printf(COLCTX("B> "));
        fflush(stdout);

        int n_bytes = read(client_fd, enc_input, BUFFER_SIZE);
        if (n_bytes <= 0) {
            break;
        } else if (n_bytes < BUFFER_SIZE) {
            buffer[n_bytes] = enc_input[n_bytes] = '\0';
        }

        int n_decrypted;
        EVP_DecryptUpdate(dec_ctx, (unsigned char *)buffer, &n_decrypted,
                          (unsigned char *)enc_input, n_bytes);
        printf("%.*s\n", BUFFER_SIZE, buffer);
    }

    printf(COLCTX("\n------------- End -------------\n"));

    EVP_CIPHER_CTX_free(dec_ctx);
    close(client_fd);
    close(server_fd);
    return 0;
}
