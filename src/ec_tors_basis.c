#include "ec_tors_basis.h"
#include "ec_mont.h"

void tors_basis_init(struct tors_basis *tb) {
    point_init(&tb->P);
    point_init(&tb->Q);
    point_init(&tb->PQd);
    mpz_init(tb->n);
}

void tors_basis_clear(struct tors_basis *tb) {
    point_clear(&tb->P);
    point_clear(&tb->Q);
    point_clear(&tb->PQd);
    mpz_clear(tb->n);
}

void tors_basis_get_subgroup(struct tors_basis *RS, mpz_t n,
                             const struct tors_basis *PQ, const fp2_t A24p,
                             const fp2_t C24) {
    // Temporarily set the order as [N/n] to not introduce additional variable
    mpz_divexact(RS->n, PQ->n, n);
    // Multiply all points in the subbasis by [N/n]
    xLADDER(RS->P, PQ->P, RS->n, A24p, C24);
    xLADDER(RS->Q, PQ->Q, RS->n, A24p, C24);
    xLADDER(RS->PQd, PQ->PQd, RS->n, A24p, C24);
    // Set the proper order of the subgroup torsion basis
    mpz_set(RS->n, n);
}
