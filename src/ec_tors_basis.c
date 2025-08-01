#include "ec_tors_basis.h"
#include "ec_mont.h"

void tors_basis_init(struct tors_basis *tb) {
    point_init(&tb->P);
    point_init(&tb->Q);
    point_init(&tb->PQd);
    pprod_init(&tb->n);
}

void tors_basis_clear(struct tors_basis *tb) {
    point_clear(&tb->P);
    point_clear(&tb->Q);
    point_clear(&tb->PQd);
    pprod_clear(&tb->n);
}

/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q)
 * of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, pprod_t n,
                             const struct tors_basis *PQ, const fp2_t A24p,
                             const fp2_t C24) {
    // Temporarily set the order as [N/n]
    mpz_divexact(PQsub->n->value, PQ->n->value, n->value);
    // Multiply all points in the subbasis by [N/n]
    xLADDER(PQsub->P, PQ->P, PQsub->n->value, A24p, C24);
    xLADDER(PQsub->Q, PQ->Q, PQsub->n->value, A24p, C24);
    xLADDER(PQsub->PQd, PQ->PQd, PQsub->n->value, A24p, C24);
    // Set the proper order of the subgroup torsion basis
    pprod_set(PQsub->n, n);
}
