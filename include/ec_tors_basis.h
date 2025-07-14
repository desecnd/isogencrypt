#pragma once

#include "ec_point_xz.h"
#include "pprod.h"

// Torsion Basis for elliptic curve 
struct tors_basis {
    point_t P, Q, PQd;
    pprod_t n;
};

void tors_basis_init(struct tors_basis *tb);
void tors_basis_clear(struct tors_basis *tb);

/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q) of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, pprod_t n, const struct tors_basis *PQ, const fp2_t A24p, const fp2_t C24);