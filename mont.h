#ifndef MONT_H
#define MONT_H

#include "fp2.h"

// Point in projective coordinates X-Z
struct point_xz { 
    fp2_t X; 
    fp2_t Z; 
};

typedef struct point_xz* point_t;

void point_set(point_t R, const point_t P);

void point_init(point_t *P);

void point_clear(point_t *P);

// Calculate the double of the point
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24);

// Calculate R: multiple [2^e]P of point P
void xDBLe(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24, const int e);

// Calculate the codomain (A24+, C24) of the 2-degree isogeny from given kernel K of order 2
void isog2_codomain(const point_t K, fp2_t A24_plus, fp2_t C24);


#endif