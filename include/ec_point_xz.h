#pragma once

#include "fp2.h"

/* 
 * @brief Transform (A : C) = (a : 1) into (A' : C') = ((a+ 2)/4 : 1) = (a + 2 : 4) suitable for xDBL.
 * Function is argument safe for A24p = A and C24 = C.
 */
void A24p_from_A(fp2_t A24p, fp2_t C24, const fp2_t A, const fp2_t C);

/* 
 * @brief Transform (A' : C') = (a + 2 : 4) into (A : C) = (a : 1). 
 * Function is argument safe for A = A24p and C = C24
 */
void A_from_A24p(fp2_t A, fp2_t C, const fp2_t A24p, const fp2_t C24);


// Point in projective coordinates X-Z
struct point_xz { 
    fp2_t X; 
    fp2_t Z; 
};

typedef struct point_xz* point_t;

void point_set(point_t R, const point_t P);

void point_init(point_t *P);

void point_clear(point_t *P);

// Set: X(P) = x, and Z(P) = 1
void point_set_str_x(point_t P, const char *x); 

// Set: X(P) = x, and Z(P) = 1
void point_set_fp2_x(point_t P, fp2_t x); 

/*
 * @brief Print the normalized (Z = 1) x-coordinate of the point P = (X : Z)
 */
void point_printx_normalized(point_t P, const char* name);

/*
 * @brief Print the x-coordinate of the point P = (X : Z) with assertion that Z == 1
 */
void point_printx(point_t P, const char* name);

// Normalize: P = (X : Z) -> (X' : 1) with X' = X/Z
void point_normalize_coords(point_t P);

/*
 * @brief Return 1 if point is normalized (P->Z == 1), 0 if P->Z != 1
 */
int point_is_normalized(point_t P);

/* 
 * @brief Compare the x-only coordinate of point by normalization of P and initalization of the const char* str.
 */
int point_equal_str_x(point_t P, const char* str);
