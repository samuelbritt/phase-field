#ifndef LINALG_H
#define LINALG_H

/* Return x in b for Ux=b where U is upper triangular */
void back_sub(float** U, float* b, int dim):

/* Return x in b for Lx=b where L is lower triangular */
void forward_sub(float** L, float* b, int dim):

/* alters A so that lower triangular part of A becomes C  and the upper part
 * becomes C^T for the decomposition A = C C^T */
void cholesky_decomp(float** A, int dim);

/* solves Ax = b by cholesky decomposition. Returns x in b*/
void solve_cholesky(float **A, float *b, int dim);

#endif
