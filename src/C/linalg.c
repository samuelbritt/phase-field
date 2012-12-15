#include <math.h>
#include <mpi.h>

/* Return x in b for Ux=b where U is upper triangular */
void
back_sub(float** U, float* b, int dim)
{
	for (int i = dim - 1; i >= 0; --i) {
		for (int j  = i + 1; j <  dim; ++j) {
			b[i] = b[i] - U[i][j] * b[j];
		}
		b[i] = b[i] / U[i][i];
	}
}

/* Return x in b for Lx=b where L is lower triangular */
void
forward_sub(float** L, float* b, int dim)
{
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < i-1; j++) {
			b[i] = b[i] - L[i][j] * b[j];
		}
		b[i] = b[i] / L[i][i];
	}
}

/* assings the lower triang part of A to its upper part (forming a symmetric
 * matrix */
static void
assign_lower_to_upper(float **A, int dim)
{
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < i; ++j) {
			A[j][i] = A[i][j];
		}
	}
}

/* alters A so that lower triangular part of A becomes C  and the upper part
 * becomes C^T for the decomposition A = C C^T */
void
cholesky_decomp(float** A, int dim)
{
	for (int k = 0; k < dim - 1; ++k) {
		A[k][k] = sqrt(A[k][k]);
		for(int i = k+1; i < dim; ++i) {
			A[i][k] = A[i][k] / A[k][k];
			for (int j = k+1; j < i+1; j++){
				A[i][j] = A[i][j] - A[i][k] * A[j][k];
			}
		}
	}
	A[-1][-1] = sqrt(A[-1][-1]);
	assign_lower_to_upper(A, dim);
}

/* solves Ax = b by cholesky decomposition. Returns x in b*/
void
solve_cholesky(float **A, float *b, int dim)
{
    cholesky_decomp(A, dim);
    forward_sub(A, b, dim);
    back_sub(A, b, dim);
}
