/*
 * 1-D Phase field solidification simulation of a pure, isotropic substance.
 */

#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <clapack.h>
#include <cblas.h>
#include <mpi.h>

#define PHI_SOLID =  0
#define PHI_LIQUID = 1

/* For Nickel, from Bragard, Karma, and Lee et al. */

/* Latent heat of melting */
static float L = 2.311e3;	// J cm^-3

/* heat capacity */
static float cp = 5.313;	// J cm^-3 K^-1

/* Melting Temp */
static float Tm = 1726.;	// K

/* capillary (interfacial) energy */
/* static float gamma0 = 0.326;	// J m^-2 */

/* kinetic coefficient; controls rate at which atoms attach to ijk surface */
/* static float mu_100 = 0.52;	// m s^-1 */
/* static float mu_110 = 0.40;	// m s^-1 */

/* Diffusivities (liquid, solid) */
static float D = 0.1;		// cm^2 s^-1

/* Interfacial mobility */
static float M = -1e6;

/* capillary anisotropy; set to 0 for isotropic interfacial energy */
static float eps_c = 0.018;

/* isotropic energy gradient coefficient \epsilon^2 */
static float eps2 = 1.25e-13;

/* Additional activation bump */
static float Q = 1;

/* for isothermal */
static float u_iso = -0.3;

/* Normalized temperature, called u(T) in Bragard et al. */
static float
u(float T)
{
	float T_norm = L / cp;	// K
	return (T - Tm) / T_norm;
}

/* double well function, called g(phi) in my notes */
static float
g(float phi)
{
	return  phi * phi * (1 - phi) * (1 - phi);
}

/* interpolation function, called p(phi) in my notes */
static float
p(float phi)
{
	float phi2 = phi * phi;
	return phi2 * phi * (6 * phi2 - 15 * phi + 10);
}

/* second term on right hand side of dphi/dt, called h(phi, u) */
/* description */
static float
h(float phi, float u)
{
	float kappa = L * L / (cp * Tm);
	return 2 * M * phi * (15 * kappa * u * phi * (phi - 1) * (phi - 1)
			      - Q * (1 - phi) * (1 - 2 * phi));
}

/* returns (via paramater reference) the boundary conditions */
static void
set_boundary(float *phi_0, float *phi_n) {
	*phi_0 = 0;
	*phi_n = 1;
}

/* solves Ax=b for nxn A */
static void
solve(float **A, float *x, float *b, size_t n) {

	/* int clapack_sgesv(const enum CBLAS_ORDER Order, const int N,
	 *                   const int NRHS, float *A, const int lda, int *ipiv,
	 *                   float *B, const int ldb); */
	int *ipiv = malloc(n * sizeof(*ipiv));
	int info = clapack_sgesv(CblasRowMajor,
				 n,		// order of A
				 1,		// num RHS
				 A[0],		// A
				 n, 		// leading dim A
				 ipiv,		// pivot indices
				 b,		// b
				 n);		// leading dim of b
	if (info != 0) {
		fprintf(stderr, "Error! LAPACK returned %d\n", info);
		exit(1);
	}

	for (int i = 0; i < n; ++i) {
		x[i] = b[i];
	}

}

/* Initializes the coefficient matrix A */
static void
init_A(float **A, size_t mat_size, float mu, float eps_2)
{
	for (int i = 0; i < mat_size; ++i) {
		for (int j = 0; j < mat_size; ++j) {
			A[i][j] = 0;
		}
	}

	float mu_eps2 = mu * eps2;

	/* the full tridiagonal portion */
	for (int i = 1; i < mat_size - 1; ++i) {
		A[i][i] = 1 + 2 * mu_eps2;
		A[i][i - 1] = - mu_eps2;
		A[i][i + 1] = - mu_eps2;
	}

	/*  the last bits */
	A[0][0] = 1 + 2 * mu_eps2;
	A[0][1] = - mu_eps2;
	A[mat_size - 1][mat_size - 2] = - mu_eps2;
	A[mat_size - 1][mat_size - 1] = 1 + 2 * mu_eps2;
}

/* main integration function. Updates phi(x, t) to phi(x, t+dt) */
static void
update_phi(float *x, float *phi, size_t size, float t, float mu, float eps2,
	   float u) {
	/* Our matrices here are of size `size - 2`, because the first and
	 * last values of phi are taken care of by the boundary conditions */
	/* size_t mat_size = size - 2; */

	float mu_eps2 = mu * eps2;

	size_t mat_size = size - 2;
	float** A = create_contiguous_array_2D(mat_size);
	init_A(A, mat_size, mu, eps2);

	float *b = create_init_array(mat_size, 0);

	/* at the boundary */
	float phi_0, phi_n;
	set_boundary(&phi_0, &phi_n);
	b[0] = mu_eps2 * phi_0;
	b[mat_size - 1] = mu_eps2 * phi_n;

	for (int i = 0; i < mat_size; ++i) {
		b[i] = phi[i+1] + h(phi[i+1], u);
	}

	solve(A, &phi[1], b, mat_size);
	phi[0] = phi_0;
	phi[size-1] = phi_n;
}


/* initializes the order parameter phi(x) */
static void
init_phi_1D(float* phi, float* x, size_t size)
{
	for (int i = 0; i < size; ++i) {
		phi[i] = 0.5 * (1 + tanh(3 * x[i]));
	}
}

/* description */
static void
output_data(float *x, float *phi, size_t size, float t, char *dir, char *tag)
{
	char fname[1024];
	snprintf(fname, 1024, "%s/%s_t%05.2g.dat", dir, tag, t);
	FILE *f = fopen(fname, "w");
	fprintf(f, "# Output data, time:\n");
	fprintf(f, "%g\n", t);
	for (int i = 0; i < size; ++i) {
		fprintf(f, "%g\t%g\n", x[i], phi[i]);
	}
	fclose(f);
}

static void mkdir(char *dir)
{
	char cmd[1024];
	sprintf(cmd, "mkdir -p %s", dir);  // not exactly portable... but easy
	system(cmd);
}

int main(int argc, char const *argv[])
{
	float x0 = -1.;
	float x_max = 10.;
	float dx = 0.05;
	size_t x_size;
	float* x_all = create_range(x0, x_max, dx, &x_size);

	float t0 = 0.;
	float t_max = 10.;
	float dt = 0.02;

	int plot_count = 10;
	int plot_every = (int) t_max / dt / plot_count;
	char *tag = "pure_1D";
	char *dir = "pure_1D";
	mkdir(dir);

	float* phi = create_empty_array(x_size);
	init_phi_1D(phi, x_all, x_size);

	float mu = dt * M / (dx * dx);

	float t = t0;
	int i = 0;
	while (t < t_max) {
		if (i % plot_every == 0) {
			output_data(x_all, phi, x_size, t, dir, tag);
		}
		i++;
		t += dt;
		update_phi(x_all, phi, x_size, t, mu, eps2, u_iso);
	}
}
