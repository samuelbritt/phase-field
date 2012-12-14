#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

/* creates and returns an empty 2D array. The memory is contiguous to play
 * nice with LAPACK (row major order) */
float **
create_contiguous_array_2D(size_t count)
{
	float** arr = malloc(count * sizeof(*arr));
	float* data = malloc(count * count * sizeof(*data));
	for (int i = 0; i < count; ++i) {
		arr[i] = &data[i * count];
	}
	return arr;
}

/* creates and returns an empty array */
float *
create_empty_array(size_t count)
{
	return malloc(count * sizeof(float));
}

/* creates and returns an array initialized to `val` */
float *
create_init_array(size_t count, float val)
{
	float *arr = create_empty_array(count);
	for (int i = 0; i < count; ++i) {
		arr[i] = val;
	}
	return arr;
}

/* creates a range from start to stop with spacing dx, and returns it. Returns
 * the size of the array in the `size` parameter */
float *
create_range(float start, float stop, float dx, size_t *size)
{
	*size = (int) (stop - start) / dx;
	float *out = create_empty_array(*size);

	float x = start;
	for (int i = 0; i < *size; ++i) {
		out[i] = x;
		x += dx;
	}
	return out;
}

void
destroy(float *arr)
{
	free(arr);
}

void
destroy_2D(float **arr, size_t count)
{
	for (int i = 0; i < count; ++i) {
		free(arr[i]);
	}
	free(arr);
}

void
print_array(float *arr, size_t size, FILE *f, int max)
{
	f = f ? f : stdout;

	int target = max < 0 ? size : max;
	fprintf(f, "[ ");
	for (int i = 0; i < target; ++i) {
		fprintf(f, "%6.2g, ", arr[i]);
	}
	fprintf(f, "]\n");
}
