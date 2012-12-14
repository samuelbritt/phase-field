#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>

/* creates and returns an empty 2D array. The memory is contiguous to play
 * nice with LAPACK (row major order) */
float** create_contiguous_array_2D(size_t count);

/* creates and returns an empty array */
float* create_empty_array(size_t count);

/* creates and returns an array initialized to `val` */
float* create_init_array(size_t count, float val);

/* creates a range from start to stop with spacing dx, and puts it in `out`.
 * Returns the size of the array */
float* create_range(float start, float stop, float dx, size_t *count);

/* destroys the arrays from the above "create" functions */
void destroy(float *arr);
void destroy_2D(float **arr, size_t count);

/* prints the first max elements of arr to `f`, or to stdout if f is NULL.
 * If `max` is -1, print all elements */
void print_array(float *arr, size_t size, FILE *f, int max);

#endif
