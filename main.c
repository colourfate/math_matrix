#include <stdio.h>
#include "math_matrix.h"

int main(void)
{
	struct matrix_t	m1, m2, A;
	float a[2][3] = {{1, 2, 3}, {4, 5, 6}};
	// 这个在栈中，不用free
	MATRIX_INIT(m2, a, 2, 3);
	matrix_t_malloc(&m1, 2, 3);
	matrix_t_malloc(&A, 2, 3);

	matrix_t_set(&m1, 1);
	matrix_t_plus(&A, &m1, &m2, 1);

	matrix_t_show("A", &A);

	matrix_t_free(&m1);
	matrix_t_free(&A);

	return 0;
}

