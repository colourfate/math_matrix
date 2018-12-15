#include <stdio.h>
#include "math_matrix.h"

int main(void)
{
	struct matrix_t	m1, m2, A;
	matrix_t_malloc(&m1, 2, 3);
	matrix_t_malloc(&m2, 2, 3);
	matrix_t_malloc(&A, 2, 3);

	matrix_t_set(&m1, 1);
	matrix_t_set(&m2, 2);
	matrix_t_plus(&A, &m1, &m2, 1);

	matrix_t_show("A", &A);

	return 0;
}

