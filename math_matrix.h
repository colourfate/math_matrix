#ifndef __MATRIX_H
#define __MATRIX_H
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)
#define LIMIT(m,a,b) MIN(MAX(m,MIN(a,b)),MAX(a,b))
#include <stdint.h>
#include <stdlib.h>

typedef uint8_t u8;
typedef uint16_t u16;

// 下面是基于结构体的矩阵运算
struct matrix_t{
	float *m;
	u8 row;
	u8 column;
};
#define MATRIX_INIT(a,b,c,d) 		\
			a.m=b;					\
			a.row=c;				\
			a.column=d
int8_t matrix_t_T(struct matrix_t *, const struct matrix_t *);
void matrix_t_show(char *name, const struct matrix_t *M);
int8_t matrix_t_plus(struct matrix_t *, const struct matrix_t *, 
			const struct matrix_t *, int8_t mode);
int8_t matrix_t_mul(struct matrix_t *, const struct matrix_t *, 
			const struct matrix_t *, int8_t mode);
int8_t matrix_t_inv(struct matrix_t *, const struct matrix_t *);
int8_t matrix_t_copy(struct matrix_t *, const struct matrix_t *);
int8_t matrix_t_eye(struct matrix_t *);
int8_t matrix_t_k(struct matrix_t *, float k, const struct matrix_t *);
int8_t matrix_t_concat(struct matrix_t *, const struct matrix_t *,
				const struct matrix_t *, u8 mode);
int8_t matrix_t_transport(struct matrix_t *, const struct matrix_t *, 
				u8 x1, u8 x2, u8 y1, u8 y2);
int8_t matrix_t_conv(struct matrix_t *, const struct matrix_t *,
				const struct matrix_t *);
void matrix_t_zero(struct matrix_t *A);
void matrix_t_malloc(struct matrix_t *A, u8 row, u8 column);
void matrix_t_free(struct matrix_t *A);
int8_t matrix_t_vector_mul(struct matrix_t *A, struct matrix_t *B, u8 n1, u8 n2,
					struct matrix_t *C, u8 n3, u8 n4, int8_t mode);
void matrix_t_set(struct matrix_t *A, float val);
int8_t matrix_t_vector_transport(struct matrix_t *A, struct matrix_t *B, u8 n1, u8 n2);
int8_t matrix_t_vector_T(struct matrix_t *A);
#endif
