all:
	gcc main.c math_matrix.c -o matrix

git:
	git add main.c math_matrix* Makefile
	git status

