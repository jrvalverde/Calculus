/*
* * * * * * * * * * * * * * * * * * * * * * * * *
*						*
*		    MATRIZ.H			*
*   Include file with the definitions needed	*
*   to use the functions of MATRIX.C		*
*						*
* * * * * * * * * * * * * * * * * * * * * * * * *
*/

#define M_SINGULAR	(-2)

#define matfree(x)	free(x)

typedef double *matriz;

extern double *matalloc();

extern double *mat2alloc();

extern minit();

extern massign();

extern mtransp();

extern mident();

extern msum();

extern msubst();

extern mmult();

extern mscprod();

extern int minv();

/*
* * * * * * * * * * * * * * * * * * * * * * *
*					    *
*		FIN DE MATRIZ.H		    *
*					    *
* * * * * * * * * * * * * * * * * * * * * * *
*/
