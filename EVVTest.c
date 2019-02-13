/*
 * EVVTest.c
 *
 *	Test Eigenvalue/Eigenvector routines (by J. R. Valverde, 2008).
 *
 *	Compile with 'cc -o evvt EVVTest.c EVV.c -lm Mat.c' and
 *	run 'evvt'. See source for check of validity.
 *
 *	José Ramón Valverde Carrillo
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "portable.h"
#include "Mat.h"


main(argc, argv)
int argc;
char *argv[];
  {
    char ch;
    int i, j, k, l, n;
    double  a[3][3],
    	    b[3][3],
    	    c[4][4],
	    d[4][4],
	    e[4][4];
    double *pa, *pb, *pc,
    	   *pmat1, *pmat2, *pmat3;
    double temp;
    double max;


    /* The problem matrix */
    a[0][0] = 2.; a[0][1] =-1.; a[0][2] = 1.;
    a[1][0] =-1.; a[1][1] = 2.; a[1][2] = 1.;
    a[2][0] = 1.; a[2][1] =-1.; a[2][2] = 2.;
    
    pa = calloc(3, sizeof(double)); /* coefficient vector */
    pa[0] = 1.0; pa[1] = 0.0; pa[2] = 0.0;

    printf("Original matrix and eigenvector guess\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf("\n");
    max = 0.0;
    i = EVVLeverrierFaddeev(3, a, pa, b);
    
    printf("LF: Final matrix, coefficients and inverse\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]", pa[i]);
	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    minv(b, a, 3);
    printf("Inverse matrix\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    i = EVVSuccApprox(3, a, pa, &max, 0.00001, 100);
    printf("SA: Final matrix, eigenvector and maximum eigenvalue\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf(" Max = %g\n\n", max); 
    free(pa);
    
    /* The problem matrix */
    c[0][0] = 8.; c[0][1] =-1.; c[0][2] = 3.; c[0][3] =-1.;
    c[1][0] =-1.; c[1][1] = 6.; c[1][2] = 2.; c[1][3] = 0.;
    c[2][0] = 3.; c[2][1] = 2.; c[2][2] = 9.; c[2][3] = 1.;
    c[3][0] =-1.; c[3][1] = 0.; c[3][2] = 1.; c[3][3] = 7.;
    
    pc = calloc(4, sizeof(double)); /* coefficient vector */
    pc[0] = 1.0; pc[1] = 0.0; pc[2] = 0.0; pc[3] = 0.0;

    printf("Original matrix and eigenvector guess\n");
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++) {
	    printf("\t%+3.3g", c[i][j]);
	}
	printf("\t[%+3.3g]\n", pc[i]);
    }
    printf("\n");
    max = 0.0;
    i = EVVLeverrierFaddeev(4, c, pc, d);
    
    printf("LF: Final matrix, coefficients and inverse\n");
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++) {
	    printf("\t%+3.3g", c[i][j]);
	}
	printf("\t[%+3.3g]", pc[i]);
	for (j = 0; j < 4; j++)
	    printf("\t%+3.3g", d[i][j]);
	printf("\n");
    }
    minv(d, c, 4);
    printf("Inverse matrix\n");
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++)
	    printf("\t%+3.3g", d[i][j]);
	printf("\n");
    }
    mmult(e, c, d, 4, 4, 4);
    printf("matrix by its inverse\n");
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++)
	    printf("\t%+3.3f", e[i][j]);
	printf("\n");
    }
    EVVSuccApprox(4, c, pc, &max, 0.00001, 100);
    printf("SA: Final matrix, eigenvector and maximum eigenvalue\n");
    for (i = 0; i < 4; i++) {
    	for (j = 0; j < 4; j++) {
	    printf("\t%+3.3f", c[i][j]);
	}
	printf("\t[%+3.3f]\n", pc[i]);
    }
    printf(" Max = %f\n\n", max); 
  
    free(pc);
    
    /* The problem matrix */
    a[0][0] = 4.; a[0][1] = 1.; a[0][2] =-1.;
    a[1][0] = 2.; a[1][1] = 3.; a[1][2] =-1.;
    a[2][0] =-2.; a[2][1] = 1.; a[2][2] = 5.;
    
    pa = calloc(3, sizeof(double)); /* coefficient vector */
    pa[0] = 1.0; pa[1] = 0.0; pa[2] = 0.0;

    printf("Original matrix and eigenvector guess\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf("\n");
    max = 0.0;
    i = EVVLeverrierFaddeev(3, a, pa, b);
    
    printf("LF: Final matrix, coefficients and inverse\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]", pa[i]);
	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    minv(b, a, 3);
    printf("Inverse matrix\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    i = EVVSuccApprox(3, a, pa, &max, 0.00001, 100);
    
    printf("SA: Final matrix, eigenvector and maximum eigenvalue\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf(" Max = %g\n\n", max); 
    free(pa);
    

    /* The problem matrix */
    a[0][0] =-2.; a[0][1] = 1.; a[0][2] =-1.;
    a[1][0] = 2.; a[1][1] =-3.; a[1][2] =-1.;
    a[2][0] =-2.; a[2][1] = 1.; a[2][2] =-1.;
    
    pa = calloc(3, sizeof(double)); /* coefficient vector */
    pa[0] = 1.0; pa[1] = 0.0; pa[2] = 0.0;

    printf("Original matrix and eigenvector guess\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf("\n");
    max = 0.0;
    i = EVVLeverrierFaddeev(3, a, pa, b);
    
    printf("LF: Final matrix, coefficients and inverse\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]", pa[i]);
	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    minv(b, a, 3);
    printf("Inverse matrix\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++)
	    printf("\t%+3.3g", b[i][j]);
	printf("\n");
    }
    
    i = EVVSuccApprox(3, a, pa, &max, 0.00001, 100);
    
    printf("SA: Final matrix, eigenvector and maximum eigenvalue\n");
    for (i = 0; i < 3; i++) {
    	for (j = 0; j < 3; j++) {
	    printf("\t%+3.3g", a[i][j]);
	}
	printf("\t[%+3.3g]\n", pa[i]);
    }
    printf(" Max = %g\n\n", max); 
    free(pa);

/*
    Should yield
    
Original matrix and eigenvector guess
	 +2	 -1	 +1	[ +1]
	 -1	 +2	 +1	[ +0]
	 +1	 -1	 +2	[ +0]

LF: Final matrix, coefficients and inverse
	 +2	 -1	 +1	[ +6]	+0.833	+0.167	-0.5
	 -1	 +2	 +1	[-11]	+0.5	+0.5	-0.5
	 +1	 -1	 +2	[ +6]	-0.167	+0.167	+0.5
Inverse matrix
	+0.833	+0.167	-0.5
	+0.5	+0.5	-0.5
	-0.167	+0.167	+0.5
SA: Final matrix, eigenvector and maximum eigenvalue
	 +2	 -1	 +1	[ +1]
	 -1	 +2	 +1	[-1.14e-05]
	 +1	 -1	 +2	[ +1]
 Max = 3.00002

Original matrix and eigenvector guess
	 +8	 -1	 +3	 -1	[ +1]
	 -1	 +6	 +2	 +0	[ +0]
	 +3	 +2	 +9	 +1	[ +0]
	 -1	 +0	 +1	 +7	[ +0]

LF: Final matrix, coefficients and inverse
	 +8	 -1	 +3	 -1	[+30]	    	+0.161	+0.0496	-0.0683	+0.0327
	 -1	 +6	 +2	 +0	[-319]	    	+0.0496	+0.196	-0.0617	+0.0159
	 +3	 +2	 +9	 +1	[+1.41e+03]	-0.0683	-0.0617	+0.151	-0.0313
	 -1	 +0	 +1	 +7	[-2.14e+03]	+0.0327	+0.0159	-0.0313	+0.152
Inverse matrix
	+0.161	+0.0496	-0.0683	+0.0327
	+0.0496	+0.196	-0.0617	+0.0159
	-0.0683	-0.0617	+0.151	-0.0313
	+0.0327	+0.0159	-0.0313	+0.152
matrix by its inverse
	+1.000	+0.000	+0.000	-0.000
	+0.000	+1.000	+0.000	+0.000
	-0.000	-0.000	+1.000	-0.000
	-0.000	-0.000	+0.000	+1.000
SA: Final matrix, eigenvector and maximum eigenvalue
	+8.000	-1.000	+3.000	-1.000	[+0.735]
	-1.000	+6.000	+2.000	+0.000	[+0.222]
	+3.000	+2.000	+9.000	+1.000	[+1.000]
	-1.000	+0.000	+1.000	+7.000	[+0.056]
 Max = 11.704324

Original matrix and eigenvector guess
	 +4	 +1	 -1	[ +1]
	 +2	 +3	 -1	[ +0]
	 -2	 +1	 +5	[ +0]

LF: Final matrix, coefficients and inverse
	 +4	 +1	 -1	[+12]	+0.333	-0.125	+0.0417
	 +2	 +3	 -1	[-44]	-0.167	+0.375	+0.0417
	 -2	 +1	 +5	[+48]	+0.167	-0.125	+0.208
Inverse matrix
	+0.333	-0.125	+0.0417
	-0.167	+0.375	+0.0417
	+0.167	-0.125	+0.208
SA: Final matrix, eigenvector and maximum eigenvalue
	 +4	 +1	 -1	[ -1]
	 +2	 +3	 -1	[ -1]
	 -2	 +1	 +5	[ +1]
 Max = 5.99998

Original matrix and eigenvector guess
	 -2	 +1	 -1	[ +1]
	 +2	 -3	 -1	[ +0]
	 -2	 +1	 -1	[ +0]

LF: Final matrix, coefficients and inverse
	 -2	 +1	 -1	[ -6]	+inf	+nan	-inf
	 +2	 -3	 -1	[ -8]	+inf	+nan	-inf
	 -2	 +1	 -1	[ +0]	-inf	+nan	+inf
Inverse matrix
	-0.75	-0.25	 -0
	-0.5	-0.5	 -0
	 -1	 +0	 +1
SA: Final matrix, eigenvector and maximum eigenvalue
	 -2	 +1	 -1	[ -1]
	 +2	 -3	 -1	[ +1]
	 -2	 +1	 -1	[ -1]
 Max = -3.99999


*/
    

}

