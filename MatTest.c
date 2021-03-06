/*
 * MatTest.c
 *
 *	Test matrix routines (by J. R. Valverde, 1988).
 *
 *	Compile with 'cc -o matt MatTest.c Mat.c -lm' and
 *	run 'matt -test_number'. See source for check of validity.
 *
 *	Jos� Ram�n Valverde Carrillo
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
    	    c[3][3];
    double *pa, *pb, *pc,
    	   *pmat1, *pmat2, *pmat3;
    double temp;

    if (argc != 2) {
	printf("Usage: mattest -# (where # is a number between 1 and 0)\n");
	exit(1);
    }

    switch (argv[1][1]) {
    case '1':
    pa = mat2alloc(3, 3);
     
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        {
          pa[(i * 3) + j] = 7.77;
         }

    for (i = 0; i < 9; i++)
      {
        printf("%1.2f\t", pa[i]);
        
      }
    /* this is an add-on for coherence with program termination */
    minit(c, 3, 3, 0.0);
    pmat1 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 1.1);
/*
	Must yield
7.77    7.77    7.77    7.77    7.77    7.77    7.77    7.77    7.77
1.10    1.10    1.10            0.00    0.00    0.00
1.10    1.10    1.10            0.00    0.00    0.00
1.10    1.10    1.10            0.00    0.00    0.00
*/
    break;

    case '2':
    /* 1. Using pointers.	*/
    pmat1 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 1.87);

    /* 2. Using matrices	*/
    minit(c, 3, 2, 1.1);
/*
	Should yield
1.87    1.87  1.87              1.10    1.10    1.10
1.87    1.87  1.87              1.10    1.10    1.10
1.87    1.87  1.87              0.00    0.00    0.00
*/
    break;

    case '3':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 1.1);
    minit(pmat2, 3, 3, 2.2);
    /* pmat1 = pmat2	*/
    massign(pmat1, pmat2, 3, 3);
     
    minit(c, 3, 3, 1.1);
    minit(b, 3, 3, 2.2);
    /* c = b */
    massign(c, b, 3, 3);
/*
	Must give
2.20    2.20  2.20              2.20    2.20    2.20
2.20    2.20  2.20              2.20    2.20    2.20
2.20    2.20  2.20              2.20    2.20    2.20
*/
    break;

    case '4':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 0.0);
    minit(pmat2, 3, 3, 0.0);

    /* first assign values to source matrices */
    temp = 0.0;
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        {
          a[i][j] = cos(temp);
          pmat2[i * 3 + j] = sin(temp);
          temp += 0.7;
        }
    /* Transpos'em to pmat1 and c */
    mtransp(pmat1, pmat2, 3, 3);
    mtransp(c, a, 3, 3);
    /* Show source matrices  */
    for (i = 0; i < 3; i++)
      {
        puts("");
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", pmat2[i * 3 + j]); 
        printf("\t");
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", a[i][j]);
      }
    puts("");
/*
	Should give
0.00    0.64  0.99              1.00    0.76    0.17
0.86    0.33  -0.35             -0.50   -0.94   -0.94
-0.87   -0.98 -0.63             -0.49   0.19    0.78

0.00    0.86  -0.87             1.00    -0.50   -0.49
0.64    0.33  -0.98             0.76    -0.94   0.19
0.99    -0.35 -0.63             0.17    -0.94   0.78
*/
    break;

    case '5':
    pmat1 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 7.81);
    mident(pmat1, 3);
    minit(c, 3, 3, 1.87);
    mident(c, 3);
/*
	Output should be
1.00    0.00  0.00              1.00    0.00    0.00
0.00    1.00  0.00              0.00    1.00    0.00
0.00    0.00  1.00              0.00    0.00    1.00
*/
    break;

    case '6':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    pmat3 = mat2alloc(3, 3);
    minit(pmat3, 3, 3, 1.45);
    minit(pmat2, 3, 3, 2.47);
    minit(pmat1, 3, 3, 1.00);
    msum(pmat1, pmat2, pmat3, 3, 3);
    
    minit(a, 3, 3, 1.45);
    minit(b, 3, 3, 2.47);
    minit(c, 3, 3, 1.00);
    pa = (double *) a;
    pb = (double *) b;
    pc = (double *) c;
    msum(c, a, b, 3, 3);
/*
	Correct output is
3.92    3.92  3.92              3.92    3.92    3.92
3.92    3.92  3.92              3.92    3.92    3.92
3.92    3.92  3.92              3.92    3.92    3.92
*/
    break;
    
    case '7':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    pmat3 = mat2alloc(3, 3);
    minit(pmat1, 3, 3, 1.00);
    minit(pmat2, 3, 3, 2.47);
    minit(pmat3, 3, 3, 1.45);
    msubst(pmat1, pmat2, pmat3, 3, 3);

    minit(a, 3, 3, 1.45);
    minit(b, 3, 3, 2.47);
    minit(c, 3, 3, 1.00);
    pa = (double *) a;
    pb = (double *) b;
    pc = (double *) c;
    msubst(c, b, a, 3, 3);
/*
	The result should be
1.02    1.02  1.02              1.02    1.02    1.02
1.02    1.02  1.02              1.02    1.02    1.02
1.02    1.02  1.02              1.02    1.02    1.02
*/
    break;
    
    case '8':
    /* Verification with a matrix and its inverse:
    The result should be the identity matrix	*/
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    pmat3 = mat2alloc(3, 3);
    minit(c, 3, 3, 0.0);
    pmat2[0] = 1.0; pmat2[1] = 2.0; pmat2[2] = 3.0;
    pmat2[3] = 1.0; pmat2[4] = 3.0; pmat2[5] = 3.0;
    pmat2[6] = 1.0; pmat2[7] = 2.0; pmat2[8] = 4.0;

    pmat1[0] = 6.0; pmat1[1] = -2.0; pmat1[2] = -3.0;
    pmat1[3] = -1.0; pmat1[4] = 1.0; pmat1[5] = 0.0;
    pmat1[6] = -1.0; pmat1[7] = 0.0; pmat1[8] = 1.0;
    
    mmult(pmat3, pmat2, pmat1, 3, 3, 3);
/*
    A more general test: pmat1[5][2] x pmat2[2][2] => pmat3[5][2]
    pmat1 = mat2alloc(5, 2);
    pmat2 = mat2alloc(2, 2);
    pmat3 = mat2alloc(5, 2);
    pmat1[0 + 0] = 1.0; pmat1[0 + 1] = 3.0;
    pmat1[2 + 0] = 4.0; pmat1[2 + 1] = 2.0;
    pmat1[4 + 0] = 1.0; pmat1[4 + 1] = 1.0;
    pmat1[6 + 0] = 6.0; pmat1[6 + 1] = 4.0;
    pmat1[8 + 0] = 3.0; pmat1[8 + 1] = 2.0;
    
    pmat2[0 + 0] = 1.0; pmat2[0 + 1] = 2.0;
    pmat2[2 + 0] = 3.0; pmat2[2 + 1] = 4.0;
    
    mmult(pmat3, pmat1, pmat2, 5, 2, 2);
*/
    for (i = 0; i < 3; i++)
      {
        puts("");
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", pmat3[i * 3 + j]); 
      }
    puts("");
    
    minit(a, 3, 3, 1.5);
    minit(b, 3, 3, 2.0);
    minit(c, 3, 3, 0.0);
    pa = (double *) a;
    pb = (double *) b;
    pc = (double *) c;
    mmult(pc, pa, pb, 3, 3, 3);
/*
	Results are
1.00    0.00  0.00
0.00    1.00  0.00
0.00    0.00  1.00

6.00    -2.00 -3.00             9.00    9.00    9.00
-1.00   1.00  0.00              9.00    9.00    9.00
-1.00   0.00  1.00              9.00    9.00    9.00

*/
    break;

    case '9':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    minit(pmat2, 3, 3, 1.5e-2);
    mscprod(pmat1, pmat2, 3, 3, 3.0);

    minit(a, 3, 3, 1.5);
    minit(c, 3, 3, 0.0);
    pa = (double *) a;
    pc = (double *) c;
    mscprod(c, a, 3, 3, 2.0);
/*
	Gives
0.04    0.04  0.04              3.00    3.00    3.00
0.04    0.04  0.04              3.00    3.00    3.00
0.04    0.04  0.04              3.00    3.00    3.00
*/
    break;

    case '0':
    pmat1 = mat2alloc(3, 3);
    pmat2 = mat2alloc(3, 3);
    pmat3 = mat2alloc(3, 3);
    minit(c, 3, 3, 0.0);
    pmat2[0] = 1.0; pmat2[1] = 2.0; pmat2[2] = 3.0;
    pmat2[3] = 1.0; pmat2[4] = 3.0; pmat2[5] = 3.0;
    pmat2[6] = 1.0; pmat2[7] = 2.0; pmat2[8] = 4.0;
    
    /* see its inverse upwards, in DEBUG == 8	*/
    
    minv(pmat1, pmat2, 3);  /* OK */
    mmult(pmat3, pmat2, pmat1, 3, 3, 3);
    
    for (i = 0; i < 3; i++)
      { /* print pmat3 */
        puts("");
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", pmat3[i * 3 + j]);
      }
    puts("\nthis was pmat3");
/*
	Must give
1.00    0.00  0.00
0.00    1.00  0.00
0.00    0.00  1.00
this was pmat3

6.00    -2.00 -3.00             0.00    0.00    0.00
-1.00   1.00  0.00              0.00    0.00    0.00
-1.00   0.00  1.00              0.00    0.00    0.00
*/
    break;
    default:
	printf("Usage: mattest -# (where # is a number between 1 and 0)\n");
	exit(1);
    break;
    }

    /* show pmat1 and c */
    for (i = 0; i < 3; i++)
      {
        puts("");
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", pmat1[i * 3 + j]);
        printf("\t"); 
        for (j = 0; j < 3; j++)
          printf("%1.2f\t", c[i][j]);
      }   
    puts("\n");
    
    /* to end	*/
    puts("\nPress [c]");
    while ((ch = getchar()) != 'c');

  }



