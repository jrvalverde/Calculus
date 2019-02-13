/*
 * MatTest.c
 *
 *	Test matrix routines (by J. R. Valverde, 1988).
 *
 *	Compile with 'cc -o matt MatTest.c Mat.c -lm' and
 *	run 'matt -test_number'. See source for check of validity.
 *
 *	José Ramón Valverde Carrillo
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "portable.h"
#include "matrix.h"


main(argc, argv)
int argc;
char *argv[];
  {
    char ch;
    int i, j, k, l, n;
    real  **a,
    	  **b,
    	  **c;
    real *bv;
    matrix A, B, C,
    	    AT, BT;
    real temp;
    int status;
    
    printf("\n\nStep 1: matrix allocation and access\n------------------------------------\n");

    status = mat_alloc(&A, 3, 3); 
    if (status != SUCCESS) {
    	printf("Error: can't allocate A\n");
	exit(1);
    }

    /* direct access to the values list */     
    a = A->values;
    for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
        {
          a[i][j] = 7.77;
         }

    for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
      {
        printf("% 1.2f\t", a[i][j]);
        
      }
    printf("\n\nThat should have been\n%s\n",
"7.77    7.77    7.77    7.77    7.77    7.77    7.77    7.77    7.77\n\n");
    mat_free(A);
    
    /* Access as a unidimensional array */
    mat_alloc(&B, 3, 3);
    mat_init(B, 0.0);
    bv = B->values[1];
    for (i = 1; i <= 3; i++) {
    	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", b[i * j]);
	printf("\n");
    }
    printf("\n");
    mat_free(B);

    /* Access using mat_element() */
    mat_alloc(&C, 3, 3);
    mat_init(C, 1.1);
    for(i = 1; i <= 3; i++) {
    	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", mat_element(C, i, j));
	printf("\n");
    }
    printf("\n\nMust yield\n%s\n%s\n%s\n\n%s\n%s\n%s\n\n",
"0.00	0.00	0.00 ",
"0.00	0.00	0.00 ",  
"0.00	0.00	0.00 ",  

"1.10    1.10    1.10 ",        
"1.10    1.10    1.10 ",        
"1.10    1.10    1.10 ");
    mat_free(C);
    
    printf("\n\nStep 2: initialization\n----------------------\n\n");
        
    mat_alloc(&A, 3, 3);
    mat_init(A, 1.87);
    
    /* using matrix information routines */
    for (i = 1; i <= mat_rows(A); i++) {
    	real *row;
	row = mat_row(A, i);
	for (j = 1; j <= mat_cols(A); j++) {
	    printf("% 1.2f\t", mat_element(A, i, j));
	}
    	printf("\n");
    }
    printf("\n");
    printf("Should have been:\n\
1.87	1.87	1.87\n\
1.87	1.87	1.87\n\
1.87	1.87	1.87\n\n");

/*
    mat_alloc(&C, 3, 3);
    mat_fill(c, 1, 3, 1, 2, 1.1);
*/

    mat_alloc(&A, 3, 3);
    mat_alloc(&B, 3, 3);
    mat_init(A, 1.1);
    mat_init(B, 2.2);
    /* A = B	*/
    mat_assign(A, B);
     
    mat_alloc(&C, 3, 3);
    mat_init(C, 1.1);
    /* B = C */
    mat_assign(B, C);
    
    a = mat_values(A);
    b = mat_values(B);
    for (i = 1; i <= 3; i++) {
    	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", a[i][j]);
	printf("\t");
	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", b[i][j]);
    	printf("\n");
    }
    printf("\n\n\
	Must give\n\
2.20    2.20  2.20              1.10    1.10    1.10\n\
2.20    2.20  2.20              1.10    1.10    1.10\n\
2.20    2.20  2.20              1.10    1.10    1.10\n\n");
    mat_free(A);
    mat_free(B);
    mat_free(C);

    printf("\n\nStep 3: transpose\n-----------------\n\n");
    
    mat_alloc(&A, 3, 3);
    mat_alloc(&B, 3, 3);
    mat_init(A, 0.0);
    mat_init(B, 0.0);

    /* first assign values to source matrices */
    temp = 0.0;
    a = mat_values(A);
    for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++)
        {
          a[i][j] = cos(temp);
          mat_set(B, i, j, sin(temp));
          temp += 0.7;
        }
	
    /* Show source matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", a[i][j]); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(B, i, j));
      }
    puts("");
    /* Transpos'em to AT and BT */
    mat_alloc(&AT, 3, 3);
    mat_alloc(&BT, 3, 3);
    mat_init(AT, 0.0);
    mat_init(BT, 0.0);
    
    mat_transpose(AT, A);
    mat_transpose(BT, B);
    /* Show transpose matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(AT, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(BT, i, j));
      }
    puts("");

    printf("\nShould give\n\
 1.00	 0.76	 0.17 		   0.00    0.64   0.99 \n\
-0.50	-0.94	-0.94		   0.86    0.33  -0.35\n\
-0.49	 0.19	 0.78 		  -0.87   -0.98  -0.63\n\
\n\
 1.00	-0.50	-0.49		   0.00    0.86  -0.87\n\
 0.76	-0.94	 0.19 		   0.64    0.33  -0.98\n\
 0.17	-0.94	 0.78 		   0.99   -0.35  -0.63\n\n");
    mat_free(A); mat_free(AT);
    mat_free(B); mat_free(BT);


    printf("\n\nStep 4: Identity matrix\n-----------------------\n");
    mat_alloc(&A, 3, 3);
    mat_init(A, 7.81);
    mat_identity(A);
    mat_alloc(&C, 3, 3);
    mat_init(C, 1.87);
    mat_identity(C);
    /* Show matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(A, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(B, i, j));
      }
    puts("");

    printf("\n\nOutput should be\n\
 1.00    0.00  0.00              1.00    0.00    0.00\n\
 0.00    1.00  0.00              0.00    1.00    0.00\n\
 0.00    0.00  1.00              0.00    0.00    1.00\n\n");
    mat_free(A);
    mat_free(C);


    printf("\n\nStep 5: Arithmetic\n------------------\n");
    mat_alloc(&A, 3, 3);
    mat_alloc(&B, 3, 3);
    mat_alloc(&C, 3, 3);

    mat_init(A, 1.45);
    mat_init(B, 2.47);
    mat_init(C, 1.00);

    mat_sum(C, A, B);
    for(i = 1; i <= 3; i++) {
    	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", mat_element(C, i, j));
	printf("\n");
    }
    printf("\n\nCorrect output is\n\
 3.92    3.92    3.92\n\
 3.92    3.92    3.92\n\
 3.92    3.92    3.92\n\n");

    mat_init(A, 1.00);
    mat_init(B, 2.47);
    mat_init(C, 1.45);
    mat_substract(A, B, C);
    for(i = 1; i <= 3; i++) {
    	for (j = 1; j <= 3; j++)
	    printf("% 1.2f\t", mat_element(A, i, j));
	printf("\n");
    }
    printf("\n\nCorrect output is\n\
 1.02    1.02    1.02\n\
 1.02    1.02    1.02\n\
 1.02    1.02    1.02\n\n");


    /* Verification with a matrix and its inverse:
    The result should be the identity matrix	*/
    mat_init(C, 0.0);
    
    a = mat_values(A);
    a[1][1] = 1.0; a[1][2] = 2.0; a[1][3] = 3.0;
    a[2][1] = 1.0; a[2][2] = 3.0; a[2][3] = 3.0;
    a[3][1] = 1.0; a[3][2] = 2.0; a[3][3] = 4.0;

    b = mat_values(B);
    b[1][1] = 6.0;  b[1][2] = -2.0; b[1][3] = -3.0;
    b[2][1] = -1.0; b[2][2] = 1.0;  b[2][3] = 0.0;
    b[3][1] = -1.0; b[3][2] = 0.0;  b[3][3] = 1.0;
    
    mat_multiply(C, A, B);
    /* Show matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(A, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(B, i, j));
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(C, i, j));
      }
    puts("");
    printf("\nThat must have been\n\
 1.00    2.00    3.00            6.00   -2.00   -3.00            1.00    0.00   0.00\n\
 1.00    3.00    3.00           -1.00    1.00    0.00            0.00    1.00   0.00\n\
 1.00    2.00    4.00           -1.00    0.00    1.00            0.00    0.00   1.00\n\n");
    mat_free(A); mat_free(B); mat_free(C);


/*
    A more general test: A[5][2] x B[2][2] => C[5][2]
*/
    mat_alloc(&A, 5, 2);
    mat_alloc(&B, 2, 2);
    mat_alloc(&C, 5, 2);
    a = mat_values(A);
    a[1][1] = 1.0; a[1][2] = 3.0;
    a[2][1] = 4.0; a[2][2] = 2.0;
    a[3][1] = 1.0; a[3][2] = 1.0;
    a[4][1] = 6.0; a[4][2] = 4.0;
    a[5][1] = 3.0; a[5][2] = 2.0;
    
    b= mat_values(B);
    b[1][1] = 1.0; b[1][2] = 2.0;
    b[2][1] = 3.0; b[2][2] = 4.0;
    
    mat_multiply(C, A, B);
    for (i = 1; i <= 5; i++) 
    	printf("% 2.2f\t% 2.2f\n", mat_element(A, i, 1), mat_element(A, i, 2));
    printf("\n");
    for (i = 1; i <= 2; i++) 
	printf("% 2.2f\t% 2.2f\n", mat_element(B, i, 1), mat_element(B, i, 2));
     printf("\n");
    for (i = 1; i <= 5; i++) 
    	printf("% 2.2f\t% 2.2f\n", mat_element(C, i, 1), mat_element(C, i, 2));
    printf("\n\n");

    mat_free(A); mat_free(B), mat_free(C);
    mat_alloc(&A, 3, 3);
    mat_alloc(&B, 3, 3);
    mat_alloc(&C, 3, 3);
    mat_init(A, 1.5);
    mat_init(B, 2.0);
    mat_init(C, 0.0);
    mat_multiply(C, A, B);
    /* Show matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(A, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(B, i, j));
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(C, i, j));
      }
    puts("");
    printf("\nResults are\n\
 1.50    1.50    1.50            2.00    2.00    2.00            9.00    9.00   9.00\n\
 1.50    1.50    1.50            2.00    2.00    2.00            9.00    9.00   9.00\n\
 1.50    1.50    1.50            2.00    2.00    2.00            9.00    9.00   9.00\n\n");


    mat_init(B, 1.5e-2);
    mat_scalar_product(A, B, 3.0);

    mat_init(B, 1.5);
    mat_init(C, 0.0);
    mat_scalar_product(C, B, 2.0);
    /* Show matrices  */
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(A, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(C, i, j));
      }
    printf("\n\nGives\n\
 0.04    0.04  0.04              3.00    3.00    3.00\n\
 0.04    0.04  0.04              3.00    3.00    3.00\n\
 0.04    0.04  0.04              3.00    3.00    3.00\n\n");


    mat_init(C, 0.0);
    a = mat_values(A);
    a[1][1] = 1.0; a[1][2] = 2.0; a[1][3] = 3.0;
    a[2][1] = 1.0; a[2][2] = 3.0; a[2][3] = 3.0;
    a[3][1] = 1.0; a[3][2] = 2.0; a[3][3] = 4.0;
    
    /* see its inverse upwards	*/
    
    mat_invert(B, A);  /* OK */
    mat_multiply(C, A, B);
    
    for (i = 1; i <= 3; i++)
      {
        puts("");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(A, i, j)); 
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(B, i, j));
        printf("\t");
        for (j = 1; j <= 3; j++)
          printf("% 1.2f\t", mat_element(C, i, j));
      }
    printf("\nMust give:\n\
 1.00    2.00    3.00            6.00   -2.00   -3.00            1.00    0.00   0.00\n\
 1.00    3.00    3.00           -1.00    1.00    0.00            0.00    1.00   0.00\n\
 1.00    2.00    4.00           -1.00    0.00    1.00            0.00    0.00   1.00\n\n");
    
}
