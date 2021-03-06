/**
 *  @file   mat_housekeep.c
 *
 *  @brief  Routines for creating/destroying 2D matrices.
 *
 *  This module defines the routines used for creating and destroying
 *  bidimensional matrices used for 2D matrix calculus. This is the 
 *  starting point before doing any worthy work.
 *
 *  In order to manipulate matrices using the libmatrix library, you
 *  need to housekeep them with the routines provided in this module.
 *  The library uses a special implementation in order to allow for
 *  easy and efficient matrix manipulation: to achieve this, special,
 *  additional information is added to the matrix data type. For this
 *  reason you should always create and destroy 2D-matrices using the
 *  routines provided.
 *
 *  An interesting side effect is that by using these routines and the
 *  ones provided in the matrix_access module you do not need to be
 *  aware of the internals of the implementation: not only means this
 *  less mnemonic effort, it also gives the implementer more freedom
 *  to adapt the implementation in the future without any need to
 *  affect your programs, thus increasing maintainability.
 *
 *  Right now, two implementations are provided, one is more straightforward
 *  and easier to understand and maintain, while the other should help
 *  yield some speed increases in some very frequent operations. You
 *  can select which one is used at compile time by defining 
 *  MAT_OPTIMIZE.
 *
 *  @pre    matrix.h
 *
 *  @see    matrix.h to learn about the definition of a 
 *	    2D-matrix and the additional housekeeping information used to
 *	    better understand this module.
 *
 *  @see    mat_access.c once you know how to create/destroy bi-dimensional 
 *	    matrices, you may procceed to learn more about how to obtain
 *	    various bits of information about your matrices and access its
 *	    data.
 *
 *  @see    mat_test.c to view some examples on using this code and
 *	    obtain debugging information.
 *
 *  @author Jos� Ram�n Valverde Carrillo    (jrvalverde@acm.org)
 *
 *  @version	3.0
 *
 *  @date   23 - february - 2004    v3.0
 *
 *  @date   11 - february - 2004    v2.0
 *
 *  @date    1 - october - 1988     Last modification of v1.0
 *
 *	COPYRIGHT:
 *	    � YoEgo.	Since I have no cash, I can't
 *	register this (nor do I believe I should). So
 *	this module is left in the PUBLIC DOMAIN.
 *	    It is furthermore forbidden its use for
 *	commercial purposes unless I get a share on
 *	the profits.
 *	    I say.
 *						YoEgo.
 *
 * $Id: mat_housekeep.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_housekeep.c,v $
 * Revision 1.1  2004/03/05 18:50:14  jr
 * Initial revision
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/nan.h>

#include "portable.h"
#include "matrix.h"

/*----------------------------------------------------------------*/
/*			MEMORY MANAGEMENT			  */
/*----------------------------------------------------------------*/

/** @defgroup matrix_housekeeping 
 * @{ 
 */

/**
 *  @fn public status mat_alloc(matrix *pmat, int rows, int cols);
 *
 *  @brief Allocate a new, uninitialized 2-D matrix
 *
 *  Allocates memory space for a bidimensional matrix
 *  of 1..rows x 1..cols elements
 *
 *  @param pmat     address of the matrix to allocate
 *  @param rows     the number of rows in the 2D-matrix
 *  @param cols     the number of columns in the 2D-matrix
 *
 *  @return	    success if all went well, an error code otherwise
 *
 *  @note   if MAT_PARANOIA has been defined, parameter checking is performed
 */

public status mat_alloc(matrix *pmat, int rows, int cols)
{
    auto int i, j;
    auto real **row_vector;
    auto real *this_row;
    matrix mat;

    *pmat = NULL;
#   ifdef MAT_PARANOIA
    /* Check dimensions for validity */
    if ((rows <= 0) || (cols <= 0)) {
	return MAT_BOUNDSCHECK;
    }
#   endif

    /* allocate structure */
    mat = (matrix) malloc(sizeof(struct matrix));
    if (mat == NULL)
	return MAT_NOMEMORY;

    /* allocate array of vectors */
    if ((row_vector = malloc(rows * sizeof(real *))) == NULL) {
	free(mat);
	return MAT_NOMEMORY;
    }
    /* make one-offset */
    row_vector--;
    mat->values = row_vector;

#   ifndef MAT_OPTIMIZE
    /* allocate each and every vector */
    for (i = 1; i <= rows; i++) {
	if ((this_row = malloc(cols * sizeof(real))) == NULL) {
	    /* roll back and free any allocated vectors */
	    for (j = i - 1; j >= 1; j--) {
		/* this shouldn't execute if initially i == 1 */
		this_row = mat->values[j];
		/* make row zero-offset */
		this_row++;
		free(this_row);
	    }
	    row_vector = mat->values;
	    row_vector++;
	    free(row_vector);
	    free(mat);
	    return MAT_NOMEMORY;
	}
	/* make row one-offset */
	this_row--;
	mat->values[i] = this_row;
    }
#   else
    /* trick: we allocate all rows at once */
    if ((this_row = malloc(rows * cols * sizeof(real))) == NULL) {
	free(mat->values);
	free(mat);
	return MAT_NOMEMORY;
    }
    /* make zero-offset */
    this_row--;
    for (i = 1; i <= rows; i++) {
	/* compute address of each row */
	/* this should be much faster than allocating each separately */
	/* and will allow for other speed-ups */
	mat->values[i] = this_row;
	this_row += cols;
    }
#   endif

    /* fill in matrix dimensions */
    mat->rows = rows;
    mat->cols = cols;

    *pmat = mat;
    return SUCCESS;
}

/**
 *  @fn public void mat_free(matrix mat)
 *
 *  @brief free a 2D-matrix
 *
 *  @param mat  the matrix to be freed, must have been allocated by mat_alloc()
 *
 *  @return SUCCESS if all went well
 */
public status mat_free(matrix mat)
{
    int i, maxrow;
    real *row;
    real **rows_vector;

    if (mat == NULL)
	return SUCCESS;

    /* NOTE: mat should have been allocated through mat_alloc */
    /* if not, or if mat->rows is corrupted, this will bomb out */
    maxrow = mat->rows;
#   ifndef MAT_OPTIMIZE
    for (i = 1; i <= maxrow; i++) {
	row = mat->values[i];
	row++;
	free(row);
    }
#   else
    /* we allocated the entire matrix as a single array */
    row = mat->values[1];
    row++;
    free(row);
#   endif
    rows_vector = mat->values;
    rows_vector++;
    free(rows_vector);
    free(mat);
    return SUCCESS;
}

/** @} */
