#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "math_util.h"

#define TRRFN(root, suffix)    root##suffix
#define TRFN(root, suffix)     TRRFN(root, suffix)
#define TFN(suffix)            TRFN(TNAME, suffix)

/**
 * Perform element-wise addition of two vectors and store the result in a third vector.
 *
 * @param[in]    a       A vector to add
 * @param[in]    b       A vector to add
 * @param[in]    len     Length of the vectors to add
 * @param[out]   r       Vector where the result is stored
 * @pre length(a) >= len
 * @pre length(b) >= len
 * @pre length(r) >= len
 */
static inline void TFN(s_add) (const TNAME *a, const TNAME *b, int len, TNAME *r)
{
    for (int i = 0; i < len; i++) {
        r[i] = a[i] + b[i];
    }
}

/**
 * Perform element-wise subtraction of two vectors and store the result in a third vector.
 *
 * @param[in]    a       The vector to subtract from
 * @param[in]    b       The vector with the amount to subtract
 * @param[in]    len     Length of the vectors to subtract
 * @param[out]   r       Vector where the result is stored
 * @pre length(a) >= len
 * @pre length(b) >= len
 * @pre length(r) >= len
 */
static inline void TFN(s_subtract) (const TNAME *a, const TNAME *b, int len, TNAME *r)
{
    for (int i = 0; i < len; i++) {
        r[i] = a[i] - b[i];
    }
}

/**
 * Multiply matrix A by matrix B and store the result in a third matrix.
 *
 * Matrix should be in row-major order, allocated in a single packed array.
 * (This is compatible with matd.)
 *
 * @param[in]    A           A matrix to multiply
 * @param[in]    Arows       Number of rows in A
 * @param[in]    Acols       Number of columns in A
 * @param[in]    B           A matrix to multiply
 * @param[in]    Brows       Number of rows in B
 * @param[in]    Bcols       Number of columns in B
 * @param[out]   R           Matrix in which to store result of A * B
 * @param[in]    Rrows       Number of rows in R
 * @param[in]    Rcols       Number of columns in R
 * @pre Acols == Brows
 * @pre Arows == Rrows
 * @pre Bcols == Rcols
 * @pre A, B are matrices stored in row-major order.
 */
static inline void TFN(s_mat_AB) (const TNAME *A,
                                  int          Arows,
                                  int          Acols,
                                  const TNAME *B,
                                  int          Brows,
                                  int          Bcols,
                                  TNAME *      R,
                                  int          Rrows,
                                  int          Rcols)
{
    assert(Acols == Brows);
    assert(Rrows == Arows);
    assert(Bcols == Rcols);

    for (int Rrow = 0; Rrow < Rrows; Rrow++) {
        for (int Rcol = 0; Rcol < Rcols; Rcol++) {
            TNAME acc = 0;
            for (int i = 0; i < Acols; i++) {
                acc += A[Rrow * Acols + i] * B[i * Bcols + Rcol];
            }
            R[Rrow * Rcols + Rcol] = acc;
        }
    }
}

/**
 * Multiply the transpose of matrix A by matrix B and store the result in a third matrix.
 *
 * Matrix should be in row-major order, allocated in a single packed array.
 * (This is compatible with matd.)
 *
 * @param[in]    A           A matrix to multiply
 * @param[in]    Arows       Number of rows in A
 * @param[in]    Acols       Number of columns in A
 * @param[in]    B           A matrix to multiply
 * @param[in]    Brows       Number of rows in B
 * @param[in]    Bcols       Number of columns in B
 * @param[out]   R           Matrix in which to store result of A^t * B
 * @param[in]    Rrows       Number of rows in R
 * @param[in]    Rcols       Number of columns in R
 * @pre Arows == Brows
 * @pre Acols == Rrows
 * @pre Bcols == Rcols
 * @pre  A, B are matrices stored in row-major order.
 */
static inline void TFN(s_mat_AtB) (const TNAME *A,
                                   int          Arows,
                                   int          Acols,
                                   const TNAME *B,
                                   int          Brows,
                                   int          Bcols,
                                   TNAME *      R,
                                   int          Rrows,
                                   int          Rcols)
{
    assert(Arows == Brows);
    assert(Rrows == Acols);
    assert(Bcols == Rcols);

    for (int Rrow = 0; Rrow < Rrows; Rrow++) {
        for (int Rcol = 0; Rcol < Rcols; Rcol++) {
            TNAME acc = 0;
            for (int i = 0; i < Arows; i++) {
                acc += A[i * Acols + Rrow] * B[i * Bcols + Rcol];
            }
            R[Rrow * Rcols + Rcol] = acc;
        }
    }
}

static inline void TFN(s_mat_AtBt) (const TNAME *A,
                                    int          Arows,
                                    int          Acols,
                                    const TNAME *B,
                                    int          Brows,
                                    int          Bcols,
                                    TNAME *      R,
                                    int          Rrows,
                                    int          Rcols)
{
    assert(Arows == Brows);
    assert(Rrows == Acols);
    assert(Rcols == Brows);

    for (int Rrow = 0; Rrow < Rrows; Rrow++) {
        for (int Rcol = 0; Rcol < Rcols; Rcol++) {
            TNAME acc = 0;
            // column Rrow of A, row Rcol of B
            for (int i = 0; i < Arows; i++) {
                acc += A[i * Acols + Rrow] * B[Rcol * Bcols + i];
            }
            R[Rrow * Rcols + Rcol] = acc;
        }
    }
}

/**
 * Multiply matrix A by the transpose of matrix B and store the result in a third matrix.
 *
 * Matrix should be in row-major order, allocated in a single packed array.
 * (This is compatible with matd.)
 *
 * @param[in]    A           A matrix to multiply
 * @param[in]    Arows       Number of rows in A
 * @param[in]    Acols       Number of columns in A
 * @param[in]    B           A matrix to multiply
 * @param[in]    Brows       Number of rows in B
 * @param[in]    Bcols       Number of columns in B
 * @param[out]   R           Matrix in which to store result of A * B^t
 * @param[in]    Rrows       Number of rows in R
 * @param[in]    Rcols       Number of columns in R
 * @pre Acols == Bcols
 * @pre Arows == Rows
 * @pre Brows == Rcols
 * @pre  A, B are matrices stored in row-major order.
 */
static inline void TFN(s_mat_ABt) (const TNAME *A,
                                   int          Arows,
                                   int          Acols,
                                   const TNAME *B,
                                   int          Brows,
                                   int          Bcols,
                                   TNAME *      R,
                                   int          Rrows,
                                   int          Rcols)
{
    assert(Acols == Bcols);
    assert(Rrows == Arows);
    assert(Brows == Rcols);

    for (int Rrow = 0; Rrow < Rrows; Rrow++) {
        for (int Rcol = 0; Rcol < Rcols; Rcol++) {
            TNAME acc = 0;
            for (int i = 0; i < Acols; i++) {
                acc += A[Rrow * Acols + i] * B[Rcol * Bcols + i];
            }
            R[Rrow * Rcols + Rcol] = acc;
        }
    }
}

/**
 * Compute a matrix inverse, with pivoting, but few special
 * considerations for ill-conditioned matrices.
 **/
int TFN(s_mat_inv) (const TNAME *A, int nrows, int ncols, TNAME *X, int nrows2, int ncols2);

/** Computes the lower triangular factor of symmetric square matrix
 * A. Returns 1 if matrix is SPD. (Return value of zero means
 * ill-conditioned.)
 **/
int TFN(s_mat_chol_upper)(const TNAME *A, int nrows, int ncols, TNAME *X, int nrows2, int ncols2);
                         int TFN(s_mat_chol_lower)(const TNAME *A, int nrows, int ncols, TNAME *X, int nrows2, int ncols2);



#undef TRRFN
#undef TRFN
#undef TFN
