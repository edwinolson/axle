#include <assert.h>
#include <math.h>
#include <string.h>

// Should not be compiled directly... see doubles.c & floats.c
//
// Functions are put the .c when they are too big to be inlined.

#define TRRFN(root, suffix)    root##suffix
#define TRFN(root, suffix)     TRRFN(root, suffix)
#define TFN(suffix)            TRFN(TNAME, suffix)

int TFN(s_mat_inv) (const TNAME *A, int nrows, int ncols, TNAME *X, int nrows2, int ncols2)
{
    const int dim = nrows;

    assert(dim == ncols);
    assert(dim == nrows2);
    assert(dim == ncols2);

    int piv[dim];
    int pivsign  = 1;
    int singular = 0;

    TNAME LU[dim * dim];
    memcpy(LU, A, dim * dim * sizeof(TNAME));

    for (int i = 0; i < dim; i++) {
        piv[i] = i;
    }

    for (int j = 0; j < dim; j++) {
        for (int i = 0; i < dim; i++) {
            int kmax = i < j ? i : j;

            TNAME acc = 0;
            for (int k = 0; k < kmax; k++) {
                acc += LU[dim * i + k] * LU[dim * k + j];
            }

            LU[dim * i + j] -= acc;
        }

        // find best pivot
        int p = j;
        for (int i = j + 1; i < dim; i++) {
            if (fabs(LU[dim * i + j]) > fabs(LU[dim * p + j])) {
                p = i;
            }
        }

        // swap rows p and j
        if (p != j) {
            TNAME tmp[dim];
            memcpy(tmp, &LU[dim * p], dim * sizeof(TNAME));
            memcpy(&LU[dim * p], &LU[dim * j], dim * sizeof(TNAME));
            memcpy(&LU[dim * j], tmp, dim * sizeof(TNAME));
            int k   = piv[p];
            piv[p]  = piv[j];
            piv[j]  = k;
            pivsign = -pivsign;
        }

        TNAME LUjj = LU[dim * j + j];

        if (fabs(LUjj) < 1.0E-8) {
            singular = 1;
        }

        if (j < dim && j < dim && LUjj != 0) {
            LUjj = 1.0 / LUjj;
            for (int i = j + 1; i < dim; i++) {
                LU[dim * i + j] *= LUjj;
            }
        }
    }

    if (!singular) {
        // fill X with permuted version of identity matrix
        memset(X, 0, sizeof(TNAME) * dim * dim);

        for (int i = 0; i < dim; i++) {
            X[dim * i + piv[i]] = 1;
        }

        // solve Ly = b
        for (int k = 0; k < dim; k++) {
            for (int i = k + 1; i < dim; i++) {
                TNAME LUik = -LU[dim * i + k];
                for (int t = 0; t < dim; t++) {
                    X[dim * i + t] += X[dim * k + t] * LUik;
                }
            }
        }

        // solve Ux = y
        for (int k = dim - 1; k >= 0; k--) {
            TNAME LUkk = 1.0 / LU[dim * k + k];
            for (int t = 0; t < dim; t++) {
                X[dim * k + t] *= LUkk;
            }

            for (int i = 0; i < k; i++) {
                TNAME LUik = -LU[dim * i + k];
                for (int t = 0; t < dim; t++) {
                    X[dim * i + t] += X[dim * k + t] * LUik;
                }
            }
        }
    }

    return singular;
}

int TFN(s_mat_chol_upper)(const TNAME *A, int nrows, int ncols, TNAME *U, int nrows2, int ncols2)
{
    assert(nrows == ncols);
    assert(nrows == nrows2);
    assert(ncols == ncols2);

    int N = nrows;

    // copy upper triangular
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            U[i*N + j] = (j < i) ? 0 : A[i*N + j];
        }
    }

    int is_spd = 1;    // (A->nrows == A->ncols);

    for (int i = 0; i < N; i++) {
        double d = U[i*N + i];
        is_spd &= (d > 0);

        if (d < 1E-8) {
            d = 1E-8;
        }
        d = 1.0 / sqrt(d);

        for (int j = i; j < N; j++) {
            U[i*N + j] *= d;
        }

        for (int j = i + 1; j < N; j++) {
            double s = U[i*N + j];

            if (s == 0) {
                continue;
            }

            for (int k = j; k < N; k++) {
                U[j*N + k] -= U[i*N + k] * s;
            }
        }
    }

    return is_spd;
}

int TFN(s_mat_chol_lower)(const TNAME *A, int nrows, int ncols, TNAME *U, int nrows2, int ncols2)
{
    if (0) {
        TNAME tmp[nrows2*ncols2];

        int res = TFN(s_mat_chol_upper)(A, nrows, ncols, tmp, nrows2, ncols2);
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < nrows; j++) {
                U[i*nrows + j] = tmp[j*nrows + i];
            }
        }
        return res;
    }

    assert(nrows == ncols);
    assert(nrows == nrows2);
    assert(ncols == ncols2);

    int N = nrows;

    // copy lower triangular
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            U[i*N + j] = (j > i) ? 0 : A[i*N + j];
        }
    }

    int is_spd = 1;    // (A->nrows == A->ncols);

    for (int i = 0; i < N; i++) {
        double d = U[i*N + i];
        is_spd &= (d > 0);

        if (d < 1E-8) {
            d = 1E-8;
        }
        d = 1.0 / sqrt(d);

        for (int j = i; j < N; j++) {
            U[j*N + i] *= d;
        }

        for (int j = i + 1; j < N; j++) {
            double s = U[j*N + i];

            if (s == 0) {
                continue;
            }

            for (int k = j; k < N; k++) {
                U[k*N + j] -= U[k*N + i] * s;
            }
        }
    }

    return is_spd;
}

#undef TRRFN
#undef TRFN
#undef TFN
