#include <math.h>
#include <stdio.h>
#include "doubles.h"

/**
   This demo illustrates use of AXLE in the context of trajectory
   interpolation. An initial and final state are hard-coded at (0,0)
   and (20,10); an initial solution is computed using linear
   interpolation; then AXLE is used to refine the estimate.

   With 2 dimensional unary constraints, the output will be that the
   vehicle travels on a constant heading between the two waypoints.


   --nstates #nstates   How long should the state vector be?

   --niter #niter       How many times should we iterate the solution?

   --unary2             Use 2 dimensional unary constraints (like the paper)

   --unary3             Use 3 dimensional unary constraints, with heading.

   --help               Show help

   --quiet              Don't show the output (useful for benchmarking)

Examples:

// with no heading constraints (unary2), heading is inferred from the end
// points (0.463 radians) and a straight-line trajectory results.

$ ./axle-demo --nstates 10 --niter 10 --unary2
ebolson@ubuntu-20:~/papers/axle/code/axle$ ./axle-demo --nstates 10 --niter 10 --unary2
  0       -0.000000        0.000000        2.484520        0.463648
  1        2.222222        1.111111        2.484520        0.463648
  2        4.444444        2.222222        2.484520        0.463648
  3        6.666667        3.333333        2.484520        0.463648
  4        8.888889        4.444444        2.484520        0.463648
  5       11.111111        5.555556        2.484520        0.463648
  6       13.333333        6.666667        2.484520        0.463648
  7       15.555556        7.777778        2.484520        0.463648
  8       17.777778        8.888889        2.484520        0.463648
  9       20.000000       10.000000        2.484520        0.463648

// with heading constraints (unary3), you can see the vehicle
// performing an "S" curve in column 4 and a higher overall speed
// required since more distance is covered.

$ ./axle-demo --nstates 10 --niter 10 --unary3
  0       -0.000116        0.000240        2.500637        0.003275
  1        2.488870        0.032431        2.512195        0.270637
  2        4.897987        0.728056        2.528550        0.472081
  3        7.138337        1.901892        2.544354        0.606093
  4        9.217854        3.375313        2.556050        0.673047
  5       11.204859        4.992683        2.561883        0.673488
  6       13.195720        6.614572        2.561841        0.607271
  7       15.287889        8.100433        2.557660        0.473547
  8       17.552457        9.290845        2.552889        0.271617
  9       20.000116        9.999760        2.552889        0.002689

**/
int main(int argc, char *argv[])
{
    int unary2 = 1;
    int unary3 = 0;
    int N = 10;
    int niter = 10;
    int help = 0;
    int quiet = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--nstates")) {
            N = atoi(argv[i+1]);
            i++;
            continue;
        }
        if (!strcmp(argv[i], "--niter")) {
            niter = atoi(argv[i+1]);
            i++;
            continue;
        }
        if (!strcmp(argv[i], "--unary2")) {
            unary2 = 1;
            unary3 = 0;
            continue;
        }
        if (!strcmp(argv[i], "--unary3")) {
            unary2 = 0;
            unary3 = 1;
            continue;
        }
        if (!strcmp(argv[i], "--help")) {
            help = 1;
            continue;
        }
        if (!strcmp(argv[i], "--quiet")) {
            quiet = 1;
            continue;
        }
    }

    if (help) {
        printf("usage: %s [--nstates #nstates] [--niter #niterations] [--unary2] [--unary3] [--help] [--quiet]\n",
               argv[0]);
        exit(1);
    }

     // state vector elements [x, y, v, theta]. For simplicity in this
     // example, we assume states are 1 second apart.
    double Xs[4*N][4];

    // create an initial state by specifying start/end conditions and
    // then linearly interpolating between them (ignoring kino-dynamic
    // constraints).
    double X0[] = { 0, 0, 0, 0 };
    double XN[] = { 20, 10, 0, 0 };

    for (int i = 0; i < N; i++) {
        double alpha = 1.0 * i / N;
        for (int j = 0; j < 4; j++) {
            Xs[i][j] = (1 - alpha)*X0[j] + alpha*XN[j];
        }
    }

    // allocate matrix storage.
    double As[2*N - 1][16];
    double Es[N][4];

    for (int iter = 0; iter < niter; iter++) {

        // build the matrices from the factors
        memset(As, 0, sizeof(As));
        memset(Es, 0, sizeof(Es));

        for (int i = 0; i < N; i++) {

            // The companion paper shows 2DOF unary constraints, i.e.,
            // without a heading observation. That makes the solution
            // to this problem uninteresting-- just drive straight
            // along the linear interpolation.

            if (i == 0 || i == N-1) {

                if (unary2) {
                    double r[2];
                    if (i == 0) {
                        r[0] = 0 - Xs[i][0];
                        r[1] = 0 - Xs[i][1];
                    } else {
                        r[0] = XN[0] - Xs[i][0];
                        r[1] = XN[1] - Xs[i][1];
                    }

                    double J[] = { -1, 0, 0, 0,
                                   0, -1, 0, 0 };

                    double W[] = { 1, 0,
                                   0, 1 };

                    double JtW[8];
                    doubles_mat_AtB(J, 2, 4,    W, 2, 2,   JtW, 4, 2);

                    double JtWJ[16];
                    doubles_mat_AB(JtW, 4, 2,   J, 2, 4,   JtWJ, 4, 4);

                    double JtWr[4];
                    doubles_mat_AB(JtW, 4, 2,   r, 2, 1,   JtWr, 4, 1);

                    doubles_add     (As[2*i], JtWJ, 16,  As[2*i]);
                    doubles_subtract(Es[i],   JtWr,  4,   Es[i]);
                }

                if (unary3) {
                    double r[3];
                    if (i == 0) {
                        r[0] = X0[0] - Xs[i][0];
                        r[1] = X0[1] - Xs[i][1];
                        r[2] = mod2pi(X0[3] - Xs[i][3]); // heading
                    } else {
                        r[0] = XN[0] - Xs[i][0];
                        r[1] = XN[1] - Xs[i][1];
                        r[2] = mod2pi(XN[3] - Xs[i][3]);
                    }

                    double J[] = { -1, 0, 0, 0,
                                   0, -1, 0, 0,
                                   0,  0, 0, -1};

                    double W[] = { 100, 0, 0,
                                   0, 100, 0,
                                   0, 0, 100};

                    double JtW[8];
                    doubles_mat_AtB(J, 3, 4,    W, 3, 3,   JtW, 4, 3);

                    double JtWJ[16];
                    doubles_mat_AB(JtW, 4, 3,   J, 3, 4,   JtWJ, 4, 4);

                    double JtWr[4];
                    doubles_mat_AB(JtW, 4, 3,   r, 3, 1,   JtWr, 4, 1);

                    doubles_add     (As[2*i], JtWJ, 16,  As[2*i]);
                    doubles_subtract(Es[i],   JtWr,  4,   Es[i]);
                }
            }

            if (i + 1 < N) {
                // do binary/kinematic links

                // assume states are spaced one second apart.
                double dt = 1;

                double r[] = { Xs[i+1][0] - (Xs[i][0] + Xs[i][2] * dt * cos(Xs[i][3])),
                               Xs[i+1][1] - (Xs[i][1] + Xs[i][2] * dt * sin(Xs[i][3])),
                               Xs[i+1][2] - Xs[i][2],
                               mod2pi(Xs[i+1][3] - Xs[i][3]) };

                double Ja[] = { -1, 0,   -dt*cos(Xs[i][3]),    Xs[i][2]*dt*sin(Xs[i][3]),
                                0, -1,   -dt*sin(Xs[i][3]),   -Xs[i][2]*dt*cos(Xs[i][3]),
                                0, 0, -1 / dt, 0,
                                0, 0, 0, -1 / dt };

                double Jb[] = { 1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1 / dt, 0,
                                0, 0, 0, 1 / dt };

                double W[] = { 1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, 1 };

                double JatW[16];
                doubles_mat_AtB(Ja, 4, 4,    W, 4, 4,    JatW, 4, 4);

                double JbtW[16];
                doubles_mat_AtB(Jb, 4, 4,    W, 4, 4,    JbtW, 4, 4);

                double JatWJa[16];
                doubles_mat_AB(JatW, 4, 4,   Ja, 4, 4,    JatWJa, 4, 4);
                doubles_add(As[2*i], JatWJa, 16, As[2*i]);

                double JatWJb[16];
                doubles_mat_AB(JatW, 4, 4,   Jb, 4, 4,    JatWJb, 4, 4);
                doubles_add(As[2*i+1], JatWJb, 16, As[2*i + 1]);

                double JbtWJb[16];
                doubles_mat_AB(JbtW, 4, 4,   Jb, 4, 4,    JbtWJb, 4, 4);
                doubles_add(As[2*(i+1)], JbtWJb, 16, As[2*(i + 1)]);

                assert(2*(i+1) < 2*N - 1);

                double JatWr[4];
                doubles_mat_AB(JatW, 4, 4,   r,  4, 1,    JatWr, 4, 1);
                doubles_subtract(Es[i], JatWr, 4, Es[i]);

                double JbtWr[4];
                doubles_mat_AB(JbtW, 4, 4,   r,  4, 1,    JbtWr, 4, 1);
                doubles_subtract(Es[i+1], JbtWr, 4, Es[i+1]);
            }
        }

        // tikahnov
        if (0) {
            for (int i = 0; i < 2 * N - 1; i+=2) {
                for (int j = 0; j < 4; j++) {
                    As[i][4*j+j] += 1;
                }
            }
        }

        // solve for L
        double Ls[2*N - 1][16];

        // inverses of the even values of Ls[]
        double LsInvs[2*N - 1][16];

        // along the way, compute the inverses of any even #'d L_i
        for (int i = 0; i < 2 * N - 1; i++) {
            if (i == 0) {
                doubles_mat_chol_lower(As[i], 4, 4,    Ls[i], 4, 4);
                doubles_mat_inv(Ls[i], 4, 4,    LsInvs[i], 4, 4);

            } else if (i & 1) {

                doubles_mat_AtBt(As[i], 4, 4,     LsInvs[i-1], 4, 4,    Ls[i], 4, 4);

            } else {
                double LLt[16];
                doubles_mat_ABt(Ls[i-1], 4, 4,   Ls[i-1], 4, 4,  LLt, 4, 4);
                double tmp[16];
                doubles_subtract(As[i],   LLt,  16,  tmp);
                doubles_mat_chol_lower(tmp, 4, 4,    Ls[i], 4, 4);

                doubles_mat_inv(Ls[i], 4, 4,    LsInvs[i], 4, 4);
            }
        }

        /*
        printf("Ls\n");
        for (int i = 0; i < 2*N - 1; i++) {
            doubles_mat_print(Ls[i], 4, 4, "%15.5f");
            printf("\n");
        }
        */

        double Us[N][4];
        for (int i = 0; i < N; i++) {
            if (i == 0) {
                doubles_mat_AB(LsInvs[2*i], 4, 4,   Es[i],  4, 1,   Us[i],  4, 1);
            } else {
                double LU[4];
                doubles_mat_AB(Ls[2*i-1], 4, 4,  Us[i-1], 4, 1,    LU, 4, 1);
                double tmp[4];
                doubles_subtract(Es[i], LU, 4, tmp);

                doubles_mat_AB(LsInvs[2*i], 4, 4,   tmp, 4, 1,   Us[i], 4, 1);
            }
        }

        double dXs[N][4];
        for (int i = N-1; i >= 0; i--) {

            if (i == N - 1) {
                doubles_mat_AtB(LsInvs[2*i], 4, 4,    Us[i], 4, 1,    dXs[i], 4, 1);
            } else {
                double LX[4];
                doubles_mat_AtB(Ls[2*i+1], 4, 4,     dXs[i+1], 4, 1,   LX, 4, 1);

                double tmp[4];
                doubles_subtract(Us[i],     LX, 4,     tmp);

                doubles_mat_AtB(LsInvs[2*i], 4, 4,   tmp, 4, 1,   dXs[i], 4, 1);
            }
        }


        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 4; j++) {
                Xs[i][j] += dXs[i][j];
            }
        }
    }


    if (!quiet) {
        for (int i = 0; i < N; i++) {
            printf("%3d %15f %15f %15f %15f\n", i, Xs[i][0], Xs[i][1], Xs[i][2], Xs[i][3]);
        }
    }


    return 0;
}
