#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include 

int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    f[0] = 3*A*y[6]*sinbeta-9*A*y[5]*cosalpha+g/a;
    f[1] = y[5]*(omega^2-g/a)+3*A*cosalpha*(1-3*y[4]);
    f[2] = y[6]*(omega^2-mb/i-omega*y[3]-omega^2i3/i1)+A*a^2*M/i1*sinbeta*(3*y[4]-1);
    f[3] = omega*y[6]*y[2];
    f[4] = y[0];
    f[5] = y[1];
    f[6] = y[2];
    f[7] = y[3];
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy,
    double dfdt[], void *params)
{
    (void)(t); /* avoid unused parameter warning */
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

int main (int argc, char** argv)
{
    double mu = 10;
    gsl_odeiv2_system sys = {func, jac, 2, &mu};

    gsl_odeiv2_driver * d =
        gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                    1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 100.0;
    double y[2] = { 1.0, 0.0 };

    for (i = 1; i <= 100; i++)
        {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS)
            {
            printf ("error, return value=%d\n", status);
            break;
            }

        printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
        }

    gsl_odeiv2_driver_free (d);
    return 0;
}