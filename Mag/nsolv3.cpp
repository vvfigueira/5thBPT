#include <stdio.h>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include "TAxis.h"
#include "TPad.h"
#include "TF1.h"
#include "vector"
#include <TLegend.h>

struct p_type {
    double A;
    double B;
    double C;
    double D;
    double E;
    double F;
};

int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t);

    struct p_type *npar = (struct p_type *)params;

    double A = npar->A;
    double B = npar->B;
    double C = npar->C;
    double D = npar->D;
    double E = npar->E;
    double F = npar->F;

    f[0] = y[5]*TMath::Power(TMath::Sin(y[6]),2)*y[2]*y[2] // \ddot{r} = r*{\dot\phi}^2*\sin^2\theta
                +y[5]*y[1]*y[1] // +r*{\dot\theta}^2
                -TMath::Cos(y[6])*B // -B*\cos\theta
                -3*A*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)
                    *TMath::Sin(y[9]-y[7])-TMath::Sin(y[8])*TMath::Sin(y[9]-F*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t))
                    /(TMath::Power(y[5],4)); // -3*A*(3*\sin^2\theta\sin\theta'\cos\qty(\phi-F*t)\sin\qty(\phi'-\phi)
                    // -\sin\theta'\sin\qty(\phi'-F*t)
                    // +3*\sin\theta\cos\theta\cos\theta'\cos\qty(\phi- Ft))/r^4
    f[1] = -2*y[1]*y[0]/y[5] // \ddot{\theta} = -2*\dot{\theta}*\dot{r}/r
                +y[2]*y[2]*TMath::Sin(2*y[6])/2 // +{\dot\phi}^2\sin\qty(2\theta)/2
                +B*TMath::Sin(y[6])/y[5] // +B\sin\theta/r
                +3*A*TMath::Cos(y[7]-F*t)*(TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Sin(y[9]-y[7])
                    +TMath::Cos(2*y[6])*TMath::Cos(y[8]))/(TMath::Power(y[5],5));
                    // +3*A*\cos\qty(\phi-F*t)*
                    // \qty(\sin\qty(2\theta)\sin\theta'\sin\qty(\phi'-\phi)+\cos\qty(2\theta)\cos\theta')/r^5
    f[2] = -2*y[2]*y[0]/y[5] // \ddot{\phi} = -2*\dot{\phi}*\dot{r}/r
                -2*y[2]*y[1]*TMath::Cos(y[6])/TMath::Sin(y[6]) // -2*\dot{\phi}\dot{\theta}*\cot\theta
                -3*A*(TMath::Sin(y[8])*TMath::Cos(2*y[7]-y[9]-F*t)
                    +TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t)/TMath::Sin(y[6]))/(TMath::Power(y[5],5)); 
                    // -3*A*\qty(\sin\theta'\cos\qty(2\phi-\phi'-Ft)
                    //+\cot\theta\cos\theta'\sin\qty(\phi-Ft))/r^5
    f[3] = -C*TMath::Sin(y[8])*y[4] // \ddot{\theta'}
                +TMath::Sin(2*y[8])*y[4]*y[4]/2
                -E*TMath::Sin(y[8])
                +(D/(TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Sin(y[9]-F*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t));
    f[4] = C*y[3]/(TMath::Sin(y[8]))
                -2*y[4]*TMath::Cos(y[8])*y[3]/TMath::Sin(y[8])
                +(D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-F*t));
    f[5] = y[0]; // Distância r 
    f[6] = y[1]; // Theta
    f[7] = y[2]; // Phi
    f[8] = y[3]; // Theta'
    f[9] = y[4]; // Phi'

    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy,
    double dfdt[], void *params)
{
    (void)(t);

    struct p_type *npar = (struct p_type *)params;

    double A = npar->A;
    double B = npar->B;
    double C = npar->C;
    double D = npar->D;
    double E = npar->E;
    double F = npar->F;

    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 10, 10);
    gsl_matrix * mat = &dfdy_mat.matrix;
    gsl_matrix_set (mat, 0, 0, 0.0);
    gsl_matrix_set (mat, 0, 1, 2*y[5]*y[1]);
    gsl_matrix_set (mat, 0, 2, 2*y[5]*TMath::Power(TMath::Sin(y[6]),2)*y[2]);
    gsl_matrix_set (mat, 0, 3, 0.0);
    gsl_matrix_set (mat, 0, 4, 0.0);
    gsl_matrix_set (mat, 0, 5, 12*A*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])
                        *TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])-TMath::Sin(y[8])*TMath::Sin(y[9]-F*t)
                        +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t))
                        /(TMath::Power(y[5],5)));
    gsl_matrix_set (mat, 0, 6, y[5]*TMath::Sin(2*y[6])*y[2]*y[2]+TMath::Sin(y[6])*B
                        -3*A*(3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)
                        *TMath::Sin(y[9]-y[7])
                        +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t))
                        /(TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 7, -3*A*(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-y[9]-F*t)
                        -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t))
                        /(TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 8, -3*A*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)
                        *TMath::Sin(y[9]-y[7])-TMath::Cos(y[8])*TMath::Sin(y[9]-F*t)
                        -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t))
                        /(TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 9, -3*A*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)
                        *TMath::Cos(y[9]-y[7])-TMath::Sin(y[8])*TMath::Cos(y[9]-F*t))
                        /(TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 1, 0, -2*y[1]/y[5]);
    gsl_matrix_set (mat, 1, 1, -2*y[0]/y[5]);
    gsl_matrix_set (mat, 1, 2, TMath::Sin(2*y[6])*y[2]); 
    gsl_matrix_set (mat, 1, 3, 0.0);
    gsl_matrix_set (mat, 1, 4, 0.0);
    gsl_matrix_set (mat, 1, 5, -TMath::Sin(y[6])*B/(y[5]*y[5])+2*y[1]*y[0]/(y[5]*y[5])
                    -5*(A/(TMath::Power(y[5],6)))
                    *(6*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t))); 
    gsl_matrix_set (mat, 1, 6, 2*TMath::Cos(2*y[6])*y[2]*y[2]/2
                    +(A/(TMath::Power(y[5],5)))
                    *(6*TMath::Cos(2*y[6])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    -6*TMath::Sin(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t))
                    +TMath::Cos(y[6])*B/y[5]);
    gsl_matrix_set (mat, 1, 7, (A/(TMath::Power(y[5],5)))
                    *(-3*TMath::Sin(2*y[6])*TMath::Cos(2*y[7]-y[9]-F*t)
                    -3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t)));
    gsl_matrix_set (mat, 1, 8, (A/(TMath::Power(y[5],5)))
                    *(-3*TMath::Cos(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)));
    gsl_matrix_set (mat, 1, 9, (A/(TMath::Power(y[5],5)))
                    *(3*TMath::Sin(2*y[6])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])));
    gsl_matrix_set (mat, 2, 0, -2*y[2]/y[5]);
    gsl_matrix_set (mat, 2, 1, -2*y[2]*TMath::Cos(y[6])/TMath::Sin(y[6]));
    gsl_matrix_set (mat, 2, 2, -2*y[0]/y[5]
                    -2*TMath::Cos(y[6])*y[1]/TMath::Sin(y[6]));
    gsl_matrix_set (mat, 2, 3, 0.0);
    gsl_matrix_set (mat, 2, 4, 0.0);
    gsl_matrix_set (mat, 2, 5, 2*y[2]*y[0]/(y[5]*y[5])
                    -5*(A/(TMath::Power(y[5],6)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-F*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t)));
    gsl_matrix_set (mat, 2, 6, 2*y[2]*y[1]/(TMath::Sin(y[6])*TMath::Sin(y[6]))
                    +(A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(2*y[7]-F*t-y[9])
                    -3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t))
                    -(A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),3)))*TMath::Cos(y[6])
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-F*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-F*t)));
    gsl_matrix_set (mat, 2, 7, (A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(6*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-F*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)));
    gsl_matrix_set (mat, 2, 8, (A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(2*y[7]-F*t-y[9])
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Sin(y[7]-F*t)));
    gsl_matrix_set (mat, 2, 9, (A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-F*t-y[9])));
    gsl_matrix_set (mat, 3, 0, 0.0);
    gsl_matrix_set (mat, 3, 1, 0.0);
    gsl_matrix_set (mat, 3, 2, 0.0);
    gsl_matrix_set (mat, 3, 3, 0.0);
    gsl_matrix_set (mat, 3, 4, -TMath::Sin(y[8])*C
                    +TMath::Sin(2*y[8])*y[4]);
    gsl_matrix_set (mat, 3, 5, 
                -3*(D/(TMath::Power(y[5],4)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Sin(y[9]-F*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)));
    gsl_matrix_set (mat, 3, 6, (D/(TMath::Power(y[5],3)))
                    *(3*TMath::Sin(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    -3*TMath::Cos(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)));
    gsl_matrix_set (mat, 3, 7, (D/(TMath::Power(y[5],3)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(2*y[7]-y[9]-F*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Sin(y[7]-F*t)));
    gsl_matrix_set (mat, 3, 8, -C*TMath::Cos(y[8])*y[4]
                    +TMath::Cos(2*y[8])*y[4]*y[4]
                    -E*TMath::Cos(y[8])
                    +(D/(TMath::Power(y[5],3)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    +TMath::Sin(y[8])*TMath::Sin(y[9]-F*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)));
    gsl_matrix_set (mat, 3, 9, (D/(TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Cos(y[9]-F*t)));
    gsl_matrix_set (mat, 4, 0, 0.0);
    gsl_matrix_set (mat, 4, 1, 0.0);
    gsl_matrix_set (mat, 4, 2, 0.0);
    gsl_matrix_set (mat, 4, 3, C/(TMath::Sin(y[8]))
                    -2*y[4]*TMath::Cos(y[8])/TMath::Sin(y[8]));
    gsl_matrix_set (mat, 4, 4, -2*TMath::Cos(y[8])*y[3]/TMath::Sin(y[8]));
    gsl_matrix_set (mat, 4, 5, -3*(D/(TMath::Power(y[5],4)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-F*t)));
    gsl_matrix_set (mat, 4, 6, (D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-F*t)));
    gsl_matrix_set (mat, 4, 7, (D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-y[9]-F*t)));
    gsl_matrix_set (mat, 4, 8, -C*y[3]*TMath::Cos(y[8])/(TMath::Sin(y[8])*TMath::Sin(y[8]))
                    +2*y[4]*y[3]/(TMath::Sin(y[8])*TMath::Sin(y[8]))
                    -(D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),3)))*TMath::Cos(y[8])
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-F*t))
                    +(D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Cos(y[9]-F*t)));
    gsl_matrix_set (mat, 4, 9, (D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    +TMath::Sin(y[8])*TMath::Sin(y[9]-F*t)));
    gsl_matrix_set (mat, 5, 0, 1.0);
    gsl_matrix_set (mat, 5, 1, 0.0);
    gsl_matrix_set (mat, 5, 2, 0.0);
    gsl_matrix_set (mat, 5, 3, 0.0);
    gsl_matrix_set (mat, 5, 4, 0.0);
    gsl_matrix_set (mat, 5, 5, 0.0);
    gsl_matrix_set (mat, 5, 6, 0.0);
    gsl_matrix_set (mat, 5, 7, 0.0);
    gsl_matrix_set (mat, 5, 8, 0.0);
    gsl_matrix_set (mat, 5, 9, 0.0);
    gsl_matrix_set (mat, 6, 0, 0.0);
    gsl_matrix_set (mat, 6, 1, 1.0);
    gsl_matrix_set (mat, 6, 2, 0.0);
    gsl_matrix_set (mat, 6, 3, 0.0);
    gsl_matrix_set (mat, 6, 4, 0.0);
    gsl_matrix_set (mat, 6, 5, 0.0);
    gsl_matrix_set (mat, 6, 6, 0.0);
    gsl_matrix_set (mat, 6, 7, 0.0);
    gsl_matrix_set (mat, 6, 8, 0.0);
    gsl_matrix_set (mat, 6, 9, 0.0);
    gsl_matrix_set (mat, 7, 0, 0.0);
    gsl_matrix_set (mat, 7, 1, 0.0);
    gsl_matrix_set (mat, 7, 2, 1.0);
    gsl_matrix_set (mat, 7, 3, 0.0);
    gsl_matrix_set (mat, 7, 4, 0.0);
    gsl_matrix_set (mat, 7, 5, 0.0);
    gsl_matrix_set (mat, 7, 6, 0.0);
    gsl_matrix_set (mat, 7, 7, 0.0);
    gsl_matrix_set (mat, 7, 8, 0.0);
    gsl_matrix_set (mat, 7, 9, 0.0);
    gsl_matrix_set (mat, 8, 0, 0.0);
    gsl_matrix_set (mat, 8, 1, 0.0);
    gsl_matrix_set (mat, 8, 2, 0.0);
    gsl_matrix_set (mat, 8, 3, 1.0);
    gsl_matrix_set (mat, 8, 4, 0.0);
    gsl_matrix_set (mat, 8, 5, 0.0);
    gsl_matrix_set (mat, 8, 6, 0.0);
    gsl_matrix_set (mat, 8, 7, 0.0);
    gsl_matrix_set (mat, 8, 8, 0.0);
    gsl_matrix_set (mat, 8, 9, 0.0);
    gsl_matrix_set (mat, 9, 0, 0.0);
    gsl_matrix_set (mat, 9, 1, 0.0);
    gsl_matrix_set (mat, 9, 2, 0.0);
    gsl_matrix_set (mat, 9, 3, 0.0);
    gsl_matrix_set (mat, 9, 4, 1.0);
    gsl_matrix_set (mat, 9, 5, 0.0);
    gsl_matrix_set (mat, 9, 6, 0.0);
    gsl_matrix_set (mat, 9, 7, 0.0);
    gsl_matrix_set (mat, 9, 8, 0.0);
    gsl_matrix_set (mat, 9, 9, 0.0);

    dfdt[0] =   -3*A*(3*F*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(y[7]-F*t)
                    *TMath::Sin(y[9]-y[7])+F*TMath::Sin(y[8])*TMath::Cos(y[9]-F*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*F*TMath::Sin(y[7]-F*t))
                    /(TMath::Power(y[5],4));
    dfdt[1] =   (A/(TMath::Power(y[5],5)))
                    *(6*TMath::Sin(y[6])*TMath::Cos(y[6])*F*TMath::Sin(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*F*TMath::Sin(y[7]-F*t));
    dfdt[2] =   (A/(TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*F*TMath::Sin(2*y[7]-F*t-y[9])
                    +F*3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-F*t));
    dfdt[3] =   (D/(TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*F*TMath::Sin(y[7]-F*t)*TMath::Sin(y[9]-y[7])
                    -F*TMath::Cos(y[8])*TMath::Cos(y[9]-F*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*F*TMath::Sin(y[7]-F*t));
    dfdt[4] =   (D/(TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*F*TMath::Sin(y[7]-F*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*F*TMath::Sin(y[9]-F*t));
    dfdt[5] = 0.0;
    dfdt[6] = 0.0;
    dfdt[7] = 0.0;
    dfdt[8] = 0.0;
    dfdt[9] = 0.0;

    return GSL_SUCCESS;
}

int main (int argc, char** argv)
{

    double t0 = 0.0, t1 = 5, O0 = 1.0, O1 = 1.0, y[10], ti;
    int divtemp = 1000, divO = 1, contpt, status;

    std::string nome;
    char * nomef;

    std::ofstream ofs;
    ofs.open("Data.tsv", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    std::ofstream eFile ("Data.tsv" ,std::ofstream::app);
    
    struct p_type parametros = {1.0, // A
                                1.0, // B
                                1.0, // C
                                1.0, // D
                                0.0, // E
                                1.0, // F
    }; 

    const double initval[10] = {0.0, // r ponto
                                0.0, // theta ponto
                                0.0, // phi ponto
                                0.0, // theta' ponto
                                0.0, // phi' ponto
                                0.015, // r
                                3.14, // theta
                                0.0, // phi
                                0.0001, // theta' 
                                0.0 // phi'
    };

    gsl_odeiv2_system sys = {func, jac, 10, &parametros};   

    std::vector<TGraph*> EpGvec, TGvec, TpGvec;

    // TGraph * teste = new TGraph();

    TCanvas *EpCv = new TCanvas();

    EpCv->SetTickx();
    EpCv->SetTicky();
    EpCv->SetGridx();
    EpCv->SetGridy();

    TCanvas *TCv = new TCanvas();

    TCv->SetTickx();
    TCv->SetTicky();
    TCv->SetGridx();
    TCv->SetGridy();

    TCanvas *TpCv = new TCanvas();

    TpCv->SetTickx();
    TpCv->SetTicky();
    TpCv->SetGridx();
    TpCv->SetGridy();

    TMultiGraph *EpMg = new TMultiGraph();
    TMultiGraph *TMg = new TMultiGraph();
    TMultiGraph *TpMg = new TMultiGraph();

    TLegend *legend = new TLegend(0.6,0.65,0.88,0.85);
    legend->SetTextFont(132);
    legend->SetFillStyle(4000);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    

    for(int k = 0; k < divO; k++)
    {

        t0 = 0.0, t1 = 5;
        
        EpGvec.push_back(new TGraph);
        TGvec.push_back(new TGraph);
        TpGvec.push_back(new TGraph);

        printf("Iteracao %i\n", k+1);

        contpt = 0;

        for (int l = 0; l< 10; l++){y[l]=initval[l];};

        parametros.E = O0 + (double)(k*O1/divO);

        gsl_odeiv2_driver * d=
            gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-6, 1e-6, 0.0);

        for (int i = 0; i < divtemp; i++)
        {
            ti = i * t1 / (double )divtemp;

            status = gsl_odeiv2_driver_apply (d, &t0, ti, y);

            if (status != GSL_SUCCESS )
            {
                printf ("error, return value = %d\n", status);
                break;
            }

            printf ("%.5e %.5e %.5e %.5e\n", t0, y[5], y[6], y[8]);

            eFile << t0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t"
                << y[4] << "\t" << y[5] << "\t" << y[6] << "\t" << y[7] << "\t" << y[8] << "\t" << y[9] << "\n";

            contpt = EpGvec[k]->GetN();

            EpGvec[k]->SetPoint(contpt, t0, y[5]);
            // teste->SetPoint(contpt, t0, 9.8*t0*t0/2);
            TGvec[k]->SetPoint(contpt, t0, y[6]);
            TpGvec[k]->SetPoint(contpt, t0, y[8]);

        }

        EpGvec[k]->SetMarkerStyle(kFullCircle);
        EpGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        EpGvec[k]->SetMarkerSize(0.65);

        // teste->SetMarkerStyle(kFullCircle);
        // teste->SetMarkerColor(kRed);
        // teste->SetMarkerSize(0.65);

        TGvec[k]->SetMarkerStyle(kFullCircle);
        TGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        TGvec[k]->SetMarkerSize(0.65);

        TpGvec[k]->SetMarkerStyle(kFullCircle);
        TpGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        TpGvec[k]->SetMarkerSize(0.65);

        EpMg->Add(EpGvec[k]);
        // EpMg->Add(teste);

        TMg->Add(TGvec[k]);

        TpMg->Add(TpGvec[k]);

        nome = "Omega ";
        nome = nome + (int)parametros.A;
        nomef = &nome[0];

        // legend->AddEntry(EpGvec[k], nomef, "p");
    
        gsl_odeiv2_driver_free (d);

    }

    EpMg->GetXaxis()->SetLabelFont(132);
    EpMg->GetXaxis()->SetTitleFont(132);
    EpMg->GetYaxis()->SetLabelFont(132);
    EpMg->GetYaxis()->SetTitleFont(132);
    EpMg->GetYaxis()->SetLabelSize(0.035);
    EpMg->GetYaxis()->SetTitleSize(0.035);
    EpMg->GetXaxis()->SetLabelSize(0.035);
    EpMg->GetXaxis()->SetTitleSize(0.035);

    EpMg->GetXaxis()->SetLimits(0, t1);
    EpMg->GetYaxis()->SetMaxDigits(2);
    EpMg->GetXaxis()->SetTitle("Tempo #bf{[s]}");
    EpMg->GetYaxis()->SetTitle("Distancia #bf{[m]}");
    EpCv->cd();
    EpMg->Draw("AP");
    legend->SetY1(0.45);
    legend->SetY2(0.85);
    legend->SetX1(0.67);
    legend->SetX2(0.97);
    legend->Draw();

    EpCv->Print("Epsilon.pdf");

    TMg->GetXaxis()->SetLabelFont(132);
    TMg->GetXaxis()->SetTitleFont(132);
    TMg->GetYaxis()->SetLabelFont(132);
    TMg->GetYaxis()->SetTitleFont(132);
    TMg->GetYaxis()->SetLabelSize(0.035);
    TMg->GetYaxis()->SetTitleSize(0.035);
    TMg->GetXaxis()->SetLabelSize(0.035);
    TMg->GetXaxis()->SetTitleSize(0.035);

    TMg->GetXaxis()->SetLimits(0, t1);
    TMg->GetYaxis()->SetMaxDigits(2);
    TMg->GetXaxis()->SetTitle("Tempo #bf{[s]}");
    TMg->GetYaxis()->SetTitle("Theta #bf{[rad]}");
    TCv->cd();
    TMg->Draw("AP");
    legend->Draw();

    TCv->Print("Theta.pdf");

    TpMg->GetXaxis()->SetLabelFont(132);
    TpMg->GetXaxis()->SetTitleFont(132);
    TpMg->GetYaxis()->SetLabelFont(132);
    TpMg->GetYaxis()->SetTitleFont(132);
    TpMg->GetYaxis()->SetLabelSize(0.035);
    TpMg->GetYaxis()->SetTitleSize(0.035);
    TpMg->GetXaxis()->SetLabelSize(0.035);
    TpMg->GetXaxis()->SetTitleSize(0.035);

    TpMg->GetXaxis()->SetLimits(0, t1);
    TpMg->GetYaxis()->SetMaxDigits(2);
    TpMg->GetXaxis()->SetTitle("Tempo #bf{[s]}");
    TpMg->GetYaxis()->SetTitle("Theta' #bf{[rad]}");
    TpCv->cd();
    TpMg->Draw("AP");
    legend->Draw();

    TpCv->Print("Theta'.pdf");

    while (!EpGvec.empty())
    {
        TGraph* f = EpGvec.back();
        EpGvec.pop_back();
        delete f;
    }

    while (!TGvec.empty())
    {
        TGraph* f = TGvec.back();
        TGvec.pop_back();
        delete f;
    }

    while (!TpGvec.empty())
    {
        TGraph* f = TpGvec.back();
        TpGvec.pop_back();
        delete f;
    }

    delete EpMg, TMg, TpMg, EpCv, TCv, TpCv;

    return 0;
}