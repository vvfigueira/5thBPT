
// BIBLIOTECAS

// // BIBLIOTECAS DO C

#include <stdio.h>
#include <fstream>
#include <vector>

// // BILBIOTECAS DO GSL --- NECESSÁRIAS PARA RESOLVER AS EDO'S

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

// // BILIOTECAS DO ROOT --- PARA GERAR GRÁFICOS COM OS DADOS OBTIDOS.
// // NÃO É NECESSÁRIA PARA RODAR O CÓDIGO, BASTA COMENTAR AS LINHAS ONDE ESTA É UTILIZADA

#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TMath.h>
#include <TLegend.h>
#include <TAxis.h>

// DEFINIÇÃO DOS PARÂMETROS RELEVANTES AO SISTEMA DE EQUAÇÕES

struct p_type {
    double omega;
    double g;
    double mu;
    double M;
    double i1;
    double i3;
    double m;
    double mp;
    double mb;
    double Omega;
};

// DEFINIÇÃO DAS EQUAÇÕES DIFERENCIAIS

int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t);

    struct p_type *npar = (struct p_type *)params;

    // DEFINIÇÃO DOS PARÂMETROS

    double omega = npar->omega;
    double g = npar->g;
    double mu = npar->mu;
    double M = npar->M;
    double i1 = npar->i1;
    double i3 = npar->i3;
    double m = npar->m;
    double mp = npar->mp;
    double mb = npar->mb;
    double Omega = npar->Omega;

    // AS VARIÁVEIS SÃO DEFINIDAS COMO :
    // y[0] = \dot{r}
    // y[1] = \dot{\theta}
    // y[2] = \dot{\phi}
    // y[3] = \dot{\theta'}
    // y[4] = \dot{\phi'}
    // y[5] = r
    // y[6] = \theta
    // y[7] = \phi
    // y[8] = \theta'
    // y[9] = \phi'


    // \ddot{r} = ...

    f[0] = y[5]*TMath::Power(TMath::Sin(y[6]),2)*y[2]*y[2]
                +y[5]*y[1]*y[1]
                -3*mu*m*mp*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)
                    *TMath::Sin(y[9]-y[7])-TMath::Sin(y[8])*TMath::Sin(y[9]-omega*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))
                    /(4*TMath::Pi()*M*TMath::Power(y[5],4))
                -TMath::Cos(y[6])*g;

    // \ddot{\theta} = ...

    f[1] = -2*y[1]*y[0]/y[5]
                +TMath::Sin(2*y[6])*y[2]*y[2]/2
                +(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(6*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))
                +TMath::Sin(y[6])*g/y[5];

    // \ddot{\phi} = ...
    
    f[2] = -2*y[2]*y[0]/y[5]
                -2*y[2]*TMath::Cos(y[6])*y[1]/TMath::Sin(y[6])
                +(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-omega*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t));

    // \ddot{\theta'} = ...
    
    f[3] = -i3*TMath::Sin(y[8])*y[4]*Omega/i1
                +TMath::Sin(2*y[8])*y[4]*y[4]/2
                -mb*TMath::Sin(y[8])/i1
                +(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Sin(y[9]-omega*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t));

    // \ddot{\phi'} = ...

    f[4] = i3*Omega*y[3]/(i1*TMath::Sin(y[8]))
                -2*y[4]*TMath::Cos(y[8])*y[3]/TMath::Sin(y[8])
                +(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t));

    // NECESSÁRIO POIS O SISTEMA DE N EQUAÇÕES DE SEGUNDA ORDEM 
    // FOI ESCRITO COMO 2N EQUAÇÕES DE PRIMEIRA ORDEM
    
    f[5] = y[0]; 
    f[6] = y[1]; 
    f[7] = y[2]; 
    f[8] = y[3]; 
    f[9] = y[4]; 

    // f[1] = 0; f[2] = 0; f[6] = 0; f[7] = 0;

    return GSL_SUCCESS;
}

// JACOBIANO DO SISTEMA

int jac (double t, const double y[], double *dfdy,
    double dfdt[], void *params)
{
    (void)(t);

    struct p_type *npar = (struct p_type *)params;

    // DEFINIÇÃO DOS PARÂMETROS

    double omega = npar->omega;
    double g = npar->g;
    double mu = npar->mu;
    double M = npar->M;
    double i1 = npar->i1;
    double i3 = npar->i3;
    double m = npar->m;
    double mp = npar->mp;
    double mb = npar->mb;
    double Omega = npar->Omega;

    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 10, 10);
    gsl_matrix * mat = &dfdy_mat.matrix;

    // PARA CADA ELEMENTO DE MATRIZ, (mat, 'n', 'm', p) => p = DERIVADA PARCIAL DE f[n]
    // EM RELAÇÃO A VARIÁVEL y[m]
    gsl_matrix_set (mat, 0, 0, 0.0);
    gsl_matrix_set (mat, 0, 1, 2*y[5]*y[1]);
    gsl_matrix_set (mat, 0, 2, 2*y[5]*TMath::Power(TMath::Sin(y[6]),2)*y[2]);
    gsl_matrix_set (mat, 0, 3, 0.0);
    gsl_matrix_set (mat, 0, 4, 0.0);
    gsl_matrix_set (mat, 0, 5, 12*mu*m*mp*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])
                        *TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])-TMath::Sin(y[8])*TMath::Sin(y[9]-omega*t)
                        +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))
                        /(4*TMath::Pi()*M*TMath::Power(y[5],5)));
    gsl_matrix_set (mat, 0, 6, y[5]*TMath::Sin(2*y[6])*y[2]*y[2]+TMath::Sin(y[6])*g
                        -3*mu*m*mp*(3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)
                        *TMath::Sin(y[9]-y[7])
                        +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))
                        /(4*TMath::Pi()*M*TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 7, -3*mu*m*mp*(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-y[9]-omega*t)
                        -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t))
                        /(4*TMath::Pi()*M*TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 8, -3*mu*m*mp*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)
                        *TMath::Sin(y[9]-y[7])-TMath::Cos(y[8])*TMath::Sin(y[9]-omega*t)
                        -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t))
                        /(4*TMath::Pi()*M*TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 0, 9, -3*mu*m*mp*(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)
                        *TMath::Cos(y[9]-y[7])-TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t))
                        /(4*TMath::Pi()*M*TMath::Power(y[5],4)));
    gsl_matrix_set (mat, 1, 0, -2*y[1]/y[5]);
    gsl_matrix_set (mat, 1, 1, -2*y[0]/y[5]);
    gsl_matrix_set (mat, 1, 2, TMath::Sin(2*y[6])*y[2]); 
    gsl_matrix_set (mat, 1, 3, 0.0);
    gsl_matrix_set (mat, 1, 4, 0.0);
    gsl_matrix_set (mat, 1, 5, -TMath::Sin(y[6])*g/(y[5]*y[5])+2*y[1]*y[0]/(y[5]*y[5])
                    -5*(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],6)))
                    *(6*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))); 
    gsl_matrix_set (mat, 1, 6, 2*TMath::Cos(2*y[6])*y[2]*y[2]/2
                    +(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(6*TMath::Cos(2*y[6])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    -6*TMath::Sin(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t))
                    +TMath::Cos(y[6])*g/y[5]);
    gsl_matrix_set (mat, 1, 7, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(-3*TMath::Sin(2*y[6])*TMath::Cos(2*y[7]-y[9]-omega*t)
                    -3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t)));
    gsl_matrix_set (mat, 1, 8, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(-3*TMath::Cos(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)));
    gsl_matrix_set (mat, 1, 9, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(3*TMath::Sin(2*y[6])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])));
    gsl_matrix_set (mat, 2, 0, -2*y[2]/y[5]);
    gsl_matrix_set (mat, 2, 1, -2*y[2]*TMath::Cos(y[6])/TMath::Sin(y[6]));
    gsl_matrix_set (mat, 2, 2, -2*y[0]/y[5]
                    -2*TMath::Cos(y[6])*y[1]/TMath::Sin(y[6]));
    gsl_matrix_set (mat, 2, 3, 0.0);
    gsl_matrix_set (mat, 2, 4, 0.0);
    gsl_matrix_set (mat, 2, 5, 2*y[2]*y[0]/(y[5]*y[5])
                    -5*(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],6)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-omega*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t)));
    gsl_matrix_set (mat, 2, 6, 2*y[2]*y[1]/(TMath::Sin(y[6])*TMath::Sin(y[6]))
                    +(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(2*y[7]-omega*t-y[9])
                    -3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t))
                    -(mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),3)))*TMath::Cos(y[6])
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(2*y[7]-omega*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Sin(y[7]-omega*t)));
    gsl_matrix_set (mat, 2, 7, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(6*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-omega*t-y[9])
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)));
    gsl_matrix_set (mat, 2, 8, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(2*y[7]-omega*t-y[9])
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Sin(y[7]-omega*t)));
    gsl_matrix_set (mat, 2, 9, (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-omega*t-y[9])));
    gsl_matrix_set (mat, 3, 0, 0.0);
    gsl_matrix_set (mat, 3, 1, 0.0);
    gsl_matrix_set (mat, 3, 2, 0.0);
    gsl_matrix_set (mat, 3, 3, 0.0);
    gsl_matrix_set (mat, 3, 4, -i3*TMath::Sin(y[8])*Omega/i1
                    +TMath::Sin(2*y[8])*y[4]);
    gsl_matrix_set (mat, 3, 5, 
                -3*(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],4)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Sin(y[9]-omega*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)));
    gsl_matrix_set (mat, 3, 6, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(3*TMath::Sin(2*y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    -3*TMath::Cos(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)));
    gsl_matrix_set (mat, 3, 7, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(2*y[7]-y[9]-omega*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*TMath::Sin(y[7]-omega*t)));
    gsl_matrix_set (mat, 3, 8, -i3*TMath::Cos(y[8])*y[4]*Omega/i1
                    +TMath::Cos(2*y[8])*y[4]*y[4]
                    -mb*TMath::Cos(y[8])/i1
                    +(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    +TMath::Sin(y[8])*TMath::Sin(y[9]-omega*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)));
    gsl_matrix_set (mat, 3, 9, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Cos(y[9]-omega*t)));
    gsl_matrix_set (mat, 4, 0, 0.0);
    gsl_matrix_set (mat, 4, 1, 0.0);
    gsl_matrix_set (mat, 4, 2, 0.0);
    gsl_matrix_set (mat, 4, 3, i3*Omega/(i1*TMath::Sin(y[8]))
                    -2*y[4]*TMath::Cos(y[8])/TMath::Sin(y[8]));
    gsl_matrix_set (mat, 4, 4, -2*TMath::Cos(y[8])*y[3]/TMath::Sin(y[8]));
    gsl_matrix_set (mat, 4, 5, -3*(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],4)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t)));
    gsl_matrix_set (mat, 4, 6, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Sin(2*y[6])*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t)));
    gsl_matrix_set (mat, 4, 7, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(2*y[7]-y[9]-omega*t)));
    gsl_matrix_set (mat, 4, 8, -i3*Omega*y[3]*TMath::Cos(y[8])/(i1*TMath::Sin(y[8])*TMath::Sin(y[8]))
                    +2*y[4]*y[3]/(TMath::Sin(y[8])*TMath::Sin(y[8]))
                    -(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),3)))*TMath::Cos(y[8])
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t))
                    +(mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Cos(y[8])*TMath::Cos(y[9]-omega*t)));
    gsl_matrix_set (mat, 4, 9, (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Cos(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    +TMath::Sin(y[8])*TMath::Sin(y[9]-omega*t)));
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

    // DERIVADAS TEMPORAIS, CADA TERMO dfdt['n'] == DERIVADA PARCIAL EM RELAÇÃO AO
    // TEMPO DA EQUAÇÃO f[n]

    dfdt[0] =   -3*mu*m*mp*(3*omega*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*TMath::Sin(y[7]-omega*t)
                    *TMath::Sin(y[9]-y[7])+omega*TMath::Sin(y[8])*TMath::Cos(y[9]-omega*t)
                    +3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*omega*TMath::Sin(y[7]-omega*t))
                    /(4*TMath::Pi()*M*TMath::Power(y[5],4));
    dfdt[1] =   (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)))
                    *(6*TMath::Sin(y[6])*TMath::Cos(y[6])*omega*TMath::Sin(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    +3*TMath::Cos(2*y[6])*TMath::Cos(y[8])*omega*TMath::Sin(y[7]-omega*t));
    dfdt[2] =   (mu*m*mp/(4*TMath::Pi()*M*TMath::Power(y[5],5)*TMath::Power(TMath::Sin(y[6]),2)))
                    *(-3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*omega*TMath::Sin(2*y[7]-omega*t-y[9])
                    +omega*3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Cos(y[8])*TMath::Cos(y[7]-omega*t));
    dfdt[3] =   (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Cos(y[8])*omega*TMath::Sin(y[7]-omega*t)*TMath::Sin(y[9]-y[7])
                    -omega*TMath::Cos(y[8])*TMath::Cos(y[9]-omega*t)
                    -3*TMath::Sin(y[6])*TMath::Cos(y[6])*TMath::Sin(y[8])*omega*TMath::Sin(y[7]-omega*t));
    dfdt[4] =   (mu*m*mp/(4*TMath::Pi()*i1*TMath::Power(y[5],3)*TMath::Power(TMath::Sin(y[8]),2)))
                    *(3*TMath::Power(TMath::Sin(y[6]),2)*TMath::Sin(y[8])*omega*TMath::Sin(y[7]-omega*t)*TMath::Cos(y[9]-y[7])
                    -TMath::Sin(y[8])*omega*TMath::Sin(y[9]-omega*t));
    dfdt[5] = 0.0;
    dfdt[6] = 0.0;
    dfdt[7] = 0.0;
    dfdt[8] = 0.0;
    dfdt[9] = 0.0;

    return GSL_SUCCESS;
}

// AQUI COMEÇA O PROGRAMA

int main (int argc, char** argv)
{

    double t0 = 0.0, t1 = 5, y[10], ti; // TEMPO INICIAL, TEMPO FINAL, ARRAY DAS VARIÁVEIS, TEMPO DE CADA PASSO
    int divtemp = 1000, status; // NUMMERO DE DIVISÕES ENTRE T. INICIAL E FINAL, STATUS INTERNO DO GSL
    
    int divO = 1, contpt; // NÚMERO DE ITERAÇÕES À SEREM RODADAS ALTERANDO ALGUM PARÂMETRO, 
                        // VARIÁVEL INTERNA DO ROOT
    double O0 = 0.0, O1 = 10000.0*2*TMath::Pi()/60; // VALOR DO PARÂMETRO INICIAL, VALOR FINAL

    std::string nome; // DEFINIÇÕES PARA ARQUIVO .TSV
    std::string arq[5];
    char * nomef, *arq_f[5];

    std::ofstream ofs;
    ofs.open("Data.tsv", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    std::ofstream eFile ("Data.tsv" ,std::ofstream::app); // PRINTA ARQUIVO COM VALORES
    
    // VALOR DOS PARÂMETROS

    struct p_type parametros = {5805*2*TMath::Pi()/60, // omega
                                9.8, // g
                                4e-7 *TMath::Pi(), // mu0
                                7.5e-3, // M
                                7.5e-3 * 0.01 * 0.01/6, // i1
                                7.5e-3 * 0.01 * 0.01/6, // i3
                                0.2375, // m
                                0.2375, // mp
                                1.0, // m'B
                                684.0 // Omega
    }; 

    // VALOR INICIAL DE CADA VARIÁVEL

    const double initval[10] = {0.0, // r ponto
                                0.0, // theta ponto
                                5805*2*TMath::Pi()/60, // phi ponto
                                0.0, // theta' ponto
                                5805*2*TMath::Pi()/60, // phi' ponto
                                0.02, // r
                                3.073, // theta
                                0.0, // phi
                                4.442, // theta' 
                                0.0 // phi'
    };

    gsl_odeiv2_system sys = {func, jac, 10, &parametros};   

    // PODE SER COMENTADO SE NÃO FOR UTILIZAR ROOT PARA PLOTAR
    
    // ##########*********************************

    std::vector<TGraph*> EpGvec, TGvec, TpGvec;

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
    
    // *********************************##########

    for(int k = 0; k < divO; k++)
    {

        t0 = 0.0, t1 = 5; 

        // ##########*********************************
        
        EpGvec.push_back(new TGraph);
        TGvec.push_back(new TGraph);
        TpGvec.push_back(new TGraph);

        // *********************************##########

        printf("Iteracao %i\n", k+1);

        contpt = 0;

        for (int l = 0; l< 10; l++){y[l]=initval[l];};

        //parametros.omega = O0 + k*O1/divO; // CASO FOREM SIMULADAS MAIS DE UMA ITERAÇÃO É NECESSÁRIO

        arq[0] = "Raio_";
        arq[0] = arq[0] + (int)parametros.Omega;
        arq_f[0] = &arq[0][0];

        arq[1] = "Phi_";
        arq[1] = arq[1] + (int)parametros.Omega;
        arq_f[1] = &arq[1][0];

        arq[2] = "Theta_";
        arq[2] = arq[2] + (int)parametros.Omega;
        arq_f[2] = &arq[2][0];

        arq[3] = "Phip_";
        arq[3] = arq[3] + (int)parametros.Omega;
        arq_f[3] = &arq[3][0];

        arq[4] = "Thetap_";
        arq[4] = arq[4] + (int)parametros.Omega;
        arq_f[4] = &arq[4][0];

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

            // printf ("%.5e %.5e %.5e %.5e\n", t0, y[5], y[6], y[8]);

            // PRINTA DADOS NO ARQUIVO

            eFile << t0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t"
                << y[4] << "\t" << y[5] << "\t" << y[6] << "\t" << y[7] << "\t" << y[8] << "\t" << y[9] << "\n";

            // ##########*********************************

            contpt = EpGvec[k]->GetN();

            EpGvec[k]->SetPoint(contpt, t0, y[5]);
            TGvec[k]->SetPoint(contpt, t0, y[6]);
            TpGvec[k]->SetPoint(contpt, t0, y[8]);

            // *********************************##########

        }

        // ##########*********************************

        EpGvec[k]->SetMarkerStyle(kFullCircle);
        EpGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        EpGvec[k]->SetMarkerSize(0.65);

        TGvec[k]->SetMarkerStyle(kFullCircle);
        TGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        TGvec[k]->SetMarkerSize(0.65);

        TpGvec[k]->SetMarkerStyle(kFullCircle);
        TpGvec[k]->SetMarkerColor(((k+1)%10) ? k+1 : 1 );
        TpGvec[k]->SetMarkerSize(0.65);

        EpMg->Add(EpGvec[k]);

        TMg->Add(TGvec[k]);

        TpMg->Add(TpGvec[k]);

        // *********************************##########

        nome = "Omega ";
        nome = nome + (int)parametros.Omega;
        nomef = &nome[0];

        // legend->AddEntry(EpGvec[k], nomef, "p");
    
        gsl_odeiv2_driver_free (d);

    }

    // ##########*********************************

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

    // *********************************##########

    return 0;
}