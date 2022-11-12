#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <TCanvas.h>
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TAxis.h"
#include "TPad.h"
#include "TF1.h"

const double omega = 7300*2*M_PI/60;
const double g = 9.8;
const double a = 0.02;
const double M = 7.5e-3;
const double i1 = M * 0.01* 0.01/6;
const double i3 = i1;
const double A = 235;
const double mb = 1000;
const double Omega = 0 ;

int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t);

    f[0] = y[1]*y[1]/2-3*A*(3*y[6]*TMath::Cos(y[7]-omega*t)-y[8]*TMath::Sin(y[9]-omega*t))+g/a;
    f[1] = y[6]*y[2]*y[2]/(1+2*y[5])+3*A*(1-3*y[5])*TMath::Cos(y[7]-omega*t)/(1+2*y[5])-g*y[6]/(a*(1+2*y[5]))-2*y[0]*y[1]/(1+2*y[5]);
    f[2] = -3*A*(1-3*y[5])*TMath::Sin(y[7]-omega*t)/y[6]-2*y[2]*y[1]/y[6];
    f[3] = y[4]*y[4]*y[8]-i3*y[8]*y[3]*y[4]/i1-i3*y[8]*y[8]*y[4]/i1-A*a*a*M*(1-3*y[5])*TMath::Sin(y[9]-omega*t)/i1-mb*y[8];
    f[4] = A*a*a*M*(3*y[5]-1)*TMath::Cos(y[9]-omega*t)/(y[8]*(i1+i3/2))-2*i1*y[4]*y[3]/(y[8]*(i1+i3/2))-i3*(Omega-y[4])/(2*(i1+i3/2));
    f[5] = y[0];
    f[6] = y[1];
    f[7] = y[2];
    f[8] = y[3];
    f[9] = y[4];

    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy,
    double dfdt[], void *params)
{
    (void)(t);

    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 10, 10);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, y[1]);
    gsl_matrix_set (m, 0, 2, 0.0);
    gsl_matrix_set (m, 0, 3, 0.0);
    gsl_matrix_set (m, 0, 4, 0.0);
    gsl_matrix_set (m, 0, 5, 0.0);
    gsl_matrix_set (m, 0, 6, -9*A*TMath::Cos(y[7]-omega*t));
    gsl_matrix_set (m, 0, 7, 9*A*y[6]*TMath::Sin(y[7]-omega*t));
    gsl_matrix_set (m, 0, 8, 3*A*TMath::Sin(y[9]-omega*t));
    gsl_matrix_set (m, 0, 9, 3*A*y[8]*TMath::Cos(y[9]-omega*t));
    gsl_matrix_set (m, 1, 0, 0.0);
    gsl_matrix_set (m, 1, 1, 0.0);
    gsl_matrix_set (m, 1, 2, 2*y[6]*y[2]);
    gsl_matrix_set (m, 1, 3, 0.0);
    gsl_matrix_set (m, 1, 4, 0.0);
    gsl_matrix_set (m, 1, 5, -9*A*TMath::Cos(y[7]-omega*t));
    gsl_matrix_set (m, 1, 6, -g/a);
    gsl_matrix_set (m, 1, 7, -3*A*(1-3*y[5])*TMath::Sin(y[7]-omega*t));
    gsl_matrix_set (m, 1, 8, 0.0);
    gsl_matrix_set (m, 1, 9, 0.0);
    gsl_matrix_set (m, 2, 0, 0.0);
    gsl_matrix_set (m, 2, 1, -2*y[2]/y[6]);
    gsl_matrix_set (m, 2, 2, -2*y[1]/y[6]);
    gsl_matrix_set (m, 2, 3, 0.0);
    gsl_matrix_set (m, 2, 4, 0.0);
    gsl_matrix_set (m, 2, 5, 9*A*TMath::Sin(y[7]-omega*t)/y[6]);
    gsl_matrix_set (m, 2, 6, 3*A*(1-3*y[5])*TMath::Sin(y[7]-omega*t)/(y[6]*y[6])+2*y[2]*y[1]/(y[6]*y[6]));
    gsl_matrix_set (m, 2, 7, 0.0);
    gsl_matrix_set (m, 2, 8, 0.0);
    gsl_matrix_set (m, 2, 9, 0.0);
    gsl_matrix_set (m, 3, 0, 0.0);
    gsl_matrix_set (m, 3, 1, 0.0);
    gsl_matrix_set (m, 3, 2, 0.0);
    gsl_matrix_set (m, 3, 3, 0.0);
    gsl_matrix_set (m, 3, 4, 0.0);
    gsl_matrix_set (m, 3, 5, 0.0);
    gsl_matrix_set (m, 3, 6, 0.0);
    gsl_matrix_set (m, 3, 7, 0.0);
    gsl_matrix_set (m, 3, 8, 0.0);
    gsl_matrix_set (m, 3, 9, 0.0);
    gsl_matrix_set (m, 4, 0, 0.0);
    gsl_matrix_set (m, 4, 1, 0.0);
    gsl_matrix_set (m, 4, 2, 0.0);
    gsl_matrix_set (m, 4, 3, 0.0);
    gsl_matrix_set (m, 4, 4, 0.0);
    gsl_matrix_set (m, 4, 5, 0.0);
    gsl_matrix_set (m, 4, 6, 0.0);
    gsl_matrix_set (m, 4, 7, 0.0);
    gsl_matrix_set (m, 4, 8, 0.0);
    gsl_matrix_set (m, 4, 9, 0.0);
    gsl_matrix_set (m, 5, 0, 1.0);
    gsl_matrix_set (m, 5, 1, 0.0);
    gsl_matrix_set (m, 5, 2, 0.0);
    gsl_matrix_set (m, 5, 3, 0.0);
    gsl_matrix_set (m, 5, 4, 0.0);
    gsl_matrix_set (m, 5, 5, 0.0);
    gsl_matrix_set (m, 5, 6, 0.0);
    gsl_matrix_set (m, 5, 7, 0.0);
    gsl_matrix_set (m, 5, 8, 0.0);
    gsl_matrix_set (m, 5, 9, 0.0);
    gsl_matrix_set (m, 6, 0, 0.0);
    gsl_matrix_set (m, 6, 1, 1.0);
    gsl_matrix_set (m, 6, 2, 0.0);
    gsl_matrix_set (m, 6, 3, 0.0);
    gsl_matrix_set (m, 6, 4, 0.0);
    gsl_matrix_set (m, 6, 5, 0.0);
    gsl_matrix_set (m, 6, 6, 0.0);
    gsl_matrix_set (m, 6, 7, 0.0);
    gsl_matrix_set (m, 6, 8, 0.0);
    gsl_matrix_set (m, 6, 9, 0.0);
    gsl_matrix_set (m, 7, 0, 0.0);
    gsl_matrix_set (m, 7, 1, 0.0);
    gsl_matrix_set (m, 7, 2, 1.0);
    gsl_matrix_set (m, 7, 3, 0.0);
    gsl_matrix_set (m, 7, 4, 0.0);
    gsl_matrix_set (m, 7, 5, 0.0);
    gsl_matrix_set (m, 7, 6, 0.0);
    gsl_matrix_set (m, 7, 7, 0.0);
    gsl_matrix_set (m, 7, 8, 0.0);
    gsl_matrix_set (m, 7, 9, 0.0);
    gsl_matrix_set (m, 8, 0, 0.0);
    gsl_matrix_set (m, 8, 1, 0.0);
    gsl_matrix_set (m, 8, 2, 0.0);
    gsl_matrix_set (m, 8, 3, 1.0);
    gsl_matrix_set (m, 8, 4, 0.0);
    gsl_matrix_set (m, 8, 5, 0.0);
    gsl_matrix_set (m, 8, 6, 0.0);
    gsl_matrix_set (m, 8, 7, 0.0);
    gsl_matrix_set (m, 8, 8, 0.0);
    gsl_matrix_set (m, 8, 9, 0.0);
    gsl_matrix_set (m, 9, 0, 0.0);
    gsl_matrix_set (m, 9, 1, 0.0);
    gsl_matrix_set (m, 9, 2, 0.0);
    gsl_matrix_set (m, 9, 3, 0.0);
    gsl_matrix_set (m, 9, 4, 1.0);
    gsl_matrix_set (m, 9, 5, 0.0);
    gsl_matrix_set (m, 9, 6, 0.0);
    gsl_matrix_set (m, 9, 7, 0.0);
    gsl_matrix_set (m, 9, 8, 0.0);
    gsl_matrix_set (m, 9, 9, 0.0);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;
    dfdt[4] = 0.0;
    dfdt[5] = 0.0;
    dfdt[6] = 0.0;
    dfdt[7] = 0.0;
    dfdt[8] = 0.0;
    dfdt[9] = 0.0;

    return GSL_SUCCESS;
}

int main (int argc, char** argv)
{

    TCanvas *Princ = new TCanvas();

    Princ->SetTickx();
    Princ->SetTicky();
    Princ->SetGridx();
    Princ->SetGridy();

    TGraph *Graph = new TGraph();

    int contpt = 0;

    gsl_odeiv2_system sys = {func, jac, 10, 0};

    gsl_odeiv2_driver * d =
        gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    
    double t = 0.0, t1 = 10.0;

    double y[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, TMath::Pi()/100, 0.0, TMath::Pi()/10, 0.0};

    printf ("Tempo  Epsilon  Theta  Theta'\n");
    printf ("%.5e %.5e %.5e %.5e\n", t, y[5], y[6], y[8]);

    for (int i = 1; i <= 1000; i++)
        {
        double ti = i * t1 / 1000.0;

        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS)
            {
            printf ("error, return value = %d\n", status);
            break;
            }
        printf ("%.5e %.5e %.5e %.5e\n", t, y[5], y[6], y[8]);

        contpt = Graph->GetN();

        Graph->SetPoint(contpt, t, y[5]);

        }

    Graph->SetMarkerStyle(kFullCircle);
    Graph->SetMarkerColor(kBlack);
    Graph->GetXaxis()->SetLabelFont(132);
    Graph->GetXaxis()->SetTitleFont(132);
    Graph->GetYaxis()->SetLabelFont(132);
    Graph->GetYaxis()->SetTitleFont(132);
    Graph->GetYaxis()->SetLabelSize(0.035);
    Graph->GetYaxis()->SetTitleSize(0.035);
    Graph->GetXaxis()->SetLabelSize(0.035);
    Graph->GetXaxis()->SetTitleSize(0.035);
    Graph->SetMarkerSize(0.7);

    Graph->GetXaxis()->SetLimits(0, t1);
    Graph->GetYaxis()->SetMaxDigits(2);
    Princ->cd();
    Graph->Draw("AP");

    Princ->Print("gp.pdf");

    gsl_odeiv2_driver_free (d);
    return 0;
}