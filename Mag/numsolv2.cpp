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
const double sinbeta = 0.1;
const double cosalpha = 0.5;
const double g = 9.8;
const double a = 0.02;
const double M = 7.5e-3;
const double i1 = M * 0.01* 0.01/6;
const double i3 = i1;
const double A = 235;
const double mb = 1;
// bool vdd =true;


int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t); /* avoid unused parameter warning */
    f[0] = g/a + 3*A*sinbeta*TMath::Sin(y[3]-omega*t);
    f[1] = A*M*a*a/(i1*sinbeta) * (3*y[2]-1)*TMath::Cos(y[3]-omega*t);
    f[2] = y[0];
    f[3] = y[1];
    // std::cout << "\n\n FUNCIONANDO \n\n";
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy,
    double dfdt[], void *params)
{
    (void)(t); /* avoid unused parameter warning */
    // std::cout << "\n\n TESTE\n\n";
    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 4, 4);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 0.0);
    gsl_matrix_set (m, 0, 2, 0.0);
    gsl_matrix_set (m, 0, 3, 3*A*sinbeta*TMath::Cos(y[3]-omega*t));
    gsl_matrix_set (m, 1, 0, 0.0);
    gsl_matrix_set (m, 1, 1, 0.0);
    gsl_matrix_set (m, 1, 2, 3*A*M*a*a/(i1*sinbeta)*TMath::Cos(y[3]-omega*t));
    gsl_matrix_set (m, 1, 3, -A*M*a*a/(i1*sinbeta) * (3*y[2]-1)*TMath::Sin(y[3]-omega*t));
    gsl_matrix_set (m, 2, 0, 1.0);
    gsl_matrix_set (m, 2, 1, 0.0);
    gsl_matrix_set (m, 2, 2, 0.0);
    gsl_matrix_set (m, 2, 3, 0.0);
    gsl_matrix_set (m, 3, 0, 0.0);
    gsl_matrix_set (m, 3, 1, 1);
    gsl_matrix_set (m, 3, 2, 0.0);
    gsl_matrix_set (m, 3, 3, 0.0);
    // if (vdd)
    // {
    //     std::cout << "\n\nPRINTADO MATRIX JACOBIANA\n\n";
    //     for(int i = 0; i< 8; i++){
    //         for(int j = 0;j<8;j++){
    //             std::cout << gsl_matrix_get (m,i,j) << "\t";
    //         }
    //         std::cout << "\n";
    //     }
    //     std::cout << "\n\n";
    //     vdd = false;
    // }
    dfdt[0] = -omega*3*A*sinbeta*TMath::Cos(y[3]-omega*t);
    dfdt[1] = omega*A*M*a*a/(i1*sinbeta) * (3*y[2]-1)*TMath::Sin(y[3]-omega*t);
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;
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

    gsl_odeiv2_system sys = {func, jac, 4, 0};

    gsl_odeiv2_driver * d =
        gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                    1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 10.0;
    double y[4] = { 0.0, omega/10, 0.0, 0.0};

    for (i = 1; i <= 100; i++)
        {
        double ti = i * t1 / 100.0;
        // std::cout << "\nChamando driver\n";
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        // std::cout << "\nDriver Executado\n";
        if (status != GSL_SUCCESS)
            {
            printf ("error, return value=%d\n", status);
            break;
            }

        printf ("%.5e %.5e %.5e\n", t, y[2], y[3]);

        contpt = Graph->GetN();

        Graph->SetPoint(contpt, t, y[2]);

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