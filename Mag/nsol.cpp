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
    double omega;
    double g;
    double a;
    double M;
    double i1;
    double i3;
    double A;
    double mb;
    double Omega;
};

int func (double t, const double y[], double f[],
    void *params)
{
    (void)(t);

    struct p_type *npar = (struct p_type *)params;

    double omega = npar->omega;
    double g = npar->g;
    double a = npar->a;
    double M = npar->M;
    double i1 = npar->i1;
    double i3 = npar->i3;
    double A = npar->A;
    double mb = npar->mb;
    double Omega = npar->Omega;

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

    struct p_type *npar = (struct p_type *)params;

    double omega = npar->omega;
    double g = npar->g;
    double a = npar->a;
    double M = npar->M;
    double i1 = npar->i1;
    double i3 = npar->i3;
    double A = npar->A;
    double mb = npar->mb;
    double Omega = npar->Omega;

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
    gsl_matrix_set (m, 2, 7, -3*A*(1-3*y[5])*TMath::Cos(y[7]-omega*t)/y[6]);
    gsl_matrix_set (m, 2, 8, 0.0);
    gsl_matrix_set (m, 2, 9, 0.0);
    gsl_matrix_set (m, 3, 0, 0.0);
    gsl_matrix_set (m, 3, 1, 0.0);
    gsl_matrix_set (m, 3, 2, 0.0);
    gsl_matrix_set (m, 3, 3, -i3*y[4]*y[8]/i1);
    gsl_matrix_set (m, 3, 4, 2*y[4]*y[8]-i3*y[3]*y[8]-2*i3*y[8]*y[4]/i1);
    gsl_matrix_set (m, 3, 5, 3*A*M*a*a*TMath::Sin(y[9]-omega*t)/i1);
    gsl_matrix_set (m, 3, 6, 0.0);
    gsl_matrix_set (m, 3, 7, 0.0);
    gsl_matrix_set (m, 3, 8, y[4]*y[4]-i3*y[4]*y[4]/i1-i3*y[4]*y[3]/i1-mb);
    gsl_matrix_set (m, 3, 9, -A*M*a*a*(1-3*y[5])*TMath::Cos(y[9]-omega*t)/i1);
    gsl_matrix_set (m, 4, 0, 0.0);
    gsl_matrix_set (m, 4, 1, 0.0);
    gsl_matrix_set (m, 4, 2, 0.0);
    gsl_matrix_set (m, 4, 3, -2*i1*y[4]/(y[8]*(i1+i3/2)));
    gsl_matrix_set (m, 4, 4, -2*i1*y[3]/(y[8]*(i1+i3/2))+i3/(2*i1+i3));
    gsl_matrix_set (m, 4, 5, A*M*a*a*3*TMath::Cos(y[9]-omega*t)/(y[8]*(i1+i3/2)));
    gsl_matrix_set (m, 4, 6, 0.0);
    gsl_matrix_set (m, 4, 7, 0.0);
    gsl_matrix_set (m, 4, 8, -A*M*a*a*(3*y[5]-1)*TMath::Cos(y[9]-omega*t)/(y[8]*y[8]*(i1+i3/2))+2*i1*y[4]*y[3]/(y[8]*y[8]*(i1+i3/2)));
    gsl_matrix_set (m, 4, 9, -A*M*a*a*(3*y[5]-1)*TMath::Sin(y[9]-omega*t)/(y[8]*(i1+i3/2)));
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

    dfdt[0] = -omega*9*A*y[6]*TMath::Sin(y[7]-omega*t)-3*omega*A*y[8]*TMath::Cos(y[9]-omega*t);
    dfdt[1] = -3*omega*A*(3*y[5]-1)*TMath::Sin(y[7]-omega*t);
    dfdt[2] = -3*omega*A*(1-3*y[5])*TMath::Cos(y[7]-omega*t)/y[6];
    dfdt[3] = A*M*a*a*(1-3*y[5])*TMath::Cos(y[9]-omega*t)*omega/i1;
    dfdt[4] = A*M*a*a*(3*y[5]-1)*TMath::Sin(y[9]-omega*t)*omega/(y[8]*(i1+i3/2));
    dfdt[5] = 0.0;
    dfdt[6] = 0.0;
    dfdt[7] = 0.0;
    dfdt[8] = 0.0;
    dfdt[9] = 0.0;

    return GSL_SUCCESS;
}

int main (int argc, char** argv)
{

    double t0 = 0.0, t1 = 10.0, O0 = 0.0, O1 = 100.0, y[10], ti;
    int divtemp = 1000, divO = 1, contpt, status;

    std::string nome;
    char * nomef;

    std::ofstream ofs;
    ofs.open("Data.tsv", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    std::ofstream eFile ("Data.tsv" ,std::ofstream::app);
    
    struct p_type parametros = {7300*2*TMath::Pi()/60, // omega
                                    9.8, // g
                                    0.02, // a
                                    7.5e-3, // M
                                    7.5e-3 * 0.01 * 0.01/6, // i1
                                    7.5e-3 * 0.01 * 0.01/6, // i3
                                    235, // A
                                    0.0, // m'B
                                    0.0 // Omega
    }; 

    const double initval[10] = {0.0, // epsilon ponto
                        0.0, // theta ponto
                        0.0, //parametros.omega/5, // phi ponto
                        0.0, // theta' ponto
                        0.0, //parametros.omega/5, // phi' ponto
                        0.0, // epsilon
                        TMath::Pi()/100, // theta
                        0.0, // phi
                        TMath::Pi()/10, // theta' 
                        0.0 // phi'
    };

    gsl_odeiv2_system sys = {func, jac, 10, &parametros};   

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
    

    for(int k = 0; k < divO; k++)
    {

        t0 = 0.0, t1 = 10.0;
        
        EpGvec.push_back(new TGraph);
        TGvec.push_back(new TGraph);
        TpGvec.push_back(new TGraph);

        printf("Iteracao %i\n", k+1);

        contpt = 0;

        for (int l = 0; l< 10; l++){y[l]=initval[l];};

        parametros.Omega = (double)(O0+k*O1/divO);

        gsl_odeiv2_driver * d=
            gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-6, 1e-6, 0.0);

        for (int i = 0; i < divtemp; i++)
        {
            ti = i * t1 / (double )divtemp;

            status = gsl_odeiv2_driver_apply (d, &t0, ti, y);

            if (status != GSL_SUCCESS)
            {
                printf ("error, return value = %d\n", status);
                break;
            }

            // printf ("%.5e %.5e %.5e %.5e\n", t0, y[5], y[6], y[8]);

            eFile << t0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t"
                << y[4] << "\t" << y[5] << "\t" << y[6] << "\t" << y[7] << "\t" << y[8] << "\t" << y[9] << "\n";

            contpt = EpGvec[k]->GetN();

            EpGvec[k]->SetPoint(contpt, t0, y[5]);
            TGvec[k]->SetPoint(contpt, t0, y[6]);
            TpGvec[k]->SetPoint(contpt, t0, y[8]);

        }

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

        nome = "Omega ";
        nome = nome + (int)parametros.Omega;
        nomef = &nome[0];

        legend->AddEntry(EpGvec[k], nomef, "p");
    
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
    EpMg->GetYaxis()->SetTitle("Epsilon #bf{[m]}");
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