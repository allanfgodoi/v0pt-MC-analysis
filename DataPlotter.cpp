#include <iostream>
#include <type_traits>
#include <iostream>
#include <fstream>
#include <cmath>
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TPad.h"

// TGraph creator function
TGraph* create_TGraph(int nPoints, const float* x, const float* y, const char* title, float xmin, int xmax, float ymin, float ymax, int style, int color, float size){
    TGraph* g = new TGraph(nPoints, x, y);
    g->SetTitle(title);
    g->GetXaxis()->SetLimits(xmin, xmax); // X axis range
    g->GetXaxis()->CenterTitle(); // Center axis label
    g->GetYaxis()->CenterTitle(); // Center axis label
    g->SetMinimum(ymin); // Y axis range
    g->SetMaximum(ymax); // Y axis range
    g->SetMarkerStyle(style);
    g->SetMarkerColor(color);
    g->SetMarkerSize(size);
    return g;
}

void customize_TGraph(TGraph *g, const char* title, float xmin, float xmax, float ymin, float ymax, int style, int color, float size){
    g->SetTitle(title);
    g->GetXaxis()->SetLimits(xmin, xmax); // X axis range
    g->GetXaxis()->CenterTitle(); // Center axis label
    g->GetYaxis()->CenterTitle(); // Center axis label
    g->SetMinimum(ymin); // Y axis range
    g->SetMaximum(ymax); // Y axis range
    g->SetMarkerStyle(style);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerSize(size);
}

void customize_TGraphAsymmErrors(TGraphAsymmErrors *g, const char* title, float xmin, float xmax, float ymin, float ymax, int style, int color, float size){
    g->SetTitle(title);
    g->GetXaxis()->SetLimits(xmin, xmax); // X axis range
    g->GetXaxis()->CenterTitle(); // Center axis label
    g->GetYaxis()->CenterTitle(); // Center axis label
    g->SetMinimum(ymin); // Y axis range
    g->SetMaximum(ymax); // Y axis range
    g->SetMarkerStyle(style);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerSize(size);
    g->SetFillColorAlpha(color, 0.3);
    g->SetFillStyle(1001);
}

void customize_TGraphErrors(TGraphErrors *g, const char* title, float xmin, float xmax, float ymin, float ymax, int style, int color, float size){
    g->SetTitle(title);
    g->GetXaxis()->SetLimits(xmin, xmax); // X axis range
    g->GetXaxis()->CenterTitle(); // Center axis label
    g->GetYaxis()->CenterTitle(); // Center axis label
    g->SetMinimum(ymin); // Y axis range
    g->SetMaximum(ymax); // Y axis range
    g->SetMarkerStyle(style);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerSize(size);
    g->SetFillColorAlpha(color, 0.3);
    g->SetFillStyle(1001);
}

pair<TGraphAsymmErrors*, TGraphAsymmErrors*> csvToTGraph(TString fipflath){
    ifstream f(fipflath);

    TGraphAsymmErrors *g_stat = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_syst = new TGraphAsymmErrors();

    double x, y, ystat_low, ystat_high, ysyst_low, ysyst_high, xsyst_low, xsyst_high;
    char c;
    int i = 0;

    while (f >> x >> c >> y >> c >> ystat_high >> c >> ystat_low >> c >> ysyst_high >> c >> ysyst_low){
        g_stat->SetPoint(i, x, y);
        g_syst->SetPoint(i, x, y);

        g_stat->SetPointError(i, 0, 0, abs(ystat_low), abs(ystat_high));
        if (fipflath.BeginsWith("./Data/Figures_ATLAS/v0_")) g_syst->SetPointError(i, 0.3, 0.3, abs(ysyst_low), abs(ysyst_high));
        else g_syst->SetPointError(i, (x-(x/1.05)), ((x*1.05)-x), abs(ysyst_low), abs(ysyst_high)); // log scale in x axis

        i++;
    }

    f.close();
    return {g_stat, g_syst};
}

void DoPlotMain(){
    
    TFile *f = TFile::Open("./data_v0pt_MC.root", "READ");

    // v0(pT)
    TGraph *gr_v0pt_pion = (TGraph*)f->Get("v0pt_ptref_pion");
    TGraph *gr_v0pt_kaon = (TGraph*)f->Get("v0pt_ptref_kaon"); 
    TGraph *gr_v0pt_proton = (TGraph*)f->Get("v0pt_ptref_proton"); 

    customize_TGraph(gr_v0pt_pion, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.3, 0.9, 20, kRed, 0.8);
    customize_TGraph(gr_v0pt_kaon, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.3, 0.9, 21, kBlue, 0.8);
    customize_TGraph(gr_v0pt_proton, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.46, 22, kGreen+2, 0.8);

    auto c = new TCanvas("c", "c", 1100, 500);
    c->Divide(2, 1);

    // v0(pT) legend
    auto legend_v0pt_text = new TLegend(0.025, 0.98, 0.5, 0.82);
    legend_v0pt_text->SetTextSize(0.055);
    legend_v0pt_text->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    legend_v0pt_text->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV,   #eta_{gap} = 1", "");
    legend_v0pt_text->SetBorderSize(0);
    legend_v0pt_text->SetFillStyle(0);

    auto legend_v0pt_cents = new TLegend(0.103, 0.81, 0.5, 0.61);
    legend_v0pt_cents->SetTextSize(0.034);
    legend_v0pt_cents->AddEntry(gr_v0pt_pion, "#pi^{#pm}", "pfl");
    legend_v0pt_cents->AddEntry(gr_v0pt_kaon, "K^{#pm}", "pfl");
    legend_v0pt_cents->AddEntry(gr_v0pt_proton, "p,#bar{p}", "pfl");
    legend_v0pt_cents->SetBorderSize(0);
    legend_v0pt_cents->SetFillStyle(0);

    auto legend_v0pt_label = new TLegend(0.78, 0.9, 0.93, 0.93);
    legend_v0pt_label->SetTextSize(0.055);
    legend_v0pt_label->AddEntry((TObject*)0, "(a)", "");
    legend_v0pt_label->SetBorderSize(0);
    legend_v0pt_label->SetFillStyle(0);

    // Drawing v0(pT) plot
    c->cd(1);
    gr_v0pt_pion->Draw("AP");
    gr_v0pt_kaon->Draw("P SAME");
    gr_v0pt_proton->Draw("P SAME");
    legend_v0pt_text->Draw();
    legend_v0pt_cents->Draw();
    legend_v0pt_label->Draw();
    gPad->SetLogx();
    gPad->SetTopMargin(0.01);

    // sv0(pT)
    TGraph *gr_sv0pt_pion = (TGraph*)f->Get("sv0pt_ptref_pion");
    TGraph *gr_sv0pt_kaon = (TGraph*)f->Get("sv0pt_ptref_kaon"); 
    TGraph *gr_sv0pt_proton = (TGraph*)f->Get("sv0pt_ptref_proton"); 

    customize_TGraph(gr_sv0pt_pion, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 15.0, -15.0, 55.0, 20, kRed, 0.8);
    customize_TGraph(gr_sv0pt_kaon, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 10.0, -15.0, 55.0, 21, kBlue, 0.8);
    customize_TGraph(gr_sv0pt_proton, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 10.0, -15.0, 55.0, 22, kGreen+2, 0.8);

    // sv0(pT) legend
    auto legend_sv0pt_text = new TLegend(0.025, 0.98, 0.5, 0.82);
    legend_sv0pt_text->SetTextSize(0.055);
    legend_sv0pt_text->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    legend_sv0pt_text->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV,   #eta_{gap} = 1", "");
    legend_sv0pt_text->SetBorderSize(0);
    legend_sv0pt_text->SetFillStyle(0);

    auto legend_sv0pt_cents = new TLegend(0.103, 0.81, 0.5, 0.61);
    legend_sv0pt_cents->SetTextSize(0.034);
    legend_sv0pt_cents->AddEntry(gr_sv0pt_pion, "#pi^{#pm}", "pfl");
    legend_sv0pt_cents->AddEntry(gr_sv0pt_kaon, "K^{#pm}", "pfl");
    legend_sv0pt_cents->AddEntry(gr_sv0pt_proton, "p,#bar{p}", "pfl");
    legend_sv0pt_cents->SetBorderSize(0);
    legend_sv0pt_cents->SetFillStyle(0);

    auto legend_sv0pt_label = new TLegend(0.78, 0.9, 0.93, 0.93);
    legend_sv0pt_label->SetTextSize(0.055);
    legend_sv0pt_label->AddEntry((TObject*)0, "(a)", "");
    legend_sv0pt_label->SetBorderSize(0);
    legend_sv0pt_label->SetFillStyle(0);

    // Drawing v0(pT) plot
    c->cd(2);
    gr_sv0pt_pion->Draw("AP");
    gr_sv0pt_kaon->Draw("P SAME");
    gr_sv0pt_proton->Draw("P SAME");
    legend_sv0pt_text->Draw();
    legend_sv0pt_cents->Draw();
    legend_sv0pt_label->Draw();
    gPad->SetLogx();
    gPad->SetTopMargin(0.01);

    // Saving canvas as pdf
    c->Update();
    c->SaveAs("./PlotMain.pdf");
    delete c;
}