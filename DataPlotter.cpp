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

void customize_TGraph(TGraphErrors *g, const char* title, float xmin, float xmax, float ymin, float ymax, int style, int color, float size){
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

void DoPlotMain(TString Filename, TString Savename){
    
    TFile *f = TFile::Open(Filename, "READ");

    // v0(pT)
    TGraph *gr_v0pt_pion = (TGraph*)f->Get("v0pt_ptref_1_pion");
    TGraph *gr_v0pt_kaon = (TGraph*)f->Get("v0pt_ptref_1_kaon"); 
    TGraph *gr_v0pt_proton = (TGraph*)f->Get("v0pt_ptref_1_proton");
    TGraph *gr_v0pt_all = (TGraph*)f->Get("v0pt_ptref_1_all"); 

    customize_TGraph(gr_v0pt_pion, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.9, 21, kRed, 0.8);
    customize_TGraph(gr_v0pt_kaon, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.9, 22, kBlue, 0.8);
    customize_TGraph(gr_v0pt_proton, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.9, 34, kGreen+2, 0.8);
    customize_TGraph(gr_v0pt_all, "; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.9, 20, kBlack, 0.8);

    auto c = new TCanvas("c", "c", 1100, 500);
    c->Divide(2, 1);

    // v0(pT) legend
    auto legend_v0pt_text = new TLegend(0.025, 0.98, 0.5, 0.82);
    legend_v0pt_text->SetTextSize(0.055);
    legend_v0pt_text->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    legend_v0pt_text->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV,   #eta_{gap} = 1", "");
    legend_v0pt_text->SetBorderSize(0);
    legend_v0pt_text->SetFillStyle(0);

    auto legend_v0pt_cents = new TLegend(0.103, 0.81, 0.5, 0.61);
    legend_v0pt_cents->SetTextSize(0.034);
    legend_v0pt_cents->AddEntry(gr_v0pt_pion, "#pi^{#pm}", "pfl");
    legend_v0pt_cents->AddEntry(gr_v0pt_kaon, "K^{#pm}", "pfl");
    legend_v0pt_cents->AddEntry(gr_v0pt_proton, "p,#bar{p}", "pfl");
    legend_v0pt_cents->AddEntry(gr_v0pt_all, "#pi^{#pm},K^{#pm},p,#bar{p}", "pfl");
    legend_v0pt_cents->SetBorderSize(0);
    legend_v0pt_cents->SetFillStyle(0);

    auto legend_v0pt_label = new TLegend(0.78, 0.9, 0.93, 0.93);
    legend_v0pt_label->SetTextSize(0.055);
    legend_v0pt_label->AddEntry((TObject*)0, "(a)", "");
    legend_v0pt_label->SetBorderSize(0);
    legend_v0pt_label->SetFillStyle(0);

    // Drawing v0(pT) plot
    c->cd(1);
    gr_v0pt_all->Draw("AP");
    gr_v0pt_pion->Draw("P SAME");
    gr_v0pt_kaon->Draw("P SAME");
    gr_v0pt_proton->Draw("P SAME");
    legend_v0pt_text->Draw();
    legend_v0pt_cents->Draw();
    legend_v0pt_label->Draw();
    gPad->SetLogx();
    gPad->SetTopMargin(0.01);

    // sv0(pT)
    TGraph *gr_sv0pt_pion = (TGraph*)f->Get("sv0pt_ptref_1_pion");
    TGraph *gr_sv0pt_kaon = (TGraph*)f->Get("sv0pt_ptref_1_kaon"); 
    TGraph *gr_sv0pt_proton = (TGraph*)f->Get("sv0pt_ptref_1_proton");
    TGraph *gr_sv0pt_all = (TGraph*)f->Get("sv0pt_ptref_1_all");

    customize_TGraph(gr_sv0pt_pion, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 15.0, 0.1, 10.0, 21, kRed, 0.8);
    customize_TGraph(gr_sv0pt_kaon, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 10.0, 0.1, 10.0, 22, kBlue, 0.8);
    customize_TGraph(gr_sv0pt_proton, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 10.0, 0.1, 10.0, 34, kGreen+2, 0.8);
    customize_TGraph(gr_sv0pt_all, "; p_{T} [GeV]; v_{0}(p_{T})/v_{0}", 0.0, 10.0, 0.1, 10.0, 20, kBlack, 0.8);

    // sv0(pT) legend
    auto legend_sv0pt_text = new TLegend(0.025, 0.98, 0.5, 0.82);
    legend_sv0pt_text->SetTextSize(0.055);
    legend_sv0pt_text->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    legend_sv0pt_text->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV,   #eta_{gap} = 1", "");
    legend_sv0pt_text->SetBorderSize(0);
    legend_sv0pt_text->SetFillStyle(0);

    auto legend_sv0pt_cents = new TLegend(0.103, 0.81, 0.5, 0.61);
    legend_sv0pt_cents->SetTextSize(0.034);
    legend_sv0pt_cents->AddEntry(gr_sv0pt_pion, "#pi^{#pm}", "pfl");
    legend_sv0pt_cents->AddEntry(gr_sv0pt_kaon, "K^{#pm}", "pfl");
    legend_sv0pt_cents->AddEntry(gr_sv0pt_proton, "p,#bar{p}", "pfl");
    legend_sv0pt_cents->AddEntry(gr_sv0pt_all, "#pi^{#pm},K^{#pm},p,#bar{p}", "pfl");
    legend_sv0pt_cents->SetBorderSize(0);
    legend_sv0pt_cents->SetFillStyle(0);

    auto legend_sv0pt_label = new TLegend(0.78, 0.9, 0.93, 0.93);
    legend_sv0pt_label->SetTextSize(0.055);
    legend_sv0pt_label->AddEntry((TObject*)0, "(b)", "");
    legend_sv0pt_label->SetBorderSize(0);
    legend_sv0pt_label->SetFillStyle(0);

    // Drawing v0(pT) plot
    c->cd(2);
    gr_sv0pt_pion->Draw("AP");
    gr_sv0pt_kaon->Draw("P SAME");
    gr_sv0pt_proton->Draw("P SAME");
    gr_sv0pt_all->Draw("P SAME");
    legend_sv0pt_text->Draw();
    legend_sv0pt_cents->Draw();
    legend_sv0pt_label->Draw();
    gPad->SetLogx();
    gPad->SetTopMargin(0.01);

    // Saving canvas as pdf
    c->Update();
    c->SaveAs(Savename);
    delete c;
}

void DoPlotRefs(TString Filename, TString Savename){

    auto c = new TCanvas("c_ptref", "c_ptref", 1100, 1000);
    c->Divide(2, 2);

    TFile *f = TFile::Open(Filename, "READ");

    // READ TGRAPHS

    // Pion
    TGraph *gr_1_pion = (TGraph*)f->Get("v0pt_ptref_1_pion");
    TGraph *gr_2_pion = (TGraph*)f->Get("v0pt_ptref_2_pion"); 
    TGraph *gr_3_pion = (TGraph*)f->Get("v0pt_ptref_3_pion"); 

    // Kaon
    TGraph *gr_1_kaon = (TGraph*)f->Get("v0pt_ptref_1_kaon");
    TGraph *gr_2_kaon = (TGraph*)f->Get("v0pt_ptref_2_kaon"); 
    TGraph *gr_3_kaon = (TGraph*)f->Get("v0pt_ptref_3_kaon");

    // Proton
    TGraph *gr_1_proton = (TGraph*)f->Get("v0pt_ptref_1_proton");
    TGraph *gr_2_proton = (TGraph*)f->Get("v0pt_ptref_2_proton"); 
    TGraph *gr_3_proton = (TGraph*)f->Get("v0pt_ptref_3_proton");

    // All
    TGraph *gr_1_all = (TGraph*)f->Get("v0pt_ptref_1_all");
    TGraph *gr_2_all = (TGraph*)f->Get("v0pt_ptref_2_all"); 
    TGraph *gr_3_all = (TGraph*)f->Get("v0pt_ptref_3_all");

    // DO PLOTS

    // Pion
    c->cd(1);
    customize_TGraph(gr_1_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 20, 2, 1.0);
    customize_TGraph(gr_2_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 21, 4, 1.0);
    customize_TGraph(gr_3_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 22, 6, 1.0);

    auto leg_title_pion = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_pion->SetTextSize(0.055);
    leg_title_pion->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    leg_title_pion->AddEntry((TObject*)0, "#pi^{#pm}     #eta_{gap} = 1", "");
    leg_title_pion->SetBorderSize(0);
    leg_title_pion->SetFillStyle(0);

    auto leg_ptref_pion = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_pion->SetTextSize(0.030);
    leg_ptref_pion->AddEntry(gr_1_pion, "0.5-2 GeV", "pfl");
    leg_ptref_pion->AddEntry(gr_2_pion, "0.5-5 GeV", "pfl");
    leg_ptref_pion->AddEntry(gr_3_pion, "1-5 GeV", "pfl");
    leg_ptref_pion->SetBorderSize(0);
    leg_ptref_pion->SetFillStyle(0);

    gr_2_pion->Draw("AP");
    gr_1_pion->Draw("P SAME");
    gr_3_pion->Draw("P SAME");
    leg_title_pion->Draw();
    leg_ptref_pion->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // Kaon
    c->cd(2);
    customize_TGraph(gr_1_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.5, 20, 2, 1.0);
    customize_TGraph(gr_2_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.5, 21, 4, 1.0);
    customize_TGraph(gr_3_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.5, 22, 6, 1.0);

    auto leg_title_kaon = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_kaon->SetTextSize(0.055);
    leg_title_kaon->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    leg_title_kaon->AddEntry((TObject*)0, "K^{#pm}     #eta_{gap} = 1", "");
    leg_title_kaon->SetBorderSize(0);
    leg_title_kaon->SetFillStyle(0);

    auto leg_ptref_kaon = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_kaon->SetTextSize(0.030);
    leg_ptref_kaon->AddEntry(gr_1_kaon, "0.5-2 GeV", "pfl");
    leg_ptref_kaon->AddEntry(gr_2_kaon, "0.5-5 GeV", "pfl");
    leg_ptref_kaon->AddEntry(gr_3_kaon, "1-5 GeV", "pfl");
    leg_ptref_kaon->SetBorderSize(0);
    leg_ptref_kaon->SetFillStyle(0);

    gr_3_kaon->Draw("AP");
    gr_1_kaon->Draw("P SAME");
    gr_2_kaon->Draw("P SAME");
    leg_title_kaon->Draw();
    leg_ptref_kaon->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // Proton
    c->cd(3);
    customize_TGraph(gr_1_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.4, 20, 2, 1.0);
    customize_TGraph(gr_2_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.4, 21, 4, 1.0);
    customize_TGraph(gr_3_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.4, 22, 6, 1.0);

    auto leg_title_proton = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_proton->SetTextSize(0.055);
    leg_title_proton->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    leg_title_proton->AddEntry((TObject*)0, "p,#bar{p}     #eta_{gap} = 1", "");
    leg_title_proton->SetBorderSize(0);
    leg_title_proton->SetFillStyle(0);

    auto leg_ptref_proton = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_proton->SetTextSize(0.030);
    leg_ptref_proton->AddEntry(gr_1_proton, "0.5-2 GeV", "pfl");
    leg_ptref_proton->AddEntry(gr_2_proton, "0.5-5 GeV", "pfl");
    leg_ptref_proton->AddEntry(gr_3_proton, "1-5 GeV", "pfl");
    leg_ptref_proton->SetBorderSize(0);
    leg_ptref_proton->SetFillStyle(0);

    gr_3_proton->Draw("AP");
    gr_1_proton->Draw("P SAME");
    gr_2_proton->Draw("P SAME");
    leg_title_proton->Draw();
    leg_ptref_proton->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // All
    c->cd(4);
    customize_TGraph(gr_1_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.75, 20, 2, 1.0);
    customize_TGraph(gr_2_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.75, 21, 4, 1.0);
    customize_TGraph(gr_3_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.75, 22, 6, 1.0);

    auto leg_title_all = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_all->SetTextSize(0.055);
    leg_title_all->AddEntry((TObject*)0, "He+He 5.36 TeV", "");
    leg_title_all->AddEntry((TObject*)0, "#pi^{#pm},K^{#pm},p,#bar{p}     #eta_{gap} = 1", "");
    leg_title_all->SetBorderSize(0);
    leg_title_all->SetFillStyle(0);

    auto leg_ptref_all = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_all->SetTextSize(0.030);
    leg_ptref_all->AddEntry(gr_1_all, "0.5-2 GeV", "pfl");
    leg_ptref_all->AddEntry(gr_2_all, "0.5-5 GeV", "pfl");
    leg_ptref_all->AddEntry(gr_3_all, "1-5 GeV", "pfl");
    leg_ptref_all->SetBorderSize(0);
    leg_ptref_all->SetFillStyle(0);

    gr_2_all->Draw("AP");
    gr_1_all->Draw("P SAME");
    gr_3_all->Draw("P SAME");
    leg_title_all->Draw();
    leg_ptref_all->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    c->Update();
    c->SaveAs(Savename);
}

void DoPlotEtas(TString Filename, TString Savename){

    auto c = new TCanvas("c_etas", "c_etas", 1100, 1000);
    c->Divide(2, 2);

    TFile *f = TFile::Open(Filename, "READ");

    // READ TGRAPHS

    // Pion
    TGraph *gr_0_pion = (TGraph*)f->Get("v0pt_eta_0_pion");
    TGraph *gr_1_pion = (TGraph*)f->Get("v0pt_ptref_1_pion");
    TGraph *gr_2_pion = (TGraph*)f->Get("v0pt_eta_2_pion"); 
    TGraph *gr_3_pion = (TGraph*)f->Get("v0pt_eta_3_pion"); 

    // Kaon
    TGraph *gr_0_kaon = (TGraph*)f->Get("v0pt_eta_0_kaon"); 
    TGraph *gr_1_kaon = (TGraph*)f->Get("v0pt_ptref_1_kaon");
    TGraph *gr_2_kaon = (TGraph*)f->Get("v0pt_eta_2_kaon"); 
    TGraph *gr_3_kaon = (TGraph*)f->Get("v0pt_eta_3_kaon");

    // Proton
    TGraph *gr_0_proton = (TGraph*)f->Get("v0pt_eta_0_proton"); 
    TGraph *gr_1_proton = (TGraph*)f->Get("v0pt_ptref_1_proton");
    TGraph *gr_2_proton = (TGraph*)f->Get("v0pt_eta_2_proton"); 
    TGraph *gr_3_proton = (TGraph*)f->Get("v0pt_eta_3_proton");

    // All
    TGraph *gr_0_all = (TGraph*)f->Get("v0pt_eta_0_all"); 
    TGraph *gr_1_all = (TGraph*)f->Get("v0pt_ptref_1_all");
    TGraph *gr_2_all = (TGraph*)f->Get("v0pt_eta_2_all"); 
    TGraph *gr_3_all = (TGraph*)f->Get("v0pt_eta_3_all");

    // DO PLOTS

    // Pion
    c->cd(1);
    customize_TGraph(gr_0_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 20, 2, 1.0);
    customize_TGraph(gr_1_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 21, 4, 1.0);
    customize_TGraph(gr_2_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 22, 6, 1.0);
    customize_TGraph(gr_3_pion, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.0, 0.75, 34, 209, 1.0);

    auto leg_title_pion = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_pion->SetTextSize(0.055);
    leg_title_pion->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    leg_title_pion->AddEntry((TObject*)0, "#pi^{#pm}     p_{T}^{ref}: 0.5-2 GeV", "");
    leg_title_pion->SetBorderSize(0);
    leg_title_pion->SetFillStyle(0);

    auto leg_ptref_pion = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_pion->SetTextSize(0.030);
    leg_ptref_pion->AddEntry(gr_0_pion, "#eta_{gap} = 0", "pfl");
    leg_ptref_pion->AddEntry(gr_1_pion, "#eta_{gap} = 1", "pfl");
    leg_ptref_pion->AddEntry(gr_2_pion, "#eta_{gap} = 2", "pfl");
    leg_ptref_pion->AddEntry(gr_3_pion, "#eta_{gap} = 3", "pfl");
    leg_ptref_pion->SetBorderSize(0);
    leg_ptref_pion->SetFillStyle(0);

    gr_1_pion->Draw("AP");
    gr_0_pion->Draw("P SAME");
    gr_3_pion->Draw("P SAME");
    gr_2_pion->Draw("P SAME");
    leg_title_pion->Draw();
    leg_ptref_pion->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // Kaon
    c->cd(2);
    customize_TGraph(gr_0_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.6, 20, 2, 1.0);
    customize_TGraph(gr_1_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.6, 21, 4, 1.0);
    customize_TGraph(gr_2_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.6, 22, 6, 1.0);
    customize_TGraph(gr_3_kaon, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.1, 0.6, 34, 209, 1.0);

    auto leg_title_kaon = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_kaon->SetTextSize(0.055);
    leg_title_kaon->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    leg_title_kaon->AddEntry((TObject*)0, "#K^{#pm}     p_{T}^{ref}: 0.5-2 GeV", "");
    leg_title_kaon->SetBorderSize(0);
    leg_title_kaon->SetFillStyle(0);

    auto leg_ptref_kaon = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_kaon->SetTextSize(0.030);
    leg_ptref_kaon->AddEntry(gr_0_kaon, "#eta_{gap} = 0", "pfl");
    leg_ptref_kaon->AddEntry(gr_1_kaon, "#eta_{gap} = 1", "pfl");
    leg_ptref_kaon->AddEntry(gr_2_kaon, "#eta_{gap} = 2", "pfl");
    leg_ptref_kaon->AddEntry(gr_3_kaon, "#eta_{gap} = 3", "pfl");
    leg_ptref_kaon->SetBorderSize(0);
    leg_ptref_kaon->SetFillStyle(0);

    gr_1_kaon->Draw("AP");
    gr_0_kaon->Draw("P SAME");
    gr_3_kaon->Draw("P SAME");
    gr_2_kaon->Draw("P SAME");
    leg_title_kaon->Draw();
    leg_ptref_kaon->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // Proton
    c->cd(3);
    customize_TGraph(gr_0_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.09, 0.4, 20, 2, 1.0);
    customize_TGraph(gr_1_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.09, 0.4, 21, 4, 1.0);
    customize_TGraph(gr_2_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.09, 0.4, 22, 6, 1.0);
    customize_TGraph(gr_3_proton, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, 0.09, 0.4, 34, 209, 1.0);

    auto leg_title_proton = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_proton->SetTextSize(0.055);
    leg_title_proton->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    leg_title_proton->AddEntry((TObject*)0, "p,#bar{p}     p_{T}^{ref}: 0.5-2 GeV", "");
    leg_title_proton->SetBorderSize(0);
    leg_title_proton->SetFillStyle(0);

    auto leg_ptref_proton = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_proton->SetTextSize(0.030);
    leg_ptref_proton->AddEntry(gr_0_proton, "#eta_{gap} = 0", "pfl");
    leg_ptref_proton->AddEntry(gr_1_proton, "#eta_{gap} = 1", "pfl");
    leg_ptref_proton->AddEntry(gr_2_proton, "#eta_{gap} = 2", "pfl");
    leg_ptref_proton->AddEntry(gr_3_proton, "#eta_{gap} = 3", "pfl");
    leg_ptref_proton->SetBorderSize(0);
    leg_ptref_proton->SetFillStyle(0);

    gr_1_proton->Draw("AP");
    gr_0_proton->Draw("P SAME");
    gr_3_proton->Draw("P SAME");
    gr_2_proton->Draw("P SAME");
    leg_title_proton->Draw();
    leg_ptref_proton->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    // All
    c->cd(4);
    customize_TGraph(gr_0_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.02, 0.75, 20, 2, 1.0);
    customize_TGraph(gr_1_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.02, 0.75, 21, 4, 1.0);
    customize_TGraph(gr_2_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.02, 0.75, 22, 6, 1.0);
    customize_TGraph(gr_3_all, " ; p_{T} [GeV]; v_{0}(p_{T})", 0.0, 10.0, -0.02, 0.75, 34, 209, 1.0);

    auto leg_title_all = new TLegend(0.04, 0.98, 0.5, 0.825);
    leg_title_all->SetTextSize(0.055);
    leg_title_all->AddEntry((TObject*)0, "Ne+Ne 5.36 TeV", "");
    leg_title_all->AddEntry((TObject*)0, "#pi^{#pm},K^{#pm},p,#bar{p}     p_{T}^{ref}: 0.5-2 GeV", "");
    leg_title_all->SetBorderSize(0);
    leg_title_all->SetFillStyle(0);

    auto leg_ptref_all = new TLegend(0.14, 0.82, 0.485, 0.68);
    leg_ptref_all->SetTextSize(0.030);
    leg_ptref_all->AddEntry(gr_0_all, "#eta_{gap} = 0", "pfl");
    leg_ptref_all->AddEntry(gr_1_all, "#eta_{gap} = 1", "pfl");
    leg_ptref_all->AddEntry(gr_2_all, "#eta_{gap} = 2", "pfl");
    leg_ptref_all->AddEntry(gr_3_all, "#eta_{gap} = 3", "pfl");
    leg_ptref_all->SetBorderSize(0);
    leg_ptref_all->SetFillStyle(0);

    gr_0_all->Draw("AP");
    gr_3_all->Draw("P SAME");
    gr_2_all->Draw("P SAME");
    gr_1_all->Draw("P SAME");
    leg_title_all->Draw();
    leg_ptref_all->Draw();
    gPad->SetLogx();
    gPad->SetLeftMargin(0.12);
    gPad->SetTopMargin(0.01);

    c->Update();
    c->SaveAs(Savename);
}
