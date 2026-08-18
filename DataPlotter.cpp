// Use the DoPlotCents function with the output file from RDF_ObsConstructor.cpp
// USAGE EXAMPLE: DoPlotCents("/eos/user/a/afloresg/MC-studies/BB/ObsData.root", "/eos/user/a/afloresg/MC-studies/BB/Plots/PlotMain-BB-0005.pdf", "0005") "0005" is the centrality range
// For the centralities, use: 0005, 0510, 1020, 2030, 3040, 4050, 5060, 6080, 80100

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

void DoPlotSystems(){
    
    const int nSystems = 9;
    const char* input_filename[nSystems] = {"HeHe_data.root", "LiLi_data.root", "BB_data.root", "OO_data.root", "NeNe_data.root", "MgMg_data.root", "ArAr_data.root", "CaCa_data.root", "KrKr_data.root"};
    TString system_name[nSystems] = {"HeHe", "LiLi", "BB", "OO", "NeNe", "MgMg", "ArAr", "CaCa", "KrKr"};

    // Configurações das Centralidades
    const int nCents = 7;
    TString cent_bins[nCents]   = {"0005", "0510", "1020", "2030", "3040", "4050", "6080"};
    TString cent_labels[nCents] = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "60-80%"};
    
    // Arrays para diferenciar visualmente cada centralidade
    int colors[nCents]  = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1};
    int markers[nCents] = {20, 21, 22, 23, 33, 34, 29};

    for (int i = 0; i < nSystems; i++){
        string input_filepath = "./Data/"; input_filepath += input_filename[i];
        TFile *f = TFile::Open(input_filepath.c_str(), "READ");
        
        // Verifica se o arquivo abriu corretamente para evitar falhas silenciosas
        if (!f || f->IsZombie()) {
            printf("Aviso: Nao foi possivel abrir o arquivo %s\n", input_filepath.c_str());
            continue;
        }

        auto c = new TCanvas("c", Form("c_%s", system_name[i].Data()), 1100, 500);
        c->Divide(2, 1);

        // Uso de TMultiGraph para lidar com múltiplos gráficos no mesmo eixo
        TMultiGraph *mg_v0pt = new TMultiGraph();
        mg_v0pt->SetTitle("v_{0}(p_{T}); p_{T}; v_{0}(p_{T})");
        
        TMultiGraph *mg_sv0pt = new TMultiGraph();
        mg_sv0pt->SetTitle("v_{0}(p_{T})/v_{0}; p_{T}; v_{0}(p_{T})/v_{0}");

        // Legenda de parâmetros do sistema (igual a anterior)
        TString legend_text = system_name[i] + " 5.36 GeV";
        auto legend_syst_params = new TLegend(0.025, 0.88, 0.5, 0.76); // Levemente ajustada para caber melhor
        legend_syst_params->SetTextSize(0.045);
        legend_syst_params->AddEntry((TObject*)0, legend_text, "");
        legend_syst_params->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV    #eta_{gap} = 1", "");
        legend_syst_params->SetBorderSize(0);
        legend_syst_params->SetFillStyle(0);

        // Nova legenda para as centralidades
        auto legend_cent = new TLegend(0.125, 0.755, 0.5, 0.385);
        legend_cent->SetBorderSize(0);
        legend_cent->SetFillStyle(0);
        legend_cent->SetTextSize(0.04);

        // Loop sobre as centralidades para buscar e formatar cada TGraph
        for (int j = 0; j < nCents; j++){
            // Formata o nome do objeto (ex: v0pt_0005_all, v0pt_0510_all)
            TString g_name_v0  = Form("v0pt_%s_all",  cent_bins[j].Data());
            TString g_name_sv0 = Form("sv0pt_%s_all", cent_bins[j].Data());

            TGraphErrors *gr_v0pt  = (TGraphErrors*)f->Get(g_name_v0);
            TGraphErrors *gr_sv0pt = (TGraphErrors*)f->Get(g_name_sv0);

            if (gr_v0pt && gr_sv0pt) {
                // Formatação visual v0pt
                gr_v0pt->SetMarkerColor(colors[j]);
                gr_v0pt->SetLineColor(colors[j]);
                gr_v0pt->SetMarkerStyle(markers[j]);
                gr_v0pt->SetMarkerSize(0.6);
                mg_v0pt->Add(gr_v0pt);


                // Formatação visual sv0pt
                gr_sv0pt->SetMarkerColor(colors[j]);
                gr_sv0pt->SetLineColor(colors[j]);
                gr_sv0pt->SetMarkerStyle(markers[j]);
                gr_sv0pt->SetMarkerSize(0.6);
                mg_sv0pt->Add(gr_sv0pt);
                
                // Adiciona na legenda (apenas precisa do primeiro, já que as cores são idênticas)
                legend_cent->AddEntry(gr_v0pt, cent_labels[j], "lp");
            } else {
                printf("  Aviso: Graficos %s ausentes no arquivo %s\n", cent_bins[j].Data(), input_filename[i]);
            }
        }

        // --- Desenho do Pad 1 (v0pt) ---
        c->cd(1);
        gPad->SetLogx();
        mg_v0pt->Draw("AP"); // "A" desenha os eixos para o MultiGraph
        mg_v0pt->GetXaxis()->SetLimits(0.485, 3.6);
        mg_v0pt->SetMinimum(-0.16);
        mg_v0pt->SetMaximum(0.55);
        mg_v0pt->GetXaxis()->SetLabelFont(42);
        //mg_v0pt->GetXaxis()->SetLabelSize(0.05);
        //mg_v0pt->GetXaxis()->SetTitleSize(0.06);
        mg_v0pt->GetXaxis()->SetTitleFont(42);
        mg_v0pt->GetXaxis()->CenterTitle(true);
        mg_v0pt->GetYaxis()->SetLabelFont(42);
        //mg_v0pt->GetYaxis()->SetLabelSize(0.05);
        //mg_v0pt->GetYaxis()->SetTitleSize(0.06);
        mg_v0pt->GetYaxis()->SetTitleOffset(1.10);
        mg_v0pt->GetYaxis()->SetTitleFont(42);
        mg_v0pt->GetYaxis()->CenterTitle(true);
        legend_syst_params->Draw();
        legend_cent->Draw();

        // --- Desenho do Pad 2 (sv0pt) ---
        c->cd(2);
        gPad->SetLogx();
        mg_sv0pt->Draw("AP");
        mg_sv0pt->GetXaxis()->SetLimits(0.485, 4.0);
        mg_sv0pt->SetMinimum(-7.0);
        mg_sv0pt->SetMaximum(30.0);
        mg_sv0pt->GetXaxis()->SetLabelFont(42);
        //mg_sv0pt->GetXaxis()->SetLabelSize(0.05);
        //mg_sv0pt->GetXaxis()->SetTitleSize(0.06);
        mg_sv0pt->GetXaxis()->SetTitleFont(42);
        mg_sv0pt->GetXaxis()->CenterTitle(true);
        mg_sv0pt->GetYaxis()->SetLabelFont(42);
        //mg_sv0pt->GetYaxis()->SetLabelSize(0.05);
        //mg_sv0pt->GetYaxis()->SetTitleSize(0.06);
        mg_sv0pt->GetYaxis()->SetTitleOffset(1.10);
        mg_sv0pt->GetYaxis()->SetTitleFont(42);
        mg_sv0pt->GetYaxis()->CenterTitle(true);
        legend_syst_params->Draw();
        legend_cent->Draw();

        // Salvar e limpar
        TString saveplot_path = "./Plots/" + system_name[i] + ".pdf";
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        //c->Range(0.04000706,-0.4080692,0.7405774,1.786744);
        c->SetFillColor(0);
        c->SetBorderMode(0);
        c->SetBorderSize(0);
        c->SetTickx(1);
        c->SetTicky(1);
        c->SetLeftMargin(0.1650854);
        c->SetRightMargin(0.0341556);
        c->SetTopMargin(0.08508404);
        c->SetBottomMargin(0.1859244);
        c->SetFrameBorderMode(0);
        c->SetFrameBorderMode(0);
        c->Update();
        c->SaveAs(saveplot_path);
        
        delete c; // Limpa o canvas
        f->Close(); // Fecha o arquivo para evitar vazamento de memória (memory leak)
    }
}

void DoPlotCentrality() {
    const int nSystems = 9;
    const char* input_filename[nSystems] = {"HeHe_data.root", "LiLi_data.root", "BB_data.root", "OO_data.root", "NeNe_data.root", "MgMg_data.root", "ArAr_data.root", "CaCa_data.root", "KrKr_data.root"};
    TString system_name[nSystems] = {"HeHe", "LiLi", "BB", "OO", "NeNe", "MgMg", "ArAr", "CaCa", "KrKr"};
    
    // Arrays para diferenciar visualmente cada sistema
    int colors[nSystems] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+1, kOrange+7, kViolet, kAzure+7};
    int markers[nSystems] = {20, 21, 22, 23, 33, 34, 29, 20, 21}; 

    // TMultiGraph resolve o problema de escala automática dos eixos para múltiplos gráficos
    TMultiGraph *mg_v0pt = new TMultiGraph();
    mg_v0pt->SetTitle("v_{0}(p_{T}); p_{T}; v_{0}(p_{T})");
    
    TMultiGraph *mg_sv0pt = new TMultiGraph();
    mg_sv0pt->SetTitle("v_{0}(p_{T})/v_{0}; p_{T}; v_{0}(p_{T})/v_{0}");

    // Legenda para os sistemas (desenhada no segundo painel ou no primeiro, onde preferir)
    auto legend_systems = new TLegend(0.125, 0.72, 0.5, 0.35);
    legend_systems->SetBorderSize(0);
    legend_systems->SetFillStyle(0);
    legend_systems->SetTextSize(0.04);

    for (int i = 0; i < nSystems; i++) {
        string input_filepath = "./Data/"; input_filepath += input_filename[i];
        TFile *f = TFile::Open(input_filepath.c_str(), "READ");
        
        // Verificação de segurança caso o arquivo não exista
        if (!f || f->IsZombie()) continue; 

        TGraphErrors *gr_v0pt = (TGraphErrors*)f->Get("v0pt_0005_all");
        TGraphErrors *gr_sv0pt = (TGraphErrors*)f->Get("sv0pt_0005_all");

        // Aplicando as cores e marcadores
        gr_v0pt->SetMarkerColor(colors[i]);
        gr_v0pt->SetLineColor(colors[i]);
        gr_v0pt->SetMarkerStyle(markers[i]);
        gr_v0pt->SetMarkerSize(0.6);
        
        gr_sv0pt->SetMarkerColor(colors[i]);
        gr_sv0pt->SetLineColor(colors[i]);
        gr_sv0pt->SetMarkerStyle(markers[i]);
        gr_sv0pt->SetMarkerSize(0.6);

        // Adicionando os gráficos aos MultiGraphs
        mg_v0pt->Add(gr_v0pt);
        mg_sv0pt->Add(gr_sv0pt);

        // Adicionando à legenda global
        legend_systems->AddEntry(gr_v0pt, system_name[i], "lp");
    }

    // Criando o Canvas FORA do loop para colocar todos no mesmo lugar
    auto c = new TCanvas("c_all", "All Collision Systems", 1100, 500);
    c->Divide(2, 1);

    // Legenda de parâmetros fixos (usada no primeiro pad)
    auto legend_syst_params = new TLegend(0.025, 0.88, 0.5, 0.76);
    legend_syst_params->SetTextSize(0.045);
    legend_syst_params->AddEntry((TObject*)0, "5.36 GeV   0-5% centrality", "");
    legend_syst_params->AddEntry((TObject*)0, "p_{T}^{ref} 0.5-2 GeV    #eta_{gap}=1", "");
    legend_syst_params->SetBorderSize(0);
    legend_syst_params->SetFillStyle(0);

    // --- PAINEL 1: v0(pT) ---
    c->cd(1);
    gPad->SetLogx();
    mg_v0pt->Draw("AP"); // "A" desenha os eixos, "P" desenha os pontos
    mg_v0pt->GetXaxis()->SetLimits(0.485, 10.0);
    //mg_v0pt->SetMaximum(0.23);
    mg_v0pt->GetXaxis()->SetLabelFont(42);
    //mg_v0pt->GetXaxis()->SetLabelSize(0.05);
    //mg_v0pt->GetXaxis()->SetTitleSize(0.06);
    mg_v0pt->GetXaxis()->SetTitleFont(42);
    mg_v0pt->GetXaxis()->CenterTitle(true);
    mg_v0pt->GetYaxis()->SetLabelFont(42);
    //mg_v0pt->GetYaxis()->SetLabelSize(0.05);
    //mg_v0pt->GetYaxis()->SetTitleSize(0.06);
    mg_v0pt->GetYaxis()->SetTitleOffset(1.10);
    mg_v0pt->GetYaxis()->SetTitleFont(42);
    mg_v0pt->GetYaxis()->CenterTitle(true);
    legend_syst_params->Draw();
    legend_systems->Draw();

    // --- PAINEL 2: v0(pT)/v0 ---
    c->cd(2);
    gPad->SetLogx();
    mg_sv0pt->Draw("AP");
    mg_sv0pt->GetXaxis()->SetLimits(0.485, 10.0);
    //mg_sv0pt->SetMaximum(14.0);
    mg_sv0pt->GetXaxis()->SetLabelFont(42);
    //mg_sv0pt->GetXaxis()->SetLabelSize(0.05);
    //mg_sv0pt->GetXaxis()->SetTitleSize(0.06);
    mg_sv0pt->GetXaxis()->SetTitleFont(42);
    mg_sv0pt->GetXaxis()->CenterTitle(true);
    mg_sv0pt->GetYaxis()->SetLabelFont(42);
    //mg_sv0pt->GetYaxis()->SetLabelSize(0.05);
    //mg_sv0pt->GetYaxis()->SetTitleSize(0.06);
    mg_sv0pt->GetYaxis()->SetTitleOffset(1.10);
    mg_sv0pt->GetYaxis()->SetTitleFont(42);
    mg_sv0pt->GetYaxis()->CenterTitle(true);
    legend_systems->Draw(); // Coloquei a legenda dos sistemas no 2º painel para não sobrepor o texto
    legend_syst_params->Draw();

    // Salvando o resultado final
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    //c->Range(0.04000706,-0.4080692,0.7405774,1.786744);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(0);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetLeftMargin(0.1650854);
    c->SetRightMargin(0.0341556);
    c->SetTopMargin(0.08508404);
    c->SetBottomMargin(0.1859244);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->Update();
    c->SaveAs("./Plots/All_Systems_Combined.pdf");
    
    // Opcional: Limpar memória do canvas se for rodar em uma macro maior
    // delete c;
}
