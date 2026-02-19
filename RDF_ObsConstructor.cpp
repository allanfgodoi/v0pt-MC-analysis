#include <cmath>
#include "TMath.h"
using namespace std;
using namespace ROOT::VecOps;

vector<RVec<float>> Transpose(vector<RVec<float>>& x){
    int rows = x.size();
    int cols = x[0].size();
    vector<RVec<float>> v(cols, RVec<float>(rows, 0.0));
    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            v[j][i] = x[i][j];
        }
    }
    return v;
}

float Bootstrap(const RVec<float>& x, int B){
    static TRandom3 rndgen;
    int N = x.size();
    float lambda = 1.0;
    vector<float> bootstrapped(B, 0.0); 
    for (int i=0; i<B; i++){
        float acc = 0.0;
        float n = 0.0;
        for (int j=0; j<N; j++){
            float rd = rndgen.Poisson(lambda);
            acc += x[j]*rd;
            n += rd;
        }
        bootstrapped[i] = acc/n;
    }
    float bootstrapped_error = TMath::StdDev(B, bootstrapped.data());
    return bootstrapped_error;
}

struct Gathered_Data{
    float mean_pt;
    float mean_pt_ref;
    vector<float> vec_dpt_ref_A;
    vector<float> vec_dpt_ref_B;
    vector<float> vec_dpt_ref_AB;
    vector<RVec<float>> vec_dn_pt_A;
    vector<RVec<float>> vec_dn_pt_B;
    vector<RVec<float>> vec_dn_pt_AB;
    vector<RVec<float>> vec_sum_dpt_dn_pt;
    vector<RVec<float>> vec_n_pt_AB;
    vector<float> pT_axis;
};

Gathered_Data DataGathering(TString Filename, float eta_gap, int nch_min, int nch_max, float ptr_min, float ptr_max, unordered_set<int> targetPID, vector<float> Xaxis_del){
    TH1::AddDirectory(kFALSE);
    ROOT::EnableImplicitMT();

    // File check
    TFile *fileCheck = TFile::Open(Filename);
    if (!fileCheck || fileCheck->IsZombie()) std::cerr << "Cannot open ROOT file: " << Filename << std::endl;
    fileCheck->Close(); delete fileCheck;

    ROOT::RDataFrame df("tree", Filename); // Create RDataFrame
    cout << "Reading data..." << endl;

    // Defining kinematic variables
    auto df_kin = df.Define("pt_raw", "sqrt(px*px + py*py)")
                    .Define("theta_raw", "Where(pt_raw>1e-5, atan2(pt_raw, pz), 0.0)") // Avoids division by zero
                    .Define("eta_raw", "-log(tan(theta_raw/2))");

    // Centrality cut
    auto df_nch = df_kin.Define("Nch", "(int)Sum(pt_raw > 0.02 && abs(eta_raw) > 3.0 && abs(eta_raw) < 5.0)"); // Counting Nch in forward region (Sum: counts the number of true entries in each row)
    TString cent_logic = TString::Format("Nch >= %d && Nch <= %d", nch_min, nch_max);
    auto df_centrality = df_nch.Filter(cent_logic.Data(), "Centrality Cut");

    // Applying PID and kinematic cuts
    TString cut_expr = "";
    if (targetPID.count(0) > 0){
        cut_expr = "1";
    } else{
        bool first = true;
        for (int pid : targetPID){
            if (!first){
                cut_expr += " || ";
            }
            cut_expr += TString::Format("abs(pid)==%d", pid);
            first = false;
        }
    }
    TString cut_logic = TString::Format("(%s) && pt_raw>0.5 && pt_raw<10.0 && abs(eta_raw)<2.4", cut_expr.Data());
    auto df_cuted = df_kin  .Define("cut", cut_logic.Data())
                            .Define("pt", "pt_raw[cut]")
                            .Define("eta", "eta_raw[cut]");

    // Separing in subevents A and B and pT-ref
    TString subA_logic = TString::Format("eta >= -2.4 && eta <= -%f/2", eta_gap);
    TString subB_logic = TString::Format("eta >= %f/2 && eta <= 2.4", eta_gap);
    TString ptref_A_logic = TString::Format("pt_A >= %f && pt_A <= %f", ptr_min, ptr_max);
    TString ptref_B_logic = TString::Format("pt_B >= %f && pt_B <= %f", ptr_min, ptr_max);
    TString ptref_AB_logic = TString::Format("pt_AB >= %f && pt_AB <= %f", ptr_min, ptr_max);
    auto df_subs = df_cuted .Define("cut_subA", subA_logic.Data())
                            .Define("cut_subB", subB_logic.Data())
                            .Define("pt_A", "pt[cut_subA]")
                            .Define("pt_B", "pt[cut_subB]")
                            .Define("pt_AB", "pt[cut_subA || cut_subB]")
                            .Define("cut_ptref_A", ptref_A_logic.Data())
                            .Define("cut_ptref_B", ptref_B_logic.Data())
                            .Define("cut_ptref_AB", ptref_AB_logic.Data())
                            .Define("pt_ref_A", "pt_A[cut_ptref_A]")
                            .Define("pt_ref_B", "pt_B[cut_ptref_B]")
                            .Define("pt_ref_AB", "pt_AB[cut_ptref_AB]");

    // Binning auxiliar lambda function
    auto binner = [Xaxis_del](const RVec<float>& x){
        int nBins = Xaxis_del.size()-1;
        RVec<float> counter(nBins, 0.0);

        for (auto val : x){
            auto it = upper_bound(Xaxis_del.begin(), Xaxis_del.end(), val); // Returns the upper bound of val's bin
            int idx = distance(Xaxis_del.begin(), it)-1; // Set the index of val's bin 
            counter[idx] += 1.0;
        }
        float N = x.size();
        if (N > 0) counter = counter/N;
        return counter;
    };

    // Constructing n_A(pT), n_B(pT) and n(pT) (using {pt} to take track's pT from each event)
    auto df_bins = df_subs  .Define("n_pt_A", binner, {"pt_A"})
                            .Define("n_pt_B", binner, {"pt_B"})
                            .Define("n_pt_AB", binner, {"pt_AB"});

    // Calculating mean pTs
    auto ptr_mean_pt_AB = df_subs.Mean("pt_AB"); float mean_pt_AB = *ptr_mean_pt_AB; // <pT>
    auto ptr_mean_pt_ref_A = df_subs.Mean("pt_ref_A"); float mean_pt_ref_A = *ptr_mean_pt_ref_A; // <pT_ref_A>
    auto ptr_mean_pt_ref_B = df_subs.Mean("pt_ref_B"); float mean_pt_ref_B = *ptr_mean_pt_ref_B; // <pT_ref_B>
    auto ptr_mean_pt_ref_AB = df_subs.Mean("pt_ref_AB"); float mean_pt_ref_AB = *ptr_mean_pt_ref_AB; // <pT_ref_AB>

    // Evaluating delta [pT]s (d[pT_ref_i] = [pT_ref_i] - <pT_ref_i>)
    auto df_dpts = df_bins  .Define("mean_evt_pt_ref_A", "Mean(pt_ref_A)")
                            .Define("mean_evt_pt_ref_B", "Mean(pt_ref_B)")
                            .Define("dpt_ref_A", [mean_pt_ref_A](double mean_evt_pt_ref){return (float)(mean_evt_pt_ref-mean_pt_ref_A);}, {"mean_evt_pt_ref_A"})
                            .Define("dpt_ref_B", [mean_pt_ref_B](double mean_evt_pt_ref){return (float)(mean_evt_pt_ref-mean_pt_ref_B);}, {"mean_evt_pt_ref_B"})
                            .Define("dpt_ref_AB", "dpt_ref_A*dpt_ref_B");                    

    // Getting the track fractions of each bin of ALL events (tranposing n_pt_i)
    auto ptr_n_pt_A = df_bins.Take<RVec<float>>("n_pt_A");
    auto ptr_n_pt_B = df_bins.Take<RVec<float>>("n_pt_B");
    auto ptr_n_pt_AB = df_bins.Take<RVec<float>>("n_pt_AB");
    vector<RVec<float>> T_n_pt_A = Transpose(*ptr_n_pt_A);
    vector<RVec<float>> T_n_pt_B = Transpose(*ptr_n_pt_B);
    vector<RVec<float>> T_n_pt_AB = Transpose(*ptr_n_pt_AB);

    // Calculating delta n(pT)s (dn(pT) = n)
    int nBins = T_n_pt_A.size();
    RVec<float> vec_mean_n_pt_A(nBins, 0.0), vec_mean_n_pt_B(nBins, 0.0);
    for (int i=0; i<nBins; i++){
        vec_mean_n_pt_A[i] = Mean(T_n_pt_A[i]);
        vec_mean_n_pt_B[i] = Mean(T_n_pt_B[i]);
    }

    auto df_final = df_dpts .Define("dn_pt_A", [vec_mean_n_pt_A](const RVec<float>& n_pt){return n_pt-vec_mean_n_pt_A;}, {"n_pt_A"})
                            .Define("dn_pt_B", [vec_mean_n_pt_B](const RVec<float>& n_pt){return n_pt-vec_mean_n_pt_B;}, {"n_pt_B"})
                            .Define("dn_pt_AB", "dn_pt_A*dn_pt_B")
                            .Define("vec_sum_dpt_dn_pt", "dn_pt_A*dpt_ref_B + dn_pt_B*dpt_ref_A"); // This is the correlation numerator

    // Setting up final struct
    Gathered_Data data;
    data.mean_pt = mean_pt_AB;
    data.mean_pt_ref = mean_pt_ref_AB;
    data.vec_dpt_ref_A = *df_final.Take<float>("dpt_ref_A");
    data.vec_dpt_ref_B = *df_final.Take<float>("dpt_ref_B");
    data.vec_dpt_ref_AB = *df_final.Take<float>("dpt_ref_AB");
    data.vec_dn_pt_A = *df_final.Take<RVec<float>>("dn_pt_A");
    data.vec_dn_pt_B = *df_final.Take<RVec<float>>("dn_pt_B");
    data.vec_sum_dpt_dn_pt = *df_final.Take<RVec<float>>("vec_sum_dpt_dn_pt");
    data.vec_n_pt_AB = *ptr_n_pt_AB;
    // Setting pT axis center
    vector<double> bins(Xaxis_del.begin(), Xaxis_del.end());
    auto ptr_profile = df_subs.Profile1D({"PpT_RDF", "Mean pT Axis", (int)bins.size()-1, bins.data()}, "pt_AB", "pt_AB");
    for(int i=1; i<=ptr_profile->GetNbinsX(); i++) data.pT_axis.push_back(ptr_profile->GetBinContent(i));
    return data;
}

void ObsConstructor(float Eta_gap, int Nch_min, int Nch_max, float pTr_Min, float pTr_Max, unordered_set<int> TargetPID, TString Name, TString Filename, TString Savename){
    
    int B = 100; // Number of Poisson bootstrap samples

    // Defining bins and plot's x axis
    vector<float> Xaxis_del = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 1.98, 2.2, 2.38, 2.98, 3.8, 4.5, 6.0, 8.0, 10.0}; // Those are the END of each bin, not the middle
    int nBins = (Xaxis_del.size()-1);

    Gathered_Data gData = DataGathering(Filename, Eta_gap, Nch_min, Nch_max, pTr_Min, pTr_Max, TargetPID, Xaxis_del);
    float mean_pt = gData.mean_pt;
    float mean_pt_ref = gData.mean_pt_ref;
    vector<float> vec_dpt_ref_A = gData.vec_dpt_ref_A;
    vector<float> vec_dpt_ref_B = gData.vec_dpt_ref_B;
    vector<float> vec_dpt_ref_AB = gData.vec_dpt_ref_AB;
    vector<RVec<float>> vec_dn_pt_A = gData.vec_dn_pt_A;
    vector<RVec<float>> vec_dn_pt_B = gData.vec_dn_pt_B;
    vector<RVec<float>> vec_sum_dpt_dn_pt = gData.vec_sum_dpt_dn_pt;
    vector<RVec<float>> vec_n_pt_AB = gData.vec_n_pt_AB;
    vector<float> pT_axis = gData.pT_axis;

    cout << "Constructing observables..." << endl;

    int nEvents = vec_dpt_ref_AB.size();

    // Calculating v0
    float sigma = sqrt(TMath::Mean(nEvents, vec_dpt_ref_AB.data()));
    float v0 = sigma/mean_pt_ref;
    // v0 uncertainty
    RVec<float> rvec_dpt_ref_AB(vec_dpt_ref_AB.data(), vec_dpt_ref_AB.size()); // Must convert this standard vector to RVec
    float unc_dpt_ref_AB = Bootstrap(rvec_dpt_ref_AB, B);
    float unc_v0 = unc_dpt_ref_AB/(2*mean_pt_ref*sqrt(unc_dpt_ref_AB)); // Uncertainty propagation formula

    cout << "v0 = " << v0 << " +- " << unc_v0 << endl;

    // Defining auxiliar vectors to calculate v0(pT)v0
    vector<RVec<float>> T_vec_sum_dpt_dn_pt = Transpose(vec_sum_dpt_dn_pt);
    vector<RVec<float>> T_vec_n_pt_AB = Transpose(vec_n_pt_AB);
    vector<float> vec_num_v0ptv0(nBins, 0.0);
    vector<float> vec_denom_v0ptv0(nBins, 0.0);
    vector<float> vec_v0ptv0(nBins, 0.0);
    vector<float> vec_unc_num_v0ptv0(nBins, 0.0);
    vector<float> vec_unc_denom_v0ptv0(nBins, 0.0);
    vector<float> vec_unc_v0ptv0(nBins, 0.0);
    // Calculating v0(pT)v0    
    for (int i=0; i<nBins; i++){
        vec_num_v0ptv0[i] = Mean(T_vec_sum_dpt_dn_pt[i]);
        vec_denom_v0ptv0[i] = Mean(T_vec_n_pt_AB[i])*mean_pt_ref;
        vec_v0ptv0[i] = 0.5*vec_num_v0ptv0[i]/vec_denom_v0ptv0[i];
        // v0(pT)v0 uncertainty
        vec_unc_num_v0ptv0[i] = Bootstrap(T_vec_sum_dpt_dn_pt[i], B);
        vec_unc_denom_v0ptv0[i] = Bootstrap(T_vec_n_pt_AB[i], B);
        vec_unc_v0ptv0[i] = (1/mean_pt_ref)*abs((vec_unc_num_v0ptv0[i]/vec_denom_v0ptv0[i])-(vec_unc_denom_v0ptv0[i]*vec_num_v0ptv0[i])); // Uncertainty propagation (tot. corr.)
    }

    // Defining auxiliar vectors to calculate v0(pT)
    vector<float> vec_v0pt(nBins, 0.0);
    vector<float> vec_unc_v0pt(nBins, 0.0);
    // Calculating v0(pT)
    for (int i=0; i<nBins; i++){
        vec_v0pt[i] = vec_v0ptv0[i]/v0;
        // v0(pT) uncertainty
        vec_unc_v0pt[i] = abs((vec_unc_v0ptv0[i]/v0)-(vec_v0ptv0[i]*unc_v0));
    }

    // Defining auxiliar vectors to calculate v0(pT)/v0
    vector<float> vec_sv0pt(nBins, 0.0);
    vector<float> vec_unc_sv0pt(nBins, 0.0);
    // Calculating v0(pT)/v0
    for (int i=0; i<nBins; i++){
        vec_sv0pt[i] = vec_v0pt[i]/v0;
        // v0(pT)/v0 uncertainty
        vec_unc_sv0pt[i] = abs((vec_unc_v0pt[i]/v0)-(vec_v0pt[i]*unc_v0));
    }

    TFile *save_file = new TFile(Savename, "UPDATE"); // Savefile

    // Defining PID labels
    if (TargetPID.size() == 1){
        if (TargetPID.count(211)) Name += "_pion";
        if (TargetPID.count(321)) Name += "_kaon";
        if (TargetPID.count(2212)) Name += "_proton";
    } else {
        Name += "_all";
    }
    TString mean_pt_name = "mean_pt_";
    mean_pt_name += Name;

    // Saving <pT>
    TParameter<float> *p_mean_pt = new TParameter<float>(mean_pt_name, mean_pt);
    p_mean_pt->Write(mean_pt_name, TObject::kOverwrite); // Does not depend on pT-ref (only eta-gap)

    vector<float> vec_zeros(nBins, 0.0);

    // Creating and writting v0(pT)v0 TGraph
    vector<float> vec_v0ptv0_plot(nBins, 0.0);
    vector<float> vec_unc_v0ptv0_plot(nBins, 0.0);
    for (int i=0; i<nBins; i++){
        vec_v0ptv0_plot[i] = 1e3*vec_v0ptv0[i];
        vec_unc_v0ptv0_plot[i] = 1e3*vec_unc_v0ptv0[i];
    }
    TGraph* gr_v0ptv0_ptref = new TGraphErrors(nBins, pT_axis.data(), vec_v0ptv0_plot.data(), vec_zeros.data(), vec_unc_v0ptv0_plot.data());
    TString v0ptv0_name = "v0ptv0_";
    v0ptv0_name += Name;
    gr_v0ptv0_ptref->SetName(v0ptv0_name);
    gr_v0ptv0_ptref->Write();

    // Creating and writting v0(pT) TGraph
    TGraphErrors* gr_v0pt = new TGraphErrors(nBins, pT_axis.data(), vec_v0pt.data(), vec_zeros.data(), vec_unc_v0pt.data());
    TString v0pt_name = "v0pt_";
    v0pt_name += Name;
    gr_v0pt->SetName(v0pt_name);
    gr_v0pt->Write();

    // Creating and writting v0(pT)/v0 TGraph
    TGraphErrors* gr_sv0pt = new TGraphErrors(nBins, pT_axis.data(), vec_sv0pt.data(), vec_zeros.data(), vec_unc_sv0pt.data());
    TString sv0pt_name = "sv0pt_";
    sv0pt_name += Name;
    gr_sv0pt->SetName(sv0pt_name);
    gr_sv0pt->Write();

    cout << "mean_pt = " << mean_pt << endl;
    save_file->Close();
}