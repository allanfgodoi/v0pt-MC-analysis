// USAGE EXAMPLE (0-5%): ObsConstructor(1.0, 105, 9999, 0.5, 2.0, {211, 321, 2212}, "0005", "/eos/cms/store/user/sdogra/ampt/inPutFiles/HeHe/*.dat", "/eos/user/a/afloresg/MC-studies/HeHe/ObsData.root")
// ObsConstructor(eta-gap, Nch-min, Nch-max, pT-ref-min, pT-ref-max, {TargetPIDs}, "centrality-label", "input-files-pattern", "output-file")

#include <cmath>
#include "TMath.h"
#include <unordered_set>
#include <fstream>
using namespace std;

vector<vector<float>> transpose(vector<vector<float>> x){
    int rows = x.size();
    int cols = x[0].size();
    vector<vector<float>> transposed(cols, vector<float>(rows, 0.0));
    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            transposed[j][i] = x[i][j];
        }
    }
    return transposed;
}

vector<float> double_vector_mean(vector<vector<float>> x, const unsigned int nBins){
    int nEvents_A = x.size();
    vector<float> vec;
    for (int i=0; i<nBins; i++){ // pT bins
        float acc_f_pt = 0;
        for (int j=0; j<nEvents_A; j++){ // Events
            acc_f_pt += x[j][i];
        }
        vec.push_back(acc_f_pt/nEvents_A);
    }
    return vec;
}

float StdPoissonBootstrap(vector<float> x, int B){
    TRandom3 rndgen;
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
    vector<float> pT_axis;
    vector<float> vec_dPt_A;
    vector<float> vec_dPt_B;
    vector<float> vec_dPt_ref_A;
    vector<float> vec_dPt_ref_B;
    vector<vector<float>> vec_n_pt_AB;
    vector<vector<float>> vec_dn_pt_A;
    vector<vector<float>> vec_dn_pt_B;
    float mean_pt_ref;
    float mean_pt;
    double total_tracks;
    vector<double> tracks_per_bin;
};

Gathered_Data DataGathering(TString Filename, float eta_gap, int nch_min, int nch_max, float ptr_min, float ptr_max, unordered_set<int> targetPID, vector<float> Xaxis_del){
    TH1::AddDirectory(kFALSE);

    // Find input files
    TString fp_pattern(Filename);
    TString dir = gSystem->DirName(fp_pattern);
    TString basename = gSystem->BaseName(fp_pattern);

    vector<TString> files;
    TSystemDirectory directory("amptDir", dir);
    TList* filelist = directory.GetListOfFiles();
    TIter next(filelist);
    TSystemFile* file;
    TRegexp re(basename, true);

    while ((file = (TSystemFile*)next())){
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.Index(re) != kNPOS){
            files.push_back(dir + "/" + fname);
        }
    }
    cout << "Found " << files.size() << " files" << endl;

    // Defining auxiliar constants
    int nBins = (Xaxis_del.size()-1);
    vector<vector<float>> Matrix_trkPt_A;
    vector<vector<float>> Matrix_trkPt_B;
    vector<vector<float>> Matrix_trkW_A;
    vector<vector<float>> Matrix_trkW_B;
    vector<vector<float>> Vec_n_pt_AB; // This vector will hold the fractions of pT of all events
    vector<vector<float>> Vec_n_pt_A; 
    vector<vector<float>> Vec_n_pt_B;
    
    TH1F *hist_all_pt_ref_AB = new TH1F("all_pt_ref", "All pT-ref", 100, 0, 100); // Create histogram to calculate <pT-ref_A> (all events)
    TH1F *hist_all_pt_ref_A = new TH1F("all_pt_ref_A", "All pT-ref from subset A", 100, 0, 100); 
    TH1F *hist_all_pt_ref_B = new TH1F("all_pt_ref_B", "All pT-ref from subset B", 100, 0, 100);
    TH1F *hist_all_pt_AB = new TH1F("all_pt", "All pT", 100, 0, 100);
    TH1F *hist_all_pt_A = new TH1F("all_pt_A", "All pT from subset A", 100, 0, 100); 
    TH1F *hist_all_pt_B = new TH1F("all_pt_B", "All pT from subset B", 100, 0, 100);
    TH1F *hist_pt_AB = new TH1F("pt", "pT", nBins, Xaxis_del.data()); // Create histogram to scale pT bins and get n(pT)
    TH1F *hist_pt_A = new TH1F("pt_A", "pT from subset A", nBins, Xaxis_del.data());
    TH1F *hist_pt_B = new TH1F("pt_B", "pT from subset B", nBins, Xaxis_del.data());
    
    TProfile *PpT = new TProfile("PpT", "pT bins mean", nBins, Xaxis_del.data());

    int fCounter = 1;
    int N_all_valid_tracks = 0;

    for (auto& fname : files){
        ifstream amptFile(fname.Data());

        // AMPT header aux variables
        int eventID, dummy, n_particles, npartP, npartT, nElaP, nInElaP, nElaT, nInElaT;
        float b, unused, x, y, z, t;

        while (amptFile >> eventID >> dummy >> n_particles >> b >> npartP >> npartT >> nElaP >> nInElaP >> nElaT >> nInElaT >> unused){
            vector<int> evt_pid;
            vector<float> evt_px, evt_py, evt_pz;

            if (eventID % 1000 == 0) cout << "File " << fCounter << " - " << "Processing event: " << eventID << endl;

            for (int i = 0; i < n_particles; i++) {
                // Holder variables for each track
                int h_pid;
                float h_px, h_py, h_pz, h_E;

                amptFile >> h_pid >> h_px >> h_py >> h_pz >> h_E >> x >> y >> z >> t;
                evt_pid.push_back(h_pid);
                evt_px.push_back(h_px);
                evt_py.push_back(h_py);
                evt_pz.push_back(h_pz);
            }

            int Nch = 0;

            // 1st track loop to find event multiplicty Nch (No filters except for centrality select criteria)
            for (int i=0; i<n_particles; i++){
                int pid = evt_pid[i];
                float px = evt_px[i];
                float py = evt_py[i];
                float pz = evt_pz[i];

                double pt = sqrt(px*px + py*py);
                if (pt <= 0.02) continue;
                double theta = 0.0;
                if (pt > 1e-5) theta = atan2(pt, pz); // Avoids division by zero
                double eta = -log(tan(theta/2.0));

                // Nch definition
                if (fabs(eta) > 3.0 && fabs(eta) < 5.0) Nch++;
            }

            if (Nch < nch_min || Nch > nch_max) continue; // Centrality selection

            vector<float> n_pt_A(nBins, 0.0); // Define vector to hold the fractions of pT in the event
            vector<float> n_pt_B(nBins, 0.0);
            vector<float> n_pt_AB(nBins, 0.0);
            vector<float> Vec_trkPt_A, Vec_trkPt_B;
            vector<float> Vec_trkW_A, Vec_trkW_B;

            float corrFac = 1.0; // Monte Carlo events. I'm using this to reuse the same code as data

            // 2nd track loop, now to calculate the quantities
            for (int i=0; i<n_particles; i++){
                int pid = evt_pid[i];
                float px = evt_px[i];
                float py = evt_py[i];
                float pz = evt_pz[i];

                // PID filter
                if (targetPID.count(abs(pid)) == 0) continue;

                double pt = sqrt(px*px + py*py);
                if (pt < 0.5 || pt > 10.0) continue; // Kinematic filters pt. 1 (same as ATLAS)
                double theta = 0.0;
                if (pt > 1e-5) theta = atan2(pt, pz); // Avoids division by zero
                double eta = -log(tan(theta/2.0));

                if (fabs(eta) > 2.4) continue; // Kinematic filter pt. 2 (same as ATLAS)

                // Subevent selection
                bool isInA = (eta >= -2.4 && eta <= -eta_gap/2.0);
                bool isInB = (eta >= eta_gap/2.0 && eta <= 2.4);
                if (isInA){
                    if (pt >= ptr_min && pt <= ptr_max){
                        Vec_trkPt_A.push_back(pt);
                        Vec_trkW_A.push_back(corrFac);
                        hist_all_pt_ref_AB->Fill(pt, corrFac);
                        hist_all_pt_ref_A->Fill(pt, corrFac);
                    }
                    hist_pt_A->Fill(pt, corrFac);
                    hist_pt_AB->Fill(pt, corrFac);
                    hist_all_pt_AB->Fill(pt, corrFac);
                    hist_all_pt_A->Fill(pt, corrFac);
                    PpT->Fill(pt, pt, corrFac);
                }
                if (isInB){
                    if (pt >= ptr_min && pt <= ptr_max){
                        Vec_trkPt_B.push_back(pt);
                        Vec_trkW_B.push_back(corrFac);
                        hist_all_pt_ref_AB->Fill(pt, corrFac);
                        hist_all_pt_ref_B->Fill(pt, corrFac);
                    }
                    hist_pt_B->Fill(pt, corrFac);
                    hist_pt_AB->Fill(pt, corrFac);
                    hist_all_pt_AB->Fill(pt, corrFac);
                    hist_all_pt_B->Fill(pt, corrFac);
                    PpT->Fill(pt, pt, corrFac);
                }
            }

            // Scaling the created hist and taking the bins content to the array
            if (hist_pt_A->Integral() != 0) hist_pt_A->Scale(1/hist_pt_A->Integral()); // Defining a normalized histogram
            if (hist_pt_B->Integral() != 0) hist_pt_B->Scale(1/hist_pt_B->Integral());
            if (hist_pt_AB->Integral() != 0) hist_pt_AB->Scale(1/hist_pt_AB->Integral());
            for (int i=0; i<nBins; i++){
                n_pt_A[i] = hist_pt_A->GetBinContent(i+1);
                n_pt_B[i] = hist_pt_B->GetBinContent(i+1);
                n_pt_AB[i] = hist_pt_AB->GetBinContent(i+1);
            }
            hist_pt_A->Reset();
            hist_pt_B->Reset();
            hist_pt_AB->Reset();
            // Setting the desired vectors
            Vec_n_pt_A.push_back(n_pt_A);
            Vec_n_pt_B.push_back(n_pt_B);
            Vec_n_pt_AB.push_back(n_pt_AB);

            Matrix_trkPt_A.push_back(Vec_trkPt_A);
            Matrix_trkPt_B.push_back(Vec_trkPt_B);
            Matrix_trkW_A.push_back(Vec_trkW_A);
            Matrix_trkW_B.push_back(Vec_trkW_B);

            Vec_trkPt_A.clear();
            Vec_trkPt_B.clear();
            Vec_trkW_A.clear();
            Vec_trkW_B.clear();
            Vec_trkPt_A.shrink_to_fit();
            Vec_trkPt_B.shrink_to_fit();
            Vec_trkW_A.shrink_to_fit();
            Vec_trkW_B.shrink_to_fit();
        }

    fCounter++;
    }

    float Mean_pt = hist_all_pt_AB->GetMean();
    float Mean_pt_A = hist_all_pt_A->GetMean();
    float Mean_pt_B = hist_all_pt_B->GetMean();
    float Mean_pt_ref = hist_all_pt_ref_AB->GetMean();
    float Mean_pt_ref_A = hist_all_pt_ref_A->GetMean();
    float Mean_pt_ref_B = hist_all_pt_ref_B->GetMean();

    vector<float> Pt_axis(nBins, 0.0);
    for (int i=0; i<nBins; i++){
        Pt_axis[i] = PpT->GetBinContent(i+1);
    }

    vector<vector<float>> Matrix_dtrkPt_ref_A; // Shape: nEvents x nTrk, but nTrk isn't fixed
    vector<vector<float>> Matrix_dtrkPt_ref_B;
    vector<vector<float>> Matrix_dtrkPt_A; 
    vector<vector<float>> Matrix_dtrkPt_B;
    vector<float> Vec_dtrkPt_ref_A;
    vector<float> Vec_dtrkPt_ref_B;
    vector<float> Vec_dtrkPt_A;
    vector<float> Vec_dtrkPt_B;
    for (int i=0; i<Matrix_trkPt_A.size(); i++){ // Events (i)
        for (int j=0; j<Matrix_trkPt_A[i].size(); j++){ // Tracks (j)
            Vec_dtrkPt_A.push_back(Matrix_trkPt_A[i][j]-Mean_pt_ref_A);
            Vec_dtrkPt_ref_A.push_back(Matrix_trkPt_A[i][j]-Mean_pt_ref_A);
        }
        Matrix_dtrkPt_ref_A.push_back(Vec_dtrkPt_ref_A);
        Matrix_dtrkPt_A.push_back(Vec_dtrkPt_A);
        Vec_dtrkPt_ref_A.clear();
        Vec_dtrkPt_A.clear();
    }
    for (int i=0; i<Matrix_trkPt_B.size(); i++){ // Events (i)
        for (int j=0; j<Matrix_trkPt_B[i].size(); j++){ // Tracks (j)
            Vec_dtrkPt_B.push_back(Matrix_trkPt_B[i][j]-Mean_pt_ref_B);
            Vec_dtrkPt_ref_B.push_back(Matrix_trkPt_B[i][j]-Mean_pt_ref_B);
        }
        Matrix_dtrkPt_ref_B.push_back(Vec_dtrkPt_ref_B);
        Matrix_dtrkPt_B.push_back(Vec_dtrkPt_B);
        Vec_dtrkPt_ref_B.clear();
        Vec_dtrkPt_B.clear();
    }

    int nValid_Events = Matrix_trkPt_A.size(); // Non-zero events. Subset A and B have the same number of events

    vector<float> Vec_dPt_ref_A(nValid_Events, 0.0);
    vector<float> Vec_dPt_ref_B(nValid_Events, 0.0);
    vector<float> Vec_dPt_A(nValid_Events, 0.0);
    vector<float> Vec_dPt_B(nValid_Events, 0.0);

    for (int i=0; i<nValid_Events; i++){
        float sum_WdPt_A = 0.0; float sum_W_A = 0.0;
        float sum_WdPt_B = 0.0; float sum_W_B = 0.0;
        float sum_WdPt_ref_A = 0.0; float sum_WdPt_ref_B = 0.0;
        for (int j=0; j<Matrix_dtrkPt_A[i].size(); j++){
            sum_WdPt_ref_A += Matrix_trkW_A[i][j]*Matrix_dtrkPt_ref_A[i][j];
            sum_WdPt_A += Matrix_trkW_A[i][j]*Matrix_dtrkPt_A[i][j];
            sum_W_A += Matrix_trkW_A[i][j];
        }
        for (int j=0; j<Matrix_dtrkPt_B[i].size(); j++){
            sum_WdPt_ref_B += Matrix_trkW_B[i][j]*Matrix_dtrkPt_ref_B[i][j];
            sum_WdPt_B += Matrix_trkW_B[i][j]*Matrix_dtrkPt_B[i][j];
            sum_W_B += Matrix_trkW_B[i][j];
        }
        if (sum_W_A == 0) sum_W_A = 1.0;
        if (sum_W_B == 0) sum_W_B = 1.0;
        Vec_dPt_ref_A[i] = sum_WdPt_ref_A/sum_W_A;
        Vec_dPt_ref_B[i] = sum_WdPt_ref_B/sum_W_B;
        Vec_dPt_A[i] = sum_WdPt_A/sum_W_A;
        Vec_dPt_B[i] = sum_WdPt_B/sum_W_B;
    }

    // Calculating delta n(pT)
    vector<float> Vec_mean_n_pt_A(nBins, 0.0);
    vector<float> Vec_mean_n_pt_B(nBins, 0.0);
    vector<vector<float>> Vec_dn_pt_A(nValid_Events, vector<float>(nBins, 0.0));
    vector<vector<float>> Vec_dn_pt_B(nValid_Events, vector<float>(nBins, 0.0));
    for (int i=0; i<nBins; i++){
        Vec_mean_n_pt_A[i] = TMath::Mean(nValid_Events, transpose(Vec_n_pt_A)[i].data());
        Vec_mean_n_pt_B[i] = TMath::Mean(nValid_Events, transpose(Vec_n_pt_B)[i].data());
    }
    for (int i=0; i<nValid_Events; i++){
        for (int j=0; j<nBins; j++){
            Vec_dn_pt_A[i][j] = Vec_n_pt_A[i][j] - Vec_mean_n_pt_A[j];
            Vec_dn_pt_B[i][j] = Vec_n_pt_B[i][j] - Vec_mean_n_pt_B[j];
        }
    }

    vector<double> Tracks_per_bin(nBins, 0.0);
    for (int i=0; i<nBins; i++){
        Tracks_per_bin[i] = PpT->GetBinEntries(i+1);
    }
    double Total_tracks = PpT->GetEntries();

    delete hist_all_pt_ref_A;
    delete hist_all_pt_ref_B;
    delete hist_all_pt_ref_AB;
    delete hist_all_pt_A;
    delete hist_all_pt_B;
    delete hist_pt_AB;
    delete hist_pt_A;
    delete hist_pt_B;

    Gathered_Data struct_data;
    struct_data.pT_axis = Pt_axis;
    struct_data.vec_dPt_A = Vec_dPt_A;
    struct_data.vec_dPt_B = Vec_dPt_B;
    struct_data.vec_dPt_ref_A = Vec_dPt_ref_A;
    struct_data.vec_dPt_ref_B = Vec_dPt_ref_B;
    struct_data.vec_n_pt_AB = Vec_n_pt_AB;
    struct_data.vec_dn_pt_A = Vec_dn_pt_A;
    struct_data.vec_dn_pt_B = Vec_dn_pt_B;
    struct_data.mean_pt_ref = Mean_pt_ref;
    struct_data.mean_pt = Mean_pt;
    struct_data.total_tracks = Total_tracks;
    struct_data.tracks_per_bin = Tracks_per_bin;
    return struct_data;
}

// Thats the function we call to construct the observable
void ObsConstructor(float Eta_gap, int Nch_min, int Nch_max, float pTr_Min, float pTr_Max, unordered_set<int> TargetPID, TString Name, TString Filename, TString Savename){
    int B = 100; // Number of Poisson bootstrap samples
    // Defining bins and plot's x axes
    vector<float> Xaxis_del = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 1.98, 2.2, 2.38, 2.98, 3.8, 4.5, 6.0, 8.0, 10.0}; // Those are the END of each bin, not the middle
    int nBins = (Xaxis_del.size()-1);
    Gathered_Data gData = DataGathering(Filename, Eta_gap, Nch_min, Nch_max, pTr_Min, pTr_Max, TargetPID, Xaxis_del);
    vector<float> pT_axis = gData.pT_axis;
    vector<float> vec_dPt_A = gData.vec_dPt_A;
    vector<float> vec_dPt_B = gData.vec_dPt_B;
    vector<float> vec_dPt_ref_A = gData.vec_dPt_ref_A;
    vector<float> vec_dPt_ref_B = gData.vec_dPt_ref_B;
    vector<vector<float>> vec_n_pt_AB = gData.vec_n_pt_AB;
    vector<vector<float>> vec_dn_pt_A = gData.vec_dn_pt_A;
    vector<vector<float>> vec_dn_pt_B = gData.vec_dn_pt_B;
    float mean_pt_ref = gData.mean_pt_ref;
    float mean_pt = gData.mean_pt;
    double total_tracks = gData.total_tracks;
    vector<double> tracks_per_bin = gData.tracks_per_bin;

    int nEvents = vec_dPt_A.size();

    cout << "Bootstrapping samples..." << endl;

    // v0
    vector<float> vec_dPt_ref_AB(nEvents, 0.0);
    for (int i=0; i<nEvents; i++){
        float h_dPt_ref_AB = vec_dPt_ref_A[i]*vec_dPt_ref_B[i];
        if (TMath::IsNaN(h_dPt_ref_AB) || !TMath::Finite(h_dPt_ref_AB)){
            vec_dPt_ref_AB[i] = 0.0;
        } else{
            vec_dPt_ref_AB[i] = h_dPt_ref_AB;
        }
    }

    float sigma = sqrt(TMath::Mean(nEvents, vec_dPt_ref_AB.data()));
    float v0 = sigma/mean_pt_ref;
    // v0 uncertainty
    float unc_dPt_ref_AB = StdPoissonBootstrap(vec_dPt_ref_AB, B);
    float unc_v0 = unc_dPt_ref_AB/(2*mean_pt_ref*sqrt(unc_dPt_ref_AB));
    cout << "v0 = " << v0 << " +- " << unc_v0 << endl;    

    // Calculating n_A(pT)*d[pT]_B + n_B(pT)*d[pT]_A
    vector<vector<float>> vec_sum_dpt_dn_pt(nEvents, vector<float>(nBins, 0.0));
    for (int i=0; i<nEvents; i++){ // Events
        for (int j=0; j<nBins; j++){ // pT bins
            vec_sum_dpt_dn_pt[i][j] = (vec_dn_pt_A[i][j]*vec_dPt_B[i])+(vec_dn_pt_B[i][j]*vec_dPt_A[i]);
        }
    }

    // Defining auxiliar vectors to calculate v0(pT)v0 and its uncertainty
    vector<float> vec_v0ptv0(nBins, 0.0);
    vector<float> vec_v0ptv0_denom(nBins, 0.0);
    vector<float> vec_mean_n_pt_AB = double_vector_mean(vec_n_pt_AB, nBins);
    vector<float> vec_unc_v0ptv0_num(nBins, 0.0);
    vector<float> vec_unc_n_pt_AB(nBins, 0.0);
    vector<float> vec_unc_v0ptv0(nBins, 0.0);

    // Calculating v0(pT)v0 and its uncertainty
    vector<float> vec_v0ptv0_num = double_vector_mean(vec_sum_dpt_dn_pt, nBins); // v0(pT) numerator
    for (int i=0; i<nBins; i++){
        // Calculating the observable v0(pT)v0
        vec_v0ptv0_denom[i] = vec_mean_n_pt_AB[i]*mean_pt_ref;
        vec_v0ptv0[i] = 0.5*vec_v0ptv0_num[i]/vec_v0ptv0_denom[i];
        // Calculating the uncertainty of v0(pT)
        vec_unc_v0ptv0_num[i] = StdPoissonBootstrap(transpose(vec_sum_dpt_dn_pt)[i], B); // Bootstraping numerator term
        vec_unc_n_pt_AB[i] = StdPoissonBootstrap(transpose(vec_n_pt_AB)[i], B); // Bootstrapping <n(pT)>
        vec_unc_v0ptv0[i] = (1/mean_pt_ref)*abs((vec_unc_v0ptv0_num[i]/vec_mean_n_pt_AB[i])-(vec_unc_n_pt_AB[i]*vec_v0ptv0_num[i])); // Uncertainty propagation (tot. corr.)
    }

    vector<float> vec_v0pt(nBins, 0.0);
    vector<float> vec_unc_v0pt(nBins, 0.0);
    // Calculating v0(pT)
    for (int i=0; i<nBins; i++){
        vec_v0pt[i] = vec_v0ptv0[i]/v0;
        vec_unc_v0pt[i] = abs((vec_unc_v0ptv0[i]/v0)-(vec_v0ptv0[i]*unc_v0));
    }

     // Calculating v0(pT) sum rules
    float sum1_v0pt = 0;
    float sum2_v0pt_left = 0;
    float acc_sum2_v0pt_right = 0;
    for (int i=0; i<nBins; i++){
        sum1_v0pt += vec_v0pt[i]*vec_mean_n_pt_AB[i];
        float pT = pT_axis[i];
        sum2_v0pt_left += pT*vec_v0pt[i]*vec_mean_n_pt_AB[i];
        acc_sum2_v0pt_right += vec_mean_n_pt_AB[i];
    }
    float sum2_v0pt_right = sigma*acc_sum2_v0pt_right;

    // Printing v0(pT) sum rules
    cout << "------------------ v0(pT) ------------------" << endl;
    cout << "Sum rule #1: " << sum1_v0pt << " = 0" << endl;
    cout << "Sum rule #2: " << sum2_v0pt_left << " = " << sum2_v0pt_right << endl << endl;

    // Calculating v0(pT)/v0
    vector<float> vec_sv0pt(nBins, 0.0);
    vector<float> vec_unc_sv0pt(nBins, 0.0);
    for (int i=0; i<nBins; i++){
        vec_sv0pt[i] = vec_v0pt[i]/v0;
        vec_unc_sv0pt[i] = abs((vec_unc_v0pt[i]/v0)-(vec_v0pt[i]*unc_v0));
    }

    // Calculating v0(pT)/v0 sum rules
    float sum1_sv0pt = 0;
    float sum2_sv0pt_left = 0;
    float acc_sum2_sv0pt_right = 0;
    for (int i=0; i<nBins; i++){
        sum1_sv0pt += vec_sv0pt[i]*vec_mean_n_pt_AB[i];
        float pT = pT_axis[i];
        sum2_sv0pt_left += pT*vec_sv0pt[i]*vec_mean_n_pt_AB[i];
        acc_sum2_sv0pt_right += vec_mean_n_pt_AB[i];
    }
    float sum2_sv0pt_right = mean_pt_ref*acc_sum2_sv0pt_right;

    cout << "----------------- v0(pT)/v0 -----------------" << endl;
    cout << "Sum rule #1: " << sum1_sv0pt << " = 0" << endl;
    cout << "Sum rule #2: " << sum2_sv0pt_left << " = " << sum2_sv0pt_right << endl;

    TFile *save_file = new TFile(Savename, "UPDATE");

    if (TargetPID.size() == 1) {
        int pid = *TargetPID.begin(); // É assim que pegamos o único elemento de um set
        if (pid == 211) Name += "_pion";
        else if (pid == 321) Name += "_kaon";
        else if (pid == 2212) Name += "_proton";
    } else {
        Name += "_all";
    }

    TString mean_pt_name = "mean_pt_";
    mean_pt_name += Name;

    TParameter<float> *p_mean_pt = new TParameter<float>(mean_pt_name, mean_pt);
    p_mean_pt->Write(mean_pt_name, TObject::kOverwrite); // Does not depend on pT-ref (only eta-gap)

    vector<float> vec_zeros(nBins, 0.0);

    // v0(pT)v0
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

    // v0(pT)
    TGraphErrors* gr_v0pt = new TGraphErrors(nBins, pT_axis.data(), vec_v0pt.data(), vec_zeros.data(), vec_unc_v0pt.data());
    TString v0pt_name = "v0pt_";
    v0pt_name += Name;
    gr_v0pt->SetName(v0pt_name);
    gr_v0pt->Write();

    // v0(pT)/v0
    TGraphErrors* gr_sv0pt = new TGraphErrors(nBins, pT_axis.data(), vec_sv0pt.data(), vec_zeros.data(), vec_unc_sv0pt.data());
    TString sv0pt_name = "sv0pt_";
    sv0pt_name += Name;
    gr_sv0pt->SetName(sv0pt_name);
    gr_sv0pt->Write();

    cout << "mean_pt = " << mean_pt << endl;
    save_file->Close();

    TString Savename_txt = Savename;
    Savename_txt.ReplaceAll(".root", ".txt");

    bool txt_file_exists = !gSystem->AccessPathName(txt_filename.Data());
    ofstream txt_file(txt_filename.Data(), ios::app);
        if (txt_file_exists){
            txt_file << "\n\n";
            txt_file << "--------------------------------------------------\n";
            txt_file << "\n";
        }
        txt_file << "Centrality: " << Name << "\n";
        txt_file << "Total tracks: " << total_tracks << "\n";
        txt_file << "Tracks per bin:\n";
        for (int i = 0; i < tracks_per_bin.size(); i++) {
            txt_file << "  Bin " << i << ": " << tracks_per_bin[i] << " tracks\n";
        }
        txt_file.close();
}