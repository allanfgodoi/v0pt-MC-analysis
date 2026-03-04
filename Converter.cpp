// Converts ampt events from .dat to .root file
// USAGE EXAMPLE: Converter("/eos/cms/store/user/sdogra/ampt/inPutFiles/BB/ampt_1M.dat", "/eos/user/a/afloresg/MC-studies/BB/ampt_1M.root");
// Do this to all .dat files and then merge the .root files with hadd

void Converter(TString filename, TString savename){

    ifstream amptFile(filename.Data());
    if (!amptFile.is_open()) std::cerr << "Cannot open AMPT file: " << filename << std::endl;
    
    TFile *f = new TFile(savename, "RECREATE");
    TTree *tree = new TTree("tree", "AMPT Events");

    // Event info
    int eventID, n_particles;
    float b;

    // Track info
    vector<int> pid;
    vector<float> px, py, pz, E;

    // Setting up branches
    tree->Branch("eventID", &eventID, "eventID/I");
    tree->Branch("n_particles", &n_particles, "n_particles/I");
    tree->Branch("b", &b, "b/F");
    tree->Branch("pid", &pid);
    tree->Branch("px", &px);
    tree->Branch("py", &py);
    tree->Branch("pz", &pz);
    tree->Branch("E", &E);

    // AMPT header aux variables
    int dummy, npartP, npartT, nElaP, nInElaP, nElaT, nInElaT;
    float unused, x, y, z, t;

    int count = 0;
    while (amptFile >> eventID >> dummy >> n_particles >> b >> npartP >> npartT >> nElaP >> nInElaP >> nElaT >> nInElaT >> unused){

        pid.clear(); px.clear(); py.clear(); pz.clear(); E.clear();
        for (int i = 0; i < n_particles; i++) {
            // Holder variables for each track
            int h_pid;
            float h_px, h_py, h_pz, h_E;

            amptFile >> h_pid >> h_px >> h_py >> h_pz >> h_E >> x >> y >> z >> t;
            pid.push_back(h_pid);
            px.push_back(h_px);
            py.push_back(h_py);
            pz.push_back(h_pz);
            E.push_back(h_E);
        }

        tree->Fill();
        count++;
        if (count % 1000 == 0) std::cout << "Converting events: " << count << std::endl;
    } 

    // Verifying file reading
    if (amptFile.eof()) {
        std::cout << "SUCCESS: The file was completely read)." << std::endl;
    } else if (amptFile.fail()) {
        std::cout << "ERROR: Loop stopped due to some error!" << std::endl;
        std::cout << "Reading stopped after event: " << count << std::endl;
        
        // Try to read the wrong line
        amptFile.clear(); // Clears error state to read the wrong line
        string wrong_line;
        getline(amptFile, wrong_line);
        std::cout << "Something wrong in the following line: " << wrong_line << std::endl;
    }

    tree->Write();
    f->Close();
}