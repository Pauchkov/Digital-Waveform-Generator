#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <vector>
using namespace std;
using namespace std::filesystem;

// Search target file in script directory first, then recursively in subfolders.
static path FindFileInScriptTree(const path& script_dir, const string& filename) {
    path direct_candidate = script_dir / filename;
    if (is_regular_file(direct_candidate)) return weakly_canonical(direct_candidate);

    for (const auto& entry : recursive_directory_iterator(script_dir, directory_options::skip_permission_denied)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().filename() == filename) return weakly_canonical(entry.path());
    }

    return path();
}

// Resolve user-provided input path against script location and local tree.
static path ResolveInputInScriptTree(const path& script_dir, const string& input_name) {
    path requested(input_name);

    if (requested.is_absolute() && is_regular_file(requested)) return weakly_canonical(requested);

    if (is_regular_file(requested)) return weakly_canonical(requested);

    path from_script_dir = script_dir / requested;
    if (is_regular_file(from_script_dir)) return weakly_canonical(from_script_dir);

    return FindFileInScriptTree(script_dir, requested.filename().string());
}

// Escape backslashes/quotes so path can be safely injected in ROOT command strings.
static string EscapeForRootString(const string& value) {
    string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
        if (c == '\\' || c == '"') escaped.push_back('\\');
        escaped.push_back(c);
    }
    return escaped;
}

// Timestamp suffix used to create unique output folders per run.
static string BuildTimestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm tm_snapshot{};
#if defined(_WIN32)
    localtime_s(&tm_snapshot, &now_time);
#else
    localtime_r(&now_time, &tm_snapshot);
#endif

    std::ostringstream ts;
    ts << std::put_time(&tm_snapshot, "%Y%m%d_%H%M%S");
    return ts.str();
}

void Digital_Waveform_Generator(string sim_input_name, bool use_mean = true) { // insert input root filename with evt_quadrant_time

   // Locate auxiliary files near the macro and resolve simulation input.

   path script_dir = weakly_canonical(absolute(path(__FILE__))).parent_path();
   path shape_file = FindFileInScriptTree(script_dir, "norm_1pe.root");
   path dc_results_file = FindFileInScriptTree(script_dir, "DC_results.txt");
   path sim_file_path = ResolveInputInScriptTree(script_dir, sim_input_name);

   if (shape_file.empty()) {
	cout << "Unable to locate norm_1pe.root near script path: " << script_dir << endl;
	exit(1);
   }
   if (dc_results_file.empty()) {
	cout << "Unable to locate DC_results.txt near script path: " << script_dir << endl;
	exit(1);
   }
   if (sim_file_path.empty()) {
	cout << "Unable to locate input file " << sim_input_name << " near script path: " << script_dir << endl;
	exit(1);
   }

	   int n_fl_n=100; // number of files 
	   int n_evt =1; //number of events per file from the simulation (I generally create many file with just 1 event) 
	   const int n_channels = 16; // number of SiPM channels to process
	   string sim_stem = sim_file_path.stem().string();
	   const string run_stamp = BuildTimestamp();
	   const string input_root_name = sim_file_path.filename().string();
	   path output_dir = script_dir / (input_root_name + "_" + run_stamp);
	   std::error_code output_dir_error;
	   create_directories(output_dir, output_dir_error);
	   if (output_dir_error) {
		cout << "Unable to create output directory " << output_dir << ": " << output_dir_error.message() << endl;
		exit(1);
	   }
	   cout << "Output directory: " << output_dir << endl;

	   // If requested, pre-average the event histogram before waveform generation.
	   path mean_output_file = output_dir / ("mean_" + sim_stem + ".root");
   path processing_file = sim_file_path;
   if (use_mean) {
	path mean_macro_file = FindFileInScriptTree(script_dir, "Mean.C");
	if (mean_macro_file.empty()) {
		cout << "Unable to locate Mean.C near script path: " << script_dir << endl;
		exit(1);
	}

	cout << "Running Mean.C for input: " << sim_file_path << endl;
	cout << "Mean output file: " << mean_output_file << endl;
	string load_mean_cmd = ".L " + mean_macro_file.string();
	gROOT->ProcessLine(load_mean_cmd.c_str());
	string run_mean_cmd = "Mean(\"" + EscapeForRootString(sim_file_path.string()) + "\",\"" + EscapeForRootString(mean_output_file.string()) + "\")";
	gROOT->ProcessLine(run_mean_cmd.c_str());
	if (!is_regular_file(mean_output_file)) {
		cout << "Mean.C did not produce output file: " << mean_output_file << endl;
		exit(1);
	}
	processing_file = mean_output_file;
   } else {
	cout << "Skipping Mean.C, using input file directly: " << sim_file_path << endl;
   }
   
    TFile *f1 = TFile::Open(shape_file.string().c_str()); // dark count measurement for the shape -- root file with the 1 photon-electron waveform (necessary to create the wf)
    if (!f1 || f1->IsZombie()) {
    	cout << "Unable to open file " << shape_file << endl;
    	exit(1);
    }
    
    ifstream inFile(dc_results_file.string());	//dark count measurements for the fit results (amplitude) -- txt file with infos about the SiPMs
    
    
    double max_ampl[2][n_channels];
    double dark_count[n_channels];
    double cross_talk[n_channels];
    if (!inFile) {
    	cout << "Unable to open file ";
    	exit(1);   // call system to stop
    }
    int J=0;
    while (J < n_channels && (inFile >> max_ampl[0][J] >> max_ampl[1][J] >> dark_count[J] >> cross_talk[J])) {
    	cout<<"j "<<J<<" "<<max_ampl[0][J]<<" "<<max_ampl[1][J]<<" "<<dark_count[J]<<" "<<cross_talk[J]<<endl;	
    	J++;
    }
    if (J != n_channels) {
    	cout << "Unexpected number of lines in " << dc_results_file << ". Expected " << n_channels << ", got " << J << endl;
    	exit(1);
    }
    
  
    gROOT->SetBatch(kTRUE);
    
    // ROOT file OUT 
    string outfile_name = (output_dir / ("out_" + sim_stem + "_wfs.root")).string();
    cout<<outfile_name<<endl;
    TFile *outfile  = TFile::Open(outfile_name.c_str(),"RECREATE");
    
	 // to store all the average waveforms for each channel   
    TMultiGraph*g_av = new TMultiGraph();
    auto* leg = new TLegend(0.5,0.6,0.9,0.9);
    leg->SetNColumns(4);

    double x_av[1024];
    double y_av[n_channels][1024]={};
    int count_evt=0; //just to check

    ////////////////////////////////////////////////
    //// creating hist for 1pe
    ////////////////////////////////////////////////
    
    vector<TGraph*> g1; 
    vector<TH1F*> h;

   
    for(int ch=0;ch<n_channels;ch++){ 
    TString g1_name(Form("norm_1pe__%02d",ch));  //1pe wf
    TGraph*g = (TGraph*)f1->Get(g1_name);
    if (!g) {
    	cout << "Unable to find graph " << g1_name << " in file " << shape_file << endl;
    	exit(1);
    }
  
    g1.push_back(g);

    int g1_n = g1[ch]->GetN();
    vector<double> tg1(g1_n);
    vector<double> yg1(g1_n);
    
    
    for (int i=0;i<g1_n;i++){ 
    	tg1[i]=g1[ch]->GetPointX(i);
    	yg1[i]=g1[ch]->GetPointY(i);
    }
     
     TH1F* h1 = new TH1F(Form("h_%02d",ch),"h",g1_n,tg1[0],tg1[g1_n-1]); //// the hist is shorter to save time
     h.push_back(h1);
     
     for(int i=1;i<g1_n;i++){
     	for(int j=0;j<int(g1[ch]->GetPointY(i)*100);j++){
     		h[ch]->Fill(g1[ch]->GetPointX(i));
    	}
     }
    
       
    // Restore original scale after integerized fill operation above.
    h[ch]->Scale(1/100.);
    outfile->WriteObject(h[ch],"1pe");	    
    }
    ////////////////////////////////////////////////
    //// ends of creating hist for 1pe
    ////////////////////////////////////////////////
    
    // Adding uncertanties on the peak high and noise
    // Use one RNG API with an explicit seed for reproducibility.
    const unsigned int rng_seed = 20260310u;
    std::mt19937 rng(rng_seed);
    cout << "RNG seed: " << rng_seed << endl;
    double noise_lvl =0.6;
    std::uniform_int_distribution<int> dc_bin_distribution(1, 1024);
    std::uniform_real_distribution<double> unit_distribution(0.0, 1.0);
    std::normal_distribution<double> noise_distribution(0.0, noise_lvl);
	

	    
    double inter[n_channels] = {};
    int bin_max[n_channels] = {}; 
    double max[n_channels]; 
    
    for(int ch=0;ch<n_channels;ch++){
		max[ch] = h[ch]->GetMaximum(); 

      for(int i=0; i<1024;i++) if(max[ch] == h[ch]->GetBinContent(i+1)) bin_max[ch] = i; //searching for the max of h copy of 1pe 
      // Integrate around the 1pe peak, used later to normalize waveform amplitude.
      for(int i=bin_max[ch]-5*1024/320;i<bin_max[ch]+9*1024/320;i++) inter[ch] += h[ch]->GetBinContent(i+1);  //calculating area h copy of 1pe 
    }

    
    vector<string> sim_files;
    sim_files.push_back(processing_file.string());
    n_fl_n = sim_files.size();

    
	for(int n_fl=0;n_fl<n_fl_n;n_fl++){
	 	cout<<n_fl<<endl;	
		string file_sim_name = sim_files[n_fl];
		cout<<"File tot: "<<file_sim_name<<endl;
		cout<<endl<<endl<<endl;
		TFile *f = TFile::Open(file_sim_name.c_str());
		if (!f || f->IsZombie()) {
			cout << "Unable to open simulation file " << file_sim_name << endl;
			exit(1);
		}
			 		     
		string sim_name = "evt_quadrant_time";
		TH3 *h_all = (TH3*)f->Get(sim_name.c_str()); //histo from the simulation
		if (!h_all) {
			cout << "Unable to find histogram " << sim_name << " in file " << file_sim_name << endl;
			exit(1);
		}
		const bool is_mean_hist = (h_all->GetNbinsX() == 1);


		int DC =0;
		for(int evt=0;evt<n_evt;evt++){
			count_evt++;
			string name_cv = "event "; 
			name_cv += to_string(n_evt*n_fl_n+evt);
				
		     	for(int ch=0;ch<n_channels;ch++) {
					double time0 = 128;  // 40*1024/320
			     	TH1F*h_sim1 = new TH1F("h sim1","h sim1",1024,0,320);
			     	h_sim1->GetXaxis()->SetTitle("Arrival time [ns]");
			     	h_sim1->GetYaxis()->SetTitle("# of photons");
			     	h_sim1->SetTitle("Output of the simulation - arrival times");
						      		
			     	for(int i=0;i<1024;i++) {
			     		double photons = h_all->GetBinContent(evt+1,ch+1,i+1);  // photons per event per channel per time-bin
			     		// If histogram already contains means, sample integer photon counts.
			     		if (is_mean_hist) {
			     			if (photons < 0.0) photons = 0.0;
			     			std::poisson_distribution<int> mean_to_counts(photons);
			     			h_sim1->AddBinContent(i+1, mean_to_counts(rng));
			     		} else {
			     			h_sim1->AddBinContent(i+1, photons);
			     		}
			     	}
			 		
				 		//Adding dark count
	  				std::poisson_distribution<int> distribution(dark_count[ch]);
					std::normal_distribution<double> pe_ampl_distribution(max_ampl[0][ch], max_ampl[1][ch]);
			
					int rnd_DC = distribution(rng);
					double CT_prob = cross_talk[ch];
					if(rnd_DC>0){
						for(int i_DC=0;i_DC<rnd_DC;i_DC++){
							DC++;
							int rnd_i = dc_bin_distribution(rng);
							h_sim1->AddBinContent(rnd_i,1);
						
						//CROSSTALK PHOTONS for dark count
				     		double rnd_CT = unit_distribution(rng);
						if( rnd_CT < CT_prob) h_sim1->AddBinContent(rnd_i,1);
						}
					}
			
				   string name_sim1 = to_string(ch); 
					name_sim1 += "_";
					name_sim1 += to_string(evt);
			   
					TH1F*h_rep2 = new TH1F("h rep2","h rep2",1024,0,320);
					h_rep2->GetXaxis()->SetTitle("Time [ns]");
					h_rep2->GetYaxis()->SetTitle("Amplitude [mV]");
					h_rep2->SetTitle("Simulated waveform for an event for a ch");
					h_rep2->SetLineWidth(2);
				
					//adding up 1pe copy hist using the photons from the simulation     
					for(int i=0;i<1024;i++){  
						for(int j=0; j<h_sim1->GetBinContent(i+1);j++) {
							TH1F* h_err = (TH1F*) h[ch]->Clone(); 
						
							// Per-photon amplitude fluctuates around measured 1pe response.
							h_err->Scale(pe_ampl_distribution(rng)/inter[ch]); //rescaling the h copy of 1pe
		
							for(int y=time0;y<g1[ch]->GetN()-i+time0;y++){
								h_rep2->AddBinContent(y+i+1,h_err->GetBinContent(y-time0+1));
							}
				 		delete h_err;
				  	}
				}
				
				// Add electronics noise after waveform synthesis.
				for (int i = 1; i <= h_rep2->GetNbinsX(); i++) h_rep2->SetBinContent(i, h_rep2->GetBinContent(i) + noise_distribution(rng)); //adding noise
				for(int i=0;i<1024;i++){
					x_av[i]=320./1024.*i;
					y_av[ch][i]+=h_rep2->GetBinContent(i+1);
				}
				
				
				string h_sim_name = "wf_ch"; 
				h_sim_name += to_string(ch).c_str();
				h_sim_name += "_"; 
				h_sim_name += to_string(n_fl);
			  	outfile->WriteObject(h_rep2,h_sim_name.c_str());	 //writing the waveform    
				             
				delete h_sim1; 
				delete h_rep2;
		 		
			} //end of the cicle channels 
			
		} 
		
		f ->Close(); //end of the cicle events

 
   } //end cicle files
   // Store average waveform per channel across processed events/files.
   for(int ch=0;ch<n_channels;ch++){
	auto* g = new TGraph(1024,x_av,y_av[ch]);
	g->SetLineColor(99-ch);
	g->SetLineWidth(2);
	g_av->Add(g);
	leg->AddEntry(g,to_string(ch).c_str());
   }
    
    string name_c_av = "Average_wfs";  
    TCanvas*c_av = new TCanvas(name_c_av.c_str(),name_c_av.c_str(),12,14,1000,800);
    g_av->Draw("ALP");
    leg->Draw("SAME");

    cout<<"Number of events converted: "<<count_evt<<endl;
    outfile->WriteObject(c_av,name_c_av.c_str());
    
    //printing a average wfs plot to check the results

    string image_name = (output_dir / ("Average_" + sim_stem + ".png")).string();
    c_av->SaveAs(image_name.c_str());

    gROOT->SetBatch(kFALSE);
    f1->Close();
    outfile->Close(); 
   }
