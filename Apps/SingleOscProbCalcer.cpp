#include "Constants/OscillatorConstants.h"
#include "OscProbCalcer/OscProbCalcerBase.h"

#include "OscProbCalcer/OscProbCalcerFactory.h"

#include <iostream>
#include <math.h>
#include <chrono>

#include "TCanvas.h"
#include "TH2.h"

#include <sys/resource.h>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << argv[0] << " InputConfig.yaml [key=value ...]" << std::endl;
    std::cerr << "  keys: sin2theta12 sin2theta23 sin2theta13 dm21 dm32 dcp rho Ye gamma0 n E0" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  std::string OscProbCalcerConfigname = argv[1];

  // CLI overrides for OscParams_Beam_wYe_wDeco. Index 6 (baseline) is loop-controlled
  // below, so it's intentionally absent from the key list.
  struct KeyIdx { const char* key; int idx; };
  static const KeyIdx KEYS[] = {
    {"sin2theta12", 0}, {"sin2theta23", 1}, {"sin2theta13", 2},
    {"dm21", 3},        {"dm32", 4},        {"dcp", 5},
    {"rho", 7},         {"Ye", 8},
    {"gamma0", 9},      {"n", 10},          {"E0", 11},
  };
  std::vector<std::pair<int,double>>            Overrides;
  std::vector<std::pair<std::string,std::string>> OverrideTags;
  for (int ai = 2; ai < argc; ++ai) {
    std::string tok = argv[ai];
    size_t eq = tok.find('=');
    if (eq == std::string::npos) {
      std::cerr << "Bad arg '" << tok << "' (expected key=value)" << std::endl;
      throw std::runtime_error("Invalid setup");
    }
    std::string k = tok.substr(0, eq);
    std::string v = tok.substr(eq + 1);
    int idx = -1;
    for (auto const& ki : KEYS) {
      if (k == ki.key) { idx = ki.idx; break; }
    }
    if (idx < 0) {
      std::cerr << "Unknown key '" << k << "'. Valid keys:";
      for (auto const& ki : KEYS) std::cerr << " " << ki.key;
      std::cerr << std::endl;
      throw std::runtime_error("Invalid setup");
    }
    Overrides.emplace_back(idx, std::stod(v));
    OverrideTags.emplace_back(k, v);
  }

  // Filename tag built from overridden keys; empty when no overrides → original names.
  std::string Tag;
  for (auto const& kv : OverrideTags) Tag += "_" + kv.first + "=" + kv.second;

  bool PrintWeights = true;

  bool useLogEnergy = true;
  float E_min = 0.1;
  float E_max = 10000.01;
  int E_Nbins = 1001; // E_Nbins equally spaced in log

  bool useLogBaseline = false;
  float BL_Min=10.0;
  float BL_Max=1510.0;

  int BL_Nbins = 1001; // Number of baselines
  //int BL_Nbins = 2; // when plottting single baseline fast (also change "baseline = ..." below)

//  std::vector<FLOAT_T> EnergyArray = logspace(E_min, E_max, E_Nbins);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,15);

  std::vector<double> EnergyArray(E_Nbins);
  std::vector<double> E_Edges(E_Nbins+1);

  if (useLogEnergy) {
  double logE_min = std::log10(E_min);
  double logE_max = std::log10(E_max);
  double dlogE = (logE_max - logE_min)/(E_Nbins-1);

  for (int i=0; i<E_Nbins; ++i) {
    EnergyArray[i] = std::pow(10, logE_min + i*dlogE);
  }

  // Compute bin edges from bin centers
  for (int i=0; i<E_Nbins-1; ++i) {
    E_Edges[i+1] = std::sqrt(EnergyArray[i]*EnergyArray[i+1]); // geometric mean
  }
  E_Edges[0] = EnergyArray[0]*EnergyArray[0]/E_Edges[1];           // first edge
  E_Edges[E_Nbins] = EnergyArray[E_Nbins-1]*EnergyArray[E_Nbins-1]/E_Edges[E_Nbins-1]; // last edge

  } else {
    // --- Linear-spaced energy ---
    double dE = (E_max - E_min)/(E_Nbins-1);

    for (int i=0; i<E_Nbins; ++i) {
      EnergyArray[i] = E_min + i*dE;
    }

    // Bin edges, half-bin outside first/last center
    double E_MinEdge = E_min - dE/2.0;
    for (int i=0; i<=E_Nbins; ++i) {
      E_Edges[i] = E_MinEdge + i*dE;
    }
  }

  std::cout << "EnergyArray.size() =" << EnergyArray.size() << std::endl;
  std::cout << "EnergyArray[0] =" << EnergyArray[0] << std::endl;

  std::vector<double> BL_Edges(BL_Nbins+1);

  if (useLogBaseline) {
    // --- Log-spaced baseline (centers + geometric mean for edges) ---
    std::vector<double> BLLogArray(BL_Nbins);
    double logBL_Min = std::log10(BL_Min);
    double logBL_Max = std::log10(BL_Max);
    double dlogBL = (logBL_Max - logBL_Min)/(BL_Nbins-1);

  for (int i=0; i<BL_Nbins; ++i) {
    BLLogArray[i] = std::pow(10, logBL_Min + i*dlogBL);
  }

  for (int i=0; i<BL_Nbins-1; ++i) {
    BL_Edges[i+1] = std::sqrt(BLLogArray[i]*BLLogArray[i+1]); // geometric mean
  }
  BL_Edges[0] = BLLogArray[0]*BLLogArray[0]/BL_Edges[1];
  BL_Edges[BL_Nbins] = BLLogArray[BL_Nbins-1]*BLLogArray[BL_Nbins-1]/BL_Edges[BL_Nbins-1];

  } else {

    double dBL = (BL_Max - BL_Min)/(BL_Nbins-1);

    // First and last edges are half a bin outside the first/last center
    double BL_MinEdge = BL_Min - dBL/2.0;
    double BL_MaxEdge = BL_Max + dBL/2.0;

    for (int i=0; i<=BL_Nbins; ++i) {
      BL_Edges[i] = BL_MinEdge + i*dBL;
    }
  }

  std::vector<FLOAT_T> OscParams_Basic = ReturnOscParams_Basic();
  std::vector<FLOAT_T> OscParams_Atm = ReturnOscParams_Atm();
  std::vector<FLOAT_T> OscParams_Beam_woYe = ReturnOscParams_Beam_woYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe = ReturnOscParams_Beam_wYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe_wDeco = ReturnOscParams_Beam_wYe_wDeco();
  for (auto const& ov : Overrides) OscParams_Beam_wYe_wDeco[ov.first] = ov.second;
  std::vector<FLOAT_T> OscParams_Beam_wYe_wLIV = ReturnOscParams_Beam_wYe_wLIV();

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::cout << "========================================================" << std::endl;
  std::cout << "Initialising " << OscProbCalcerConfigname << std::endl;

  OscProbCalcerFactory* OscProbCalcFactory = new OscProbCalcerFactory();
  OscProbCalcerBase* Calcer = OscProbCalcFactory->CreateOscProbCalcer(OscProbCalcerConfigname);
  delete OscProbCalcFactory;
  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  Calcer->SetEnergyArray(EnergyArray);
  if (!Calcer->ReturnCosineZIgnored()) {
    Calcer->SetCosineZArray(CosineZArray);
  }
  Calcer->Setup();

  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  struct rusage usage; // For memory check

  //Don't plot by default
  bool Plot = true;

  if (Plot) {
    TCanvas* Canv = new TCanvas;
    TString OutputName     = TString("Probability") + Tag.c_str() + ".pdf";
    TString OutputRootName = TString("Probability") + Tag.c_str() + ".root";
    Canv->Print(OutputName+"[");

    TFile f(OutputRootName,"RECREATE");

    auto NuFlavGreek = [](int flav){
        switch(flav) {
            case 1: return "#nu_{e}";
            case 2: return "#nu_{#mu}";
            case 3: return "#nu_{#tau}";
            default: return "#nu_{?}";
        }
    };

  // Pre-create histograms in a vector
  std::vector<TH2D*> hists;

  for (int iNuType = 0; iNuType < 2; iNuType++) {
    int NuType = (iNuType == 1) ? -1 : 1;

    for (int iGenFlav = 1; iGenFlav < 4; iGenFlav++) {
      for (int iDetFlav = 1; iDetFlav < 4; iDetFlav++) {
        int GenFlav = iGenFlav * NuType;
        int DetFlav = iDetFlav * NuType;

        TString nameT = Form("2D_OscProb_%i_%i_%i", NuType, GenFlav, DetFlav);
        TString title = TString(Form("P_{%s#rightarrow%s};Energy [GeV];Baseline [km]",
            NuFlavGreek(std::abs(GenFlav)),
            NuFlavGreek(std::abs(DetFlav))));

        TH2D* h = new TH2D(nameT, title,
            E_Nbins, E_Edges.data(),
            BL_Nbins, BL_Edges.data());
	h->SetDirectory(&f);
        hists.push_back(h);
      }
    }
  }

  auto histIndex = [](int NuType, int GenFlav, int DetFlav) {
    int nuIndex = (NuType > 0) ? 0 : 1;  // +1 → 0, -1 → 1
    int g = GenFlav - 1;                 // 1..3 → 0..2
    int d = DetFlav - 1;                 // 1..3 → 0..2
    return nuIndex * 9 + g * 3 + d;      // 2*3*3 = 18 histograms
  };

  int lastNuType = -1;
  int lastGenFlav = -1;
  int lastDetFlav = -1;
  int xBin = 1;

  for (int iBL=0;iBL<BL_Nbins;iBL++) {
  int yBin = iBL+1;

  // Set baseline
  double baseline = 0.5 * (BL_Edges[iBL] + BL_Edges[iBL+1]);
  //baseline = 295.0; // T2K
  //baseline = 1284.9; // DUNE
  OscParams_Beam_wYe_wDeco[6] = baseline;

  for (int i=0; i < OscParams_Beam_wYe_wDeco.size(); i++) {
    std::cout << "OscParams_Beam_wYe_wDeco[ " << i << " ] =" << OscParams_Beam_wYe_wDeco[i] << std::endl;
  }

  // Reweight and calculate oscillation probabilities

  // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
  // If you have some NSO calculater, then it will work providing the length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
  if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
    Calcer->Reweight(OscParams_Beam_woYe);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
    Calcer->Reweight(OscParams_Beam_wYe);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_wYe_wDeco.size()) {
    Calcer->Reweight(OscParams_Beam_wYe_wDeco); 
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_wYe_wLIV.size()) {
    Calcer->Reweight(OscParams_Beam_wYe_wLIV);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Atm.size()) {
    Calcer->Reweight(OscParams_Atm);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Basic.size()) {
    Calcer->Reweight(OscParams_Basic);
  } else {
    std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
    std::cerr << "Oscillator->ReturnNOscParams():" << Calcer->ReturnNOscParams() << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (PrintWeights) {
    std::vector<NuOscillator::OscillationProbability> OscProbs = Calcer->ReturnProbabilities();
    for (int iOscProb=0;iOscProb<(int)OscProbs.size();iOscProb++) {
      std::cout << iOscProb << " " << OscProbs[iOscProb].NuType << " " << OscProbs[iOscProb].OscChan.GeneratedFlavour << " " << OscProbs[iOscProb].OscChan.DetectedFlavour << " " << OscProbs[iOscProb].Energy << " " << OscProbs[iOscProb].CosineZ << " " << baseline << " " << OscProbs[iOscProb].Probability << std::endl;

      int idx = histIndex( OscProbs[iOscProb].NuType, OscProbs[iOscProb].OscChan.GeneratedFlavour, OscProbs[iOscProb].OscChan.DetectedFlavour );
      TH2D* Hist = hists[idx];

    // Check if the flavor combination has changed
    if (OscProbs[iOscProb].NuType != lastNuType || 
        OscProbs[iOscProb].OscChan.GeneratedFlavour != lastGenFlav || 
        OscProbs[iOscProb].OscChan.DetectedFlavour != lastDetFlav) 
    {
        xBin = 1;  // reset x-bin
        lastNuType = OscProbs[iOscProb].NuType;
        lastGenFlav = OscProbs[iOscProb].OscChan.GeneratedFlavour;
        lastDetFlav = OscProbs[iOscProb].OscChan.DetectedFlavour;
    }


     std::cout << "xBin =" << xBin << std::endl;
     std::cout << "yBin =" << yBin << std::endl;
     std::cout << "iBL =" << iBL << std::endl;
     std::cout << "EnergyArray[xBin] =" << EnergyArray[xBin-1] << std::endl;
     std::cout << "baseline =" << baseline << std::endl;

     std::cout << "lastNuType =" << lastNuType << std::endl;
     std::cout << "lastGenFlav =" << lastGenFlav << std::endl;
     std::cout << "lastDetFlav =" << lastDetFlav << std::endl;

      Hist->SetBinContent(xBin,yBin,OscProbs[iOscProb].Probability);

      // Check memory
      getrusage(RUSAGE_SELF, &usage); std::cout << "Memory (MB): " << usage.ru_maxrss/1024.0 << std::endl;
      ++xBin;
    }
  }

  } // Baseline loop

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;

  for ( auto &h : hists ) {
    Canv->SetLogx(useLogEnergy);
    Canv->SetLogy(useLogBaseline);

    h->SetStats(kFALSE);
    h->Write(); // Write to ROOT file
    h->Draw("COLZ");
    h->GetXaxis()->SetTitleOffset(1.3);
    Canv->Print(OutputName); // 2D plot to pdf
    //Canv->Write(); // 2D canvas plot to root

    int yBin = h->GetYaxis()->FindBin(295.0); // DUNE baseline is 1284.9; T2K baseline in 295.0
    TH1D* hEnergySlice = h->ProjectionX(Form("%s_Eslice", h->GetName()), yBin, yBin);
    Canv->SetLogx(useLogEnergy);
    Canv->SetLogy(false);

    // Optional: set title
    hEnergySlice->SetTitle(Form("%s at baseline = %.1f km; Energy [GeV]; %s",
    	h->GetTitle(), h->GetYaxis()->GetBinCenter(yBin), h->GetTitle()));
    hEnergySlice->GetXaxis()->SetTitleOffset(1.3);

    hEnergySlice->SetStats(kFALSE);
    hEnergySlice->Write();
    hEnergySlice->Draw();
    Canv->Print(OutputName);
    delete hEnergySlice;

    int xBin = h->GetXaxis()->FindBin(0.6);
    TH1D* hBaselineSlice = h->ProjectionY(Form("%s_Bslice", h->GetName()), xBin, xBin);
    Canv->SetLogx(useLogBaseline);
    Canv->SetLogy(false);

    // Optional: set title
    hBaselineSlice->SetTitle(Form("%s at energy = %.1f GeV; Baseline [km]; %s",
    	h->GetTitle(), h->GetXaxis()->GetBinCenter(xBin), h->GetTitle()));
    hBaselineSlice->GetXaxis()->SetTitleOffset(1.3);

    hBaselineSlice->SetStats(kFALSE);
//    hBaselineSlice->GetYaxis()->SetRangeUser(0.0, 1.0); // force y-axis from 0 to 1
    hBaselineSlice->Write();
    hBaselineSlice->Draw();
    Canv->Print(OutputName);
    delete hBaselineSlice;
  }
  Canv->Print(OutputName+"]");
  f.Close();
} // Plot

}

