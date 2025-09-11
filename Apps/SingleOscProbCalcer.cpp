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
  if (argc != 2) {
    std::cerr << argv[0] << " InputConfig.yaml" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  std::string OscProbCalcerConfigname = argv[1];

  bool PrintWeights = true;

  float E_min = 1e-1;
  float E_max = 1.0e4;
  int E_Nbins = 1001; // E_Nbins equally spaced in log

//  std::vector<FLOAT_T> EnergyArray = logspace(E_min, E_max, E_Nbins);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,15);

std::vector<double> EnergyArray(E_Nbins);
double logE_min = std::log10(E_min);
double logE_max = std::log10(E_max);
double dlogE = (logE_max - logE_min)/(E_Nbins-1);

for (int i=0; i<E_Nbins; ++i) {
    EnergyArray[i] = std::pow(10, logE_min + i*dlogE);
}

// Compute bin edges from bin centers
std::vector<double> E_edges(E_Nbins+1);
for (int i=0; i<E_Nbins-1; ++i) {
    E_edges[i+1] = std::sqrt(EnergyArray[i]*EnergyArray[i+1]); // geometric mean
}
E_edges[0] = EnergyArray[0]*EnergyArray[0]/E_edges[1];           // first edge
E_edges[E_Nbins] = EnergyArray[E_Nbins-1]*EnergyArray[E_Nbins-1]/E_edges[E_Nbins-1]; // last edge

  std::cout << "EnergyArray.size() =" << EnergyArray.size() << std::endl;
  std::cout << "EnergyArray[0] =" << EnergyArray[0] << std::endl;

  float BL_Min=100.0;
  float BL_Max=1500.0;
  int BL_num = 1001; // Number of baselines

  // Now compute the histogram edges so that bin centers = these values
  double BL_MinEdge = BL_Min - ( (BL_Max - BL_Min) / (BL_num - 1) ) / 2.0;
  double BL_MaxEdge = BL_Max + ( (BL_Max - BL_Min) / (BL_num - 1) ) / 2.0;

  //std::vector<FLOAT_T> BaselineArray = linspace(BL_Min,BL_Max,BL_num);

//  std::cout << "BaselineArray.size() =" << BaselineArray.size() << std::endl;
//  std::cout << "BaselineArray[0] =" << BaselineArray[0] << std::endl;

  std::vector<FLOAT_T> OscParams_Basic = ReturnOscParams_Basic();
  std::vector<FLOAT_T> OscParams_Atm = ReturnOscParams_Atm();
  std::vector<FLOAT_T> OscParams_Beam_woYe = ReturnOscParams_Beam_woYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe = ReturnOscParams_Beam_wYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe_wDeco = ReturnOscParams_Beam_wYe_wDeco();
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



struct rusage usage;

  //Don't plot by default
  bool Plot = true;

  if (Plot) {
    TCanvas* Canv = new TCanvas;
    TString OutputName = "Probability.pdf";
    Canv->Print(OutputName+"[");
    Canv->SetLogx(true);

    TFile f("Probability.root","RECREATE");

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
//for (int NuType = 0; NuType < 2; ++NuType) {
//    for (int GenFlav = 0; GenFlav < 3; ++GenFlav) {
//        for (int DetFlav = 0; DetFlav < 3; ++DetFlav) {


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
//                               EnergyArray.size()-1, EnergyArray.data(),
                   E_Nbins, E_edges.data(),
                   BL_num, BL_MinEdge, BL_MaxEdge);
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

  for (int iBL=0;iBL<BL_num;iBL++) {
//    FLOAT_T BL = iBL*10;

  int yBin = iBL+1;

  double baseline = BL_MinEdge + (iBL + 0.5) * ( (BL_MaxEdge - BL_MinEdge) / BL_num );
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
getrusage(RUSAGE_SELF, &usage); std::cout << "Memory (MB): " << usage.ru_maxrss/1024.0 << std::endl;
      ++xBin;
    }
  }

  } // Baseline loop

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;

  for ( auto &h : hists ) {
    h->SetStats(kFALSE);
    h->Write(); // Write to ROOT file
    h->Draw("COLZ");
    h->GetXaxis()->SetTitleOffset(1.3);
    Canv->Print("Probability.pdf"); // 2D plot to pdf
    //Canv->Write(); // 2D canvas plot to root

    int yBin = h->GetYaxis()->FindBin(150.9); // DUNE baseline is 1284.9; T2K baseline in 295.0
    TH1D* hEnergySlice = h->ProjectionX(Form("%s_Eslice", h->GetName()), yBin, yBin);

    // Optional: set title
    hEnergySlice->SetTitle(Form("%s at baseline = %.1f km; Energy [GeV]; %s",
    	h->GetTitle(), h->GetYaxis()->GetBinCenter(yBin), h->GetTitle()));
    hEnergySlice->GetXaxis()->SetTitleOffset(1.3);

    hEnergySlice->SetStats(kFALSE);
    hEnergySlice->Write();
    hEnergySlice->Draw();
    Canv->Print("Probability.pdf");
    delete hEnergySlice;

    int xBin = h->GetXaxis()->FindBin(2.0);
    TH1D* hBaselineSlice = h->ProjectionY(Form("%s_Bslice", h->GetName()), xBin, xBin);

    // Optional: set title
    hBaselineSlice->SetTitle(Form("%s at energy = %.1f GeV; Baseline [km]; %s",
    	h->GetTitle(), h->GetXaxis()->GetBinCenter(xBin), h->GetTitle()));
    hBaselineSlice->GetXaxis()->SetTitleOffset(1.3);

    hBaselineSlice->SetStats(kFALSE);
    hBaselineSlice->Write();
    hBaselineSlice->Draw();
    Canv->Print("Probability.pdf");
    delete hBaselineSlice;
  }
  Canv->Print(OutputName+"]");
  f.Close();
} // Plot

}
