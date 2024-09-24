#ifndef __OSCILLATOR_OSCPROB_H__
#define __OSCILLATOR_OSCPROB_H__

#include "OscProbCalcerBase.h"

/**
 * @file OscProbCalcer_OscProb.h
 *
 * @class OscProbCalcerOscProb
 *
 * @brief Oscillation calculation engine for linear propagation in NuFAST.
 */
class OscProbCalcerOscProb : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param ConfigName_ Name of config used to setup the OscProbCalcerOscProb() instance
   * @param Instance_ Which entry of the OscProbCalcerSetup config block should be read in the case where there are multiple OscProbCalcers to be initialised
   */
  OscProbCalcerOscProb(std::string ConfigName_="", int Instance_=0);

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerOscProb();

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup NuFAST specific variables
   */  
  void SetupPropagator() override;
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in NuFAST. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) override;

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param OscChanIndex The index in #fOscillationChannels to return the pointer for 
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1) override;
  
  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because NuFAST is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPRODH, kNOscParams};
  
  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};
  
};

#endif
