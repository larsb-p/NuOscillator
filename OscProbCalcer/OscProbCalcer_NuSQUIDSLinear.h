#ifndef __OSCILLATOR_NUSQUIDSLINEAR_H__
#define __OSCILLATOR_NUSQUIDSLINEAR_H__

#include "OscProbCalcerBase.h"

#include "nuSQuIDS/nuSQuIDS.h"
#include "examples/Decoherence/nuSQUIDSDecoh.h"

/**
 * @file OscProbCalcer_NuSQUIDSLinear.h
 *
 * @class OscProbCalcerNuSQUIDSLinear
 *
 * @brief Oscillation calculation engine for linear propagation in NuSQUIDS.
 */
class OscProbCalcerNuSQUIDSLinear : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerNuSQUIDSLinear() instance
   */
  OscProbCalcerNuSQUIDSLinear(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerNuSQUIDSLinear(std::string ConfigName_) : OscProbCalcerNuSQUIDSLinear(YAML::LoadFile(ConfigName_)) {}
  
  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerNuSQUIDSLinear();

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup NuSQUIDS specific variables
   */  
  void SetupPropagator() override;
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in NuSQUIDS. This function both calculates and stores
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
   * This is implementation specific because because NuSQUIDS is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
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
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kELECDENS, kNOscParams};
  
  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  nusquids::nuSQUIDSDecoh nus;
  nusquids::nuSQUIDSDecoh nubars;

  squids::Const units;

  double zenith_angle;
  double integration_step;
  double nus_rel_error;
  double nus_abs_error;
  double nubars_rel_error;
  double nubars_abs_error;
  double nus_gamma_strength;
  double nus_gamma_energy_dependence;
  double nus_gamma_energy_scale;
  double nubars_gamma_strength;
  double nubars_gamma_energy_dependence;
  double nubars_gamma_energy_scale;
};

#endif
