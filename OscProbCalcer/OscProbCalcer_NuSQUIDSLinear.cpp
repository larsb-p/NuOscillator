#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

OscProbCalcerNuSQUIDSLinear::OscProbCalcerNuSQUIDSLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  //ZenithAngle
  if (!Config_["OscProbCalcerSetup"]["ZenithAngle"]) {
    std::cerr << "Expected to find a 'ZenithAngle' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  zenith_angle = Config_["OscProbCalcerSetup"]["ZenithAngle"].as<double>();

  //IntegrationStep
  if (!Config_["OscProbCalcerSetup"]["IntegrationStep"]) {
    std::cerr << "Expected to find a 'IntegrationStep' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  integration_step = Config_["OscProbCalcerSetup"]["IntegrationStep"].as<double>();

  //Errors
  //NusRelativeError
  if (!Config_["OscProbCalcerSetup"]["NusRelativeError"]) {
    std::cerr << "Expected to find a 'NusRelativeError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

//  prem_model = Config_["OscProbCalcerSetup"]["NusRelativeError"].as<std::string>();
  nus_rel_error = Config_["OscProbCalcerSetup"]["NusRelativeError"].as<double>();

  //NusAbsoluteError
  if (!Config_["OscProbCalcerSetup"]["NusAbsoluteError"]) {
    std::cerr << "Expected to find a 'NusAbsoluteError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_abs_error = Config_["OscProbCalcerSetup"]["NusAbsoluteError"].as<double>();

  //NubarsRelativeError
  if (!Config_["OscProbCalcerSetup"]["NubarsRelativeError"]) {
    std::cerr << "Expected to find a 'NubarsRelativeError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_rel_error = Config_["OscProbCalcerSetup"]["NubarsRelativeError"].as<double>();

  //NubarsAbsoluteError
  if (!Config_["OscProbCalcerSetup"]["NubarsAbsoluteError"]) {
    std::cerr << "Expected to find a 'NubarsAbsoluteError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_abs_error = Config_["OscProbCalcerSetup"]["NubarsAbsoluteError"].as<double>();

  //Decoherence setup
  //NusGammaStrength
//  if (!Config_["OscProbCalcerSetup"]["NusGammaStrength"]) {
//    std::cerr << "Expected to find a 'NusGammaStrength' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
//    throw;
//  }

//  nus_gamma_strength = Config_["OscProbCalcerSetup"]["NusGammaStrength"].as<double>();

  //NusGammaEnergyDependence
//  if (!Config_["OscProbCalcerSetup"]["NusGammaEnergyDependence"]) {
//    std::cerr << "Expected to find a 'NusGammaEnergyDependence' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
//    throw;
//  }

//  nus_gamma_energy_dependence = Config_["OscProbCalcerSetup"]["NusGammaEnergyDependence"].as<double>();

  //NusGammaEnergyScale
  if (!Config_["OscProbCalcerSetup"]["NusGammaEnergyScale"]) {
    std::cerr << "Expected to find a 'NusGammaEnergyScale' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_gamma_energy_scale = Config_["OscProbCalcerSetup"]["NusGammaEnergyScale"].as<double>();

  //Nubars
  //NubarsGammaStrength
//  if (!Config_["OscProbCalcerSetup"]["NubarsGammaStrength"]) {
//    std::cerr << "Expected to find a 'NubarsGammaStrength' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
//    throw;
//  }

//  nubars_gamma_strength = Config_["OscProbCalcerSetup"]["NusGammaStrength"].as<double>();

  //NubarsGammaEnergyDependence
//  if (!Config_["OscProbCalcerSetup"]["NubarsGammaEnergyDependence"]) {
//    std::cerr << "Expected to find a 'NubarsGammaEnergyDependence' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
//    throw;
//  }

//  nubars_gamma_energy_dependence = Config_["OscProbCalcerSetup"]["NubarsGammaEnergyDependence"].as<double>();

  //NubarsGammaEnergyScale
  if (!Config_["OscProbCalcerSetup"]["NubarsGammaEnergyScale"]) {
    std::cerr << "Expected to find a 'NubarsGammaEnergyScale' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_gamma_energy_scale = Config_["OscProbCalcerSetup"]["NubarsGammaEnergyScale"].as<double>();

  //BSMModel
  if (!Config_["OscProbCalcerSetup"]["BSMModel"]) {
    std::cerr << "Expected to find a 'BSMModel' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  bsm_model = Config_["OscProbCalcerSetup"]["BSMModel"].as<std::string>();

  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  std::string OscMatrix = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fOscType = PMNS_StrToInt(OscMatrix);
  fNOscParams = GetNOscParams(fOscType);

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcerNuSQUIDSLinear::~OscProbCalcerNuSQUIDSLinear() {
}

void OscProbCalcerNuSQUIDSLinear::SetupPropagator() {

  nusquids::marray<double,1> E_range{fNEnergyPoints};
  for(int i=0; i < fNEnergyPoints; i++){
    E_range[i] = fEnergyArray[i]*units.GeV;
  }

  std::string OscProbCalcerNuSQUIDSLinear::set_bsm_model(std::string BSMModel) {
    if ( BSMModel == "Decoherence" ) {
      return kDecoherence;
    }
  }

  //Here we define the trajectory that the particle follows and the object for more examples
  // of how construct a track and object look body_track example.
  //zenith angle, neutrinos crossing the earth
  double phi = acos(zenith_angle);

  double eps_mutau = 1.0e-2;

  fBSMModel = BSM_StrToInt(bsm_model);

  switch() {
    case kSM:
      nus_pmns = new nusquids::nuSQUIDS(E_range, NuOscillator::kTau, nusquids::neutrino, false);
      nubars_pmns = new nusquids::nuSQUIDS(E_range, NuOscillator::kTau, nusquids::neutrino, false);
      nus_base = nus_pmns;

      //      SetupPropNuSQUIDSSM(fBSMModel);
//      CalcProbNuSQUIDSSM(fBSMModel);
      break;

    case kDecoherence:
      nus_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::neutrino, true);
      nubars_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::antineutrino, false); // anti-neutrinos
      nus_base = nus_decoh;

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus_decoh->Set_h_max( integration_step*units.km );
  nubars_decoh->Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus_decoh->Set_rel_error(nus_rel_error);
  nus_decoh->Set_abs_error(nus_abs_error);
  nubars_decoh->Set_rel_error(nubars_rel_error);
  nubars_decoh->Set_abs_error(nubars_abs_error);

  //Set the decoherence model and parameters for neutrinos
  nus_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV); // reference (default) value: nus_gamma_strength = 9.48e-18 
  nus_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence); // reference (default) value: nus_gamma_energy_dependence = 2
  nus_decoh->Set_DecoherenceGammaEnergyScale(nus_gamma_energy_scale*units.GeV); // reference (default) value: nus_gamma_energy_dependence = 1.0

  //Set the decoherence model and parameters for anti-neutrinos
  nubars_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV);
  nubars_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence);
  nubars_decoh->Set_DecoherenceGammaEnergyScale(nubars_gamma_energy_scale*units.GeV);

      break;

    case kLIV:
      nus_LV = new nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau, nusquids::neutrino, false);
      nubars_LV = new nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau, nusquids::neutrino, false);
      nus_base = nus_decoh;

      break;

    case kNSI:
      nus_NSI = new nuSQUIDSNSI(eps_mutau, E_range, NuOscillator::kTau, nusquids::neutrino, true);
      nubars_NSI = new nuSQUIDSNSI(eps_mutau, E_range, NuOscillator::kTau, nusquids::neutrino, true);
      nus_base = nus_NSI;
      break;
  }

  nus_base->Set_h_max(integration_step*units.km );

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus_base->Set_h_max( integration_step*units.km );
  nubars_base->Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus_base->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars_base->Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus_base->Set_rel_error(nus_rel_error);
  nus_base->Set_abs_error(nus_abs_error);
  nubars_base->Set_rel_error(nubars_rel_error);
  nubars_base->Set_abs_error(nubars_abs_error);

  switch() {
    case kSM:

    break;

    case kDecoherence:
      nus_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::neutrino, true);
      nubars_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::antineutrino, false); // anti-neutrinos
      nus_base = nus_decoh;

  //Set the decoherence model and parameters for neutrinos
  nus_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV); // reference (default) value: nus_gamma_strength = 9.48e-18 
  nus_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence); // reference (default) value: nus_gamma_energy_dependence = 2
  nus_decoh->Set_DecoherenceGammaEnergyScale(nus_gamma_energy_scale*units.GeV); // reference (default) value: nus_gamma_energy_dependence = 1.0

  //Set the decoherence model and parameters for anti-neutrinos
  nubars_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV);
  nubars_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence);
  nubars_decoh->Set_DecoherenceGammaEnergyScale(nubars_gamma_energy_scale*units.GeV);

      break;

    case kLIV:

      break;

    case kNSI:



      break;
  }

//  nus_NSI = new nuSQUIDSNSI(eps_mutau, E_range, NuOscillator::kTau, nusquids::neutrino, true);

//  nus_LV = new nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau, nusquids::neutrino, false);

//  nus_pmns = new nusquids::nuSQUIDS(E_range, NuOscillator::kTau, nusquids::neutrino, false);

//  nus_base = nus_decoh;
/*
  nus_base->Set_h_max(integration_step*units.km );
  //  nus_base->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, nus_gamma_strength*units.eV);


  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus_decoh->Set_h_max( integration_step*units.km );
  nubars_decoh->Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus_decoh->Set_rel_error(nus_rel_error);
  nus_decoh->Set_abs_error(nus_abs_error);
  nubars_decoh->Set_rel_error(nubars_rel_error);
  nubars_decoh->Set_abs_error(nubars_abs_error);

  //Set the decoherence model and parameters for neutrinos
  nus_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV); // reference (default) value: nus_gamma_strength = 9.48e-18 
  nus_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence); // reference (default) value: nus_gamma_energy_dependence = 2
  nus_decoh->Set_DecoherenceGammaEnergyScale(nus_gamma_energy_scale*units.GeV); // reference (default) value: nus_gamma_energy_dependence = 1.0

  //Set the decoherence model and parameters for anti-neutrinos
  nubars_decoh->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, kGammaStrength*units.eV);
  nubars_decoh->Set_DecoherenceGammaEnergyDependence(kGammaEnergyDependence);
  nubars_decoh->Set_DecoherenceGammaEnergyScale(nubars_gamma_energy_scale*units.GeV);
*/
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  // Set mixing angles and masses for neutrinos
  nus_base->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nus_base->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nus_base->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nus_base->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nus_base->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nus_base->Set_CPPhase(0,2,OscParams[kDCP]);

  // Set mixing angles and masses for anti-neutrinos
  nubars_base->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nubars_base->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nubars_base->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nubars_base->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nubars_base->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nubars_base->Set_CPPhase(0,2,OscParams[kDCP]);

  const double layer_2 = OscParams[kPATHL]*units.km;
  std::shared_ptr<nusquids::ConstantDensity> constdens_env1 = std::make_shared<nusquids::ConstantDensity>(OscParams[kDENS],OscParams[kELECDENS]); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<nusquids::ConstantDensity::Track> track_env1 = std::make_shared<nusquids::ConstantDensity::Track>(layer_2);

  // Set energy density for neutrinos
  nus_base->Set_Body(constdens_env1);
  nus_base->Set_Track(track_env1);

  // Set energy density for anti-neutrinos
  nubars_base->Set_Body(constdens_env1);
  nubars_base->Set_Track(track_env1);

  // Construct the initial state
  // E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus_decoh->GetERange();
  // Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),static_cast<size_t>(NuOscillator::kTau)};

  // Index counter to have a handle on where neutrino oscillation probs are stored in array fWeightArray
  int index_counter = 0;

  // Loop over all neutrino flavors, set the initial state, propagate the neutrinos and store osc probs in fWeightarray in
  // order nu_e->nu_e,nu_e->nu_mu, nu_e->nu_tau,
  //       nu_mu->nu_e, nu_mu->nu_mu, nu_mu->nu_tau,
  //       nu_tau->nu_e, nu_tau->nu_mu, nu_tau->nu_tau,
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
      // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
      for ( int i = 0 ; i < inistate.extent(0); i++){
          for ( int k = 0; k < inistate.extent(1); k ++){
            inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
          }
      }

      //Set the initial state in nuSQuIDS object
      nus_base->Set_initial_state(inistate,nusquids::flavor);

    switch() {
      case kSM:
      nus_pmns->EvolveState();
//      SetupPropNuSQUIDSSM(fBSMModel);
//      CalcProbNuSQUIDSSM(fBSMModel);
      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_pmns->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
      break;

      case kDecoherence:
//      CalcProbNuSQUIDSDecoherence

      //Propagate the neutrinos in the earth for the path defined in path
      nus_decoh->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }

      break;

      case kLIV:

      nus_LIV->Set_ProgressBar(true);
      nus_LIV->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_LIV->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }

      break;

      case kNSI:

      nus_NSI->Set_ProgressBar(true);
      nus_NSI->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_NSI->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }

      break;
  }


      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
//      for(int fl=0; fl<NuOscillator::kTau; fl++){
//        for(int i = 0; i < fNEnergyPoints; i++) {
//          fWeightArray[index_counter] = nus_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
//          index_counter++;
//        }
//      }
  }

  // Now the same for anti-neutrinos:
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
      // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
      for ( int i = 0 ; i < inistate.extent(0); i++){
          for ( int k = 0; k < inistate.extent(1); k ++){
            inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
          }
      }

      //Set the initial state in nuSQuIDS object
      nubars_base->Set_initial_state(inistate,nusquids::flavor);

      switch() {
      case kSM:
//      SetupPropNuSQUIDSSM(fBSMModel);
//      CalcProbNuSQUIDSSM(fBSMModel);

      nubars_pmns->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_pmns->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
	
      break;

      case kDecoherence:
//      CalcProbNuSQUIDSDecoherence
      //Propagate the neutrinos in the earth for the path defined in path
      nubars_decoh->EvolveState();
      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
      break;

      case kLIV:

      nubars_LIV->EvolveState();
      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_LIV->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
      break;

      case kNSI:

      nubars_NSI->EvolveState();
      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_NSI->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
      break;
  }

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
//      for(int fl=0; fl<NuOscillator::kTau; fl++){
//        for(int i = 0; i < fNEnergyPoints; i++) {
//          fWeightArray[index_counter] = nubars_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
//          index_counter++;
//        }
//      }
  }
}

void OscProbCalcerNuSQUIDSLinear::CalcProbNuSQUIDSSM() {

}

int OscProbCalcerNuSQUIDSLinear::PMNS_StrToInt(std::string PMNSType) {
  return -1;
}

int OscProbCalcerNuSQUIDSLinear::BSMModel_StrToInt(std::string BSMModel) {
  if (BSMModel = "SM") {
    return kSM;
  }
  else if (BSMModel = "Decoherence") {
    return kDecoherence;
  }
  else if (BSMModel = "LIV") {
    return kLIV;
  }
  else if (BSMModel = "NSI") {
    return kNSI;
  }
  else {
    std::cerr << "Invalid BSM model provided:" << BSMModel << std::endl;
  }
  return -1;
}

int OscProbCalcerNuSQUIDSLinear::GetNOscParams(int OscType) {
  return kNOscParams;
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
