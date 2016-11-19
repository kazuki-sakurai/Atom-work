///
/// @file  ATLAS_CONF_2016_093.cc
/// @brief Implementation of ATLAS_CONF_2016_093 analysis
/// @author Kazuki
/// @date created 11/14/2016
/// @date last revision 11/14/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

  class ATLAS_CONF_2016_093 : public Analysis {
  public:

    ATLAS_CONF_2016_093()
      : Analysis("ATLAS_CONF_2016_093") {
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

      FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

      IsoElectron ele( Range(PT > 10. & abseta < 2.47) );
      ele.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      ele.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
      ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );

      IsoMuon mu( Range(PT > 10. & abseta < 2.4) );
      mu.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      mu.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );
      FastJets jets(fsbase, 
              hadRange & Range( PT > 20 & abseta < 2.8 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      HeavyFlavorJets bjets(jets);
      bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
      bjets.setCurrentWorkingPoint( 0.77 );

      const double dRprong = 0.2;
      Range tau_range = Range(PT > 20) & (Range(ETA, -2.47, 2.47) - Range(ETA, -1.37, -1.52) - Range(ETA, 1.37, 1.52));
      ParamTauFinder taus(jets, 0.4, dRprong, tau_range);
      taus.setTaggingParams( getTauTag("Tau_Ident_BDT_Medium_ATLAS") );

      // Overlap removal
      NearIsoParticle taus_clean(taus);
      taus_clean.addFilter(ele, 0.2);
      taus_clean.addFilter(mu, 0.2);

      NearIsoParticle jets_clean1(jets);
      jets_clean1.addFilter(ele, 0.2);

      NearIsoParticle ele_clean(ele);
      ele_clean.addFilter(jets_clean1, 0.4);

      NearIsoParticle jets_clean(jets_clean1);
      jets_clean.addFilter(mu, 0.4);
      jets_clean.addFilter(taus_clean, 0.4);

      MergedFinalState leptons(ele_clean, mu);

      addProjection(leptons,     "Leptons");
      addProjection(jets_clean,  "Jets");
      addProjection(bjets,       "BJets");
      addProjection(taus_clean,  "TAUJets");

      MergedFinalState met_seed(jets_clean, leptons);
      MissingMomentum met( fsbase, met_seed );
      met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
      addProjection(met, "MissingEt");


      // Projection booking section -- do not edit/remove this comment
      /// @todo define projections (see examples and manual)
      // FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
      // addProjection(jets, "Jets");

      // Histogram booking section -- do not edit/remove this comment
      /// @todo book histograms
      // triplets of numbers correspond to HepData notation d??-x??-y??
      // bookHisto1D(1,1,1, "Meff");

      // Efficiency booking section -- do not edit/remove this comment
      /// @todo book the efficiencies
      // bookEfficiency("SR1");
      // bookEfficiency("SR2","Description of signal region 2 efficiency");
      // bookEfficiency("CR1","Control region 1", true);

      // Cuts booking section -- do not edit/remove this comment
      /// @todo book the cuts
      // bookCut("CutNJets");
      // bookCut("CutPTJ1","description goes here");
      // bookCut("CutEtaJet","this is a control region cut", true);

      // End init section -- do not edit/remove this comment
    }


    /// Perform the per-event analysis
    /// param[in]   event    the event to be analyzed
    void analyze(const Event& event) {

      const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      const Particles& taus = applyProjection<NearIsoParticle>(event, "TAUJets").particlesByPt();      
      const Particles& leps = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt();      
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& LF_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
      double MET = met.pT();

      cout << taus.size() <<"  "<< bjets.size() << endl;
      // Projection application section -- do not edit/remove this comment
      /// @todo apply projections
      // const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);

      // Analysis body section -- do not edit/remove this comment
      /// @todo apply cuts
      // if(!cut(jets.size(), CUT_GT, 2, "CutNJets")) {
      // 	vetoEvent;
      // }

      // if(!cut(jets[0].pT(), CUT_GT, 50., "CutPTJ1")) {
      // 	vetoEvent;
      // }

      // if(!cut(fabs(jets[0].eta()), CUT_OUT, make_pair(0.3,0.8), "CutEtaJet")) {
      // 	vetoEvent;
      // }

      /// @todo fill efficiencies
      // pass("SR1");

      /// @todo fill histograms 
      // fillPlot("Meff", meff);

      // End analyze section -- do not edit/remove this comment
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Histogram normalization section -- do not edit/remove this comment
      /// @todo normalize the histograms
      // scale("Mjj");

      // End finalize section -- do not edit/remove this comment
    }

    //@}

  };

  // This global object acts as a hook for the plugin system
  AtomPlugin(ATLAS_CONF_2016_093)
}
