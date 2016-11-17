///
/// @file  ATLAS_CONF_2016_098.cc
/// @brief Implementation of ATLAS_CONF_2016_098 analysis
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

  class ATLAS_CONF_2016_098 : public Analysis {
  public:

    ATLAS_CONF_2016_098()
        : Analysis("ATLAS_CONF_2016_098") {
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

      FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

      IsoElectron ele( Range(PT > 25. & abseta < 2.47) );
      ele.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      ele.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
      ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

      IsoMuon mu( Range(PT > 25. & abseta < 2.5) );
      mu.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      mu.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );
      FastJets jets(fsbase, 
              hadRange & Range( PT > 30 & abseta < 2.5 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      // Overlap removal
      NearIsoParticle jets_clean(jets);
      jets_clean.addFilter(ele, 0.2);

      HeavyFlavorJets bjets(jets_clean);
      bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
      bjets.setCurrentWorkingPoint( 0.77 );

      NearIsoParticle ele_clean(ele);
      ele_clean.addFilter(jets_clean, 0.4);

      NearIsoParticle mu_clean(mu);
      mu_clean.addFilter( bjets, 0.2);

      NearIsoParticle signal_jets(jets_clean);
      signal_jets.addFilter(mu_clean, 0.2);

      MergedFinalState leptons(ele_clean, mu_clean);

      addProjection(leptons,     "Leptons");
      addProjection(signal_jets, "Jets");
      addProjection(bjets,       "BJets");

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
      bookEfficiency("SRZ");
      bookEfficiency("SRZ_ee");
      bookEfficiency("SRZ_mm");                  
      bookEfficiency("SR-low_12-61");
      bookEfficiency("SR-low_12-81");
      bookEfficiency("SR-low_12-101");
      bookEfficiency("SR-low_81-101");
      bookEfficiency("SR-low_101-201");
      bookEfficiency("SR-low_101-301");
      bookEfficiency("SR-low_201-401");
      bookEfficiency("SR-low_301-501");
      bookEfficiency("SR-low_gt501");
      bookEfficiency("SR-medium_12-61");
      bookEfficiency("SR-medium_12-81");
      bookEfficiency("SR-medium_12-101");
      bookEfficiency("SR-medium_81-101");
      bookEfficiency("SR-medium_101-201");
      bookEfficiency("SR-medium_101-301");
      bookEfficiency("SR-medium_201-401");
      bookEfficiency("SR-medium_gt401");
      bookEfficiency("SR-high_12-61");
      bookEfficiency("SR-high_12-81");
      bookEfficiency("SR-high_12-101");
      bookEfficiency("SR-high_81-101");
      bookEfficiency("SR-high_101-201");
      bookEfficiency("SR-high_201-401");
      bookEfficiency("SR-high_gt401");

      // Cuts booking section -- do not edit/remove this comment
      /// @todo book the cuts
      bookCut("SFOS lepton pair");
      bookCut("Nj >= 2");
      bookCut("dphimin12(j,MET) > 0.4");
      bookCut("pT(lep1) > 50: on-Z");
      bookCut("MET > 225: on-Z");
      bookCut("HT_inc > 600: on-Z");
      bookCut("81 < mll < 101: on-Z");
      bookCut("MET > 200: SR-low");
      bookCut("12 < mll < 61: SR-low");
      bookCut("12 < mll < 81: SR-low");
      bookCut("12 < mll < 101: SR-low");
      bookCut("81 < mll < 101: SR-low");
      bookCut("101 < mll < 201: SR-low");
      bookCut("101 < mll < 301: SR-low");
      bookCut("201 < mll < 401: SR-low");
      bookCut("301 < mll < 501: SR-low");
      bookCut("mll > 501: SR-low");
      bookCut("MET > 400: SR-medium");
      bookCut("12 < mll < 61: SR-medium");
      bookCut("12 < mll < 81: SR-medium");
      bookCut("12 < mll < 101: SR-medium");
      bookCut("81 < mll < 101: SR-medium");
      bookCut("101 < mll < 201: SR-medium");
      bookCut("101 < mll < 301: SR-medium");
      bookCut("201 < mll < 401: SR-medium");
      bookCut("mll > 401: SR-medium");
      bookCut("MET > 700: SR-high");
      bookCut("12 < mll < 61: SR-high");
      bookCut("12 < mll < 81: SR-high");
      bookCut("12 < mll < 101: SR-high");
      bookCut("81 < mll < 101: SR-high");
      bookCut("101 < mll < 201: SR-high");
      bookCut("201 < mll < 401: SR-high");
      bookCut("mll > 401: SR-high");


      // End init section -- do not edit/remove this comment
    }


    /// Perform the per-event analysis
    /// param[in]   event    the event to be analyzed
    void analyze(const Event& event) {

      const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      const Particles& leps = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt();      
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& untagged_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
      double MET = met.pT();

      if(leps.size() < 2) vetoEvent;

      const Particle lep1 = leps[0];
      const Particle lep2 = leps[1];

      double mll = (lep1.momentum() + lep2.momentum()).mass();

      bool is_SFOS = false;
      bool is_mm   = false;      
      bool is_ee   = false;

      if( lep1.pdgId() + lep2.pdgId() == 0 ) is_SFOS = true;
      if( abs(lep1.pdgId()) == 11 && abs(lep2.pdgId()) == 11 ) is_ee = true;
      if( abs(lep1.pdgId()) == 13 && abs(lep2.pdgId()) == 13 ) is_mm = true;

      double dphimin2j = 100;
      if(jets.size() > 1){
        for(int i=0; i<2; i++){
          double dphi = deltaPhi(jets[i], met);
          if(dphi < dphimin2j) dphimin2j = dphi;        
        }
      }

      double HT = 0;
      for(int i=0; i<jets.size(); i++){
        HT += jets[i].pT();
      }
      double HT_inc = HT + lep1.pT() + lep2.pT();

      // on-Z
      if( cut( is_SFOS, "SFOS lepton pair" ) ){
        if( cut( jets.size(), CUT_GE, 2, "Nj >= 2" ) ){
          if( cut( dphimin2j, CUT_GT, 0.4, "dphimin12(j,MET) > 0.4" ) ){
            // on-Z
            if( cut( lep1.pT(), CUT_GT, 50., "pT(lep1) > 50: on-Z" ) ){
              if( cut( MET, CUT_GT, 225., "MET > 225: on-Z" ) ){
                if( cut( HT_inc, CUT_GT, 600., "HT_inc > 600: on-Z" ) ){
                  if( cut( mll, CUT_IN, make_pair(81., 101.), "81 < mll < 101: on-Z" ) ){
                    pass("SRZ");
                    if( is_ee ) pass("SRZ_ee");
                    if( is_mm ) pass("SRZ_mm");                  
                  }
                }
              }
            }
            // edge
            if( cut( MET, CUT_GT, 200., "MET > 200: SR-low" ) ){
              if( cut( mll, CUT_IN, make_pair(12.,  61.),  "12 < mll < 61: SR-low" )   ) pass("SR-low_12-61");
              if( cut( mll, CUT_IN, make_pair(12.,  81.),  "12 < mll < 81: SR-low" )   ) pass("SR-low_12-81");
              if( cut( mll, CUT_IN, make_pair(12.,  101.), "12 < mll < 101: SR-low" )  ) pass("SR-low_12-101");
              if( cut( mll, CUT_IN, make_pair(81.,  101.), "81 < mll < 101: SR-low" )  ) pass("SR-low_81-101");
              if( cut( mll, CUT_IN, make_pair(101., 201.), "101 < mll < 201: SR-low" ) ) pass("SR-low_101-201");
              if( cut( mll, CUT_IN, make_pair(101., 301.), "101 < mll < 301: SR-low" ) ) pass("SR-low_101-301");
              if( cut( mll, CUT_IN, make_pair(201., 401.), "201 < mll < 401: SR-low" ) ) pass("SR-low_201-401");
              if( cut( mll, CUT_IN, make_pair(301., 501.), "301 < mll < 501: SR-low" ) ) pass("SR-low_301-501");
              if( cut( mll, CUT_GT, 501., "mll > 501: SR-low" ) ) pass("SR-low_gt501");
            }
            if( cut( MET, CUT_GT, 400., "MET > 400: SR-medium" ) ){
              if( cut( mll, CUT_IN, make_pair(12.,  61.),  "12 < mll < 61: SR-medium" )   ) pass("SR-medium_12-61");
              if( cut( mll, CUT_IN, make_pair(12.,  81.),  "12 < mll < 81: SR-medium" )   ) pass("SR-medium_12-81");
              if( cut( mll, CUT_IN, make_pair(12.,  101.), "12 < mll < 101: SR-medium" )  ) pass("SR-medium_12-101");
              if( cut( mll, CUT_IN, make_pair(81.,  101.), "81 < mll < 101: SR-medium" )  ) pass("SR-medium_81-101");
              if( cut( mll, CUT_IN, make_pair(101., 201.), "101 < mll < 201: SR-medium" ) ) pass("SR-medium_101-201");
              if( cut( mll, CUT_IN, make_pair(101., 301.), "101 < mll < 301: SR-medium" ) ) pass("SR-medium_101-301");
              if( cut( mll, CUT_IN, make_pair(201., 401.), "201 < mll < 401: SR-medium" ) ) pass("SR-medium_201-401");
              if( cut( mll, CUT_GT, 401., "mll > 401: SR-medium" ) ) pass("SR-medium_gt401");
            }
            if( cut( MET, CUT_GT, 700., "MET > 700: SR-high" ) ){
              if( cut( mll, CUT_IN, make_pair(12.,  61.),  "12 < mll < 61: SR-high" )   ) pass("SR-high_12-61");
              if( cut( mll, CUT_IN, make_pair(12.,  81.),  "12 < mll < 81: SR-high" )   ) pass("SR-high_12-81");
              if( cut( mll, CUT_IN, make_pair(12.,  101.), "12 < mll < 101: SR-high" )  ) pass("SR-high_12-101");
              if( cut( mll, CUT_IN, make_pair(81.,  101.), "81 < mll < 101: SR-high" )  ) pass("SR-high_81-101");
              if( cut( mll, CUT_IN, make_pair(101., 201.), "101 < mll < 201: SR-high" ) ) pass("SR-high_101-201");
              if( cut( mll, CUT_IN, make_pair(201., 401.), "201 < mll < 401: SR-high" ) ) pass("SR-high_201-401");
              if( cut( mll, CUT_GT, 401., "mll > 401: SR-high" ) ) pass("SR-high_gt401");
            }
          }
        }
      }

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
  AtomPlugin(ATLAS_CONF_2016_098)
}
