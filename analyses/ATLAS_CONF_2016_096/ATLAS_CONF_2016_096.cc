///
/// @file  ATLAS_CONF_2016_096.cc
/// @brief Implementation of ATLAS_CONF_2016_096 analysis
/// @author Kazuki
/// @date created 11/14/2016
/// @date last revision 11/14/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"
#include "Rivet/Tools/RivetMT2.hh"

using namespace std;

namespace Atom {

  class ATLAS_CONF_2016_096 : public Analysis {
  public:

    ATLAS_CONF_2016_096()
        : Analysis("ATLAS_CONF_2016_096") {
    }

    /// @name Analysis methods
    //@{

    double get_mT(FourMomentum p1, FourMomentum p2){
        double mTsq = 2. * ( p1.pT()*p2.pT() - p1.px()*p2.px() - p1.py()*p2.py() );
        double mT = sqrt(mTsq);
        return mT;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

      FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

      IsoElectron ele( Range(PT > 10. & abseta < 2.47) );
      ele.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      ele.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
      ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

      IsoMuon mu( Range(PT > 10. & abseta < 2.5) );
      mu.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      mu.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );
      FastJets jets(fsbase, 
              hadRange & Range( PT > 20 & abseta < 4.5 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      // Overlap removal
      NearIsoParticle jets_clean(jets);
      jets_clean.addFilter(ele, 0.2);

      HeavyFlavorJets bjets(jets_clean);
      bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
      bjets.setCurrentWorkingPoint( 0.77 );

      MergedFinalState lep_base(ele, mu);

      NearIsoParticle leptons(lep_base);
      leptons.addFilter(jets_clean, 0.4);

      addProjection(leptons,     "Leptons");
      addProjection(jets_clean, "Jets");
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
      bookEfficiency("SR2l-A");
      bookEfficiency("SR2l-B");
      bookEfficiency("SR2l-C");
      bookEfficiency("SR3l-I");
      bookEfficiency("SR3l-H");

      // Cuts booking section -- do not edit/remove this comment
      /// @todo book the cuts
      bookCut("Nlep == 2: SR2l");
      bookCut("OS lepton: SR2l");
      bookCut("pT(lep1) > 25: SR2l");
      bookCut("pT(lep2) > 20: SR2l");
      bookCut("CLjet veto 20: SR2l SF");
      bookCut("bjet veto: SR2l SF");
      bookCut("Fjet veto: SR2l SF");
      bookCut("Z veto: SR2l SF");
      bookCut("mT2 > 90: SR2l SF");
      bookCut("mT2 > 120: SR2l SF");
      bookCut("mT2 > 150: SR2l SF");
      bookCut("CLjet veto 30: SR2l DF");
      bookCut("bjet veto: SR2l DF");
      bookCut("Fjet veto: SR2l DF");
      bookCut("mT2 > 90: SR2l DF");
      bookCut("mT2 > 120: SR2l DF");
      bookCut("mT2 > 150: SR2l DF");
      bookCut("Nlep == 3: SR3l");
      bookCut("SFOS lepton: SR3l");
      bookCut("bjet veto: SR3l");
      bookCut("mT > 110: SR3l");
      bookCut("mSFOS not in (81, 101): SR3l-I");
      bookCut("pT(lep3) > 30: SR3l-I");
      bookCut("MET > 120: SR3l-I");
      bookCut("mSFOS > 101: SR3l-H");
      bookCut("pT(lep3) > 80: SR3l-H");
      bookCut("MET > 60: SR3l-H");

      // End init section -- do not edit/remove this comment
    }


    /// Perform the per-event analysis
    /// param[in]   event    the event to be analyzed
    void analyze(const Event& event) {

      const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();      
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& Ljets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
      double MET = met.pT();
      double mZ = 91.2;

      int Nlep = leps.size();
      bool is_OS = false;
      bool is_SF = false;
      bool is_DF = false;
      if( Nlep == 2 ){
        if(leps[0].pdgId()*leps[0].pdgId() < 0) is_OS = true;
        if( abs(leps[0].pdgId()) == abs(leps[1].pdgId()) ){
          is_SF = true;
        }else{
          is_DF = true;
        }
      }

      int N_CLjets_20 = 0;
      int N_CLjets_30 = 0;      
      for(int i=0; i<Ljets.size(); i++){
        if(Ljets[i].pT() > 20 && Ljets[i].abseta() < 2.4) N_CLjets_20++;        
        if(Ljets[i].pT() > 30 && Ljets[i].abseta() < 2.4) N_CLjets_30++;                
      }

      int N_bjets = 0;
      for(int i=0; i<bjets.size(); i++){
        if(bjets[i].pT() > 20 && bjets[i].abseta() < 2.4) N_bjets++;        
      }

      int N_Fjets = 0;
      for(int i=0; i<jets.size(); i++){
        if(jets[i].pT() > 30 && jets[i].abseta() > 2.4) N_Fjets++;        
      }

      bool SR2l_A_pass = false;
      bool SR2l_B_pass = false;
      bool SR2l_C_pass = false;
      // Two-lepton
      if( cut( Nlep, CUT_EQ, 2, "Nlep == 2: SR2l" ) ){
        if( cut(is_OS, "OS lepton: SR2l") ){
          if( cut(leps[0].pT(), CUT_GT, 25, "pT(lep1) > 25: SR2l" ) ){
            if( cut(leps[1].pT(), CUT_GT, 20, "pT(lep2) > 20: SR2l" ) ){

              double minv = 0;
              double mt2 = Rivet::mT2::mT2(leps[0].momentum(), leps[1].momentum(), met, minv);

              if( is_SF ){
                if( cut(N_CLjets_20, CUT_EQ, 0, "CLjet veto 20: SR2l SF" ) ){
                  if( cut(N_bjets, CUT_EQ, 0, "bjet veto: SR2l SF" ) ){
                    if( cut(N_Fjets, CUT_EQ, 0, "Fjet veto: SR2l SF" ) ){
                      double mll = (leps[0].momentum() + leps[1].momentum()).mass();
                      if( cut( fabs(mll - mZ), CUT_GT, 10, "Z veto: SR2l SF" )){
                        if( cut( mt2, CUT_GT,  90, "mT2 > 90: SR2l SF" ))  SR2l_A_pass = true;
                        if( cut( mt2, CUT_GT, 120, "mT2 > 120: SR2l SF" )) SR2l_B_pass = true;
                        if( cut( mt2, CUT_GT, 150, "mT2 > 150: SR2l SF" )) SR2l_C_pass = true;
                      }
                    }
                  }
                }
              }

              if( is_DF ){
                if( cut(N_CLjets_30, CUT_EQ, 0, "CLjet veto 30: SR2l DF" ) ){
                  if( cut(N_bjets, CUT_EQ, 0, "bjet veto: SR2l DF" ) ){
                    if( cut(N_Fjets, CUT_EQ, 0, "Fjet veto: SR2l DF" ) ){
                      if( cut( mt2, CUT_GT,  90, "mT2 > 90: SR2l DF" ))  SR2l_A_pass = true;
                      if( cut( mt2, CUT_GT, 120, "mT2 > 120: SR2l DF" )) SR2l_B_pass = true;
                      if( cut( mt2, CUT_GT, 150, "mT2 > 150: SR2l DF" )) SR2l_C_pass = true;
                    }
                  }
                }
              }

            }
          }
        }
      }

      if(SR2l_A_pass) pass("SR2l-A");
      if(SR2l_B_pass) pass("SR2l-B");
      if(SR2l_C_pass) pass("SR2l-C");

      // Three-lepton
      if( cut( Nlep, CUT_EQ, 3, "Nlep == 3: SR3l" ) ){

        double mT = get_mT(leps[0].momentum(), met);

        bool is_SFOS = false;
        double mdm = 10000.;
        double mSFOS = -100.;
        for(int i=0; i<leps.size()-1; i++){
          for(int j=i+1; j<leps.size(); j++){
            if( leps[i].pdgId() + leps[j].pdgId() == 0 ){
              is_SFOS = true;
              double mll = (leps[i].momentum() + leps[j].momentum()).mass();
              if( fabs(mll - mZ) < mdm ){
                mdm = fabs(mll - mZ);
                mSFOS = mll;
              }
            }
          }
        }

        if( cut(is_SFOS, "SFOS lepton: SR3l" ) ){
          if( cut(N_bjets, CUT_EQ, 0, "bjet veto: SR3l" ) ){
            if( cut(mT, CUT_GT, 110, "mT > 110: SR3l" ) ){ 

              // SR3l-I
              if( cut(mSFOS, CUT_OUT, make_pair(81.2, 101.2), "mSFOS not in (81, 101): SR3l-I" )) {
                if( cut(leps[2].pT(), CUT_GT, 30, "pT(lep3) > 30: SR3l-I" ) ){
                  if( cut(MET, CUT_GT, 120., "MET > 120: SR3l-I" ) ){
                    pass("SR3l-I");
                  }
                }
              }
              // SR3l-H
              if( cut(mSFOS, CUT_GT, 101.2, "mSFOS > 101: SR3l-H" )) {
                if( cut(leps[2].pT(), CUT_GT, 80, "pT(lep3) > 80: SR3l-H" ) ){
                  if( cut(MET, CUT_GT, 60., "MET > 60: SR3l-H" ) ){
                    pass("SR3l-H");
                  }
                }
              }

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
  AtomPlugin(ATLAS_CONF_2016_096)
}
