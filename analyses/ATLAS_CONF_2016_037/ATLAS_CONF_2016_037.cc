///
/// @file  ATLAS_CONF_2016_037.cc
/// @brief Implementation of ATLAS_CONF_2016_037 analysis
/// @author Kazuki
/// @date created 08/20/2016
/// @date last revision 08/20/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_CONF_2016_037 : public Analysis {
	public:

		ATLAS_CONF_2016_037()
			: Analysis("ATLAS_CONF_2016_037") {
			//setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

      useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

      FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );
      FastJets jets(fsbase, 
              hadRange & Range( PT > 20 & abseta < 2.8 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      IsoMuon mu_base(Range(PT > 10. & abseta < 2.5));
      mu_base.addConeIso(TRACK_ISO_PT, 0.3,  0.06,  0.0, 0.01, CALO_ALL);
      mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      IsoElectron ele_base( Range(PT > 10.) & Range(abseta, 0., 2.47) - Range(abseta, 1.52, 1.37) );      
      ele_base.addConeIso(TRACK_ISO_PT, 0.2,  0.06,  0.0, 0.01, CALO_ALL);
      ele_base.addConeIso(CALO_ISO_ET, 0.2,  0.06,  0.0, 0.01, CALO_ALL);      
      ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
      ele_base.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

      // Overlap removal
      NearIsoParticle jets_clean(jets);
      jets_clean.addFilter(ele_base, 0.2);
      addProjection(jets_clean, "Jets");

      HeavyFlavorJets bjets(jets_clean, Range(abseta < 2.5));
      bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
      bjets.setCurrentWorkingPoint( 0.7 );
      addProjection(bjets, "BJets");

      MergedFinalState leptons_base(ele_base, mu_base);
      NearIsoParticle leptons_clean(leptons_base);
      leptons_clean.addFilter(jets_clean, 0.4);
      addProjection(leptons_clean, "Leptons");

      MergedFinalState met_seed(jets_clean, leptons_clean);
      MissingMomentum met( fsbase, met_seed );
      met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
      addProjection(met, "MissingEt");

      //=================================//
      bookEfficiency("SR3L1");
      bookEfficiency("SR3L2");
      bookEfficiency("SR0b1");
      bookEfficiency("SR0b2");
      bookEfficiency("SR1b");
      bookEfficiency("SR1b-GG");
      bookEfficiency("SR3b");
      bookEfficiency("SR1b-DD");
      bookEfficiency("SR3b-DD");

      bookCut(">=3lep(pT>20,20,10): SR3L");
      bookCut("no b-jet(pT>20): SR3L");
      bookCut(">=4jets(pT>40): SR3L");
      bookCut("MET>150: SR3L1");
      bookCut("MET>200: SR3L2");
      bookCut("meff>1500: SR3L2");
      bookCut(">=2SSlep(pT>20): SR0b");
      bookCut("no b-jet(pT>20): SR0b");
      bookCut(">=6jets(pT>25): SR0b1");
      bookCut("MET>150: SR0b1");
      bookCut("meff>500: SR0b1");
      bookCut(">=6jets(pT>25): SR0b2");
      bookCut("MET>150: SR0b2");
      bookCut("meff>900: SR0b2");
      bookCut(">=2SSlep(pT>20): SR1b(-GG)");
      bookCut(">=1b-jet(pT>20): SR1b(-GG)");
      bookCut(">=6jets(pT>25): SR1b");
      bookCut("MET>200: SR1b");
      bookCut("meff>650: SR1b");
      bookCut(">=6jets(pT>50): SR1b-GG");
      bookCut("meff>1800: SR1b-GG");
      bookCut(">=2SSlep(pT>20): SR3b");
      bookCut(">=3b-jet(pT>20): SR3b");
      bookCut(">=6jets(pT>25): SR3b");
      bookCut("MET>150: SR3b");
      bookCut("meff>600: SR3b");
      bookCut(">=2SSlep(q<0,pT>20): SR-DD");
      bookCut(">=1b-jet(pT>20): SR1b-DD");
      bookCut(">=4jets(pT>50): SR1b-DD");
      bookCut("meff>1200: SR1b-DD");
      bookCut(">=3b-jet(pT>20): SR3b-DD");
      bookCut(">=4jets(pT>50): SR3b-DD");
      bookCut("meff>1000: SR3b-DD");

			// End init section -- do not edit/remove this comment
		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

      const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      //const Particles& large_jets = applyProjection<FastJets>(event, "Large_Jets").particlesByPt();      
      const Particles& leps = applyProjection<NearIsoParticle>(event, "Signal_Leptons").particlesByPt();
      const Particles& leps_cand = applyProjection<MergedFinalState>(event, "Candidate_Leptons").particlesByPt();            
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& untagged_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

      double MET = met.pT();

      double meff = MET;
      for(int i=0; i<jets.size(); i++) meff += jets[i].pT();
      for(int i=0; i<leps.size(); i++) meff += leps[i].pT();

      int nj25 = 0; int nj40 = 0; int nj50 = 0;
      for(int i=0; i<jets.size(); i++){
        if( jets[i].pT() > 25. ) nj25++;
        if( jets[i].pT() > 40. ) nj40++;
        if( jets[i].pT() > 50. ) nj50++;        
      }

      bool three_lepton = false;
      if( leps.size() >= 3 && leps[1].pT() > 20. && leps[2].pT() > 10. ) three_lepton = true;

      //SR3L
      if( cut(three_lepton, ">=3lep(pT>20,20,10): SR3L") ){
        if( cut(bjets.size(), CUT_EQ, 0, "no b-jet(pT>20): SR3L") ){
          if( cut(nj40, CUT_GE, 4, ">=4jets(pT>40): SR3L") ){
            // SR3L1
            if( cut(MET, CUT_GT, 150., "MET>150: SR3L1") ){
              pass("SR3L1");
            }
            // SR3L2
            if( cut(MET, CUT_GT, 200., "MET>200: SR3L2") ){
              if( cut(meff, CUT_GT, 1500., "meff>1500: SR3L2") ){
                pass("SR3L2");
              }
            }
          }
        }
      }  

      bool SS_lepton = false;
      for(int i=0; i<leps.size()-1; i++){
        for(int j=i+1; j<leps.size(); j++){
          if( leps[j].pT() > 20. && leps[i].pdgId()*leps[j].pdgId() > 0 ) SS_lepton = true; 
        }
      }

      //SR0b
      if( cut(SS_lepton, ">=2SSlep(pT>20): SR0b") ){
        if( cut(bjets.size(), CUT_EQ, 0, "no b-jet(pT>20): SR0b") ){
          //SR0b1
          if( cut(nj25, CUT_GE, 6, ">=6jets(pT>25): SR0b1") ){
            if( cut(MET, CUT_GT, 150., "MET>150: SR0b1") ){
              if( cut(meff, CUT_GT, 500., "meff>500: SR0b1") ){
                pass("SR0b1");
              }
            }
          }
          //
          if( cut(nj40, CUT_GE, 6, ">=6jets(pT>25): SR0b2") ){
            if( cut(MET, CUT_GT, 150., "MET>150: SR0b2") ){
              if( cut(meff, CUT_GT, 900., "meff>900: SR0b2") ){
                pass("SR0b2");
              }
            }
          }
        }
      }            

      //SR1b
      if( cut(SS_lepton, ">=2SSlep(pT>20): SR1b(-GG)") ){
        if( cut(bjets.size(), CUT_GE, 1, ">=1b-jet(pT>20): SR1b(-GG)") ){
          //SR1b
          if( cut(nj25, CUT_GE, 6, ">=6jets(pT>25): SR1b") ){
            if( cut(MET, CUT_GT, 200., "MET>200: SR1b") ){
              if( cut(meff, CUT_GT, 650., "meff>650: SR1b") ){
                pass("SR1b");
              }
            }
          }
          //SR1b-GG
          if( cut(nj50, CUT_GE, 6, ">=6jets(pT>50): SR1b-GG") ){
            if( cut(meff, CUT_GT, 1800., "meff>1800: SR1b-GG") ){
              pass("SR1b-GG");
            }
          }          
        }
      }
      //SR3b
      if( cut(SS_lepton, ">=2SSlep(pT>20): SR3b") ){
        if( cut(bjets.size(), CUT_GE, 3, ">=3b-jet(pT>20): SR3b") ){
          if( cut(nj25, CUT_GE, 6, ">=6jets(pT>25): SR3b") ){
            if( cut(MET, CUT_GT, 150., "MET>150: SR3b") ){
              if( cut(meff, CUT_GT, 600., "meff>600: SR3b") ){
                pass("SR3b");
              }
            }
          }
        }
      }

      //SR-DD
      bool DD_lepton = false;
      for(int i=0; i<leps.size()-1; i++){
        for(int j=i+1; j<leps.size(); j++){
          if( leps[j].pT() > 20. && leps[i].pdgId() > 0 && leps[j].pdgId() > 0 ) DD_lepton = true; 
        }
      }

      if( cut(DD_lepton, ">=2SSlep(q<0,pT>20): SR-DD") ){
        //1b
        if( cut(bjets.size(), CUT_GE, 1, ">=1b-jet(pT>20): SR1b-DD") ){
          if( cut(nj50, CUT_GE, 4, ">=4jets(pT>50): SR1b-DD") ){
            if( cut(meff, CUT_GT, 1200., "meff>1200: SR1b-DD") ){
              pass("SR1b-DD");
            }
          }
        }
        //3b
        if( cut(bjets.size(), CUT_GE, 3, ">=3b-jet(pT>20): SR3b-DD") ){
          if( cut(nj50, CUT_GE, 4, ">=4jets(pT>50): SR3b-DD") ){
            if( cut(meff, CUT_GT, 1000., "meff>1000: SR3b-DD") ){
              pass("SR3b-DD");
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
	AtomPlugin(ATLAS_CONF_2016_037)
}
