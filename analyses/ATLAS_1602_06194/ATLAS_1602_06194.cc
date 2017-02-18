///
/// @file  ATLAS_1602_06194.cc
/// @brief Implementation of ATLAS_1602_06194 analysis
/// @author Kazuki
/// @date created 06/18/2016
/// @date last revision 06/18/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1602_06194 : public Analysis {
	public:

		ATLAS_1602_06194()
			: Analysis("ATLAS_1602_06194") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2015" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            
            IsoElectron ele_base( Range(PT > 10. & abseta < 2.47) );
            ele_base.addIso(TRACK_ISO_PT, 0.3,  0.01,  0.0, 0.01, CALO_ALL);
            ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele_base.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );

            IsoMuon mu_base(Range(PT > 10. & abseta < 2.5));
            mu_base.addIso(TRACK_ISO_PT, 0.3,  0.01,  0.0, 0.01, CALO_ALL);
            mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range( PT > 20 & abseta < 4.5 ), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // Overlap removal
            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(ele_base, 0.2);

            NearIsoParticle signal_jets(jets_clean, Range(PT > 40. & abseta < 2.8));
            addProjection(signal_jets, "Jets");

            HeavyFlavorJets bjets(jets_clean, Range(PT > 40. & abseta < 2.5));
            bjets.setTaggingEfficiency( *getBJetEff("BJet_Ident_MV1_ATLAS") );
            bjets.setCurrentWorkingPoint( 0.70 );
            addProjection(bjets, "BJets");

            //addProjection(Sphericity(signal_jets), "Sphericity");
            //addProjection(Sphericity(fsbase), "Sphericity");

            MergedFinalState base_leptons(ele_base, mu_base);
            NearIsoParticle leptons_clean(base_leptons);
            leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(leptons_clean, "Leptons");

            MergedFinalState met_seed(jets_clean, base_leptons);
            MissingMomentum met( fsbase, met_seed );
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

            /// @todo book the efficiencies

            bookEfficiency("8j50");
            bookEfficiency("8j50-1b");
            bookEfficiency("8j50-2b");
            bookEfficiency("9j50");
            bookEfficiency("9j50-1b");
            bookEfficiency("9j50-2b");
            bookEfficiency("10j50");
            bookEfficiency("10j50-1b");
            bookEfficiency("10j50-2b");

            bookEfficiency("7j80");
            bookEfficiency("7j80-1b");
            bookEfficiency("7j80-2b");
            bookEfficiency("8j80");
            bookEfficiency("8j80-1b");
            bookEfficiency("8j80-2b");

			// Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

            bookCut("Lepton Veto");
            bookCut("Trigger: n50");
            bookCut("n50 >= 8: 8j50");
            bookCut("MET/sqrt(HT) > 4: 8j50");
            bookCut("nb >= 1: 8j50");
            bookCut("nb >= 2: 8j50");
            bookCut("n50 >= 9: 9j50");
            bookCut("MET/sqrt(HT) > 4: 9j50");
            bookCut("nb >= 1: 9j50");
            bookCut("nb >= 2: 9j50");
            bookCut("n50 >= 10: 10j50");
            bookCut("MET/sqrt(HT) > 4: 10j50");
            bookCut("nb >= 1: 10j50");
            bookCut("nb >= 2: 10j50");
            bookCut("Trigger: n80");
            bookCut("n80 >= 8: 7j80");
            bookCut("MET/sqrt(HT) > 4: 7j80");
            bookCut("nb >= 1: 7j80");
            bookCut("nb >= 2: 7j80");
            bookCut("n80 >= 8: 8j80");
            bookCut("MET/sqrt(HT) > 4: 8j80");
            bookCut("nb >= 1: 8j80");
            bookCut("nb >= 2: 8j80");

			// End init section -- do not edit/remove this comment
		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

            double MET = met.pT();
            int nb = bjets.size();
            int nlep = leps.size();

	        if( !cut(nlep, CUT_EQ, 0, "Lepton Veto") ) vetoEvent;

            int n50_trigger = 0; int n80_trigger = 0;
            int n50 = 0; int n80 = 0;
            double HT = 0.;
            for(int i=0; i<jets.size(); i++){
            	//cout << i <<"  "<< jets[i].pT() << endl;
            	if( jets[i].pT() > 50. && jets[i].abseta() < 2.0 ) n50++;
            	if( jets[i].pT() > 80. && jets[i].abseta() < 2.0 ) n80++;
            	if( jets[i].pT() > 45. && jets[i].abseta() < 2.4 ) n50_trigger++;
            	if( jets[i].pT() > 70. && jets[i].abseta() < 2.4 ) n80_trigger++;
            	if( jets[i].pT() > 40. && jets[i].abseta() < 2.8 ) HT += jets[i].pT();
            }

            //cout << jets.size() <<"  "<< n50_trigger << endl;

            //n50
            if( cut(n50_trigger, CUT_GE, 6, "Trigger: n50") ){
            	//8j50
	            if( cut(n50, CUT_GE, 8, "n50 >= 8: 8j50") ){
		            if( cut(MET/sqrt(HT), CUT_GT, 4., "MET/sqrt(HT) > 4: 8j50") ){
		            	pass("8j50");
			            if( cut(nb, CUT_GE, 1, "nb >= 1: 8j50") ) pass("8j50-1b");
			            if( cut(nb, CUT_GE, 2, "nb >= 2: 8j50") ) pass("8j50-2b");				            
		            }
            	}

            	//9j50
	            if( cut(n50, CUT_GE, 9, "n50 >= 9: 9j50") ){
		            if( cut(MET/sqrt(HT), CUT_GT, 4., "MET/sqrt(HT) > 4: 9j50") ){
		            	pass("9j50");
			            if( cut(nb, CUT_GE, 1, "nb >= 1: 9j50") ) pass("9j50-1b");
			            if( cut(nb, CUT_GE, 2, "nb >= 2: 9j50") ) pass("9j50-2b");				            
		            }
            	}

            	//10j50
	            if( cut(n50, CUT_GE, 10, "n50 >= 10: 10j50") ){
		            if( cut(MET/sqrt(HT), CUT_GT, 4., "MET/sqrt(HT) > 4: 10j50") ){
		            	pass("10j50");
			            if( cut(nb, CUT_GE, 1, "nb >= 1: 10j50") ) pass("10j50-1b");
			            if( cut(nb, CUT_GE, 2, "nb >= 2: 10j50") ) pass("10j50-2b");				            
		            }
            	}
            }  


            //n80
            if( cut(n80_trigger, CUT_GE, 5, "Trigger: n80") ){
            	//7j80
	            if( cut(n80, CUT_GE, 7, "n80 >= 8: 7j80") ){
		            if( cut(MET/sqrt(HT), CUT_GT, 4., "MET/sqrt(HT) > 4: 7j80") ){
		            	pass("7j80");
			            if( cut(nb, CUT_GE, 1, "nb >= 1: 7j80") ) pass("7j80-1b");
			            if( cut(nb, CUT_GE, 2, "nb >= 2: 7j80") ) pass("7j80-2b");				            
		            }
            	}

            	//8j80
	            if( cut(n80, CUT_GE, 8, "n80 >= 8: 8j80") ){
		            if( cut(MET/sqrt(HT), CUT_GT, 4., "MET/sqrt(HT) > 4: 8j80") ){
		            	pass("8j80");
			            if( cut(nb, CUT_GE, 1, "nb >= 1: 8j80") ) pass("8j80-1b");
			            if( cut(nb, CUT_GE, 2, "nb >= 2: 8j80") ) pass("8j80-2b");				            
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
	AtomPlugin(ATLAS_1602_06194)
}
