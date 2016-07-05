///
/// @file  ATLAS_1603_09222.cc
/// @brief Implementation of ATLAS_1603_09222 analysis
/// @author Seng Pei Liew
/// @date created 06/30/2016
/// @date last revision 06/30/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1603_09222 : public Analysis {
	public:

		ATLAS_1603_09222()
			: Analysis("ATLAS_1603_09222") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );
            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            // electrons
            IsoElectron ele( Range(PT > 25.) & Range(abseta, 0.0, 2.47) - Range(abseta, 1.37, 1.52) );
            ele.addIso(CALO_ISO_PT,  0.3,  0.18,  0.0, 0.01, CALO_ALL);                        
            ele.addIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );
            addProjection(ele, "Electrons");

            // muons
            IsoMuon mu( Range(PT > 25.) & Range(abseta < 2.4) );
            mu.addIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);                        
            mu.addIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );
            addProjection(mu, "Muons");

            bookHisto1D("Mee", 40, 70, 110, "Mee", "Mee", "Arbitrary");
            bookHisto1D("Mmm", 40, 70, 110, "Mmm", "Mmm", "Arbitrary");
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

			// Projection application section -- do not edit/remove this comment
            /// @todo apply projections
			// const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);
			const Particles& mu = applyProjection<IsoMuon>(event, "Muons").particlesByPt();
			const Particles& ele = applyProjection<IsoMuon>(event, "Electrons").particlesByPt();

			if( mu.size() == 2 ){
				
	                double minv = (mu[0].momentum() + mu[1].momentum()).mass();                
    	            fillPlot("Mmm", minv);
				
			}

			if( ele.size() == 2 ){
				
	                double minv = (ele[0].momentum() + ele[1].momentum()).mass();                
    	            fillPlot("Mee", minv);
				
			}
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
			scale("Mmm");
			scale("Mee");
			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1603_09222)
}
