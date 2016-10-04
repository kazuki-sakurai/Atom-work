///
/// @file  CMS_1402_0923.cc
/// @brief Implementation of CMS_1402_0923 analysis
/// @author Kazuki
/// @date created 06/25/2016
/// @date last revision 06/25/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class CMS_1402_0923 : public Analysis {
	public:

		CMS_1402_0923()
			: Analysis("CMS_1402_0923") {
			//setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );
            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            // electrons
            IsoElectron ele( Range(PT > 25.) & Range(abseta, 0.0, 2.5) - Range(abseta, 1.44, 1.57) );
            ele.addConeIso(CALO_ISO_PT,  0.4,  0.15,  0.0, 0.01, CALO_ALL);                        
            ele.addConeIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_Ceffective_CMS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_mediumWP_CMS" ) );
            addProjection(ele, "Electrons");

            // muons
            IsoMuon mu( Range(PT > 25.) & Range(abseta < 2.1) );
            mu.addConeIso(CALO_ISO_PT,  0.4,  0.15,  0.0, 0.01, CALO_ALL);                        
            mu.addConeIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_SIDRA_CMS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_PF_CMS" ) );
            addProjection(mu, "Muons");

            bookHisto1D("Mee", 60, 60, 120, "Mee", "Mee", "Arbitrary");
            bookHisto1D("Mmm", 60, 60, 120, "Mmm", "Mmm", "Arbitrary");

			// End init section -- do not edit/remove this comment
		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& eles = applyProjection<IsoElectron>(event, "Electrons").particlesByPt();
            const Particles& mus  = applyProjection<IsoMuon>(event, "Muons").particlesByPt();

            if( eles.size() == 2){
                double minv = (eles[0].momentum() + eles[1].momentum()).mass();                
                fillPlot("Mee", minv);
                //cout << minv << endl;
            }            

            if( mus.size() == 2){
                double minv = (mus[0].momentum() + mus[1].momentum()).mass();                
                fillPlot("Mmm", minv);
            }            

            // fillPlot("Meff", meff);


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
            scaleToLumi("Mee");
            scaleToLumi("Mmm");

			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(CMS_1402_0923)
}
