///
/// @file  ATLAS_1407_3935.cc
/// @brief Implementation of ATLAS_1407_3935 analysis
/// @author Kazuki
/// @date created 06/30/2016
/// @date last revision 06/30/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1407_3935 : public Analysis {
	public:

		ATLAS_1407_3935()
			: Analysis("ATLAS_1407_3935") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );

            // muons
            IsoMuon mu( Range(PT > 20.) & Range(abseta < 2.7) );
            mu.addIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);                        
            mu.addIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_PF_CMS" ) );
            addProjection(mu, "Muons");

            bookHisto1D("Mmm", 40, 80, 100, "Mmm", "Mmm", "Arbitrary");

			// End init section -- do not edit/remove this comment
		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

			const Particles& mus = applyProjection<IsoMuon>(event, "Muons").particlesByPt();

			if( mus.size() == 2 ){
				if( mus[0].pT() > 25 && mus[1].pT() > 20){
	                double minv = (mus[0].momentum() + mus[1].momentum()).mass();                
    	            fillPlot("Mmm", minv);
				}
			}

			// End analyze section -- do not edit/remove this comment
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			// Histogram normalization section -- do not edit/remove this comment
			/// @todo normalize the histograms
            scale("Mmm");

			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1407_3935)
}
