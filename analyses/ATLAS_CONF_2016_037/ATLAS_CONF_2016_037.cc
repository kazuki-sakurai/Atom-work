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
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

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
	AtomPlugin(ATLAS_CONF_2016_037)
}
