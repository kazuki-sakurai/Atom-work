///
/// @file  ATLAS_1408_7084.cc
/// @brief Implementation of ATLAS_1408_7084 analysis
/// @author Seng Pei Liew
/// @date created 08/19/2016
/// @date last revision 08/19/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1408_7084 : public Analysis {
	public:

		ATLAS_1408_7084()
			: Analysis("ATLAS_1408_7084") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		double get_thrust(FourMomentum p1, FourMomentum p2){
			FourMomentum pplus = p1+p2;
			FourMomentum pminus = p1-p2;
			Vector3 pp = Vector3(pplus.px(),pplus.py(),0);
			Vector3 pm = Vector3(pminus.px(),pminus.py(),0);
			double cross =  abs(pp.x()*pm.y()-pp.y()*pm.x());
			double pmsqrt = sqrt(pm.x()*pm.x()+pm.y()*pm.y());
			return cross/pmsqrt;
		}

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );
			FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

	//photon

            IsoPhoton gam( Range(PT > 25.) & Range(abseta, 0.0, 2.37) - Range(abseta, 1.37, 1.56) );
          //  gam.addIso(CALO_ISO_ET,  0.4,  0.0, 6, 0.0, CALO_ALL);                        
          //  gam.addIso(TRACK_ISO_PT, 0.2,  0.0,  2.6, 0.0, CALO_ALL, 1.0);
            gam.addIso(CALO_ISO_ET,  0.4,  0.0, 6, 0.01, CALO_ALL);                        
            gam.addIso(TRACK_ISO_PT, 0.2,  0.0,  2.6, 0.01, CALO_ALL, 1.0);
            gam.setSmearingParams( getPhotonSim( "Photon_Smear_2013_ATLAS" ) );
            gam.setEfficiencyParams( getPhotonEff( "Photon_Ident_Tight_2011_ATLAS" ) );

        //    addProjection(Thrust(gam),"Thrust");
            addProjection(gam,"Photon");

            bookHisto1D("central_highpt", 67, 108, 141.5, "Mgg", "Mgg", "Arbitrary");
            bookHisto1D("forward_lowpt", 67, 108, 141.5, "Mgg", "Mgg", "Arbitrary");

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
			const Particles& photon = applyProjection<IsoPhoton>(event, "Photon").particlesByPt();
			// const Thrust& thrust = applyProjection<Thrust>(event, "Thrust");
			//if(photon.size()> 0) 				cout << "found photon" << "\n";      
			if(photon.size() == 2){
			//	cout << "found photon" << "\n";
			    double minv = (photon[0].momentum() + photon[1].momentum()).mass();   
			    if (photon[0].Et()> 0.35*minv && photon[1].Et()> 0.25*minv) {
			    	double ptt = get_thrust(photon[0].momentum(),photon[1].momentum());
			    	if (ptt > 70 && photon[0].momentum().abseta() < 0.95 && photon[1].momentum().abseta() < 0.95){
			    	fillPlot("central_highpt", minv);      
			    }
			    if (ptt <= 70 && photon[0].momentum().abseta() >= 0.95 && photon[1].momentum().abseta() >= 0.95){
			    	fillPlot("forward_lowpt", minv);      
			    }
			    }
			}
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
			scale("central_highpt",1);
			scale("forward_lowpt",1);
			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1408_7084)
}
