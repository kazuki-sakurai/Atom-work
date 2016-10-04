///
/// @file  ATLAS_1407_5063.cc
/// @brief Implementation of ATLAS_1407_5063 analysis
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

	class ATLAS_1407_5063 : public Analysis {
	public:

		ATLAS_1407_5063()
			: Analysis("ATLAS_1407_5063") {
			//setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );
            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            // electrons
            IsoElectron ele0( Range(PT > 27.) & Range(abseta < 2.47) );
            ele0.addConeIso(CALO_ISO_PT,  0.3,  0.18,  0.0, 0.01, CALO_ALL);                        
            ele0.addConeIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL, 0.4);
            ele0.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele0.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            NearIsoParticle ele1(ele0, Range(abseta, 1.63, 1.74));
            NearIsoParticle ele2(ele0, Range(abseta, 2.3, 2.4));
            NearIsoParticle ele3(ele0, Range(abseta, 1.0, 2.4));

            addProjection(ele0, "Electron_inc");            
            addProjection(ele1, "Electron1");
            addProjection(ele2, "Electron2");
            addProjection(ele3, "Electron3");

            bookHisto1D("Mee", 40, 80, 100, "Mee", "Mee", "Arbitrary");
            bookHisto1D("Mee_eta", 40, 80, 100, "Mee", "Mee", "Arbitrary");
            bookHisto1D("Mee_eta_inc", 40, 80, 100, "Mee", "Mee", "Arbitrary");

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& ele_inc = applyProjection<IsoElectron>(event, "Electron_inc").particlesByPt();
            const Particles& ele1 = applyProjection<NearIsoParticle>(event, "Electron1").particlesByPt();
            const Particles& ele2 = applyProjection<NearIsoParticle>(event, "Electron2").particlesByPt();
            const Particles& ele3 = applyProjection<NearIsoParticle>(event, "Electron3").particlesByPt();

            if( ele_inc.size() == 2 ){
                double minv = (ele_inc[0].momentum() + ele_inc[1].momentum()).mass();                
                fillPlot("Mee", minv);
                //cout << minv << endl;
            }            

            if( ele1.size() == 1 && ele2.size() == 1){
                double minv = (ele1[0].momentum() + ele2[0].momentum()).mass();                
                fillPlot("Mee_eta", minv);
                //cout << minv << endl;
            }            

            if( ele3.size() == 2 ){
                double minv = (ele3[0].momentum() + ele3[1].momentum()).mass();                
                fillPlot("Mee_eta_inc", minv);
                //cout << minv << endl;
            }            


			// End analyze section -- do not edit/remove this comment
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			// Histogram normalization section -- do not edit/remove this comment
			/// @todo normalize the histograms
            scaleToLumi("Mee");
            scaleToLumi("Mee_eta");
            scaleToLumi("Mee_eta_inc");

			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1407_5063)
}
