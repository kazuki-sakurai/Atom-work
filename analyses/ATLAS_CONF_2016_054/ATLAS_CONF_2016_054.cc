///
/// @file  ATLAS_CONF_2016_054.cc
/// @brief Implementation of ATLAS_CONF_2016_054 analysis
/// @author Seng Pei Liew
/// @date created 09/08/2016
/// @date last revision 09/08/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"
//#include "Rivet/Projections/Sphericity.hh"

using namespace std;

namespace Atom {

	class ATLAS_CONF_2016_054 : public Analysis {
	public:

		ATLAS_CONF_2016_054()
			: Analysis("ATLAS_CONF_2016_054") {
			setNeedsCrossSection(true);
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

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );

      FastJets jets_pre(fsbase, 
              hadRange & Range( PT > 20 & abseta < 4.5 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets_pre.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets_pre.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      HeavyFlavorJets bjets_pre(jets_pre, Range(PT > 30. & abseta < 2.5));
      bjets_pre.setTaggingEfficiency( *getBJetEff("BJet_Ident_MV1_ATLAS") );
      bjets_pre.setCurrentWorkingPoint( 0.77 );
      

      IsoMuon mu_base(Range(PT > 6. & abseta < 2.5));
      mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      IsoElectron ele_base( Range(PT > 7. & abseta < 2.5) );
      ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );

      IsoMuon mu(Range(PT > 6. & abseta < 2.5));
      mu.addIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.0, CALO_ALL,1.);
      mu.addIso(CALO_ISO_PT, 0.2,  0.1,  0.0, 0.1, CALO_ALL);
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      IsoElectron ele( Range(PT > 7. & abseta < 2.5) );
      ele.addIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.0, CALO_ALL,1.);
      ele.addIso(CALO_ISO_PT, 0.2,  0.1,  0.0, 0.1, CALO_ALL);
      ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
			
            // Overlap removal

      NearIsoParticle ele_pre(ele_base);
      ele_pre.addFilter(jets_pre,0.2);

      NearIsoParticle ele_sig(ele);
      ele_sig.addFilter(jets_pre,0.2);

      NearIsoParticle bjets(bjets_pre);
      bjets.addFilter(mu_base,0.2);
      addProjection(bjets, "BJets");

      NearIsoParticle mu_pre(mu_base);
      mu_pre.addFilter(jets_pre,0.2);
      mu_pre.addFilter(ele_base,0.01);

      MergedFinalState leps_pre(ele_base,mu_base);
      
      NearIsoParticle mu_sig(mu);
      mu_sig.addFilter(jets_pre,0.2);
      mu_sig.addFilter(ele_base,0.01);

      MergedFinalState leps_sig(ele_sig,mu_sig);
      addProjection(leps_sig,"Signal_Leptons");

      NearIsoParticle jets_sig(jets_pre,Range(PT > 30. & abseta < 2.8));
      addProjection(jets_sig,"Jets");

      MergedFinalState met_seed(jets_pre,leps_pre);
      MissingMomentum met( fsbase, met_seed );
      met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
      addProjection(met, "MissingEt");

      // implement aplanarity

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
//================================// 
			bookEfficiency("GG2J");
			bookEfficiency("GG6J-bulk");
			bookEfficiency("GG6J-high-mass");
			bookEfficiency("GG4J-high-x");
			bookEfficiency("GG4J-low-x");
			bookEfficiency("GG4J-low-x-b-veto");
			bookEfficiency("SS4J-x=1/2");
			bookEfficiency("SS5J-x=1/2");
			bookEfficiency("SS4J-low-x");
			bookEfficiency("SS5J-high-x");
     	//================================// 
			// End init section -- do not edit/remove this comment
                  bookHisto1D("GG2JMET", 6, 250, 670, "GG2JMET", "MET", "Arbitrary");
                  bookHisto1D("GG6J-bulkMT", 6, 65, 545, "GG6J-bulkMT", "MT", "Arbitrary");
                  bookHisto1D("GG6J-high-massaplan", 6, 0, 0.12, "GG6J-high-massaplan", "APLAN", "Arbitrary");
                  bookHisto1D("GG4J-high-xMT", 4, 175, 575, "GG4J-high-xMT", "MT", "Arbitrary");
                  bookHisto1D("GG4J-low-xMT", 4, 25, 425, "GG4J-low-xMT", "MT", "Arbitrary");
                  bookHisto1D("GG4J-low-x-b-vetoaplan", 6, 0, 0.09, "GG4J-low-x-b-vetoaplan", "APLAN", "Arbitrary");
                  bookHisto1D("SS4J-x=1/2aplan", 4, 0, 0.16, "SS4J-x=1/2aplan", "APLAN", "Arbitrary");
                  bookHisto1D("SS5J-x=1/2MET", 6, 200, 500, "SS5J-x=1/2MET", "MET", "Arbitrary");
                  bookHisto1D("SS4J-low-xmt", 7, 50, 400, "SS4J-low-xmt", "MT", "Arbitrary");
                  bookHisto1D("SS5J-high-xaplan", 3, 0, 0.09, "SS5J-high-xaplan", "APLAN", "Arbitrary");


		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

			// temporary
 	  const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
     // const Particles& large_jets = applyProjection<FastJets>(event, "Large_Jets").particlesByPt();      
      const Particles& leps = applyProjection<NearIsoParticle>(event, "Signal_Leptons").particlesByPt();
  //    const Particles& leps_cand = applyProjection<MergedFinalState>(event, "Candidate_Leptons").particlesByPt();            
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& untagged_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

      double MET = met.pT(); 
      int Njet = jets.size(); 
      int Nb = bjets.size();
      int Nlep = leps.size();
   //   int Nlep_cand = leps_cand.size();      
  //    int N_J = large_jets.size();
      // temporary

      		if (Nlep != 1) vetoEvent; 
      		
      		double mT = -100;
     		mT = get_mT(leps[0].momentum(), met);

     		double meff_inc = MET+leps[0].pT();
     		for(int i = 0; i < Njet; i++){
				meff_inc += jets[i].pT();
      } 
		//redifine aplanarity once bug is fixed
			double aplan_jet = 1;
			double aplan_lep = 1;

		//GG2J
      if(Njet >1 and leps[0].pT() <= 35.){
      	if(jets[0].pT()>200 and jets[1].pT()>30 and mT >100 and MET >460 and MET/meff_inc > 0.35)
      				pass("GG2J");

      }

      //GG6J
	if(Njet >5 and leps[0].pT() > 35. and jets[0].pT()>125 and jets[1].pT()>30 and jets[2].pT()>30 and jets[3].pT()>30 and jets[4].pT()>30 and jets[5].pT()>30 and mT>225 and MET > 250 and aplan_jet > 0.04){
      	
      		if(meff_inc > 1000 and MET/meff_inc > 0.2){
	pass("GG6J-bulk");
      		}
      		if(meff_inc > 2000 and MET/meff_inc > 0.1){
	pass("GG6J-high-mass");
      		}
      	
      }
      	
      //GG4j
      	if(Njet > 3 and jets[0].pT()> 100){
      		if(  jets[1].pT()> 100 and jets[2].pT()>100 and jets[3].pT()>100 and mT > 125 and MET > 250 and meff_inc > 2000){
      			//aplanarity
      			if(Nb==0 and aplan_jet > 0.03){
      				pass("GG4J-low-x-b-veto");
      			}
      		if( aplan_jet > 0.06){
      				pass("GG4J-low-x");
      			}
      		}
      		if (leps[0].pT() > 35 and jets[0].pT()> 400 and jets[1].pT() > 30 and jets[2].pT() > 30 and jets[2].pT() < 100 and mT > 475 and MET > 250 and meff_inc > 1600 and MET/meff_inc > 0.3){
      	pass("GG4J-high-x");
                  }
      	}
      				
      	//SS4j
      	if(Nb==0 and Njet > 3 and jets[0].pT()> 50 and leps[0].pT() > 35){
      		if(  jets[1].pT()> 50 and jets[2].pT()>50 and jets[3].pT()>50 and mT > 175 and MET > 300 and meff_inc > 1200 and aplan_lep > 0.08){
      			      		pass("SS4J-x=1/2");
      		}
      		if( jets[0].pT()> 250 and jets[1].pT()> 250 and jets[2].pT()> 30 and jets[3].pT()>30 and mT > 150 and mT < 400 and MET > 250  and aplan_lep > 0.03){
      			      		pass("SS4J-low-x");
      		}
      	}

      	//SS5j

      	if(Nb==0 and Njet > 4 and jets[0].pT()> 30 and jets[1].pT()> 30 and jets[2].pT()> 30 and jets[3].pT()> 30 and jets[4].pT()> 30 and leps[0].pT() > 35){
      		if(jets[0].pT()> 50 and jets[1].pT()> 50 and jets[2].pT()>50 and jets[3].pT()>50 and jets[4].pT()>50 and mT > 175 and MET > 300 and MET/meff_inc > 0.2 ){
      			      		pass("SS5J-x=1/2");
      		}
      		if(  mT > 400 and MET > 400  and aplan_lep > 0.03){
      			      		pass("SS5J-high-x");
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

			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_CONF_2016_054)
}
