///
/// @file  ATLAS_CONF_2016_052.cc
/// @brief Implementation of ATLAS_CONF_2016_052 analysis
/// @author Kazuki
/// @date created 08/10/2016
/// @date last revision 08/10/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_CONF_2016_052 : public Analysis {
	public:

		ATLAS_CONF_2016_052()
			: Analysis("ATLAS_CONF_2016_052") {
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
      FastJets jets(fsbase, 
              hadRange & Range( PT > 20 & abseta < 2.8 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      HeavyFlavorJets bjets(jets, Range(PT > 20. & abseta < 2.5));
      bjets.setTaggingEfficiency( *getBJetEff("BJet_Ident_MV1_ATLAS") );
      bjets.setCurrentWorkingPoint( 0.77 );
      addProjection(bjets, "BJets");

      IsoMuon mu(Range(PT > 20. & abseta < 2.5));
      mu.addIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );
      
      IsoElectron ele_base( Range(PT > 20. & abseta < 2.47) );
      ele_base.addIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
      ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
			NearIsoParticle ele_cand(ele_base);
      ele_cand.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );
			NearIsoParticle ele_sig(ele_base);
      ele_sig.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

      // Overlap removal
      NearIsoParticle jets_clean(jets);
      jets_clean.addFilter(ele_base, 0.2);

      MergedFinalState leptons_sig(ele_sig, mu);
      NearIsoParticle leptons_clean(leptons_sig);
      leptons_clean.addFilter(jets_clean, 0.4);
      addProjection(leptons_clean, "Signal_Leptons");

      MergedFinalState leptons_cand(ele_cand, mu);
      addProjection(leptons_cand, "Candidate_Leptons");

      FastJets large_jets(jets_clean, 
              hadRange & Range( PT > 100 & abseta < 2.0 ), 
              muDetRange, FastJets::ANTIKT, 0.8 );
      addProjection(large_jets, "Large_Jets");

			NearIsoParticle jets_sig(jets_clean, Range( PT > 30 )) ;
      addProjection(jets_sig, "Jets");

      MergedFinalState met_seed(jets_clean, leptons_cand);
      MissingMomentum met( fsbase, met_seed );
      met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
      addProjection(met, "MissingEt");

     	//================================// 
			bookEfficiency("SR-Gbb-A");
			bookEfficiency("SR-Gbb-B");
			bookEfficiency("SR-Gtt-0L-A");
			bookEfficiency("SR-Gtt-0L-A");
			bookEfficiency("SR-Gtt-1L-A");
			bookEfficiency("SR-Gtt-1L-B");
			bookEfficiency("SR-Gtt-1L-C");
     	//================================// 

      bookCut("Njet >= 4: SR-Gbb(base)");
      bookCut("Lepton Veto: SR-Gbb(base)");
      bookCut("dphimin4j > 0.4: SR-Gbb(base)");
      bookCut("pTj(1,4) > 70: SR-Gbb-A");
      bookCut("Nb >= 3: SR-Gbb-A");
      bookCut("MET > 450: SR-Gbb-A");
      bookCut("meff_4j > 1900: SR-Gbb-A");
      bookCut("pTj(1,4) > 30: SR-Gbb-B");
      bookCut("Nb >= 4: SR-Gbb-B");
      bookCut("MET > 300: SR-Gbb-B");
      bookCut("meff_4j > 1000: SR-Gbb-B");
      bookCut("Njet >= 8: SR-Gtt-0L(base)");
      bookCut("Nb >= 3: SR-Gtt-0L(base)");
      bookCut("Lepton Veto: SR-Gtt-0L(base)");
      bookCut("dphimin4j > 0.4: SR-Gtt-0L(base)");
      bookCut("mT_minb > 80: SR-Gtt-0L(base)");
      bookCut("MET > 400: SR-Gtt-0L(base)");
      bookCut("meff_inc > 2000: SR-Gtt-0L-A");
      bookCut("M_J > 200: SR-Gtt-0L-A");
      bookCut("meff_inc > 1250: SR-Gtt-0L-A");
      bookCut("Nlep(sig) >= 1: SR-Gtt-1L(base)");
      bookCut("Njet >= 6: SR-Gtt-1L(base)");
      bookCut("Nb >= 3: SR-Gtt-1L(A-B)");
      bookCut("mT > 200: SR-Gtt-1L(A-B)");
      bookCut("mT_minb > 120: SR-Gtt-1L(A-B)");
      bookCut("MET > 200: SR-Gtt-1L-A");
      bookCut("meff_inc > 2000: SR-Gtt-1L-A");
      bookCut("M_J > 200: SR-Gtt-1L-A");
      bookCut("MET > 350: SR-Gtt-1L-B");
      bookCut("meff_inc > 1500: SR-Gtt-1L-B");
      bookCut("M_J > 150: SR-Gtt-1L-B");
      bookCut("Nb >= 4: SR-Gtt-1L-C");
      bookCut("mT > 150: SR-Gtt-1L-C");
      bookCut("mT_minb > 80: SR-Gtt-1L-C");
      bookCut("MET > 200: SR-Gtt-1L-C");
      bookCut("meff_inc > 500: SR-Gtt-1L-C");

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

      const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      const Particles& large_jets = applyProjection<FastJets>(event, "Large_Jets").particlesByPt();      
      const Particles& leps = applyProjection<NearIsoParticle>(event, "Signal_Leptons").particlesByPt();
      const Particles& leps_cand = applyProjection<MergedFinalState>(event, "Candidate_Leptons").particlesByPt();            
      const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& untagged_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

      double MET = met.pT(); 
      int Njet = jets.size(); 
      int Nb = bjets.size();
      int Nlep = leps.size();
      int Nlep_cand = leps_cand.size();      
      int N_J = large_jets.size();

      double meff_inc = MET; 
      double meff_4j = MET;
      for(int i = 0; i < Njet; i++){
				meff_inc += jets[i].pT();
				if(i < 4) meff_4j += jets[i].pT();
      } 
      for(int i=0; i<Nlep; i++){
      	meff_inc += leps[i].pT();
      }

      double mT = -100;
      if( Nlep > 0 ) mT = get_mT(leps[0].momentum(), met);

      double mT_minb = 10000;
      if(Nb < 1) mT_minb = -100;
	    for(int i=0; i<min(Nb, 3); i++){
	    	double dm = get_mT(bjets[i].momentum(), met);
	    	if(dm < mT_minb) mT_minb = dm;
	    }

	    double M_J = 0;
	    for(int i=0; i<min(N_J, 3); i++){
	    	M_J += large_jets[i].momentum().mass();
			}	

			double dphimin4j = -100;
			if(Njet >= 4){
				dphimin4j = 100;
		    for(int i=0; i<4; i++){
		    	double dphi = deltaPhi(jets[i], met);
		    	if(dphi < dphimin4j) dphimin4j = dphi;		    
		    }
		  }

		  int Nj70 = 0; int Nj30 = 0;
		  for(int i=0; i<Njet; i++){
		  	double pT = jets[i].pT();
		  	if(pT > 30) Nj30++;
		  	if(pT > 70) Nj70++;		  	
		  }

			//======== Gbb =========//
			if( cut( Njet, CUT_GE, 4, "Njet >= 4: SR-Gbb(base)" ) ){
				if( cut( Nlep_cand, CUT_EQ, 0, "Lepton Veto: SR-Gbb(base)" ) ){
					if( cut( dphimin4j, CUT_GT, 0.4, "dphimin4j > 0.4: SR-Gbb(base)" ) ){
						// Gbb-A
						if( cut( Nj70, CUT_GE, 4, "pTj(1,4) > 70: SR-Gbb-A" ) ){
							if( cut( Nb, CUT_GE, 3, "Nb >= 3: SR-Gbb-A" ) ){
								if( cut( MET, CUT_GT, 450, "MET > 450: SR-Gbb-A" ) ){
									if( cut( meff_4j, CUT_GT, 1900, "meff_4j > 1900: SR-Gbb-A" ) ){
										pass("SR-Gbb-A");
									}
								}
							}
						}
						// Gbb-A
						if( cut( Nj30, CUT_GE, 4, "pTj(1,4) > 30: SR-Gbb-B" ) ){
							if( cut( Nb, CUT_GE, 4, "Nb >= 4: SR-Gbb-B" ) ){
								if( cut( MET, CUT_GT, 300, "MET > 300: SR-Gbb-B" ) ){
									if( cut( meff_4j, CUT_GT, 1000, "meff_4j > 1000: SR-Gbb-B" ) ){
										pass("SR-Gbb-B");
									}
								}
							}
						}
					}
				}
			}

			//======== Gtt-0L =========//
			if( cut( Njet, CUT_GE, 8, "Njet >= 8: SR-Gtt-0L(base)" ) ){
				if( cut( Nb, CUT_GE, 3, "Nb >= 3: SR-Gtt-0L(base)" ) ){
					if( cut( Nlep_cand, CUT_EQ, 0, "Lepton Veto: SR-Gtt-0L(base)" ) ){
						if( cut( dphimin4j, CUT_GT, 0.4, "dphimin4j > 0.4: SR-Gtt-0L(base)" ) ){
							if( cut( mT_minb, CUT_GT, 80, "mT_minb > 80: SR-Gtt-0L(base)" ) ){
								if( cut( MET, CUT_GT, 400, "MET > 400: SR-Gtt-0L(base)" ) ){
									// Gtt-0L-A
									if( cut( meff_inc, CUT_GT, 2000, "meff_inc > 2000: SR-Gtt-0L-A" ) ){
										if( cut( M_J, CUT_GT, 200, "M_J > 200: SR-Gtt-0L-A" ) ){
											pass("SR-Gtt-0L-A");
										}
									}
									// Gtt-0L-B
									if( cut( meff_inc, CUT_GT, 1250, "meff_inc > 1250: SR-Gtt-0L-A" ) ){
										pass("SR-Gtt-0L-A");
									}
								}
							}
						}
					}
				}
			}

			//======== Gtt-1L =========//
			if( cut( Nlep, CUT_GE, 1, "Nlep(sig) >= 1: SR-Gtt-1L(base)" ) ){
				if( cut( Njet, CUT_GE, 6, "Njet >= 6: SR-Gtt-1L(base)" ) ){
					// A-B
					if( cut( Nb, CUT_GE, 3, "Nb >= 3: SR-Gtt-1L(A-B)" ) ){
						if( cut( mT, CUT_GT, 200, "mT > 200: SR-Gtt-1L(A-B)" ) ){
							if( cut( mT_minb, CUT_GT, 120, "mT_minb > 120: SR-Gtt-1L(A-B)" ) ){
								// Gtt-1L-A
								if( cut( MET, CUT_GT, 200, "MET > 200: SR-Gtt-1L-A" ) ){
									if( cut( meff_inc, CUT_GT, 2000, "meff_inc > 2000: SR-Gtt-1L-A" ) ){
										if( cut( M_J, CUT_GT, 200, "M_J > 200: SR-Gtt-1L-A" ) ){
											pass("SR-Gtt-1L-A");
										}
									}
								}
								// Gtt-1L-B
								if( cut( MET, CUT_GT, 350, "MET > 350: SR-Gtt-1L-B" ) ){
									if( cut( meff_inc, CUT_GT, 1500, "meff_inc > 1500: SR-Gtt-1L-B" ) ){
										if( cut( M_J, CUT_GT, 150, "M_J > 150: SR-Gtt-1L-B" ) ){
											pass("SR-Gtt-1L-B");
										}
									}
								}
							}
						}
					}
					// Gtt-1L-C
					if( cut( Nb, CUT_GE, 4, "Nb >= 4: SR-Gtt-1L-C" ) ){
						if( cut( mT, CUT_GT, 150, "mT > 150: SR-Gtt-1L-C" ) ){
							if( cut( mT_minb, CUT_GT, 80, "mT_minb > 80: SR-Gtt-1L-C" ) ){
								if( cut( MET, CUT_GT, 200, "MET > 200: SR-Gtt-1L-C" ) ){
									if( cut( meff_inc, CUT_GT, 500, "meff_inc > 500: SR-Gtt-1L-C" ) ){
										pass("SR-Gtt-1L-C");
									}
								}
							}
						}
					}
				}
			}

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
	AtomPlugin(ATLAS_CONF_2016_052)
}
