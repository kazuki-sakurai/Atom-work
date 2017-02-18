///
/// @file  ATLAS_1407_0600.cc
/// @brief Implementation of ATLAS_1501_07110 analysis
/// @author Kazuki Sakurai
/// @date created 04/15/2015
/// @date last revision 04/15/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1407_0600 : public Analysis {
	public:

		ATLAS_1407_0600()
			: Analysis("ATLAS_1407_0600") {
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

            useDetector( "ATLAS_CMS_all" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            
            IsoElectron base_ele( Range(PT, 20., 8000.) & Range(ETA, -2.47, 2.47) );
            base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoElectron ele( Range(PT, 20., 8000.) & Range(ETA, -2.47, 2.47) );
            //                         cone  frac  abs  inner 
            ele.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            ele.addIso(CALO_ISO_ET,  0.3,  0.18,  0.0, 0.01, CALO_ALL);
            ele.setVariableThreshold(0.0);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon base_mu(Range(PT, 10., 8000.) & Range(ETA, -2.5, 2.5));
            base_mu.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            base_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon mu( Range(PT, 20., 8000.) & Range(ETA, -2.5, 2.5) );
            //                         cone  frac  abs  inner 
            mu.addIso(TRACK_ISO_PT, 0.3,  0.12,  0.0, 0.01, CALO_ALL);
            mu.addIso(CALO_ISO_ET,  0.3,  0.23,  0.0, 0.01, CALO_ALL);
            mu.setVariableThreshold(0.0);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets base_jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -4.5, 4.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            base_jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            base_jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // Overlap removal

            NearIsoParticle base_jets_clean(base_jets);
            base_jets_clean.addFilter(base_ele, 0.2);

            NearIsoParticle jets_clean(base_jets_clean, Range(ETA, -2.8, 2.8));
            addProjection(jets_clean, "Jets");

            Range bjrange = Range(PT, 20.0, 8000.0) & Range(ETA, -2.5, 2.5);
            HeavyFlavorJets bjets(jets_clean, bjrange);
            bjets.setEfficiencyParams( getBJetEff("BJet_Ident_MV1_ATLAS") );
            addProjection(bjets, "BJets");

            MergedFinalState base_leptons(base_ele, base_mu);
            NearIsoParticle base_leptons_clean(base_leptons);
            base_leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(base_leptons_clean, "Base_Leptons");

            MergedFinalState leptons(ele, mu);
            NearIsoParticle leptons_clean(leptons);
            leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(leptons_clean, "Leptons");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(base_jets_clean, base_leptons_clean);            
            MissingMomentum met( fsbase, met_seed );            
            //MissingMomentum met( fsbase );                        
            //met.setSmearingParams( FinalState::SELECTED, &dp.metEff( "Jet_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

            /// @todo book the efficiencies
            bookEfficiency("0l-4j-A");
            bookEfficiency("0l-4j-B");
            bookEfficiency("0l-4j-C");
            bookEfficiency("0l-7j-A");
            bookEfficiency("0l-7j-B");
            bookEfficiency("0l-7j-C");
            bookEfficiency("1l-6j-A");
            bookEfficiency("1l-6j-B");
            bookEfficiency("1l-6j-C");

            /// @todo book the cuts
            bookCut(">= 4 jets: preselection");
            bookCut(">= 4 jets (pT > 30): preselection");
            bookCut(">= 1 jet (pT > 90): preselection");
            bookCut("MET > 150: preselection");
            bookCut(">= 3 bjets: preselection");
            bookCut(">= 3 bjets (pT > 30): preselection");
            bookCut("lepton veto: 0-lep preselection");
            bookCut("dPhi_4jmin > 0.5: 0-lep preselection");
            bookCut("MET/meff_4j > 0.2: 0-lep preselection");
            bookCut(">= 4 jets (pT > 50): 0l-4j(A, B)");
            bookCut(">= 3 bjets (pT > 50): 0l-4j(A, B)");
            bookCut("MET > 250: 0l-4j-A");
            bookCut("meff_4j > 1300: 0l-4j-A");
            bookCut("MET > 350: 0l-4j-B");
            bookCut("meff_4j > 1100: 0l-4j-B");
            bookCut("anti-b leading jet: 0l-4j-C");
            bookCut("MET > 400: 0l-4j-C");
            bookCut("meff_4j > 1000: 0l-4j-C");
            bookCut("MET/sqrt(HT_4j) > 16: 0l-4j-C");
            bookCut(">= 7 jets (pT > 30) : 0l-7j preselection");
            bookCut("MET > 200: 0l-7j-A");
            bookCut("meff_incl > 1000: 0l-7j-A");
            bookCut("MET > 350: 0l-7j-B");
            bookCut("meff_incl > 1000: 0l-7j-B");
            bookCut("MET > 250: 0l-7j-C");
            bookCut("meff_incl > 1500: 0l-7j-C");
            bookCut(">= 1 lep: 1-lep preselection");
            bookCut(">= 1 lep (pT > 25): 1-lep preselection");
            bookCut(">= 6 jets (pT > 30): 1-lep preselection");
            bookCut("MET > 175: 1l-6j-A");
            bookCut("mT > 140: 1l-6j-A");
            bookCut("meff_incl_lep > 700: 1l-6j-A");
            bookCut("MET > 225: 1l-6j-B");
            bookCut("mT > 140: 1l-6j-B");
            bookCut("meff_incl_lep > 800: 1l-6j-B");
            bookCut("MET > 275: 1l-6j-C");
            bookCut("mT > 160: 1l-6j-C");
            bookCut("meff_incl_lep > 900: 1l-6j-C");

            //bookCut("CutEtaJet","this is a control region cut", true);

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt(&event);
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt(&event);
            const Particles& base_leps = applyProjection<NearIsoParticle>(event, "Base_Leptons").particlesByPt(&event);
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(&event); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

            double MET = met.pT();
            // ********************************************* //

            // preselection
            if( !cut( jets.size(), CUT_GE, 4, ">= 4 jets: preselection" )) vetoEvent;            
            if( !cut( jets[3].pT(), CUT_GT, 30., ">= 4 jets (pT > 30): preselection" )) vetoEvent;
            if( !cut( jets[0].pT(), CUT_GT, 90., ">= 1 jet (pT > 90): preselection" )) vetoEvent;
            if( !cut( MET, CUT_GT, 150., "MET > 150: preselection" )) vetoEvent;
            if( !cut( bjets.size(), CUT_GE, 3, ">= 3 bjets: preselection" )) vetoEvent;            
            if( !cut( bjets[2].pT(), CUT_GT, 30., ">= 3 bjets (pT > 30): preselection" )) vetoEvent;

            // variable definition

            double meff_incl = MET;
            double meff_4j = MET;
            double Nj30 = 0;          
            for( int i=0; i<jets.size(); i++){
                if( jets[i].pT() > 30. ){ 
                    meff_incl += jets[i].pT();
                    Nj30++;
                }
                if( i < 4 ) meff_4j += jets[i].pT();
            }

            double meff_incl_lep = meff_incl;
            for( int i=0; i<leps.size(); i++){
                if( leps[i].pT() > 20. ) meff_incl_lep += leps[i].pT();
            } 

            double delPhi_4jmin = 1000.;
            for( int i=0; i<4; i++){
                double dPhi = deltaPhi( jets[i].momentum(), met );
                if( delPhi_4jmin < dPhi) delPhi_4jmin = dPhi;
            }

            double HT_4j = meff_4j - MET;

            double mT = 0;
            if( leps.size() > 0 ){
                mT = get_mT(leps[0].momentum(), met);
            }

            //=============================================//

            // 0-lepton
            if( cut( base_leps.size(), CUT_EQ, 0, "lepton veto: 0-lep preselection" )){
                if( cut(delPhi_4jmin, CUT_GT, 0.5, "dPhi_4jmin > 0.5: 0-lep preselection") ){
                    if( cut(MET/meff_4j, CUT_GT, 0.2, "MET/meff_4j > 0.2: 0-lep preselection") ){

                        // 4j
                        if( cut( jets[3].pT(), CUT_GT, 50., ">= 4 jets (pT > 50): 0l-4j(A, B)") ){
                            if( cut( bjets[2].pT(), CUT_GT, 50., ">= 3 bjets (pT > 50): 0l-4j(A, B)") ){
                                if( cut( MET, CUT_GT, 250., "MET > 250: 0l-4j-A") ){
                                    if( cut( meff_4j, CUT_GT, 1300., "meff_4j > 1300: 0l-4j-A") ){
                                        pass("0l-4j-A");
                                    }
                                }
                                if( cut( MET, CUT_GT, 350., "MET > 350: 0l-4j-B") ){
                                    if( cut( meff_4j, CUT_GT, 1100., "meff_4j > 1100: 0l-4j-B") ){
                                        pass("0l-4j-B");
                                    }
                                }
                            }
                        }

                        if(untagged_jets.size() > 0){
                            if( cut( untagged_jets[0].pT(), CUT_GT, jets[1].pT(), "anti-b leading jet: 0l-4j-C") ){
                                if( cut(MET, CUT_GT, 400., "MET > 400: 0l-4j-C") ){
                                    if( cut(meff_4j, CUT_GT, 1000., "meff_4j > 1000: 0l-4j-C") ){
                                        if( cut(MET/sqrt(HT_4j), CUT_GT, 16., "MET/sqrt(HT_4j) > 16: 0l-4j-C") ){
                                            pass("0l-4j-C");   
                                        }
                                    }
                                }
                            }
                        }

                        // 7j
                        if( cut(Nj30, CUT_GE, 7, ">= 7 jets (pT > 30) : 0l-7j preselection") ){

                            if( cut(MET, CUT_GT, 200., "MET > 200: 0l-7j-A") ){
                                if( cut(meff_incl, CUT_GT, 1000., "meff_incl > 1000: 0l-7j-A") ){
                                    pass("0l-7j-A");
                                }
                            
                            }
                            if( cut(MET, CUT_GT, 350., "MET > 350: 0l-7j-B") ){
                                if( cut(meff_incl, CUT_GT, 1000., "meff_incl > 1000: 0l-7j-B") ){
                                    pass("0l-7j-B");
                                }
                            
                            }
                            if( cut(MET, CUT_GT, 250., "MET > 250: 0l-7j-C") ){
                                if( cut(meff_incl, CUT_GT, 1500., "meff_incl > 1500: 0l-7j-C") ){
                                    pass("0l-7j-C");
                                }
                            
                            }

                        }
                    }
                }
            }


            // 1-lepton
            if( !cut( leps.size(), CUT_GE, 1, ">= 1 lep: 1-lep preselection" )) vetoEvent;
            if( !cut( leps[0].pT(), CUT_GT, 25, ">= 1 lep (pT > 25): 1-lep preselection" )) vetoEvent;
            if( !cut( Nj30, CUT_GE, 6, ">= 6 jets (pT > 30): 1-lep preselection") ) vetoEvent;

            if( cut( MET, CUT_GT, 175., "MET > 175: 1l-6j-A") ){
                if( cut( mT, CUT_GT, 140., "mT > 140: 1l-6j-A") ){
                    if( cut( meff_incl_lep, CUT_GT, 700., "meff_incl_lep > 700: 1l-6j-A") ){
                        pass("1l-6j-A");
                    }
                }
            }
            if( cut( MET, CUT_GT, 225., "MET > 225: 1l-6j-B") ){
                if( cut( mT, CUT_GT, 140., "mT > 140: 1l-6j-B") ){
                    if( cut( meff_incl_lep, CUT_GT, 800., "meff_incl_lep > 800: 1l-6j-B") ){
                        pass("1l-6j-B");
                    }
                }
            }
            if( cut( MET, CUT_GT, 275., "MET > 275: 1l-6j-C") ){
                if( cut( mT, CUT_GT, 160., "mT > 160: 1l-6j-C") ){
                    if( cut( meff_incl_lep, CUT_GT, 900., "meff_incl_lep > 900: 1l-6j-C") ){
                        pass("1l-6j-C");
                    }
                }
            }

		}

		/// Normalise histograms etc., after the run
		void finalize() {
			/// @todo normalize the histograms
			// scale("Mjj");
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1407_0600)
}
