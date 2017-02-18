///
/// @file  1405_7875.cc
/// @brief Implementation of 1405_7875 analysis
/// @author Kazuki
/// @date created 04/15/2015
/// @date last revision 04/15/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1405_7875 : public Analysis {
	public:

		ATLAS_1405_7875()
			: Analysis("ATLAS_1405_7875") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

        int nev;
        double cut1, cut2, cut3, cut4, cut5;

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            
            IsoElectron ele_base( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            ele_base.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele_base.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoMuon mu_base(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            mu_base.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -4.5, 4.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // Overlap removal
            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(ele_base, 0.2);

            MergedFinalState base_leptons(ele_base, mu_base);
            NearIsoParticle leptons_clean(base_leptons);
            leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(leptons_clean, "Leptons");

            NearIsoParticle signal_jets(jets_clean, Range(ETA, -2.8, 2.8));
            addProjection(signal_jets, "Jets");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets_clean, base_leptons);            
            MissingMomentum met( fsbase, met_seed );            
            //MissingMomentum met( fsbase );                        
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Const_PlaceHolder" ) );
            met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Interpolation_PlaceHolder" ) );            
            addProjection(met, "MissingEt");

            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

            /// @todo book the efficiencies

            bookEfficiency("2jl");
            bookEfficiency("2jm");
            bookEfficiency("2jt");
            bookEfficiency("2jW");
            bookEfficiency("4jW");
            bookEfficiency("3j");
            bookEfficiency("4jl-");
            bookEfficiency("4jl");
            bookEfficiency("4jm");
            bookEfficiency("4jt");
            bookEfficiency("5j");
            bookEfficiency("6jl");
            bookEfficiency("6jm");
            bookEfficiency("6jt");
            bookEfficiency("6jt+");


            /// @todo book the cuts
            bookCut("lepton veto");
            bookCut("MET > 160");
            bookCut("pTj1 > 130");
            bookCut("pTj2 > 60");
            bookCut("pTj3 > 60");
            bookCut("pTj4 > 60");
            bookCut("pTj5 > 60");
            bookCut("pTj6 > 60");
            bookCut("dPhiMin_123 > 0.4: 2j");
            bookCut("MET/sqrt(HT) > 8: 2jl");
            bookCut("meff_inc > 800: 2jl");
            bookCut("MET/sqrt(HT) > 15: 2j(m,t)");
            bookCut("meff_inc > 1200: 2jm");
            bookCut("meff_inc > 1600: 2jt");
            bookCut("2(W->j): 2jW");
            bookCut("MET/meff_Nj > 0.25: 2jW");
            bookCut("meff_inc > 1800: 2jW");
            bookCut("pTj4 > 40: 4jW");
            bookCut("dPhiMin_123 > 0.4: 4jW");
            bookCut("dPhiMin_all > 0.2: 4jW");
            bookCut("W->j: 4jW");
            bookCut("W->jj: 4jW");
            bookCut("MET/meff_Nj > 0.35: 4jW");
            bookCut("meff_inc > 1100: 4jW");
            bookCut("dPhiMin_123 > 0.4: 3j");
            bookCut("MET/meff_Nj > 0.3: 3j");
            bookCut("meff_inc > 2200: 3j");
            bookCut("dPhiMin_123 > 0.4: 4j");
            bookCut("dPhiMin_all > 0.2: 4j");
            bookCut("MET/sqrt(HT) > 10: 4j(l, l-)");
            bookCut("meff_inc > 700: 4jl-");
            bookCut("meff_inc > 1000: 4jl");
            bookCut("MET/meff_Nj > 0.4: 4jm");
            bookCut("meff_inc > 1300: 4jm");
            bookCut("MET/meff_Nj > 0.25: 4jt");
            bookCut("meff_inc > 2200: 4jt");
            bookCut("dPhiMin_123 > 0.4: 5j");
            bookCut("dPhiMin_all > 0.2: 5j");
            bookCut("MET/meff_Nj > 0.2: 5j");
            bookCut("meff_inc > 1200: 5j");
            bookCut("dPhiMin_123 > 0.4: 6j");
            bookCut("dPhiMin_all > 0.2: 6j");
            bookCut("MET/meff_Nj > 0.2: 6j(l,m)");
            bookCut("meff_inc > 900: 6jl");
            bookCut("meff_inc > 1200: 6jm");
            bookCut("MET/meff_Nj > 0.2: 6jt");
            bookCut("meff_inc > 1500: 6jt");
            bookCut("MET/meff_Nj > 0.15: 6jt+");
            bookCut("meff_inc > 1700: 6jt+");
            //bookCut("CutEtaJet","this is a control region cut", true);

            nev = 0;
            cut1 = 0;
            cut2 = 0;
            cut3 = 0;
            cut4 = 0;
            cut5 = 0;

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            nev++;
            //cout << "Event " << nev << endl;

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();
            //const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            //const Particles& bjets = bjproj.getTaggedJets(&event); 
            //const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

            //cout << jets.size() <<"  "<< leps.size() << endl;

            double MET = met.pT();

            int Njet = jets.size();

            double meff_Nj[6] = {};
            for(int i = 0; i < min(Njet, 6); i++){
                meff_Nj[i] = MET;                
                for(size_t j = 0; j < i + 1; j++){
                    meff_Nj[i] += jets[j].pT();
                }
            } 

            double H_T = 0; 
            for(int i = 0; i < Njet; i++){
                if( jets[i].pT() > 40. ){
                    H_T += jets[i].pT();
                }
            } 
            double meff_inc = H_T + MET; 

            double MET_over_meff_Nj[6] = {};
            for(int i = 0; i < min(Njet, 6); i++){
                MET_over_meff_Nj[i] = MET/meff_Nj[i];
            }

            double dPhiMin_123 = 1000;
            for(int i=0; i < min(Njet, 3); i++){
                if(jets[i].pT() > 40.){
                    double dPhi = deltaPhi(jets[i].momentum(), met);
                    if( dPhi < dPhiMin_123 ) dPhiMin_123 = dPhi;
                }
            }

            double dPhiMin_all = 1000;
            for(int i=0; i < Njet; i++){
                if(jets[i].pT() > 40.){
                    double dPhi = deltaPhi(jets[i].momentum(), met);
                    if( dPhi < dPhiMin_all ) dPhiMin_all = dPhi;
                }
            }

            int N_Wj_candi = 0;
            for(int i=0; i<Njet; i++){
                if( 60 < jets[i].mass() && jets[i].mass() < 100) 
                    N_Wj_candi++;
            }

            int N_W2j_candi = 0;
            // define this!

            //#######################################


            // ********************************************* //

            if( !cut( leps.size(), CUT_EQ, 0, "lepton veto" ) ) vetoEvent;
            if( !cut( MET, CUT_GT, 160., "MET > 160" ) ) vetoEvent;

            if( jets.size() < 2 ) vetoEvent;
            if( !cut( jets[0].pT(), CUT_GT, 130., "pTj1 > 130" ) ) vetoEvent;
            if( !cut( jets[1].pT(), CUT_GT, 60., "pTj2 > 60" ) ) vetoEvent;

            // 2j SRs
            if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 2j" ) ){

                if( cut( MET/sqrt(H_T), CUT_GT, 8., "MET/sqrt(HT) > 8: 2jl" ) ){
                    if( cut( meff_inc, CUT_GT, 800., "meff_inc > 800: 2jl" ) ){
                        pass("2jl");
                    }
                }
                
                if( cut( MET/sqrt(H_T), CUT_GT, 15., "MET/sqrt(HT) > 15: 2j(m,t)" ) ){

                    if( cut( meff_inc, CUT_GT, 1200., "meff_inc > 1200: 2jm" ) ){
                        pass("2jm");
                    }
                    if( cut( meff_inc, CUT_GT, 1600., "meff_inc > 1600: 2jt" ) ){
                        pass("2jt");
                    }                    
                }

                if( cut( N_Wj_candi, CUT_GE, 2, "2(W->j): 2jW" ) ){
                    if( cut( MET/meff_Nj[1], CUT_GT, 0.25, "MET/meff_Nj > 0.25: 2jW" ) ){
                        if( !cut( meff_inc, CUT_GT, 1800., "meff_inc > 1800: 2jW" ) ){
                            pass("2jW");
                        }
                    }
                }

            }

            // 4jW
            if( jets.size() > 3){
                if( cut( jets[3].pT(), CUT_GT, 40., "pTj4 > 40: 4jW" ) ){
                    if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 4jW" ) ){
                        if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 4jW" ) ){
                            if( cut( N_Wj_candi, CUT_GE, 1, "W->j: 4jW" ) ){
                                if( cut( N_W2j_candi, CUT_GE, 1, "W->jj: 4jW" ) ){                                        
                                    if( cut( MET/meff_Nj[3], CUT_GT, 0.35, "MET/meff_Nj > 0.35: 4jW" ) ){
                                        if( cut( meff_inc, CUT_GT, 1100., "meff_inc > 1100: 4jW" ) ){
                                            pass("4jW");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // 3j 
            if( jets.size() < 3) vetoEvent;
            if( !cut( jets[2].pT(), CUT_GT, 60., "pTj3 > 60" ) ) vetoEvent;
                    
            if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 3j" ) ){
                if( cut( MET/meff_Nj[2], CUT_GT, 0.3, "MET/meff_Nj > 0.3: 3j" ) ){                    
                    if( cut( meff_inc, CUT_GT, 2200., "meff_inc > 2200: 3j" ) ){
                        pass("3j");
                    }
                }                
            }

            // 4j SRs
            if( jets.size() < 4) vetoEvent;
            if( !cut( jets[3].pT(), CUT_GT, 60., "pTj4 > 60" ) ) vetoEvent;

            if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 4j" ) ){
                if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 4j" ) ){

                    if( cut( MET/sqrt(H_T), CUT_GT, 10., "MET/sqrt(HT) > 10: 4j(l, l-)" ) ){
                        if( cut( meff_inc, CUT_GT, 700., "meff_inc > 700: 4jl-" ) ){
                            pass("4jl-");
                        }                    
                        if( cut( meff_inc, CUT_GT, 1000., "meff_inc > 1000: 4jl" ) ){
                            pass("4jl");
                        }                    
                    }

                   if( cut( MET/meff_Nj[3], CUT_GT, 0.4, "MET/meff_Nj > 0.4: 4jm" ) ){                    
                        if( cut( meff_inc, CUT_GT, 1300., "meff_inc > 1300: 4jm" ) ){
                            pass("4jm");
                        }
                    }
                   if( cut( MET/meff_Nj[3], CUT_GT, 0.25, "MET/meff_Nj > 0.25: 4jt" ) ){
                        if( cut( meff_inc, CUT_GT, 2200., "meff_inc > 2200: 4jt" ) ){
                            pass("4jt");
                        }
                    }

                }
            }             

            // 5j 
            if( jets.size() < 5) vetoEvent;
            if( !cut( jets[4].pT(), CUT_GT, 60., "pTj5 > 60" ) ) vetoEvent;

            if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 5j" ) ){
                if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 5j" ) ){
                   if( cut( MET/meff_Nj[4], CUT_GT, 0.2, "MET/meff_Nj > 0.2: 5j" ) ){                    
                        if( cut( meff_inc, CUT_GT, 1200., "meff_inc > 1200: 5j" ) ){
                            pass("5j");
                        }
                    }
                }
            }

            // 6j SRs 
            if( jets.size() < 6) vetoEvent;
            if( !cut( jets[5].pT(), CUT_GT, 60., "pTj6 > 60" ) ) vetoEvent;

            if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 6j" ) ){
                if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 6j" ) ){

                   if( cut( MET/meff_Nj[5], CUT_GT, 0.2, "MET/meff_Nj > 0.2: 6j(l,m)" ) ){
                        if( cut( meff_inc, CUT_GT, 900., "meff_inc > 900: 6jl" ) ){
                            pass("6jl");
                        }
                        if( cut( meff_inc, CUT_GT, 1200., "meff_inc > 1200: 6jm" ) ){
                            pass("6jm");
                        }                        
                    }

                   if( cut( MET/meff_Nj[5], CUT_GT, 0.25, "MET/meff_Nj > 0.2: 6jt" ) ){
                        if( cut( meff_inc, CUT_GT, 1500., "meff_inc > 1500: 6jt" ) ){
                            pass("6jt");
                        }
                    }
                   if( cut( MET/meff_Nj[5], CUT_GT, 0.15, "MET/meff_Nj > 0.15: 6jt+" ) ){
                        if( cut( meff_inc, CUT_GT, 1700., "meff_inc > 1700: 6jt+" ) ){
                            pass("6jt+");
                        }
                    }

                }
            }

		}

		/// Normalise histograms etc., after the run
		void finalize() {

            cout << "nev " << nev << endl;
            cout << "cut1 " << cut1/nev << endl;
            cout << "cut2 " << cut2/nev << endl;
            cout << "cut3 " << cut3/nev << endl;
            cout << "cut4 " << cut4/nev << endl;
            cout << "cut5 " << cut5/nev << endl;

			/// @todo normalize the histograms
			// scale("Mjj");
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1405_7875)
}
