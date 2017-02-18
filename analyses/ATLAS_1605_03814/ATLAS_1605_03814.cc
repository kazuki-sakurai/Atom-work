///
/// @file  ATLAS_1605_03814.cc
/// @brief Implementation of ATLAS_1605_03814 analysis
/// @author Kazuki
/// @date created 06/18/2016
/// @date last revision 06/18/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1605_03814 : public Analysis {
	public:

		ATLAS_1605_03814()
			: Analysis("ATLAS_1605_03814") {
			//setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            
            IsoElectron ele_base( Range(PT > 10. & abseta < 2.47) );
            ele_base.addConeIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele_base.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoMuon mu_base(Range(PT > 10. & abseta < 2.5));
            mu_base.addConeIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range( PT > 20 & abseta < 2.8 ), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // Overlap removal
            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(ele_base, 0.2);

            NearIsoParticle signal_jets(jets_clean, Range(PT > 50.));
            addProjection(signal_jets, "Jets");

            //addProjection(Sphericity(signal_jets), "Sphericity");
            //addProjection(Sphericity(fsbase), "Sphericity");

            HeavyFlavorJets bjets(jets_clean, Range(PT > 50. & abseta < 2.5));
            bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
            bjets.setCurrentWorkingPoint( 0.77 );
            addProjection(bjets, "BJets");

            MergedFinalState base_leptons(ele_base, mu_base);
            NearIsoParticle leptons_clean(base_leptons);
            leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(leptons_clean, "Leptons");

            MergedFinalState met_seed(jets_clean, base_leptons);
            MissingMomentum met( fsbase, met_seed );
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

            /// @todo book the efficiencies

            bookEfficiency("2jl");
            bookEfficiency("2jm");
            bookEfficiency("2jt");
            bookEfficiency("4jt");
            bookEfficiency("5j");
            bookEfficiency("6jm");
            bookEfficiency("6jt");


            /// @todo book the cuts
            bookCut("lepton veto");
            bookCut("MET > 200");
            bookCut("pTj1 > 200");
            bookCut("Nj >= 2: 2jl");
            bookCut("dPhiMin_123 > 0.8: 2jl");
            bookCut("pTj2 > 200: 2jl");
            bookCut("MET/sqrt(HT) > 15: 2jl");
            bookCut("meff_inc > 1200: 2jl");
            bookCut("pTj1 > 300: 2jm");
            bookCut("Nj >= 2: 2jm");
            bookCut("dPhiMin_123 > 0.4: 2jm");
            bookCut("pTj2 > 50: 2jm");
            bookCut("MET/sqrt(HT) > 15: 2jm");
            bookCut("meff_inc > 1600: 2jm");
            bookCut("Nj >= 2: 2jt");
            bookCut("dPhiMin_123 > 0.8: 2jt");
            bookCut("pTj2 > 200: 2jt");
            bookCut("MET/sqrt(HT) > 20: 2jt");
            bookCut("meff_inc > 2000: 2jt");
            bookCut("Nj >= 4: 4jt");
            bookCut("dPhiMin_123 > 0.4: 4jt");
            bookCut("dPhiMin_all > 0.2: 4jt");            
            bookCut("pTj2 > 100: 4jt");
            bookCut("pTj4 > 100: 4jt");
            bookCut("Aplanarity > 0.04: 4jt");
            bookCut("MET/meff_Nj > 0.2: 4jt");
            bookCut("meff_inc > 2200: 4jt");
            bookCut("Nj >= 5: 5j");
            bookCut("dPhiMin_123 > 0.4: 5j");
            bookCut("dPhiMin_all > 0.2: 5j");
            bookCut("pTj2 > 100: 5j");
            bookCut("pTj5 > 50: 5j");
            bookCut("Aplanarity > 0.04: 5j");
            bookCut("MET/meff_Nj > 0.25: 5j");
            bookCut("meff_inc > 1600: 5j");
            bookCut("Nj >= 6: 6jm");
            bookCut("dPhiMin_123 > 0.4: 6jm");
            bookCut("dPhiMin_all > 0.2: 6jm");
            bookCut("pTj2 > 100: 6jm");
            bookCut("pTj6 > 50: 6jm");
            bookCut("Aplanarity > 0.04: 6jm");
            bookCut("MET/meff_Nj > 0.25: 6jm");
            bookCut("meff_inc > 1600: 6jm");
            bookCut("Nj >= 6: 6jt");
            bookCut("dPhiMin_123 > 0.4: 6jt");
            bookCut("dPhiMin_all > 0.2: 6jt");
            bookCut("pTj2 > 100: 6jt");
            bookCut("pTj6 > 50: 6jt");
            bookCut("Aplanarity > 0.04: 6jt");
            bookCut("MET/meff_Nj > 0.2: 6jt");
            bookCut("meff_inc > 2000: 6jt");

			// Projection booking section -- do not edit/remove this comment
            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

			// Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            //bookHisto1D(1,1,1, "meff_inc_2jl");
            // bookHisto1D(1,1,2, "meff_inc_2jm");
            // bookHisto1D(1,1,3, "meff_inc_2jt");

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
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            //const Sphericity& spher = applyProjection<Sphericity>(event, "Sphericity");

            //cout << jets.size() <<"  "<< leps.size() << endl;

            double MET = met.pT();
            int Njet = jets.size();
            //double aplanarity = spher.aplanarity();
            double aplanarity = 1;
            //cout << aplanarity << endl;

            double meff_Nj[6] = {};
            for(int i = 0; i < min(Njet, 6); i++){
                meff_Nj[i] = MET;                
                for(size_t j = 0; j < i + 1; j++){
                    meff_Nj[i] += jets[j].pT();
                }
            } 

            double H_T = 0;
            for(int i = 0; i < Njet; i++){
                H_T += jets[i].pT();
            } 
            double meff_inc = H_T + MET; 

            double dPhiMin_123 = 1000;
            for(int i=0; i < min(Njet, 3); i++){
                double dPhi = deltaPhi(jets[i].momentum(), met);
                if( dPhi < dPhiMin_123 ) dPhiMin_123 = dPhi;
            }

            double dPhiMin_all = 1000;
            for(int i=0; i < Njet; i++){
                double dPhi = deltaPhi(jets[i].momentum(), met);
                if( dPhi < dPhiMin_all ) dPhiMin_all = dPhi;
            }
            //#######################################


            // ********************************************* //

            if( !cut( leps.size(), CUT_EQ, 0, "lepton veto" ) ) vetoEvent;
            if( !cut( MET, CUT_GT, 200., "MET > 200" ) ) vetoEvent;

            if( jets.size() < 1 ) vetoEvent;
            if( !cut( jets[0].pT(), CUT_GT, 200., "pTj1 > 200" ) ) vetoEvent;

            // 2jl SRs
            if( cut( jets.size(), CUT_GE, 2, "Nj >= 2: 2jl" ) ){            
	            if( cut( dPhiMin_123, CUT_GT, 0.8, "dPhiMin_123 > 0.8: 2jl" ) ){
		            if( cut( jets[1].pT(), CUT_GT, 200., "pTj2 > 200: 2jl" ) ){
	                    if( cut( MET/sqrt(H_T), CUT_GT, 8., "MET/sqrt(HT) > 15: 2jl" ) ){
	                        //fillPlot("meff_inc_2jl", meff_inc);                        
	                        if( cut( meff_inc, CUT_GT, 1200., "meff_inc > 1200: 2jl" ) ){
	                            pass("2jl");
	                        }
	                    }
	                }
	            }
	        }

            // 2jm SRs
	        if( cut( jets[0].pT(), CUT_GT, 300., "pTj1 > 300: 2jm" ) ){
	            if( cut( jets.size(), CUT_GE, 2, "Nj >= 2: 2jm" ) ){            
	                if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 2jm" ) ){
	                	if( cut( jets[1].pT(), CUT_GT, 50., "pTj2 > 50: 2jm" ) ){	                    	
	                        if( cut( MET/sqrt(H_T), CUT_GT, 15., "MET/sqrt(HT) > 15: 2jm" ) ){
	                            //fillPlot("meff_inc_2jm", meff_inc);                                                    
	                            if( cut( meff_inc, CUT_GT, 1600., "meff_inc > 1600: 2jm" ) ){
	                                pass("2jm");
	                            }
	                        }
	                    }
	                }
	            }
	        }


            // 2jt SRs
	        if( cut( jets.size(), CUT_GE, 2, "Nj >= 2: 2jt" ) ){            
                if( cut( dPhiMin_123, CUT_GT, 0.8, "dPhiMin_123 > 0.8: 2jt" ) ){
      	            if( cut( jets[1].pT(), CUT_GT, 200., "pTj2 > 200: 2jt" ) ){                    	
                        if( cut( MET/sqrt(H_T), CUT_GT, 20., "MET/sqrt(HT) > 20: 2jt" ) ){
                            //fillPlot("meff_inc_2jm", meff_inc);                                                                                
                            if( cut( meff_inc, CUT_GT, 2000., "meff_inc > 2000: 2jt" ) ){
                                pass("2jt");
                            }
                        }
                    }
                }
            }

            // 4jt
	        if( cut( jets.size(), CUT_GE, 4, "Nj >= 4: 4jt" ) ){                        	
                if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 4jt" ) ){
                    if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 4jt" ) ){                    
	                    if( cut( jets[1].pT(), CUT_GT, 100., "pTj2 > 100: 4jt" ) ){
	                	    if( cut( jets[3].pT(), CUT_GT, 100., "pTj4 > 100: 4jt" ) ){
                            	if( cut( aplanarity, CUT_GT, 0.04, "Aplanarity > 0.04: 4jt" ) ){
                                	if( cut( MET/meff_Nj[3], CUT_GT, 0.2, "MET/meff_Nj > 0.2: 4jt" ) ){
                                    	if( cut( meff_inc, CUT_GT, 2200., "meff_inc > 2200: 4jt" ) ){
                                        	pass("4jt");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // 5j
	        if( cut( jets.size(), CUT_GE, 5, "Nj >= 5: 5j" ) ){                        	            	
                if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 5j" ) ){
                    if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 5j" ) ){
                 		if( cut( jets[1].pT(), CUT_GT, 100., "pTj2 > 100: 5j" ) ){                   	
	                 		if( cut( jets[4].pT(), CUT_GT, 50., "pTj5 > 50: 5j" ) ){
	                            if( cut( aplanarity, CUT_GT, 0.04, "Aplanarity > 0.04: 5j" ) ){                            
    	                            if( cut( MET/meff_Nj[4], CUT_GT, 0.25, "MET/meff_Nj > 0.25: 5j" ) ){
        	                            if( cut( meff_inc, CUT_GT, 1600., "meff_inc > 1600: 5j" ) ){
            	                            pass("5j");
                	                    }
                	                }
                                }
                            }
                        }
                    }
                }
            }

            // 6jm
	        if( cut( jets.size(), CUT_GE, 6, "Nj >= 6: 6jm" ) ){                        	            	
                if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 6jm" ) ){
                    if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 6jm" ) ){
		                if( cut( jets[1].pT(), CUT_GT, 100., "pTj2 > 100: 6jm" ) ){
			                if( cut( jets[5].pT(), CUT_GT,  50., "pTj6 > 50: 6jm" ) ){
	                            if( cut( aplanarity, CUT_GT, 0.04, "Aplanarity > 0.04: 6jm" ) ){                            
    	                            if( cut( MET/meff_Nj[5], CUT_GT, 0.25, "MET/meff_Nj > 0.25: 6jm" ) ){
        	                            if( cut( meff_inc, CUT_GT, 1600., "meff_inc > 1600: 6jm" ) ){
            	                            pass("6jm");
            	                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // 6jt
	        if( cut( jets.size(), CUT_GE, 6, "Nj >= 6: 6jt" ) ){                        	            	            	
                if( cut( dPhiMin_123, CUT_GT, 0.4, "dPhiMin_123 > 0.4: 6jt" ) ){
                    if( cut( dPhiMin_all, CUT_GT, 0.2, "dPhiMin_all > 0.2: 6jt" ) ){
		                if( cut( jets[1].pT(), CUT_GT, 100., "pTj2 > 100: 6jt" ) ){
			                if( cut( jets[5].pT(), CUT_GT, 50., "pTj6 > 50: 6jt" ) ){
    	                        if( cut( aplanarity, CUT_GT, 0.04, "Aplanarity > 0.04: 6jt" ) ){                            
        	                        if( cut( MET/meff_Nj[5], CUT_GT, 0.2, "MET/meff_Nj > 0.2: 6jt" ) ){
            	                        if( cut( meff_inc, CUT_GT, 2000., "meff_inc > 2000: 6jt" ) ){
                	                        pass("6jt");
                    	                }
                    	            }
                                }
                            }
                        }
                    }
                }
            }

			// End analyze section -- do not edit/remove this comment
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			// Histogram normalization section -- do not edit/remove this comment
			/// @todo normalize the histograms
			//scale("meff_inc_2jl");
            //scale("meff_inc_2jm");
            //scale("meff_inc_2jt");
            //binIntegrate("meff_inc");
			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1605_03814)
}
