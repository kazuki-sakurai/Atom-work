 ///
/// @file  ATLAS_1403_5294.cc
/// @brief Implementation of ATLAS_1403_5294 analysis
/// @author Kazuki Sakurai <kazuki.sakurai@kcl.ac.uk>
/// @date created 04/21/2015
/// @date last revision 04/21/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"
#include "Rivet/Tools/RivetMT2.hh"

using namespace std;

namespace Atom {

	class ATLAS_1403_5294 : public Analysis {
	public:

		ATLAS_1403_5294()
			: Analysis("ATLAS_1403_5294") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT, 15., 8000.) & Range(ETA, -4.9, 4.9), 
                            //hadRange & Range(PT, 20., 8000.) & Range(ETA, -4.9, 4.9),                             
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            IsoElectron ele_base_0( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            ele_base_0.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            ele_base_0.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele_base_0.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoElectron ele_0( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            ele_0.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            ele_0.addIso(CALO_ISO_ET,  0.3,  0.18,  0.0, 0.01, CALO_ALL);            
            ele_0.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele_0.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon mu_1(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            mu_1.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            mu_1.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_1.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS_tune" ) );

            IsoMuon mu_base_1(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            mu_base_1.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            mu_base_1.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_base_1.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            // Overlap removal

            NearIsoParticle ele_base_1(ele_base_0);
            ele_base_1.addFilter(ele_base_0, 0.05, NearIsoParticle::DROP_LOWER_PT, 0.01);

            NearIsoParticle ele_1(ele_0);
            ele_1.addFilter(ele_0, 0.05, NearIsoParticle::DROP_LOWER_PT, 0.01);

            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(ele_base_1, 0.2);

            NearIsoParticle C_jets(jets_clean, Range(abseta < 2.4) );            
            addProjection(C_jets, "C_Jets");
            HeavyFlavorJets bjets(jets_clean, Range(ETA, -2.4, 2.4));
            bjets.setTaggingEfficiency( *getBJetEff("BJet_Ident_MV1_ATLAS"));
            bjets.setCurrentWorkingPoint( 0.8 );
            addProjection(bjets, "BJets");

            const double dRprong = 0.2;
            ParamTauFinder taus(jets_clean, dRprong, Range(ETA, -2.5, 2.5));
            NearIsoParticle taus_clean(taus);
            taus_clean.addFilter(mu_base_1, 0.2);
            addProjection(taus_clean, "Taus");

            NearIsoParticle F_jets(jets_clean, Range(PT, 30.0, 8000.0) & (Range(ETA, -4.5, 4.5) - Range(ETA, -2.4, 2.4)) );
            addProjection(F_jets, "F_Jets");

            // for signal leptons
            NearIsoParticle ele_2(ele_1);
            ele_2.addFilter(jets_clean, 0.4);

            NearIsoParticle mu_2(mu_1);
            mu_2.addFilter(jets_clean, 0.4);

            NearIsoParticle mu_3(mu_2);
            mu_3.addFilter(mu_2, 0.05, NearIsoParticle::DROP_LOWER_PT, 0.01);

            MergedFinalState leptons(ele_2, mu_3);
            addProjection(leptons, "Leptons");

            // for lepton candidates
            NearIsoParticle ele_base_2(ele_base_1);
            ele_base_2.addFilter(jets_clean, 0.4);

            NearIsoParticle mu_base_2(mu_base_1);
            mu_base_2.addFilter(jets_clean, 0.4);

            NearIsoParticle mu_base_3(mu_base_2);
            mu_base_3.addFilter(mu_base_2, 0.05, NearIsoParticle::DROP_LOWER_PT, 0.01);

            MergedFinalState base_leptons(ele_base_2, mu_base_3);
            addProjection(base_leptons, "baseLeptons");

            MergedFinalState met_seed(jets_clean, base_leptons);
            MissingMomentum met( fsbase, met_seed );
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

            /// @todo book the efficiencies
            bookEfficiency("mT2_90:SF");
            bookEfficiency("mT2_120:SF");
            bookEfficiency("mT2_150:SF");
            bookEfficiency("WWa:SF");
            bookEfficiency("WWb:SF");
            bookEfficiency("WWc:SF");
            bookEfficiency("Zjets");
            bookEfficiency("mT2_90:DF");
            bookEfficiency("mT2_120:DF");
            bookEfficiency("mT2_150:DF");
            bookEfficiency("WWa:DF");
            bookEfficiency("WWb:DF");
            bookEfficiency("WWc:DF");

            /// @todo book the cuts
	        bookCut("= 2 signal lepsons");
            bookCut("Opposite Charge");
            bookCut("pTl1 > 35");
            bookCut("pTl2 > 20");
            bookCut("= 2 lepton candidates");
            bookCut("mll > 20");
            bookCut("tau veto");
            bookCut("SF after preselection");
            bookCut("Jet veto: SF");
            bookCut("Z veto: SF");
            bookCut("mT2 > 90: SF");
            bookCut("mT2 > 120: SF");
            bookCut("mT2 > 150: SF");
            bookCut("pTll > 80: WWa: SF");
            bookCut("METrel > 80: WWa: SF");
            bookCut("mll < 120: WWa: SF");
            bookCut("mT2 > 90: WWb: SF");
            bookCut("mll < 170: WWb: SF");
            bookCut("mT2 > 100: WWc: SF");
            bookCut(">= 2 C_LFjets: Zjets: SF");
            bookCut("No b or F jets: Zjets: SF");
            bookCut("Z window: Zjets: SF");
            bookCut("pTll > 80: Zjets: SF");
            bookCut("METrel > 80: Zjets: SF");
            bookCut("0.3 < dRll < 1.5: Zjets: SF");
            bookCut("50 < mjj < 100: Zjets: SF");
            bookCut("pTj2 > 45: Zjets: SF");
            bookCut("DF after preselection");
            bookCut("Jet veto: DF");
            bookCut("Z veto: DF");
            bookCut("mT2 > 90: DF");
            bookCut("mT2 > 120: DF");
            bookCut("mT2 > 150: DF");
            bookCut("pTll > 80: WWa: DF");
            bookCut("METrel > 80: WWa: DF");
            bookCut("mll < 120: WWa: DF");
            bookCut("mT2 > 90: WWb: DF");
            bookCut("mll < 170: WWb: DF");
            bookCut("mT2 > 100: WWc: DF");

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& F_jets = applyProjection<NearIsoParticle>(event, "F_Jets").particlesByPt(&event);
            const Particles& leps = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt(&event);            
            const Particles& leps_candi = applyProjection<MergedFinalState>(event, "baseLeptons").particlesByPt(&event);                        
            const Particles& taus = applyProjection<NearIsoParticle>(event, "Taus").particlesByPt(&event);                        
            const Particles& jets = applyProjection<NearIsoParticle>(event, "C_Jets").particlesByPt(&event);                                    
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& C_bjets = bjproj.getTaggedJets(&event); 
            const Particles& C_LFjets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();

            //cout << jets.size() <<" | "<< C_bjets.size() <<" + "<< C_LFjets.size() << " | "<< jets.size() - C_bjets.size() - C_LFjets.size() << endl;

            // cout <<"F_jets: "<<  F_jets.size() << endl;
            // cout <<"leps: "<<  leps.size() << endl;
            // cout <<"taus: "<<  taus.size() << endl;
            // cout <<"C_bjets: "<<  C_bjets.size() << endl;
            // cout <<"C_LFjets: "<<  C_LFjets.size() << endl;
            // cout <<"met: "<<  MET << endl;

            // ---------------------------------------- //
            //              preselection                //
            // ---------------------------------------- //

            if(!cut(leps.size(), CUT_EQ, 2, "= 2 signal lepsons")) vetoEvent;           
            Particle lep1 = leps[0];
            Particle lep2 = leps[1];
            if(!cut(lep1.pdgId() * lep2.pdgId() < 0, "Opposite Charge")) vetoEvent;           
            if(!cut(lep1.pT(), CUT_GT, 35., "pTl1 > 35")) vetoEvent;           
            if(!cut(lep2.pT(), CUT_GT, 20., "pTl2 > 20")) vetoEvent;           
            if(!cut(leps_candi.size(), CUT_EQ, 2, "= 2 lepton candidates")) vetoEvent;           
            double mll = (lep1.momentum() + lep2.momentum()).mass();
            if(!cut( mll, CUT_GT, 20., "mll > 20")) vetoEvent;           
            if(!cut( taus.size(), CUT_EQ, 0, "tau veto")) vetoEvent;           

            // ---------------------------------------- //
            //           Variable Calculation           //
            // ---------------------------------------- //

            double pi = 4.*atan(1.);
            double mZ = 91.1876;                        

            bool SF = false; 
            bool DF = false;
            if( abs(lep1.pdgId()) == abs(lep2.pdgId()) ){ 
                SF = true;
            }else{
                DF = true;
            }

            double pTll = (lep1.momentum() + lep2.momentum()).pT();

            // Compute $E_T^{miss,rel}$       
            double dphi1 = deltaPhi(lep1.momentum(), met);
            double dphi2 = deltaPhi(lep2.momentum(), met);                        
            double dphi_min = min(dphi1, dphi2);
            for(size_t i=0; i<C_bjets.size(); i++){
                double dphi = deltaPhi(C_bjets[i].momentum(), met);                        
                if(dphi < dphi_min) dphi_min = dphi;
            }
            for(size_t i=0; i<C_LFjets.size(); i++){
                double dphi = deltaPhi(C_LFjets[i].momentum(), met);                        
                if(dphi < dphi_min) dphi_min = dphi;
            }
            double METrel = MET;
            if( dphi_min < pi/2. ){
                METrel = MET * sin(dphi_min);
            }

            // Compute $m_{T2}$
            double minv = 0.0;
            double mt2 = Rivet::mT2::mT2(lep1.momentum(), lep2.momentum(), met,  minv );

            // ---------------------------------------- //
            //              Same Flavour                //
            // ---------------------------------------- //

            if(cut(SF, "SF after preselection")){

                int Njet = C_LFjets.size() + C_bjets.size() + F_jets.size(); 
                if(cut(Njet, CUT_EQ, 0, "Jet veto: SF")){
                    if(cut(abs(mll - mZ), CUT_GT, 10., "Z veto: SF")){

                        //--- mT2 SRs ---//
                        if(cut(mt2, CUT_GT, 90.,  "mT2 > 90: SF"))  pass("mT2_90:SF");
                        if(cut(mt2, CUT_GT, 120., "mT2 > 120: SF")) pass("mT2_120:SF");
                        if(cut(mt2, CUT_GT, 150., "mT2 > 150: SF")) pass("mT2_150:SF");

                        //--- WW SRs ---//                        
                        // WWa
                        if(cut(pTll, CUT_GT, 80., "pTll > 80: WWa: SF")){
                            if(cut(METrel, CUT_GT, 80., "METrel > 80: WWa: SF")){
                                if(cut(mll, CUT_LT, 120., "mll < 120: WWa: SF")){
                                    pass("WWa:SF");
                                }
                            }
                        }
                        // WWb
                        if(cut(mt2, CUT_GT, 90., "mT2 > 90: WWb: SF")){
                            if(cut(mll, CUT_LT, 170., "mll < 170: WWb: SF")){
                                pass("WWb:SF");
                            }
                        }
                        // WWc
                        if(cut(mt2, CUT_GT, 100., "mT2 > 100: WWc: SF")){
                            pass("WWc:SF");
                        }
                    }
                }

                //--- Zjets SR --//
                if(cut(C_LFjets.size(), CUT_GT, 1, ">= 2 C_LFjets: Zjets: SF")){
                    int Njet = C_bjets.size() + F_jets.size(); 
                    if(cut(Njet, CUT_EQ, 0, "No b or F jets: Zjets: SF")){
                        if(cut(abs(mll - mZ), CUT_LT, 10., "Z window: Zjets: SF")){
                            if(cut(pTll, CUT_GT, 80., "pTll > 80: Zjets: SF")){
                                if(cut(METrel, CUT_GT, 80., "METrel > 80: Zjets: SF")){
                                    double dRll = deltaR(lep1.momentum(), lep2.momentum());
                                    if(cut(dRll, CUT_IN, make_pair(0.3, 1.5), "0.3 < dRll < 1.5: Zjets: SF")){
                                        double mjj = (C_LFjets[0].momentum() + C_LFjets[1].momentum()).mass();
                                        if(cut(mjj, CUT_IN, make_pair(50., 100.), "50 < mjj < 100: Zjets: SF")){
                                            if(cut(C_LFjets[1].pT(), CUT_GT, 45., "pTj2 > 45: Zjets: SF")){
                                                pass("Zjets");
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }

            // ---------------------------------------- //
            //         Different Flavour                //
            // ---------------------------------------- //

            if(cut(DF, "DF after preselection")){

                int Njet = C_LFjets.size() + C_bjets.size() + F_jets.size(); 
                if(cut(Njet, CUT_EQ, 0, "Jet veto: DF")){
                    if(cut(abs(mll - mZ), CUT_GT, 10., "Z veto: DF")){

                        //--- mT2 SRs ---//
                        if(cut(mt2, CUT_GT, 90.,  "mT2 > 90: DF"))  pass("mT2_90:DF");
                        if(cut(mt2, CUT_GT, 120., "mT2 > 120: DF")) pass("mT2_120:DF");
                        if(cut(mt2, CUT_GT, 150., "mT2 > 150: DF")) pass("mT2_150:DF");

                        //--- WW SRs ---//                        
                        // WWa
                        if(cut(pTll, CUT_GT, 80., "pTll > 80: WWa: DF")){
                            if(cut(METrel, CUT_GT, 80., "METrel > 80: WWa: DF")){
                                if(cut(mll, CUT_LT, 120., "mll < 120: WWa: DF")){
                                    pass("WWa:DF");
                                }
                            }
                        }
                        // WWb
                        if(cut(mt2, CUT_GT, 90., "mT2 > 90: WWb: DF")){
                            if(cut(mll, CUT_LT, 170., "mll < 170: WWb: DF")){
                                pass("WWb:DF");
                            }
                        }
                        // WWc
                        if(cut(mt2, CUT_GT, 100., "mT2 > 100: WWc: DF")){
                            pass("WWc:DF");
                        }
                    }
                }
                
            }


			/// @todo fill histograms 
			// fillPlot("Meff", meff);
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			/// @todo normalize the histograms
			// scale("Mjj");
		}

		//@}

		detector_parameter_t dp;


	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1403_5294)
}
