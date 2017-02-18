///
/// @file  ATLAS_1403_4853.cc
/// @brief Implementation of ATLAS_1403_4853 analysis
/// @author Kazuki
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

	class ATLAS_1403_4853 : public Analysis {
	public:

		ATLAS_1403_4853()
			: Analysis("ATLAS_1403_4853") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            IsoElectron base_ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoElectron ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            ele.addIso(TRACK_ISO_PT, 0.2,  0.1,  0.0, 0.01, CALO_ALL);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon mu(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            mu.addIso(TRACK_ISO_PT, 0.2,  0.0,  1.8, 0.01, CALO_ALL);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -2.5, 2.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // Overlap removal

            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(base_ele, 0.2);
            addProjection(jets_clean, "Jets");

            MergedFinalState leptons(ele, mu);
            NearIsoParticle leptons_clean(leptons);
            leptons_clean.addFilter(jets_clean, 0.4);
            addProjection(leptons_clean, "Leptons");

            Range bjrange = Range(PT, 20.0, 8000.0) & Range(ETA, -2.5, 2.5);
            HeavyFlavorJets bjets(jets_clean, bjrange);
            addProjection(bjets, "BJets");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets_clean, leptons_clean);            
            MissingMomentum met( fsbase, met_seed );            
            //MissingMomentum met( fsbase );                        
            //met.setSmearingParams( FinalState::SELECTED, &dp.metEff( "Jet_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            /// @todo book the efficiencies
            bookEfficiency("H160");
            bookEfficiency("L90");
            bookEfficiency("L100");
            bookEfficiency("L110");
            bookEfficiency("L120");

            /// @todo book the cuts
            bookCut("= 2 leptons");
            bookCut("= 2 leptons: SF");
            bookCut("OS lepton: SF");
            bookCut("m_ll > 20: SF");
            bookCut("pT(lep1) > 25: SF");
            bookCut("pT(lep1) > 25: SF");
            bookCut("= 2 leptons: DF");
            bookCut("OS lepton: DF");
            bookCut("m_ll > 20: DF");
            bookCut("pT(lep1) > 25: DF");
            bookCut("pT(lep1) > 25: DF");

            bookCut("H160: =2 b-jets: SF");
            bookCut("H160: mT2(b-jet) > 160: SF");
            bookCut("H160: mT2 < 90: SF");
            bookCut("H160: pT(lep1) < 60: SF");
            bookCut("H160: =2 b-jets: DF");
            bookCut("H160: mT2(b-jet) > 160: DF");
            bookCut("H160: mT2 < 90: DF");
            bookCut("H160: pT(lep1) < 60: DF");

            bookCut("Z veto: SF");
            bookCut("Dphi_j > 1.0: SF");
            bookCut("Dphi_b < 1.5: SF");
            bookCut("mT2 > 90: SF");
            bookCut("mT2 > 120: SF");
            bookCut("L100: mT2 > 100: SF");
            bookCut("L100: pT(j1) > 100: SF");
            bookCut("L100: pT(j2) > 50: SF");
            bookCut("L110: mT2 > 110: SF");
            bookCut("L110: pT(j1) > 20: SF");
            bookCut("L110: pT(j2) > 20: SF");
            bookCut("Dphi_j > 1.0: DF");
            bookCut("Dphi_b < 1.5: DF");
            bookCut("mT2 > 90: DF");
            bookCut("mT2 > 120: DF");
            bookCut("L100: mT2 > 100: DF");
            bookCut("L100: pT(j1) > 100: DF");
            bookCut("L100: pT(j2) > 50: DF");
            bookCut("L110: mT2 > 110: DF");
            bookCut("L110: pT(j1) > 20: DF");
            bookCut("L110: pT(j2) > 20: DF");
           
            //bookCut("CutEtaJet","this is a control region cut", true);

            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            /// @todo apply projections
			// const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt(&event);
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt(&event);
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(&event); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

            double MET = met.pT();
            // ********************************************* //

            //=============================//
            //         preselecton         //
            //=============================//

            if( !cut( leps.size(), CUT_EQ, 2, "= 2 leptons" ) ) vetoEvent;
            Particle lep1 = leps[0];
            Particle lep2 = leps[1];
            bool SF = false; bool DF = false;
            if(cut( abs( abs(lep1.pdgId()) - abs(lep2.pdgId()) ), CUT_EQ, 0, "= 2 leptons: SF" )) SF = true;  
            if(cut( abs( abs(lep1.pdgId()) - abs(lep2.pdgId()) ), CUT_GT, 1, "= 2 leptons: DF" )) DF = true;  

            // SF preselection
            bool pass_preselec_SF = false;
            if(SF){
                if( cut(lep1.pdgId() * lep2.pdgId(), CUT_LT, 0, "OS lepton: SF" ) ){
                    double mll = (lep1.momentum() + lep2.momentum()).mass();
                    if( cut(mll, CUT_GT, 20., "m_ll > 20: SF" ) ){
                        if( cut(lep1.pT(), CUT_GT, 25., "pT(lep1) > 25: SF" ) ){
                            if( cut(lep1.pT(), CUT_GT, 25., "pT(lep1) > 25: SF" ) ){
                                pass_preselec_SF = true;
                            }
                        }
                    }
                }
            }

            // DF preselection
            bool pass_preselec_DF = false;
            if(DF){
                    if( cut(lep1.pdgId() * lep2.pdgId(), CUT_LT, 0, "OS lepton: DF" ) ){
                        double mll = (lep1.momentum() + lep2.momentum()).mass();
                        if( cut(mll, CUT_GT, 20., "m_ll > 20: DF" ) ){
                            if( cut(lep1.pT(), CUT_GT, 25., "pT(lep1) > 25: DF" ) ){
                                if( cut(lep1.pT(), CUT_GT, 25., "pT(lep1) > 25: DF" ) ){
                                    pass_preselec_SF = true;
                                }
                            }
                        }
                }
            }

            if( pass_preselec_SF == false && pass_preselec_DF == false ) vetoEvent;

            //=============================//
            //     Calculate Variables     // 
            //=============================//

            double mll = (lep1.momentum() + lep2.momentum()).mass();

            double Dphi_jet = 100000.;
            for(unsigned int i=0; i<jets.size(); i++){
                double dphi = deltaPhi(jets[i].momentum(), met);
                if(dphi < Dphi_jet) Dphi_jet = dphi;
            }        

            FourMomentum pdm = met + lep1.momentum() + lep2.momentum();
            FourMomentum pllTb = FourMomentum(pdm.pT(), pdm.px(), pdm.py(), 0.0); 
            double Dphi_b = deltaPhi( pllTb, met );

            ///compute mT2
            double minv = 0.0;
            double mt2 = Rivet::mT2::mT2(lep1.momentum(), lep2.momentum(), met,  minv);

            //=============================//
            //      Hadronic mT2 SRs       //
            //=============================//

            // SF
            bool pass_hadronic_SF = false;            
            if(SF){     
                if(cut(bjets.size(), CUT_EQ, 2, "H160: =2 b-jets: SF")){
                    ///compute mT2b
                    double minv = 0.0;
                    double mt2b = Rivet::mT2::mT2(bjets[0].momentum(), bjets[1].momentum(), pllTb,  minv);
                    if(cut(mt2b, CUT_GT, 160., "H160: mT2(b-jet) > 160: SF")){
                        if(cut(mt2, CUT_LT, 90., "H160: mT2 < 90: SF")){
                            if(cut(lep1.pT(), CUT_LT, 60., "H160: pT(lep1) < 60: SF")){
                                pass_hadronic_SF = true;
                            }
                        }
                    }
                }
            }

            // DF
            bool pass_hadronic_DF = false;            
            if(DF){     
                if(cut(bjets.size(), CUT_EQ, 2, "H160: =2 b-jets: DF")){
                    ///compute mT2b
                    double minv = 0.0;
                    double mt2b = Rivet::mT2::mT2(bjets[0].momentum(), bjets[1].momentum(), pllTb,  minv);
                    if(cut(mt2b, CUT_GT, 160., "H160: mT2(b-jet) > 160: DF")){
                        if(cut(mt2, CUT_LT, 90., "H160: mT2 < 90: DF")){
                            if(cut(lep1.pT(), CUT_LT, 60., "H160: pT(lep1) < 60: DF")){
                                pass_hadronic_DF = true;
                            }
                        }
                    }
                }
            }

            if( pass_hadronic_SF || pass_hadronic_DF ) pass("H160");

            //=============================//
            //      Leeptonic mT2 SRs      //
            //=============================//

            // SF
            bool pass_L90_SF = false;
            bool pass_L100_SF = false;
            bool pass_L110_SF = false;
            bool pass_L120_SF = false;
            if(SF){
                if( cut(mll, CUT_OUT, make_pair(71., 111.), "Z veto: SF") ){
                    if( cut(Dphi_jet, CUT_GT, 1., "Dphi_j > 1.0: SF") ) {
                        if( cut(Dphi_b, CUT_LT, 1.5, "Dphi_b < 1.5: SF") ){

                            if( cut(mt2, CUT_GT, 90., "mT2 > 90: SF") ){
                                pass_L90_SF = true;
                            } 
                            if( cut(mt2, CUT_GT, 120., "mT2 > 120: SF") ){
                                pass_L120_SF = true;
                            } 

                            if(cut(mt2, CUT_GT, 100., "L100: mT2 > 100: SF")){
                                if(cut(jets[0].pT(), CUT_GT, 100., "L100: pT(j1) > 100: SF")){
                                    if(cut(jets[1].pT(), CUT_GT, 50., "L100: pT(j2) > 50: SF")){
                                        pass_L100_SF = true;
                                    }
                                }
                            }
                            if(cut(mt2, CUT_GT, 110., "L110: mT2 > 110: SF")){
                                if(cut(jets[0].pT(), CUT_GT, 20., "L110: pT(j1) > 20: SF")){
                                    if(cut(jets[1].pT(), CUT_GT, 20., "L110: pT(j2) > 20: SF")){
                                        pass_L110_SF = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // DF
            bool pass_L90_DF = false;
            bool pass_L100_DF = false;
            bool pass_L110_DF = false;
            bool pass_L120_DF = false;
            if(DF){
                if( cut(Dphi_jet, CUT_GT, 1., "Dphi_j > 1.0: DF") ) {
                    if( cut(Dphi_b, CUT_LT, 1.5, "Dphi_b < 1.5: DF") ){

                        if( cut(mt2, CUT_GT, 90., "mT2 > 90: DF") ){
                            pass_L90_DF = true;
                        } 
                        if( cut(mt2, CUT_GT, 120., "mT2 > 120: DF") ){
                            pass_L120_DF = true;
                        } 

                        if(cut(mt2, CUT_GT, 100., "L100: mT2 > 100: DF")){
                            if(cut(jets[0].pT(), CUT_GT, 100., "L100: pT(j1) > 100: DF")){
                                if(cut(jets[1].pT(), CUT_GT, 50., "L100: pT(j2) > 50: DF")){
                                    pass_L100_DF = true;
                                }
                            }
                        }
                        if(cut(mt2, CUT_GT, 110., "L110: mT2 > 110: DF")){
                            if(cut(jets[0].pT(), CUT_GT, 20., "L110: pT(j1) > 20: DF")){
                                if(cut(jets[1].pT(), CUT_GT, 20., "L110: pT(j2) > 20: DF")){
                                    pass_L110_DF = true;
                                }
                            }
                        }
                    }
                }
            }

            if( pass_L90_SF || pass_L90_DF ) pass("L90");
            if( pass_L100_SF || pass_L100_DF ) pass("L100");
            if( pass_L110_SF || pass_L110_DF ) pass("L110");
            if( pass_L120_SF || pass_L120_DF ) pass("L120");

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
	AtomPlugin(ATLAS_1403_4853)
}
