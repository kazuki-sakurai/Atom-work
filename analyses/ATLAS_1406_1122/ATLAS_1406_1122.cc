///
/// @file  ATLAS_1406_1122.cc
/// @brief Implementation of ATLAS_1406_1122 analysis
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

	class ATLAS_1406_1122 : public Analysis {
	public:

		ATLAS_1406_1122()
			: Analysis("ATLAS_1406_1122") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

        double get_mT(FourMomentum p1, FourMomentum p2){
            double mTsq = 2. * ( p1.pT()*p2.pT() - p1.px()*p2.px() - p1.py()*p2.py() );
            double mT = sqrt(mTsq);
            return mT;
        }

        pair<int, int> get_nearest_pair( Particles p ){

            double min_measure = 1000.;       
            int i1 = -1; int i2 = -1;
            for(int i=0; i<p.size()-1; i++){
                for(int j=i+1; j<p.size(); j++){
                    double dR = deltaR( p[i], p[j] );
                    if( dR < min_measure ){
                        min_measure = dR;
                        i1 = i;
                        i2 = j;
                    }
                }
            }

            return make_pair(i1, i2);            

        }

        pair<double, double> get_mbjjs( Particles jets, Particles bjets ){

            Particles jet_array = jets;
            int i_b1, i_b2;
            FourMomentum pb1 = bjets[0].momentum();
            FourMomentum pb2 = bjets[1].momentum();                                        
            for(int i=0; i<jet_array.size(); i++){
                if( jet_array[i].momentum() == pb1 ) i_b1 = i;
                if( jet_array[i].momentum() == pb2 ) i_b2 = i;
            }    
            jet_array.erase(jet_array.begin() + max(i_b1, i_b2));
            jet_array.erase(jet_array.begin() + min(i_b1, i_b2));

            pair<int, int> wcandi = get_nearest_pair(jet_array);
            int i_w1 = wcandi.first;
            int i_w2 = wcandi.second;

            FourMomentum pW1 = jet_array[i_w1].momentum() + jet_array[i_w2].momentum();
            jet_array.erase(jet_array.begin() + max(i_w1, i_w2));
            jet_array.erase(jet_array.begin() + min(i_w1, i_w2));

            double dR1 = deltaR( pW1, bjets[0].momentum() );
            double dR2 = deltaR( pW1, bjets[1].momentum() );
            double mbjj1 = 0.;            
            int ib_remain = -1;
            if( dR1 < dR2 ){
                mbjj1 = (pW1 + bjets[0].momentum()).mass();
                ib_remain = 1;
            }else{
                mbjj1 = (pW1 + bjets[1].momentum()).mass();
                ib_remain = 0;                                            
            }

            pair<int, int> wcandi2 = get_nearest_pair(jet_array);
            i_w1 = wcandi2.first;
            i_w2 = wcandi2.second;
            FourMomentum pW2 = jet_array[i_w1].momentum() + jet_array[i_w2].momentum();
            double mbjj2 = (pW2 + bjets[ib_remain].momentum()).mass();

            return make_pair(mbjj1, mbjj2);

        }

        double get_mT_min( Particles jets, Particles bjets, FourMomentum met ){

            double min_measure = 10000.;
            int i_min = -1;
            for(int i=0; i<jets.size(); i++){
                bool is_b = false;
                for(int j=0; j<bjets.size(); j++){
                    if( jets[i].momentum() == bjets[j].momentum() ) is_b = true;
                }
                if( !is_b ){
                    double dPhi = deltaPhi( jets[i].momentum(), met );
                    if( dPhi < min_measure ){ 
                        i_min = i;
                        min_measure = dPhi;
                    }
                }
            }
            double mT_min = 0.;
            if(i_min > 0) mT_min = get_mT(jets[i_min].momentum(), met);

            return mT_min;
        }

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            
            IsoElectron base_ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );

            IsoElectron ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            //                         cone  frac  abs  inner 
            ele.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            ele.addIso(CALO_ISO_ET,  0.3,  0.18,  0.0, 0.01, CALO_ALL);
            ele.setVariableThreshold(0.0);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon base_mu(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            base_mu.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            base_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon mu( Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4) );
            //                         cone  frac  abs  inner 
            mu.addIso(TRACK_ISO_PT, 0.3,  0.12,  0.0, 0.01, CALO_ALL);
            mu.addIso(CALO_ISO_ET,  0.3,  0.12,  0.0, 0.01, CALO_ALL);
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

            MergedFinalState base_leptons(base_ele, base_mu);
            NearIsoParticle base_leptons_clean(base_leptons);
            base_leptons_clean.addFilter(base_jets_clean, 0.4);
            addProjection(base_leptons_clean, "Leptons");

            MergedFinalState leptons(ele, mu);
            NearIsoParticle leptons_clean(leptons);
            leptons_clean.addFilter(base_jets_clean, 0.4);
            addProjection(leptons_clean, "CR_Leptons");

            NearIsoParticle jets_clean(base_jets_clean, Range(PT, 35., 8000.) & Range(ETA, -2.8, 2.8));
            addProjection(jets_clean, "Jets");

            FastJets megajets08(jets_clean, FastJets::ANTIKT, 0.8 );
            addProjection(megajets08, "MegaJets08");

            FastJets megajets12(jets_clean, FastJets::ANTIKT, 1.2 );
            addProjection(megajets12, "MegaJets12");

            Range bjrange = Range(PT, 35.0, 8000.0) & Range(ETA, -2.5, 2.5);
            HeavyFlavorJets bjets(jets_clean, bjrange);
            addProjection(bjets, "BJets");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets_clean, base_leptons_clean);            
            MissingMomentum met( fsbase, met_seed );            
            MissingMomentum met( fsbase );                        
            met.setSmearingParams( FinalState::SELECTED, &dp.metEff( "Jet_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            ChargedFinalState tracks = ChargedFinalState(Range(PT, 0.5, 8000.) & Range(ETA, -2.5, 2.5));
            addProjection(tracks, "Tracks");

            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

            /// @todo book the efficiencies
            bookEfficiency("SRA1");
            bookEfficiency("SRA2");
            bookEfficiency("SRA3");
            bookEfficiency("SRA4");
            bookEfficiency("SRB1");
            bookEfficiency("SRB2");
            bookEfficiency("SRC1");
            bookEfficiency("SRC2");
            bookEfficiency("SRC3");

            /// @todo book the cuts
            bookCut("Lepton Veto");
            bookCut("MET > 150");
            bookCut("Nj >= 6: SRA");
            bookCut("pTj2 > 80: SRA");
            bookCut("dPhi_jet_met_123 > pi/5: SRA");
            bookCut("dPhi_MET_trackMET < pi/3: SRA");
            bookCut("Nb >= 2: SRA");
            bookCut("tau veto: SRA");
            bookCut("mTb_min > 175: SRA");
            bookCut("mbjj0 < 225: SRA(1-2)");
            bookCut("mbjj1 < 250: SRA(1-2)");
            bookCut("MET > 150: SRA1");
            bookCut("MET > 250: SRA2");
            bookCut("50 < mbjj0 < 250: SRA(3-4)");
            bookCut("50 < mbjj1 < 400: SRA(3-4)");
            bookCut("50 < mbjj1 < 400: SRA(3-4)");
            bookCut("mTj_min > 50: SRA(3-4)");
            bookCut("MET > 300: SRA3");
            bookCut("MET > 350: SRA4");
            bookCut("Nj = 4 or 5: SRB");
            bookCut("pTj2 > 80: SRB");
            bookCut("dPhi_jet_met_123 > pi/5: SRB");
            bookCut("dPhi_MET_trackMET < pi/3: SRB");
            bookCut("Nb >= 2: SRB");
            bookCut("mTb_min > 175: SRB");
            bookCut("mJ_R12_0 > 80: SRB");
            bookCut("mJ_R08_0 > 50: SRB");
            bookCut(">= 2 J_R12: SRB");
            bookCut("Amt < 0.5: SRB1");
            bookCut("60 < mJ_R12_1 < 200: SRB1");
            bookCut("mT_min > 175: SRB1");
            bookCut("mT_j3_min > 280 for Nj=4: SRB1");
            bookCut("MET > 325: SRB1");
            bookCut("Nj = 5: SRB2");
            bookCut("Amt > 0.5: SRB2");
            bookCut("pT(J_R12_0) > 350: SRB2");
            bookCut("140 < mJ_R12_0 < 500: SRB2");
            bookCut("70 < mJ_R08_0 < 300: SRB2");
            bookCut("mT_min > 125: SRB2");
            bookCut("MET/sqrt(HT) > 17: SRB2");
            bookCut("MET > 400: SRB2");
            bookCut("Nj = 5: SRC");
            bookCut("pTj2 > 80: SRC");
            bookCut("dPhi_jet_met_123 > pi/5: SRC");
            bookCut("dPhi_MET_trackMET < pi/3: SRC");
            bookCut("Nb >= 2: SRC");
            bookCut("tau veto: SRC");
            bookCut("dPhi_bb > 0.2Pi: SRC");
            bookCut("mTb_min > 185: SRC1");
            bookCut("mTb_max > 205: SRC1");
            bookCut("MET > 160: SRC1");
            bookCut("mTb_min > 200: SRC2");
            bookCut("mTb_max > 290: SRC2");
            bookCut("MET > 160: SRC2");
            bookCut("mTb_min > 200: SRC3");
            bookCut("mTb_max > 325: SRC3");
            bookCut("MET > 215: SRC3");

            //bookCut("CutEtaJet","this is a control region cut", true);

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt(&event);
            const Particles& megajets08 = applyProjection<FastJets>(event, "MegaJets08").particlesByPt(&event);
            const Particles& megajets12 = applyProjection<FastJets>(event, "MegaJets12").particlesByPt(&event);
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt(&event);
            //const Particles& CR_leps = applyProjection<NearIsoParticle>(event, "CR_Leptons").particlesByPt(&event);
            const Particles& tracks = applyProjection<ChargedFinalState>(event, "Tracks").particlesByPt(&event);                  
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(&event); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");

            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();
            double pi = 4.*atan(1.);

            //=============================//
            //         preselecton         //
            //=============================//

            if(!cut(leps.size(), CUT_EQ, 0, "Lepton Veto")) vetoEvent;
            if(!cut(MET, CUT_GT, 150., "MET > 150")) vetoEvent;

            //=============================//
            //    Variable calculation     //
            //=============================//

            FourMomentum metTracks(0,0,0,0);
            for(int i=0; i<tracks.size(); i++){
                metTracks -= tracks[i].momentum();
            }
            double dPhi_met_metTrack = deltaPhi(met, metTracks);

            double dPhi_jet_met_min123 = -100.;
            if(jets.size() >= 3){
                dPhi_jet_met_min123 = 100.;
                for(int i=0; i<3; i++){
                    double dPhi = deltaPhi( jets[i].momentum(), met );
                    if( dPhi < dPhi_jet_met_min123 ) dPhi_jet_met_min123 = dPhi;
                }
            }

            int ib_min = -1;
            int ib_max = -1;                        
            double min_measure =  1000.;
            double max_measure = -1000.;            
            for(int i=0; i<bjets.size(); i++){
                double dPhi = deltaPhi( bjets[i].momentum(), met );
                if(dPhi < min_measure){
                    min_measure = dPhi;
                    ib_min = i;
                }
                if(dPhi > max_measure){
                    max_measure = dPhi;
                    ib_max = i;
                }                
            }

            double mTbmin = 0.;
            double mTbmax = 0.;            
            if(bjets.size() > 0){
                mTbmin = get_mT( bjets[ib_min], met );
                mTbmax = get_mT( bjets[ib_max], met );                
            }

            //=============================//
            //             SRA             //
            //=============================//

            int N_tau_candi = 0; // Define!!

            if(cut( jets.size(), CUT_GE, 6, "Nj >= 6: SRA" )){
                if(cut( jets[1].pT(), CUT_GT, 80., "pTj2 > 80: SRA" )){
                    if(cut( dPhi_jet_met_min123, CUT_GT, pi/5., "dPhi_jet_met_123 > pi/5: SRA" )){
                        if(cut( dPhi_met_metTrack, CUT_LT, pi/3., "dPhi_MET_trackMET < pi/3: SRA" )){
                            if(cut( bjets.size(), CUT_GE, 2, "Nb >= 2: SRA" )){
                                if(cut( N_tau_candi, CUT_EQ, 0, "tau veto: SRA" )){
                                    if(cut( mTbmin, CUT_GT, 175., "mTb_min > 175: SRA" )){
                                        pair<double, double> mbjjs = get_mbjjs(jets, bjets);
                                        double mbjj0 = mbjjs.first;
                                        double mbjj1 = mbjjs.second;
                                        if(cut( mbjj0, CUT_LT, 225., "mbjj0 < 225: SRA(1-2)" )){
                                            if(cut( mbjj1, CUT_LT, 250., "mbjj1 < 250: SRA(1-2)" )){
                                                if(cut( MET, CUT_GT, 150., "MET > 150: SRA1" )) pass("SRA1");
                                                if(cut( MET, CUT_GT, 250., "MET > 250: SRA2" )) pass("SRA2");
                                            }
                                        }
                                        if(cut( mbjj0, CUT_IN, make_pair(50., 250.), "50 < mbjj0 < 250: SRA(3-4)" )){
                                            if(cut( mbjj1, CUT_IN, make_pair(50., 400.), "50 < mbjj1 < 400: SRA(3-4)" )){
                                                if(cut( mbjj1, CUT_IN, make_pair(50., 400.), "50 < mbjj1 < 400: SRA(3-4)" )){
                                                    double mTj_min = 1000.;
                                                    for(int i=0; i<jets.size(); i++){
                                                        double mT = get_mT( jets[i].momentum(), met );
                                                        if( mT < mTj_min ) mTj_min = mT;
                                                    }
                                                    if(cut( mTj_min, CUT_GT, 50., "mTj_min > 50: SRA(3-4)" )){
                                                        if(cut( MET, CUT_GT, 300., "MET > 300: SRA3" )) pass("SRA3");
                                                        if(cut( MET, CUT_GT, 350., "MET > 350: SRA4" )) pass("SRA4");
                                                    }

                                                }
                                            }
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }

            //=============================//
            //             SRB             //
            //=============================//

            if(cut( jets.size(), CUT_IN, make_pair(4, 5), "Nj = 4 or 5: SRB" )){
                if(cut( jets[1].pT(), CUT_GT, 80., "pTj2 > 80: SRB" )){
                    if(cut( dPhi_jet_met_min123, CUT_GT, pi/5., "dPhi_jet_met_123 > pi/5: SRB" )){
                        if(cut( dPhi_met_metTrack, CUT_LT, pi/3., "dPhi_MET_trackMET < pi/3: SRB" )){
                            if(cut( bjets.size(), CUT_GE, 2, "Nb >= 2: SRB" )){
                                if(cut( mTbmin, CUT_GT, 175., "mTb_min > 175: SRB" )){

                                    double mJ_R12_0 = 0;
                                    if(megajets12.size() > 0) mJ_R12_0 = megajets12[0].mass();
                                    double mJ_R08_0 = 0;
                                    if(megajets08.size() > 0) mJ_R08_0 = megajets08[0].mass();
                                    if(cut( mJ_R12_0, CUT_GT, 80., "mJ_R12_0 > 80: SRB" )){
                                        if(cut( mJ_R08_0, CUT_GT, 50., "mJ_R08_0 > 50: SRB" )){
                                            if(cut( megajets12.size(), CUT_GE, 2, ">= 2 J_R12: SRB" )){

                                                double mJ_R12_1 = megajets12[1].mass();
                                                double Amt = abs(mJ_R12_0 - mJ_R12_1)/(mJ_R12_0 + mJ_R12_1); 
                                                double mT_j3_met = 0;
                                                if(jets.size() == 4) mT_j3_met = get_mT( jets[3].momentum(), met );
                                                double mT_min = get_mT_min(jets, bjets, met); 

                                                // SRB1
                                                if(cut( Amt, CUT_LT, 0.5, "Amt < 0.5: SRB1" )){
                                                    if(cut( mJ_R12_1, CUT_IN, make_pair(60., 200.), "60 < mJ_R12_1 < 200: SRB1" )){
                                                        if(cut( mT_min, CUT_GT, 175., "mT_min > 175: SRB1" )){
                                                            if(cut( mT_j3_met, CUT_GT, 280., "mT_j3_min > 280 for Nj=4: SRB1" )){
                                                                if(cut( MET, CUT_GT, 325., "MET > 325: SRB1" )){
                                                                    pass("SRB1");
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                // SRB2
                                                if(cut( jets.size(), CUT_EQ, 5, "Nj = 5: SRB2" )){
                                                    if(cut( Amt, CUT_GT, 0.5, "Amt > 0.5: SRB2" )){
                                                        if(cut( megajets12[0].pT(), CUT_GT, 350., "pT(J_R12_0) > 350: SRB2" )){
                                                            if(cut( mJ_R12_0, CUT_IN, make_pair(140., 500.), "140 < mJ_R12_0 < 500: SRB2" )){
                                                                if(cut( mJ_R08_0, CUT_IN, make_pair(70., 300.), "70 < mJ_R08_0 < 300: SRB2" )){
                                                                    
                                                                    if(cut( mT_min, CUT_GT, 125., "mT_min > 125: SRB2" )){
                                                                        double HT = 0;
                                                                        for(int i=0; i<jets.size(); i++) HT += jets[i].pT();
                                                                        if(cut( MET/sqrt(HT), CUT_GT, 17., "MET/sqrt(HT) > 17: SRB2" )){
                                                                            if(cut( MET, CUT_GT, 400., "MET > 400: SRB2" )){
                                                                                pass("SRB2");
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //=============================//
            //             SRC             //
            //=============================//

            if(cut( jets.size(), CUT_EQ, 5, "Nj = 5: SRC" )){
                if(cut( jets[1].pT(), CUT_GT, 80., "pTj2 > 80: SRC" )){
                    if(cut( dPhi_jet_met_min123, CUT_GT, pi/5., "dPhi_jet_met_123 > pi/5: SRC" )){
                        if(cut( dPhi_met_metTrack, CUT_LT, pi/3., "dPhi_MET_trackMET < pi/3: SRC" )){
                            if(cut( bjets.size(), CUT_GE, 2, "Nb >= 2: SRC" )){
                                if(cut( N_tau_candi, CUT_EQ, 0, "tau veto: SRC" )){

                                    double dPhi_bb = deltaPhi( bjets[0], bjets[1] );
                                    if(cut( dPhi_bb, CUT_GT, 0.2*pi, "dPhi_bb > 0.2Pi: SRC" )){

                                        //SRC1
                                        if(cut( mTbmin, CUT_GT, 185., "mTb_min > 185: SRC1")){
                                            if(cut( mTbmax, CUT_GT, 205., "mTb_max > 205: SRC1")){
                                                if(cut( MET, CUT_GT, 160., "MET > 160: SRC1")){
                                                    pass("SRC1");
                                                }
                                            }
                                        }

                                        //SRC2
                                        if(cut( mTbmin, CUT_GT, 200., "mTb_min > 200: SRC2")){
                                            if(cut( mTbmax, CUT_GT, 290., "mTb_max > 290: SRC2")){
                                                if(cut( MET, CUT_GT, 160., "MET > 160: SRC2")){
                                                    pass("SRC2");
                                                }
                                            }
                                        }

                                        //SRC1
                                        if(cut( mTbmin, CUT_GT, 200., "mTb_min > 200: SRC3")){
                                            if(cut( mTbmax, CUT_GT, 325., "mTb_max > 325: SRC3")){
                                                if(cut( MET, CUT_GT, 215., "MET > 215: SRC3")){
                                                    pass("SRC3");
                                                }
                                            }
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ********************************************* //

            // if( !cut( leps.size(), CUT_EQ, 1, "Preselection: =1 lepton" )){
            //     vetoEvent;
            // }

		}

		/// Normalise histograms etc., after the run
		void finalize() {
			/// @todo normalize the histograms
			// scale("Mjj");
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1406_1122)
}
