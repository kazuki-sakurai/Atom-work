///
/// @file  ATLAS_1407_0583.cc
/// @brief Implementation of ATLAS_1407_0583 analysis
/// @author Kazuki
/// @date created 04/26/2015
/// @date last revision 04/26/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

#include "Rivet/Tools/RivetMT2.hh"
#include "Atom/Tools/MT2bl.h"


using namespace std;

namespace Atom {

	class ATLAS_1407_0583 : public Analysis {
	public:

		ATLAS_1407_0583()
			: Analysis("ATLAS_1407_0583") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

        double get_mT(FourMomentum p1, FourMomentum p2){
            double mTsq = 2. * ( p1.pT()*p2.pT() - p1.px()*p2.px() - p1.py()*p2.py() );
            double mT = sqrt(mTsq);
            return mT;
        }

        double AsymMt2(Particle lep, Particle b1, Particle b2, FourMomentum met){
            FourMomentum qlep = lep.momentum();
            FourMomentum qb1 = b1.momentum();
            FourMomentum qb2 = b2.momentum();
            double pl[4]  = {qlep.E(), qlep.px(), qlep.py(), qlep.pz()};  // El, plx, ply, plz,     (visible lepton)
            double pb1[4] = { qb1.E(),  qb1.px(),  qb1.py(),  qb1.pz()};  // Eb1, pb1x, pb1y, pb1z  (bottom on the same side as the visible lepton)
            double pb2[4] = { qb2.E(),  qb2.px(),  qb2.py(),  qb2.pz()};  // Eb2, pb2x, pb2y, pb2z  (other bottom, paired with the invisible W)
            double pmiss[3] = { 0., met.px(), met.py()};                  // <unused>, pmx, pmy     (missing pT)
            mt2bl_bisect::mt2bl mt2bl;
            mt2bl.set_momenta(pl, pb1, pb2, pmiss);
            return mt2bl.get_mt2bl();
        }

        double get_amT2(Particle lep, Particle b1, Particle b2, FourMomentum met){
            double amt2_1 = AsymMt2(lep, b1, b2, met);
            double amt2_2 = AsymMt2(lep, b2, b1, met);
            return min(amt2_1, amt2_2);
        }

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets base_jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -4.5, 4.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            base_jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            base_jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            IsoElectron base_ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );

            IsoElectron ele( Range(PT, 25., 8000.) & Range(ETA, -2.47, 2.47) );
            //                         cone  frac  abs  inner 
            ele.addIso(TRACK_ISO_PT, 0.2,  0.10,  0.0, 0.01, CALO_ALL);
            ele.setVariableThreshold(0.0);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoElectron soft_ele( Range(PT, 7., 8000.) & Range(ETA, -2.47, 2.47) );
            //                         cone  frac  abs  inner 
            soft_ele.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            soft_ele.setVariableThreshold(0.0);
            soft_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            soft_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoMuon base_mu(Range(PT, 10., 8000.) & Range(ETA, -2.4, 2.4));
            base_mu.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            base_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon mu( Range(PT, 25., 8000.) & Range(ETA, -2.4, 2.4) );
            //                         cone  frac  abs  inner 
            mu.addIso(TRACK_ISO_PT, 0.2,  0.0,  1.8, 0.01, CALO_ALL);
            mu.setVariableThreshold(0.0);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon mu_R( Range(PT, 25., 8000.) & Range(ETA, -2.4, 2.4) );
            //                         cone  frac  abs  inner 
            mu_R.addIso(TRACK_ISO_PT, 0.2,  0.12,  0.0, 0.01, CALO_ALL);
            mu_R.setVariableThreshold(0.0);
            mu_R.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu_R.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon soft_mu( Range(PT, 6., 8000.) & Range(ETA, -2.4, 2.4) );
            //                         cone  frac  abs  inner 
            soft_mu.addIso(TRACK_ISO_PT, 0.3,  0.12,  0.0, 0.01, CALO_ALL);
            soft_mu.setVariableThreshold(0.0);
            soft_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            soft_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );


            // Overlap removal

            HeavyFlavorJets bjets(base_jets, getBJetEff("BJet_Ident_MV1_ATLAS"), Range(ETA, -2.5, 2.5));
            addProjection(bjets, "BJets");

            NearIsoParticle base_ele_1(base_ele);
            base_ele_1.addFilter(bjets, 0.2);

            NearIsoParticle jets_clean(base_jets);
            jets_clean.addFilter(base_ele_1, 0.2);

            MergedFinalState base_leps(base_ele_1, base_mu);

            NearIsoParticle base_leps_clean(base_leps);
            base_leps_clean.addFilter(jets_clean, 0.4);
            addProjection(base_leps_clean, "baseLeptons");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets_clean, base_leps_clean);            
            MissingMomentum met( fsbase, met_seed );            
            //MissingMomentum met( fsbase );                        
            //met.setSmearingParams( FinalState::SELECTED, &dp.metEff( "Jet_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            FastJets tau_seeds(fsbase, 
                            hadRange & Range(PT, 15., 8000.) & Range(ETA, -2.47, 2.47), 
                            muDetRange, FastJets::ANTIKT, 0.4 );

            const double dRprong = 0.2;
            ParamTauFinder taus(tau_seeds, dRprong);

            NearIsoParticle taus_clean(taus);
            taus_clean.addFilter(base_leps_clean, 0.2);
            addProjection(taus_clean, "Taus");

            NearIsoParticle jets(jets_clean, Range(ETA, -2.5, 2.5));
            addProjection(jets, "Jets");

            NearIsoParticle jets25(jets, Range(PT, 25., 8000.));
            addProjection(jets25, "Jets25");

            MergedFinalState leptons(ele, mu);
            NearIsoParticle leps_clean(leptons);
            leps_clean.addFilter(jets_clean, 0.4);
            addProjection(leps_clean, "Leptons");

            MergedFinalState leps_R(ele, mu_R);
            NearIsoParticle leps_R_clean(leps_R);
            leps_R_clean.addFilter(jets_clean, 0.4);
            addProjection(leps_R_clean, "Leptons_R");

            MergedFinalState soft_leps(soft_ele, soft_mu);
            NearIsoParticle soft_leps_clean(soft_leps);
            soft_leps_clean.addFilter(jets_clean, 0.4);
            addProjection(soft_leps_clean, "softLeptons");


            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

            /// @todo book the efficiencies
            bookEfficiency("bCa_med");
            bookEfficiency("bCb_med1");
            bookEfficiency("bCa_low");
            bookEfficiency("bCb_high");
            bookEfficiency("bCc_diag");
            bookEfficiency("bCb_med2_<250_<120");
            bookEfficiency("bCb_med2_<250_>120");
            bookEfficiency("bCb_med2_>250_<120");
            bookEfficiency("bCb_med2_>250_>120");
            bookEfficiency("bCd_bulk_<250_<120");
            bookEfficiency("bCd_bulk_<250_>120");
            bookEfficiency("bCd_bulk_>250_<120");
            bookEfficiency("bCd_bulk_>250_>120");
            bookEfficiency("bCd_high1");
            bookEfficiency("bCd_high2");
            bookEfficiency("3body_<90_<120");
            bookEfficiency("3body_<90_>120");
            bookEfficiency("3body_>90_<120");
            bookEfficiency("3body_>90_>120");            

            /// @todo book the cuts
            bookCut("CutEtaJet");
            bookCut("MET > 300: bCa_med");
            bookCut("= 1 base lepton: bCa_med");
            bookCut("= 1 soft lepton: bCa_med");
            bookCut("pTl < 50: bCa_med");
            bookCut(">= 3 jets (pT>25): bCa_med");
            bookCut("pTj1 > 180: bCa_med");
            bookCut(">= 1 bjets: bCa_med");
            bookCut("j1 is not b-tagged: bCa_med");
            bookCut("mT > 100: bCa_med");
            bookCut("MET/meff > 0.3: bCa_med");
            bookCut("pTl1 < 25: bCa_med");
            bookCut("MET > 150: bCb_med1");
            bookCut("= 1 base lepton: bCb_med1");
            bookCut("= 1 soft lepton: bCb_med1");
            bookCut("pTl < 25: bCb_med1");
            bookCut(">= 2 jets: bCb_med1");
            bookCut("pTj2 > 60: bCb_med1");
            bookCut("j1 and j2 are b-tagged: bCb_med1");
            bookCut("dPhi_j_met_12 > 0.4: bCb_med1");
            bookCut("HT_2 < 50: bCb_med1");
            bookCut("mbb > 150: bCb_med1");
            bookCut("amT2 > 170: bCb_med1");
            bookCut("MET > 370: bCa_low");
            bookCut("= 1 base lepton: bCa_low");
            bookCut("= 1 soft lepton: bCa_low");
            bookCut("pTl < 50: bCa_low");
            bookCut(">= 2 jets (pT>25): bCa_low");
            bookCut("pTj1 > 180: bCa_low");
            bookCut(">= 1 bjets: bCa_low");
            bookCut("j1 is not b-tagged: bCa_low");
            bookCut("mT > 90: bCa_low");
            bookCut("MET/meff > 0.35: bCa_low");
            bookCut("pTl1 < 25: bCa_low");
            bookCut("MET > 250: bCb_high");
            bookCut("= 1 base lepton: bCb_high");
            bookCut("= 1 soft lepton: bCb_high");
            bookCut("pTl < 25: bCb_high");
            bookCut(">= 2 jets: bCb_high");
            bookCut("pTj2 > 60: bCb_high");
            bookCut("j1 and j2 are b-tagged: bCb_high");
            bookCut("dPhi_j_met_12 > 0.4: bCb_high");
            bookCut("mbb > 150: bCb_high");
            bookCut("amT2 > 200: bCb_high");
            bookCut("MET_trigger (MET > 80)");
            bookCut("SL_trigger");
            bookCut("Trigger");
            bookCut("= 1 base lepton: preselection");
            bookCut("= 1 lepton: preselection");
            bookCut(">= 3 jets (pT>25): preselection");
            bookCut("pTj1 > 80: bCc_diag");
            bookCut("pTj2 > 40: bCc_diag");
            bookCut("pTj3 > 30: bCc_diag");
            bookCut("no bjets with pT>25: bCc_diag");
            bookCut("|eta(lep);| < 1.2: bCc_diag");
            bookCut("dPhi_j1_met > 2: bCc_diag");
            bookCut("dPhi_j2_met > 0.8: bCc_diag");
            bookCut("dR_lep_j1 [0.8, 2.4]: bCc_diag");
            bookCut("MET > 140: bCc_diag");
            bookCut("MET/sqrt(HT); > 5: bCc_diag");
            bookCut("mT > 120: bCc_diag");
            bookCut(">= 4 jets (pT>25): bCb_med2");
            bookCut("pTj1 > 80: bCb_med2");
            bookCut("pTj2 > 60: bCb_med2");
            bookCut("pTj3 > 40: bCb_med2");
            bookCut(">= 2 bjets: bCb_med2");
            bookCut("pTb1 > 140: bCb_med2");
            bookCut("pTb2 > 75: bCb_med2");
            bookCut("dPhi_j1_met > 0.8: bCb_med2");
            bookCut("dPhi_j2_met > 0.8: bCb_med2");
            bookCut("MET > 170: bCb_med2");
            bookCut("MET/sqrt(HT); > 6: bCb_med2");
            bookCut("mT > 60: bCb_med2");
            bookCut("amT2 > 250: bCb_med2");
            bookCut("Veto on iso tracks: bCb_med2");
            bookCut("Veto on tight tau: bCb_med2");
            bookCut(">= 4 jets (pT>25): bCd_bulk");
            bookCut("pTj1 > 80: bCd_bulk");
            bookCut("pTj2 > 60: bCd_bulk");
            bookCut("pTj3 > 40: bCd_bulk");
            bookCut(">= 1 bjets: bCd_bulk");
            bookCut("dPhi_j1_met > 0.8: bCd_bulk");
            bookCut("dPhi_j2_met > 0.8: bCd_bulk");
            bookCut("MET > 150: bCd_bulk");
            bookCut("MET/sqrt(HT); > 7: bCd_bulk");
            bookCut("mT > 120: bCd_bulk");
            bookCut("amT2 > 175: bCd_bulk");
            bookCut("Veto on iso tracks: bCd_bulk");
            bookCut("Veto on tight tau: bCd_bulk");
            bookCut(">= 4 jets (pT>25): bCd_high1");
            bookCut("pTj1 > 80: bCd_high1");
            bookCut("pTj2 > 60: bCd_high1");
            bookCut("pTj3 > 40: bCd_high1");
            bookCut(">= 2 bjets: bCd_high1");
            bookCut("pTb2 > 75: bCd_high1");
            bookCut("dPhi_j1_met > 0.8: bCd_high1");
            bookCut("dPhi_j2_met > 0.8: bCd_high1");
            bookCut("MET > 150: bCd_high1");
            bookCut("MET/sqrt(HT); > 9: bCd_high1");
            bookCut("mT > 120: bCd_high1");
            bookCut("meff > 600: bCd_high1");
            bookCut("amT2 > 200: bCd_high1");
            bookCut("Veto on iso tracks: bCd_high1");
            bookCut("Veto on tight tau: bCd_high1");
            bookCut(">= 4 jets (pT>25): bCd_high2");
            bookCut("pTj1 > 80: bCd_high2");
            bookCut("pTj2 > 60: bCd_high2");
            bookCut("pTj3 > 40: bCd_high2");
            bookCut(">= 2 bjets: bCd_high2");
            bookCut("pTb1 > 170: bCd_high2");
            bookCut("pTb2 > 80: bCd_high2");
            bookCut("dPhi_j1_met > 0.8: bCd_high2");
            bookCut("dPhi_j2_met > 0.8: bCd_high2");
            bookCut("MET > 160: bCd_high2");
            bookCut("MET/sqrt(HT); > 8: bCd_high2");
            bookCut("mT > 120: bCd_high2");
            bookCut("amT2 > 250: bCd_high2");
            bookCut("Veto on iso tracks: bCd_high2");
            bookCut("Veto on tight tau: bCd_high2");
            bookCut(">= 4 jets (pT>25): 3body");
            bookCut("pTj1 > 80: 3body");
            bookCut(">= 1 bjets: 3body");
            bookCut("dPhi_j1_met > 0.2: 3body");
            bookCut("dPhi_j2_met > 0.2: 3body");
            bookCut("dR_lep_j1 > 1.25: 3body");
            bookCut("dR_lep_j2 > 2: 3body");
            bookCut("MET > 150: 3body");
            bookCut("MET/sqrt(HT); > 5: 3body");
            bookCut("mT > 120: 3body");
            bookCut("80 < amT2 < 100: 3body");
            bookCut("Veto on tight tau: 3body");

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt(&event);
            const Particles& jets25 = applyProjection<NearIsoParticle>(event, "Jets25").particlesByPt(&event);            
            const Particles& base_leps = applyProjection<NearIsoParticle>(event, "baseLeptons").particlesByPt(&event);
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt(&event);
            const Particles& soft_leps = applyProjection<NearIsoParticle>(event, "softLeptons").particlesByPt(&event);
            const Particles& leps_R = applyProjection<NearIsoParticle>(event, "Leptons_R").particlesByPt(&event);
            const Particles& taus = applyProjection<NearIsoParticle>(event, "Taus").particlesByPt(&event);
            const Particles& bjets60 = applyProjection<HeavyFlavorJets>(event, "BJets").getTaggedJets(&event);;
            const Particles& bjets70 = applyProjection<HeavyFlavorJets>(event, "BJets").getTaggedJets(&event);;
            const Particles& bjets80 = applyProjection<HeavyFlavorJets>(event, "BJets").getTaggedJets(&event);;            
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");

            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();

            //=============================//
            //     Variable Calculation    //
            //=============================//

            double HT = 0.;
            for(int i = 0; i < min(4, (int)jets.size()); i++){ 
                HT += jets[i].pT();
            }

            double dPhi_j_met_12 = -1;
            if(jets.size() > 1){
                double dPhi1 = deltaPhi( jets[0].momentum(), met );
                double dPhi2 = deltaPhi( jets[1].momentum(), met );
                dPhi_j_met_12 = min(dPhi1, dPhi2);
            }

            double HT2 = 0.;
            for(int i=2; i<jets.size(); i++){ 
                HT2 += jets[i].pT();
            }

            bool pass_tau_veto = true; // define!!
            bool pass_track_veto = true; // define!!

            //###################################################//
            //                Soft Lepton SRs                    //
            //###################################################//

            //=============================//
            //           bCa_med           //
            //=============================//
            if(cut(MET, CUT_GT, 300., "MET > 300: bCa_med")){
                if(cut(base_leps.size(), CUT_EQ, 1, "= 1 base lepton: bCa_med")){
                    if(cut(soft_leps.size(), CUT_EQ, 1, "= 1 soft lepton: bCa_med")){
                        Particle lep1 = soft_leps[0];
                        if(cut(lep1.pT(), CUT_LT, 50., "pTl < 50: bCa_med")){
                            if(cut(jets25.size(), CUT_GE, 3, ">= 3 jets (pT>25): bCa_med")){
                                if(cut(jets[0].pT(), CUT_GT, 180., "pTj1 > 180: bCa_med")){
                                    if(cut(bjets70.size(), CUT_GE, 1, ">= 1 bjets: bCa_med")){                                        
                                        if(cut(jets[0].pT() > bjets70[0].pT(), "j1 is not b-tagged: bCa_med")){
                                            double mT = get_mT( lep1.momentum(), met );                                            
                                            if(cut(mT, CUT_GT, 100., "mT > 100: bCa_med")){
                                                double meff = HT + lep1.pT() + MET;
                                                if(cut(MET/meff, CUT_GT, 0.3, "MET/meff > 0.3: bCa_med")){
                                                    if(cut(lep1.pT(), CUT_LT, 25., "pTl1 < 25: bCa_med")){
                                                        pass("bCa_med");
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
            //           bCb_med1          //
            //=============================//
            if(cut(MET, CUT_GT, 150., "MET > 150: bCb_med1")){
                if(cut(base_leps.size(), CUT_EQ, 1, "= 1 base lepton: bCb_med1")){
                    if(cut(soft_leps.size(), CUT_EQ, 1, "= 1 soft lepton: bCb_med1")){
                        Particle lep1 = soft_leps[0];
                        if(cut(lep1.pT(), CUT_LT, 25., "pTl < 25: bCb_med1")){
                            if(cut(jets.size(), CUT_GE, 2, ">= 2 jets: bCb_med1")){
                                if(cut(jets[1].pT(), CUT_GT, 60, "pTj2 > 60: bCb_med1")){
                                    if(cut(bjets60.size() > 1 && bjets60[1].pT() == jets[1].pT(), "j1 and j2 are b-tagged: bCb_med1")){
                                        Particles bjets = bjets60;
                                        if(cut(dPhi_j_met_12, CUT_GT, 0.4, "dPhi_j_met_12 > 0.4: bCb_med1")){
                                            if(cut(HT2, CUT_LT, 50., "HT_2 < 50: bCb_med1")){
                                                double mbb = (bjets[0].momentum() + bjets[1].momentum()).mass();
                                                if(cut(mbb, CUT_GT, 150., "mbb > 150: bCb_med1")){
                                                    double amt2 = get_amT2(lep1, bjets[0], bjets[1], met);
                                                    if(cut(amt2, CUT_GT, 170., "amT2 > 170: bCb_med1")){
                                                        pass("bCb_med1");
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
            //           bCa_low           //
            //=============================//
            if(cut(MET, CUT_GT, 370., "MET > 370: bCa_low")){
                if(cut(base_leps.size(), CUT_EQ, 1, "= 1 base lepton: bCa_low")){
                    if(cut(soft_leps.size(), CUT_EQ, 1, "= 1 soft lepton: bCa_low")){
                        Particle lep1 = soft_leps[0];
                        if(cut(lep1.pT(), CUT_LT, 50., "pTl < 50: bCa_low")){
                            if(cut(jets25.size(), CUT_GE, 2, ">= 2 jets (pT>25): bCa_low")){
                                if(cut(jets[0].pT(), CUT_GT, 180, "pTj1 > 180: bCa_low")){
                                    if(cut(bjets70.size(), CUT_GE, 1, ">= 1 bjets: bCa_low")){
                                        if(cut(jets[0].pT() > bjets70[0].pT(), "j1 is not b-tagged: bCa_low")){
                                            double mT = get_mT( lep1.momentum(), met );
                                            if(cut(mT, CUT_GT, 90., "mT > 90: bCa_low")){
                                                double meff = HT + lep1.pT() + MET;                                                
                                                if(cut(MET/meff, CUT_GT, 0.35, "MET/meff > 0.35: bCa_low")){
                                                    if(cut(lep1.pT(), CUT_LT, 25., "pTl1 < 25: bCa_low")){
                                                        pass("bCa_low");
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
            //          bCb_high           //
            //=============================//
            if(cut(MET, CUT_GT, 250., "MET > 250: bCb_high")){
                if(cut(base_leps.size(), CUT_EQ, 1, "= 1 base lepton: bCb_high")){
                    if(cut(soft_leps.size(), CUT_EQ, 1, "= 1 soft lepton: bCb_high")){
                        Particle lep1 = soft_leps[0];
                        if(cut(lep1.pT(), CUT_LT, 25., "pTl < 25: bCb_high")){
                            if(cut(jets.size(), CUT_GE, 2, ">= 2 jets: bCb_high")){
                                if(cut(jets[1].pT(), CUT_GT, 60, "pTj2 > 60: bCb_high")){
                                    if(cut(bjets60.size() > 1 && bjets60[1].pT() == jets[1].pT(), "j1 and j2 are b-tagged: bCb_high")){
                                        Particles bjets = bjets60;                                        
                                        if(cut(dPhi_j_met_12, CUT_GT, 0.4, "dPhi_j_met_12 > 0.4: bCb_high")){
                                            double mbb = (bjets[0].momentum() + bjets[1].momentum()).mass();
                                            if(cut(mbb, CUT_GT, 150., "mbb > 150: bCb_high")){
                                                double amt2 = get_amT2(lep1, bjets[0], bjets[1], met);
                                                if(cut(amt2, CUT_GT, 200., "amT2 > 200: bCb_high")){
                                                    pass("bCb_high");
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

            //###################################################//
            //            SRs without soft leptons               //
            //###################################################//

            //=============================//
            //          Trigger            //
            //=============================//

            bool pass_MET_trigger = false;
            if(cut(MET, CUT_GT, 80., "MET_trigger (MET > 80)")) pass_MET_trigger = true;

            bool pass_SL_trigger = false;
            if( leps_R.size() > 0 && leps_R[0].pT() > 24 ) pass_SL_trigger = true;
            if( base_leps.size() > 0 ){ 
                if( abs(base_leps[0].pdgId()) == 11 && base_leps[0].pT() > 60. ) pass_SL_trigger = true;
                if( abs(base_leps[0].pdgId()) == 13 && base_leps[0].pT() > 36. ) pass_SL_trigger = true;
            }
            cut(pass_SL_trigger, "SL_trigger");

            bool pass_trigger = false; 
            if(!cut( (pass_MET_trigger || pass_SL_trigger), "Trigger")) vetoEvent;

            //=============================//
            //        Common cuts          //
            //=============================//
            if(!cut(base_leps.size(), CUT_EQ, 1, "= 1 base lepton: preselection")) vetoEvent;
            if(!cut(leps.size(), CUT_EQ, 1, "= 1 lepton: preselection")) vetoEvent;
            Particle lep1 = leps[0];
            double mT = get_mT(lep1.momentum(), met);
            double meff = HT + lep1.pT() + MET;
            if(!cut(jets25.size(), CUT_GE, 3, ">= 3 jets (pT>25): preselection")) vetoEvent;


            //=============================//
            //          bCc_diag           //
            //=============================//
            if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: bCc_diag")){
                if(cut(jets[1].pT(), CUT_GT, 40., "pTj2 > 40: bCc_diag")){
                    if(cut(jets[2].pT(), CUT_GT, 30., "pTj3 > 30: bCc_diag")){
                        Particles bjets = bjets70;
                        bool pass_bveto = false;
                        if( bjets.size() == 0 ) pass_bveto = true;
                        if( bjets.size() > 0 && bjets[0].pT() < 25. ) pass_bveto = true;
                        if(cut(pass_bveto, "no bjets with pT>25: bCc_diag")){
                            if(cut(lep1.abseta(), CUT_LE, 1.2, "|eta(lep)| < 1.2: bCc_diag")){
                                double dPhi1 = deltaPhi(jets[0].momentum(), met);
                                if(cut(dPhi1, CUT_GT, 2., "dPhi_j1_met > 2: bCc_diag")){
                                    double dPhi2 = deltaPhi(jets[1].momentum(), met);
                                    if(cut(dPhi2, CUT_GT, 0.8, "dPhi_j2_met > 0.8: bCc_diag")){
                                        double dR_lep_j1 = deltaR(lep1, jets[0]);
                                        if(cut(dR_lep_j1, CUT_IN, make_pair(0.8, 2.4), "dR_lep_j1 [0.8, 2.4]: bCc_diag")){
                                            if(cut(MET, CUT_GT, 140., "MET > 140: bCc_diag")){
                                                if(cut(MET/sqrt(HT), CUT_GT, 5., "MET/sqrt(HT) > 5: bCc_diag")){
                                                    if(cut(mT, CUT_GT, 120., "mT > 120: bCc_diag")){
                                                        pass("bCc_diag");
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
            //          bCb_med2           //
            //=============================//
            if(cut(jets25.size(), CUT_GE, 4, ">= 4 jets (pT>25): bCb_med2")){
                if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: bCb_med2")){
                    if(cut(jets[1].pT(), CUT_GT, 60., "pTj2 > 60: bCb_med2")){
                        if(cut(jets[2].pT(), CUT_GT, 40., "pTj3 > 40: bCb_med2")){
                            if(cut(bjets80.size(), CUT_GE, 2, ">= 2 bjets: bCb_med2")){
                                Particles bjets = bjets80;
                                if(cut(bjets[0].pT(), CUT_GT, 140., "pTb1 > 140: bCb_med2")){
                                    if(cut(bjets[1].pT(), CUT_GT, 75., "pTb2 > 75: bCb_med2")){
                                        double dPhi1 = deltaPhi(jets[0].momentum(), met);
                                        if(cut(dPhi1, CUT_GT, 0.8, "dPhi_j1_met > 0.8: bCb_med2")){
                                            double dPhi2 = deltaPhi(jets[1].momentum(), met);
                                            if(cut(dPhi2, CUT_GT, 0.8, "dPhi_j2_met > 0.8: bCb_med2")){
                                                if(cut(MET, CUT_GT, 170., "MET > 170: bCb_med2")){
                                                    if(cut(MET/sqrt(HT), CUT_GT, 6., "MET/sqrt(HT) > 6: bCb_med2")){
                                                        if(cut(mT, CUT_GT, 60., "mT > 60: bCb_med2")){                                                            

                                                            double amt2 = get_amT2(lep1, bjets[1], bjets[0], met);

                                                            // only for cut-flow
                                                            if(cut(amt2, CUT_GT, 250., "amT2 > 250: bCb_med2")){
                                                                if(cut(pass_track_veto, "Veto on iso tracks: bCb_med2")){
                                                                    cut(pass_tau_veto, "Veto on tight tau: bCb_med2"); 
                                                                }
                                                            }

                                                            if(pass_track_veto && pass_tau_veto){
                                                                if( 175. < amt2 && amt2 < 250. ){
                                                                    if( 90. < mT && mT < 120. ) pass("bCb_med2_<250_<120");
                                                                    if( mT > 120. ) pass("bCb_med2_<250_>120");
                                                                }
                                                                if( amt2 > 250. ){
                                                                    if( 90. < mT && mT < 120. ) pass("bCb_med2_>250_<120");
                                                                    if( mT > 120. ) pass("bCb_med2_>250_>120");
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
            //          bCd_bulk           //
            //=============================//
            if(cut(jets25.size(), CUT_GE, 4, ">= 4 jets (pT>25): bCd_bulk")){
                if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: bCd_bulk")){
                    if(cut(jets[1].pT(), CUT_GT, 60., "pTj2 > 60: bCd_bulk")){
                        if(cut(jets[2].pT(), CUT_GT, 40., "pTj3 > 40: bCd_bulk")){
                            if(cut(bjets70.size(), CUT_GE, 1, ">= 1 bjets: bCd_bulk")){
                                Particles bjets = bjets70;
                                double dPhi1 = deltaPhi(jets[0].momentum(), met);
                                if(cut(dPhi1, CUT_GT, 0.8, "dPhi_j1_met > 0.8: bCd_bulk")){
                                    double dPhi2 = deltaPhi(jets[1].momentum(), met);
                                    if(cut(dPhi2, CUT_GT, 0.8, "dPhi_j2_met > 0.8: bCd_bulk")){
                                        if(cut(MET, CUT_GT, 150., "MET > 150: bCd_bulk")){
                                            if(cut(MET/sqrt(HT), CUT_GT, 7., "MET/sqrt(HT) > 7: bCd_bulk")){

                                                double amt2 = get_amT2(lep1, bjets[1], bjets[0], met);

                                                // only for cut-flow
                                                if(cut(mT, CUT_GT, 120., "mT > 120: bCd_bulk")){
                                                    if(cut(amt2, CUT_GT, 175., "amT2 > 175: bCd_bulk")){
                                                        if(cut(pass_track_veto, "Veto on iso tracks: bCd_bulk")){
                                                            cut(pass_tau_veto, "Veto on tight tau: bCd_bulk"); 
                                                        }
                                                    }
                                                }


                                                if(pass_track_veto && pass_tau_veto){
                                                    if( 175. < amt2 && amt2 < 250. ){
                                                        if( 90. < mT && mT < 120. ) pass("bCd_bulk_<250_<120");
                                                        if( mT > 120. ) pass("bCd_bulk_<250_>120");
                                                    }

                                                    if( amt2 > 250. ){
                                                        if( 90. < mT && mT < 120. ) pass("bCd_bulk_>250_<120");
                                                        if( mT > 120. ) pass("bCd_bulk_>250_>120");
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
            //          bCd_high1          //
            //=============================//
            if(cut(jets25.size(), CUT_GE, 4, ">= 4 jets (pT>25): bCd_high1")){
                if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: bCd_high1")){
                    if(cut(jets[1].pT(), CUT_GT, 60., "pTj2 > 60: bCd_high1")){
                        if(cut(jets[2].pT(), CUT_GT, 40., "pTj3 > 40: bCd_high1")){
                            if(cut(bjets80.size(), CUT_GE, 2, ">= 2 bjets: bCd_high1")){
                                Particles bjets = bjets80;
                                if(cut(bjets[1].pT(), CUT_GT, 75., "pTb2 > 75: bCd_high1")){
                                    double dPhi1 = deltaPhi(jets[0].momentum(), met);
                                    if(cut(dPhi1, CUT_GT, 0.8, "dPhi_j1_met > 0.8: bCd_high1")){
                                        double dPhi2 = deltaPhi(jets[1].momentum(), met);
                                        if(cut(dPhi2, CUT_GT, 0.8, "dPhi_j2_met > 0.8: bCd_high1")){
                                            if(cut(MET, CUT_GT, 150., "MET > 150: bCd_high1")){
                                                if(cut(MET/sqrt(HT), CUT_GT, 9., "MET/sqrt(HT) > 9: bCd_high1")){
                                                    if(cut(mT, CUT_GT, 120., "mT > 120: bCd_high1")){                                                            
                                                        if(cut(meff, CUT_GT, 600., "meff > 600: bCd_high1")){                                                            
                                                            double amt2 = get_amT2(lep1, bjets[1], bjets[0], met);
                                                            if(cut(amt2, CUT_GT, 200., "amT2 > 200: bCd_high1")){
                                                                if(cut(pass_track_veto, "Veto on iso tracks: bCd_high1")){
                                                                    if(cut(pass_tau_veto, "Veto on tight tau: bCd_high1")){
                                                                        pass("bCd_high1");
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
            //          bCd_high2          //
            //=============================//
            if(cut(jets25.size(), CUT_GE, 4, ">= 4 jets (pT>25): bCd_high2")){
                if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: bCd_high2")){
                    if(cut(jets[1].pT(), CUT_GT, 60., "pTj2 > 60: bCd_high2")){
                        if(cut(jets[2].pT(), CUT_GT, 40., "pTj3 > 40: bCd_high2")){
                            if(cut(bjets80.size(), CUT_GE, 2, ">= 2 bjets: bCd_high2")){
                                Particles bjets = bjets80;
                                if(cut(bjets[0].pT(), CUT_GT, 170., "pTb1 > 170: bCd_high2")){
                                    if(cut(bjets[1].pT(), CUT_GT, 80., "pTb2 > 80: bCd_high2")){
                                        double dPhi1 = deltaPhi(jets[0].momentum(), met);
                                        if(cut(dPhi1, CUT_GT, 0.8, "dPhi_j1_met > 0.8: bCd_high2")){
                                            double dPhi2 = deltaPhi(jets[1].momentum(), met);
                                            if(cut(dPhi2, CUT_GT, 0.8, "dPhi_j2_met > 0.8: bCd_high2")){
                                                if(cut(MET, CUT_GT, 160., "MET > 160: bCd_high2")){
                                                    if(cut(MET/sqrt(HT), CUT_GT, 8., "MET/sqrt(HT) > 8: bCd_high2")){
                                                        if(cut(mT, CUT_GT, 120., "mT > 120: bCd_high2")){                                                            
                                                            double amt2 = get_amT2(lep1, bjets[1], bjets[0], met);
                                                            if(cut(amt2, CUT_GT, 250., "amT2 > 250: bCd_high2")){
                                                                if(cut(pass_track_veto, "Veto on iso tracks: bCd_high2")){
                                                                    if(cut(pass_tau_veto, "Veto on tight tau: bCd_high2")){
                                                                        pass("bCd_high2");
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
            //           3body             //
            //=============================//
            if(cut(jets25.size(), CUT_GE, 4, ">= 4 jets (pT>25): 3body")){
                if(cut(jets[0].pT(), CUT_GT, 80., "pTj1 > 80: 3body")){

                    if(cut(bjets70.size(), CUT_GE, 1, ">= 1 bjets: 3body")){
                        Particles bjets = bjets70;
                        double dPhi1 = deltaPhi(jets[0].momentum(), met);
                        if(cut(dPhi1, CUT_GT, 0.2, "dPhi_j1_met > 0.2: 3body")){
                            double dPhi2 = deltaPhi(jets[1].momentum(), met);
                            if(cut(dPhi2, CUT_GT, 0.2, "dPhi_j2_met > 0.2: 3body")){
                                double dR_lep_j1 = deltaR(lep1, jets[0]);
                                if(cut(dR_lep_j1, CUT_GT, 1.25, "dR_lep_j1 > 1.25: 3body")){
                                    double dR_lep_j2 = deltaR(lep1, jets[1]);
                                    if(cut(dR_lep_j2, CUT_GT, 2., "dR_lep_j2 > 2: 3body")){
                                        if(cut(MET, CUT_GT, 150., "MET > 150: 3body")){
                                            if(cut(MET/sqrt(HT), CUT_GT, 5., "MET/sqrt(HT) > 5: 3body")){

                                                double amt2 = get_amT2(lep1, bjets[1], bjets[0], met);

                                                // only for cut-flow
                                                if(cut(mT, CUT_GT, 120., "mT > 120: 3body")){
                                                    if(cut(amt2, CUT_IN, make_pair(80., 100.), "80 < amT2 < 100: 3body")){
                                                        cut(pass_tau_veto, "Veto on tight tau: 3body"); 
                                                    }
                                                }


                                                if(pass_tau_veto){
                                                    if( 80. < amt2 && amt2 < 90. ){
                                                        if( 90. < mT && mT < 120. ) pass("3body_<90_<120");
                                                        if( mT > 120. ) pass("3body_<90_>120");
                                                    }

                                                    if( 90. < amt2 && amt2 < 100. ){
                                                        if( 90. < mT && mT < 120. ) pass("3body_>90_<120");
                                                        if( mT > 120. ) pass("3body_>90_>120");
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


		/// Normalise histograms etc., after the run
		void finalize() {
			/// @todo normalize the histograms
			// scale("Mjj");
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1407_0583)
}
