///
/// @file  ATLAS_1402_7029.cc
/// @brief Implementation of ATLAS_1402_7029 analysis
/// @author Kazuki
/// @date created 04/24/2015
/// @date last revision 04/24/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"
#include "Rivet/Tools/RivetMT2.hh"

using namespace std;

namespace Atom {

	class ATLAS_1402_7029 : public Analysis {
	public:

		ATLAS_1402_7029()
			: Analysis("ATLAS_1402_7029") {
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

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -2.5, 2.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            IsoElectron base_ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            //base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoElectron ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            ele.addIso(TRACK_ISO_PT, 0.3,  0.16,  0.0, 0.01, CALO_ALL);
            ele.addIso(CALO_ISO_PT,  0.3,  0.18,  0.0, 0.01, CALO_ALL);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon base_mu(Range(PT, 10., 8000.) & Range(ETA, -2.5, 2.5));
            //base_mu.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            base_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            IsoMuon mu(Range(PT, 10., 8000.) & Range(ETA, -2.5, 2.5));
            mu.addIso(TRACK_ISO_PT, 0.3,  0.12,  0.0, 0.01, CALO_ALL);
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            // Overlap removal

            NearIsoParticle ele_base1(base_ele);
            ele_base1.addFilter(base_ele, 0.1);

            NearIsoParticle ele_1(ele);
            ele_1.addFilter(base_ele, 0.1);

            NearIsoParticle jets_clean1(jets);
            jets_clean1.addFilter(ele_base1, 0.2);

            const double dRprong = 0.2;
            ParamTauFinder taus(jets_clean1, dRprong, Range(ETA, -2.47, 2.47));
              
            NearIsoParticle taus_clean(taus);
            taus_clean.addFilter(base_mu, 0.2);
            addProjection(taus_clean, "Taus");
            
            NearIsoParticle ele_clean(ele_1);
            ele_clean.addFilter(jets_clean1, 0.4);
            addProjection(ele_clean, "Electrons");

            NearIsoParticle mu_clean(mu);
            mu_clean.addFilter(jets_clean1, 0.4);
            addProjection(mu_clean, "Muons");

            MergedFinalState leptons(ele_clean, mu_clean);
            //MergedFinalState leptons(base_ele, base_mu);

            addProjection(leptons, "Leptons");

            MergedFinalState inc_leptons(leptons, taus_clean);
            addProjection(inc_leptons, "incLeptons");
                
            NearIsoParticle jets_clean(jets_clean1);
            jets_clean.addFilter(taus_clean, 0.2);
            addProjection(jets_clean, "Jets");

            Range bjrange = Range(PT, 20.0, 8000.0) & Range(ETA, -2.5, 2.5);
            HeavyFlavorJets bjets(jets_clean, bjrange);
            addProjection(bjets, "BJets");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState base_leps(base_ele, base_mu);
            MergedFinalState met_seed(jets_clean1, base_leps);
            MissingMomentum met( fsbase, met_seed );
            met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Const_PlaceHolder" ) );            
            addProjection(met, "MissingEt");

            /// @todo book the efficiencies
            bookEfficiency("SR0ta1");
            bookEfficiency("SR0ta2");
            bookEfficiency("SR0ta3");
            bookEfficiency("SR0ta4");
            bookEfficiency("SR0ta5");
            bookEfficiency("SR0ta6");
            bookEfficiency("SR0ta7");
            bookEfficiency("SR0ta8");
            bookEfficiency("SR0ta9");
            bookEfficiency("SR0ta10");
            bookEfficiency("SR0ta11");
            bookEfficiency("SR0ta12");
            bookEfficiency("SR0ta13");
            bookEfficiency("SR0ta14");
            bookEfficiency("SR0ta15");
            bookEfficiency("SR0ta16");
            bookEfficiency("SR0ta17");
            bookEfficiency("SR0ta18");
            bookEfficiency("SR0ta19");
            bookEfficiency("SR0ta20");
            bookEfficiency("SR1t");
            bookEfficiency("SR2ta");
            bookEfficiency("SR2tb");

            /// @todo book the cuts
            bookCut("= 3 leptons(e,mu,tau)");
            bookCut("dRmin(l1,l2,l3) > 0.3");
            bookCut("at least one e or mu");
            bookCut("pTl1 > 25: SL trigger");
            bookCut("DL trigger");
            bookCut("Trigger");
            bookCut("= 3 LF leptons: SR0ta");
            bookCut("at least one SFOS pair: SR0ta");
            bookCut("12 < mSFOS < 40: SR0ta(1-4)");
            bookCut("b veto: SR0ta(1-4)");
            bookCut("50 < MET < 90: SR0ta1");
            bookCut("mT < 80: SR0ta1");
            bookCut("MET > 90: SR0ta2");
            bookCut("mT < 80: SR0ta2");
            bookCut("50 < MET < 75: SR0ta3");
            bookCut("mT > 80: SR0ta3");
            bookCut("MET > 75: SR0ta4");
            bookCut("mT > 80: SR0ta4");
            bookCut("40 < mSFOS < 60: SR0ta(5-8)");
            bookCut("b veto: SR0ta(5-8)");
            bookCut("50 < MET < 75: SR0ta5");
            bookCut("mT < 80: SR0ta5");
            bookCut("Z veto: SR0ta5");
            bookCut("MET > 75: SR0ta6");
            bookCut("mT < 80: SR0ta6");
            bookCut("50 < MET < 135: SR0ta7");
            bookCut("mT > 80: SR0ta7");
            bookCut("MET > 135: SR0ta8");
            bookCut("mT > 80: SR0ta8");
            bookCut("60 < mSFOS < 81: SR0ta(9-12)");
            bookCut("b veto: SR0ta(9-12)");
            bookCut("50 < MET < 75: SR0ta9");
            bookCut("mT < 80: SR0ta9");
            bookCut("Z veto: SR0ta9");
            bookCut("50 < MET < 75: SR0ta10");
            bookCut("mT > 80: SR0ta10");
            bookCut("MET > 75: SR0ta11");
            bookCut("mT < 110: SR0ta11");
            bookCut("MET > 75: SR0ta12");
            bookCut("mT > 110: SR0ta12");
            bookCut("81 < mSFOS < 101: SR0ta(13-16)");
            bookCut("b veto: SR0ta(13-16)");
            bookCut("50 < MET < 90: SR0ta13");
            bookCut("mT < 110: SR0ta13");
            bookCut("Z veto: SR0ta13");
            bookCut("MET > 90: SR0ta14");
            bookCut("mT < 110: SR0ta14");
            bookCut("50 < MET < 135: SR0ta15");
            bookCut("mT > 110: SR0ta15");
            bookCut("MET > 135: SR0ta16");
            bookCut("mT > 110: SR0ta16");
            bookCut("mSFOS > 101: SR0ta(17-20)");
            bookCut("b veto: SR0ta(17-20)");
            bookCut("50 < MET < 210: SR0ta17");
            bookCut("mT < 180: SR0ta17");
            bookCut("50 < MET < 210: SR0ta18");
            bookCut("mT > 180: SR0ta18");
            bookCut("MET > 210: SR0ta19");
            bookCut("mT < 120: SR0ta19");
            bookCut("MET > 210: SR0ta20");
            bookCut("mT > 120: SR0ta20");
            bookCut("= 3 LF leptons: SR0tb");
            bookCut("no SFOS pair: SR0tb");
            bookCut("b veto: SR0tb");
            bookCut("MET > 50: SR0tb");
            bookCut("pTl3 > 20: SR0tb");
            bookCut("dRmin_OS < 1: SR0tb");
            bookCut("= 1 tau: SR1t");
            bookCut("SS lepton pair: SR1t");
            bookCut("Z->ee veto: SR1t");
            bookCut("b veto: SR1t");
            bookCut("MET > 50: SR1t");
            bookCut("pTl1 + pTl2 > 70: SR1t");
            bookCut("pTl1 + pTl2 > 70: SR1t");
            bookCut("pTl2 > 30: SR1t");
            bookCut("m(lep,tau) < 120: SR1t");
            bookCut("= 2 taus: SR2ta");
            bookCut("= 1 LF lepton: SR2ta");
            bookCut("b veto: SR2ta");
            bookCut("MET > 50: SR2ta");
            bookCut("mT2_max > 100: SR2ta");
            bookCut("= 2 OS taus: SR2tb");
            bookCut("= 1 LF lepton: SR2tb");
            bookCut("b veto: SR2tb");
            bookCut("MET > 50: SR2tb");
            bookCut("70 < m(tau,tau) < 120: SR2tb");
            bookCut("pTtau1 + pTtau2 > 110: SR2tb");


		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets          = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt(&event);
            const Particles& eles          = applyProjection<NearIsoParticle>(event, "Electrons").particlesByPt(&event);
            const Particles& mus           = applyProjection<NearIsoParticle>(event, "Muons").particlesByPt(&event);
            const Particles& leps          = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt(&event);
            const Particles& inc_leps      = applyProjection<MergedFinalState>(event, "incLeptons").particlesByPt(&event);            
            const Particles& taus          = applyProjection<NearIsoParticle>(event, "Taus").particlesByPt(&event);            
            const HeavyFlavorJets& bjproj  = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets         = bjproj.getTaggedJets(&event); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet    = applyProjection<MissingMomentum>(event, "MissingEt");

            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();

            //=============================//
            //         preselecton         //
            //=============================//

            if(!cut(inc_leps.size(), CUT_EQ, 3, "= 3 leptons(e,mu,tau)")) vetoEvent;

            double dRmin = 1000.;
            for(int i=0; i<inc_leps.size()-1; i++){
                for(int j=i+1; j<inc_leps.size(); j++){
                    double dR = deltaR(inc_leps[i], inc_leps[j]);
                    if(dR < dRmin) dRmin = dR;
                }
            }
            if(!cut(dRmin, CUT_GT, 0.3, "dRmin(l1,l2,l3) > 0.3")) vetoEvent;

            if(!cut(leps.size(), CUT_GE, 1, "at least one e or mu")) vetoEvent;

            //------------------------------//
            //           Trigger            //
            //------------------------------//

            bool pass_SL_trigger = false;
            if(cut(leps[0].pT(), CUT_GT, 25., "pTl1 > 25: SL trigger")) 
                pass_SL_trigger = true;

            bool pass_DL_trigger = false;
            if( eles.size() > 1 ){ 
                if( eles[1].pT() > 14. ) pass_DL_trigger = true;
                if( eles[0].pT() > 25. ) pass_DL_trigger = true;
            }
            if( mus.size() > 1){
                if( mus[1].pT() > 14. ) pass_DL_trigger = true;
                if( mus[0].pT() > 18. ) pass_DL_trigger = true;                
            }
            if( eles.size() > 0 && mus.size() > 0 ){
                if( eles[0].pT() > 14. ) pass_DL_trigger = true;
                if( mus[0].pT() > 18. ) pass_DL_trigger = true;
            }
            cut(pass_DL_trigger, "DL trigger");

            if(!cut( (pass_SL_trigger || pass_DL_trigger) , "Trigger")) vetoEvent;

            //-------------------------------//
            //     Variable calculations     //
            //-------------------------------//

            // find SFOS
            double min_measure = 1000000.;
            double mZ = 91.;
            int iSFOS1 = 0; 
            int iSFOS2 = 0;            
            int N_SFOS = 0;
            double mSFOS = 0;
            for(int i=0; i<leps.size()-1; i++){
                for(int j=i+1; j<leps.size(); j++){
                    if( leps[i].pdgId() + leps[j].pdgId() == 0 ){ // SFOS
                        double dm = (leps[i].momentum() + leps[j].momentum()).mass();
                        if( abs(dm - mZ) < min_measure ){
                            min_measure = abs(dm - mZ);
                            iSFOS1 = i;
                            iSFOS2 = j;
                            mSFOS = dm;
                            N_SFOS++;
                        }
                    }
                }
            }

            int i_notSFOS = 0;
            if(N_SFOS > 0){
                for(int i=0; i<leps.size(); i++){
                    if( i != iSFOS1 && i != iSFOS2 ) i_notSFOS = i;
                }
            }

            double dRminOS = 1000.;
            for(int i=0; i<leps.size()-1; i++){
                for(int j=i+1; j<leps.size(); j++){
                    if( leps[i].pdgId() * leps[i].pdgId() < 0 ){
                        double dR = deltaR(leps[i], leps[j]);
                        if(dR < dRminOS) dRminOS = dR;
                    }
                }
            }

            //------------------------------//
            //            SR0ta             //
            //------------------------------//

            if(cut( leps.size(), CUT_EQ, 3, "= 3 LF leptons: SR0ta")){
                if(cut( N_SFOS, CUT_GE, 1, "at least one SFOS pair: SR0ta")){

                    double mT = get_mT(leps[i_notSFOS].momentum(), met);
                    double m3l = (leps[0].momentum() + leps[1].momentum() + leps[2].momentum() ).mass();

                    //SR0ta(1-4)
                    if(cut( mSFOS, CUT_IN, make_pair(12., 40.), "12 < mSFOS < 40: SR0ta(1-4)")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0ta(1-4)")){
                            
                            if(cut( MET, CUT_IN, make_pair(50., 90.), "50 < MET < 90: SR0ta1")){
                                if(cut( mT, CUT_LT, 80., "mT < 80: SR0ta1")){
                                    pass("SR0ta1");
                                }
                            }
                            if(cut( MET, CUT_GT, 90., "MET > 90: SR0ta2")){
                                if(cut( mT, CUT_LT, 80., "mT < 80: SR0ta2")){
                                    pass("SR0ta2");
                                }
                            }
                            if(cut( MET, CUT_IN, make_pair(50., 75.), "50 < MET < 75: SR0ta3")){
                                if(cut( mT, CUT_GT, 80., "mT > 80: SR0ta3")){
                                    pass("SR0ta3");
                                }
                            }
                            if(cut( MET, CUT_GT, 75., "MET > 75: SR0ta4")){
                                if(cut( mT, CUT_GT, 80., "mT > 80: SR0ta4")){
                                    pass("SR0ta4");
                                }
                            }

                        }
                    }

                    //SR0ta(5-8)
                    if(cut( mSFOS, CUT_IN, make_pair(40., 60.), "40 < mSFOS < 60: SR0ta(5-8)")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0ta(5-8)")){

                            if(cut( MET, CUT_IN, make_pair(50., 75.), "50 < MET < 75: SR0ta5")){
                                if(cut( mT, CUT_LT, 80., "mT < 80: SR0ta5")){
                                    if(cut( abs(m3l - mZ), CUT_GT, 10., "Z veto: SR0ta5")){
                                        pass("SR0ta5");
                                    }
                                }
                            }
                            if(cut( MET, CUT_GT, 75., "MET > 75: SR0ta6")){
                                if(cut( mT, CUT_LT, 80., "mT < 80: SR0ta6")){
                                    pass("SR0ta6");
                                }
                            }
                            if(cut( MET, CUT_IN, make_pair(50., 135.), "50 < MET < 135: SR0ta7")){
                                if(cut( mT, CUT_GT, 80., "mT > 80: SR0ta7")){
                                    pass("SR0ta7");
                                }
                            }
                            if(cut( MET, CUT_GT, 135., "MET > 135: SR0ta8")){
                                if(cut( mT, CUT_GT, 80., "mT > 80: SR0ta8")){
                                    pass("SR0ta8");
                                }
                            }

                        }
                    }

                    //SR0ta(9-12)
                    if(cut( mSFOS, CUT_IN, make_pair(60., 81.2), "60 < mSFOS < 81: SR0ta(9-12)")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0ta(9-12)")){

                            if(cut( MET, CUT_IN, make_pair(50., 75.), "50 < MET < 75: SR0ta9")){
                                if(cut( mT, CUT_LT, 80., "mT < 80: SR0ta9")){
                                    if(cut( abs(m3l - mZ), CUT_GT, 10., "Z veto: SR0ta9")){
                                        pass("SR0ta9");
                                    }
                                }
                            }
                            if(cut( MET, CUT_IN, make_pair(50., 75.), "50 < MET < 75: SR0ta10")){
                                if(cut( mT, CUT_GT, 80., "mT > 80: SR0ta10")){
                                    pass("SR0ta9");
                                }
                            }
                            if(cut( MET, CUT_GT, 75., "MET > 75: SR0ta11")){
                                if(cut( mT, CUT_LT, 110., "mT < 110: SR0ta11")){
                                    pass("SR0ta11");
                                }
                            }
                            if(cut( MET, CUT_GT, 75., "MET > 75: SR0ta12")){
                                if(cut( mT, CUT_GT, 110., "mT > 110: SR0ta12")){
                                    pass("SR0ta12");
                                }
                            }

                        }
                    }

                    //SR0ta(13-16)
                    if(cut( mSFOS, CUT_IN, make_pair(81.2, 101.2), "81 < mSFOS < 101: SR0ta(13-16)")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0ta(13-16)")){

                            if(cut( MET, CUT_IN, make_pair(50., 90.), "50 < MET < 90: SR0ta13")){
                                if(cut( mT, CUT_LT, 110., "mT < 110: SR0ta13")){
                                    if(cut( abs(m3l - mZ), CUT_GT, 10., "Z veto: SR0ta13")){
                                        pass("SR0ta13");
                                    }
                                }
                            }
                            if(cut( MET, CUT_GT, 90., "MET > 90: SR0ta14")){
                                if(cut( mT, CUT_LT, 110., "mT < 110: SR0ta14")){
                                    pass("SR0ta14");
                                }
                            }

                            if(cut( MET, CUT_IN, make_pair(50., 135.), "50 < MET < 135: SR0ta15")){
                                if(cut( mT, CUT_GT, 110., "mT > 110: SR0ta15")){
                                    pass("SR0ta15");
                                }
                            }
                            if(cut( MET, CUT_GT, 135., "MET > 135: SR0ta16")){
                                if(cut( mT, CUT_GT, 110., "mT > 110: SR0ta16")){
                                    pass("SR0ta16");
                                }
                            }

                        }
                    }

                    //SR0ta(17-20)
                    if(cut( mSFOS, CUT_GT, 101.2, "mSFOS > 101: SR0ta(17-20)")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0ta(17-20)")){

                            if(cut( MET, CUT_IN, make_pair(50., 210.), "50 < MET < 210: SR0ta17")){
                                if(cut( mT, CUT_LT, 180., "mT < 180: SR0ta17")){
                                    pass("SR0ta17");
                                }
                            }
                            if(cut( MET, CUT_IN, make_pair(50., 210.), "50 < MET < 210: SR0ta18")){
                                if(cut( mT, CUT_GT, 180., "mT > 180: SR0ta18")){
                                    pass("SR0ta18");
                                }
                            }
                            if(cut( MET, CUT_GT, 210., "MET > 210: SR0ta19")){
                                if(cut( mT, CUT_LT, 120., "mT < 120: SR0ta19")){
                                    pass("SR0ta19");
                                }
                            }
                            if(cut( MET, CUT_GT, 210., "MET > 210: SR0ta20")){
                                if(cut( mT, CUT_GT, 120., "mT > 120: SR0ta20")){
                                    pass("SR0ta20");
                                }
                            }

                        }
                    }

                }
            }            

            //------------------------------//
            //            SR0tb             //
            //------------------------------//
            if(cut( leps.size(), CUT_EQ, 3, "= 3 LF leptons: SR0tb")){
                if(cut( N_SFOS, CUT_EQ, 0, "no SFOS pair: SR0tb")){
                    if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR0tb")){
                        if(cut( MET, CUT_GT, 50., "MET > 50: SR0tb")){
                            if(cut( leps[2].pT(), CUT_GT, 20., "pTl3 > 20: SR0tb")){
                                if(cut( dRminOS, CUT_LT, 1., "dRmin_OS < 1: SR0tb")){
                                    pass("SR0tb");
                                }
                            }
                        }
                    }
                }
            }

            //------------------------------//
            //            SR1t              //
            //------------------------------//
            if(cut( taus.size(), CUT_EQ, 1, "= 1 tau: SR1t")){
                bool pass_SS;
                if( leps.size() == 2 && leps[0].pdgId() * leps[1].pdgId() > 0 ){
                    pass_SS = true;
                }else{
                    pass_SS = false;
                }
                if(cut( pass_SS, "SS lepton pair: SR1t")){
                    double mee = 0;
                    if( leps[0].pdgId() * leps[1].pdgId() == 11 * 11 ){
                        mee = (leps[0].momentum() + leps[1].momentum()).mass();
                    }
                    if(cut( mee, CUT_OUT, make_pair(81.2, 101.2), "Z->ee veto: SR1t")){
                        if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR1t")){
                            if(cut( MET, CUT_GT, 50., "MET > 50: SR1t")){
                                if(cut( leps[0].pT() + leps[1].pT(), CUT_GT, 70., "pTl1 + pTl2 > 70: SR1t")){
                                    if(cut( leps[0].pT() + leps[1].pT(), CUT_GT, 70., "pTl1 + pTl2 > 70: SR1t")){
                                        if(cut( leps[1].pT(), CUT_GT, 30., "pTl2 > 30: SR1t")){
                                            double mtl1 = (leps[0].momentum() + taus[0].momentum()).mass();
                                            double mtl2 = (leps[1].momentum() + taus[0].momentum()).mass();
                                            double mtl = mtl1;
                                            double mh = 125.;
                                            if( abs(mh - mtl2) < abs(mh - mtl1) ) mtl = mtl2; 
                                            if(cut( mtl, CUT_LT, 120., "m(lep,tau) < 120: SR1t")){
                                                pass("SR1t");
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }

            //------------------------------//
            //           SR2ta              //
            //------------------------------//
            if(cut( taus.size(), CUT_EQ, 2, "= 2 taus: SR2ta")){
                if(cut( leps.size(), CUT_EQ, 1, "= 1 LF lepton: SR2ta")){
                    if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR2ta")){
                        if(cut( MET, CUT_GT, 50., "MET > 50: SR2ta")){

                            double mT2max = -100.;
                            for(int i=0; i<inc_leps.size()-1; i++){
                                for(int j=i+1; j<inc_leps.size(); j++){
                                    double minv = 0.0;
                                    double mt2 = Rivet::mT2::mT2(inc_leps[i].momentum(), inc_leps[j].momentum(), met, minv);
                                    if(mt2 > mT2max) mT2max = mt2;
                                }
                            }

                            if(cut( mT2max, CUT_GT, 100., "mT2_max > 100: SR2ta")){
                                pass("SR2ta");
                            }
                        }
                    }
                }
            }

            if(cut( taus.size() == 2 && taus[0].pdgId() * taus[1].pdgId() < 0, "= 2 OS taus: SR2tb")){
                if(cut( leps.size(), CUT_EQ, 1, "= 1 LF lepton: SR2tb")){
                    if(cut( bjets.size(), CUT_EQ, 0, "b veto: SR2tb")){
                        if(cut( MET, CUT_GT, 50., "MET > 50: SR2tb")){
                            double mtt = (taus[0].momentum() + taus[1].momentum()).mass();
                            if(cut( mtt, CUT_IN, make_pair(70., 120.), "70 < m(tau,tau) < 120: SR2tb")){
                                double pTtt = taus[0].pT() + taus[1].pT();
                                if(cut( pTtt, CUT_GT, 110., "pTtau1 + pTtau2 > 110: SR2tb")){
                                    pass("SRt2b");
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
	AtomPlugin(ATLAS_1402_7029)
}
