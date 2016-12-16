///
/// @file  ATLAS_CONF_2016_077.cc
/// @brief Implementation of ATLAS_CONF_2016_077 analysis
/// @author Kazuki
/// @date created 11/19/2016
/// @date last revision 11/19/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

    class ATLAS_CONF_2016_077 : public Analysis {
    public:

        ATLAS_CONF_2016_077()
            : Analysis("ATLAS_CONF_2016_077") {
       //     setNeedsCrossSection(true);
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

            IsoElectron ele( Range(PT > 10. & abseta < 2.47) );
            ele.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
            ele.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Loose_2012_ATLAS" ) );
        
            IsoMuon mu( Range(PT > 10. & abseta < 2.4) );
            mu.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
            mu.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );
        
            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                    hadRange & Range( PT > 20 & abseta < 2.8 ), 
                    muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );
            
            FastJets jets08(fsbase, 
                    hadRange & Range( PT > 20 & abseta < 2.8 ), 
                    muDetRange, FastJets::ANTIKT, 0.8 );
            jets08.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets08.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            FastJets jets12(fsbase, 
                    hadRange & Range( PT > 20 & abseta < 2.8 ), 
                    muDetRange, FastJets::ANTIKT, 1.2 );
            jets12.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets12.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );


            HeavyFlavorJets bjets(jets);
            bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
            bjets.setCurrentWorkingPoint( 0.77 );
        
            const double dRprong = 0.2;
            Range tau_range = Range(PT > 20) & (Range(ETA, -2.47, 2.47) - Range(ETA, -1.37, -1.52) - Range(ETA, 1.37, 1.52));
            ParamTauFinder taus(jets, 0.4, dRprong, tau_range);
            taus.setTaggingParams( getTauTag("Tau_Ident_BDT_Medium_ATLAS") );
        
            // Overlap removal
            NearIsoParticle taus_clean(taus);
            taus_clean.addFilter(ele, 0.2);
            taus_clean.addFilter(mu, 0.2);
        
            NearIsoParticle jets_clean1(jets);
            jets_clean1.addFilter(ele, 0.2);
        
            NearIsoParticle ele_clean(ele);
            ele_clean.addFilter(jets_clean1, 0.4);
        
            NearIsoParticle jets_clean(jets_clean1);
            jets_clean.addFilter(mu, 0.4);
            jets_clean.addFilter(taus_clean, 0.4);
            
            NearIsoParticle jets08_(jets08);
            NearIsoParticle jets12_(jets12);

            MergedFinalState leptons(ele_clean, mu);
        
            addProjection(leptons,     "Leptons");
            addProjection(jets_clean,  "Jets");
            addProjection(bjets,       "BJets");
            addProjection(taus_clean,  "TAUJets");
            addProjection(jets08_,       "Jets08");
            addProjection(jets12_,       "Jets12");
        
            MergedFinalState met_seed(jets_clean, leptons);
            MissingMomentum met( fsbase, met_seed );
            met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");
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
            bookEfficiency("SRA-TT");
            bookEfficiency("SRA-TW");
            bookEfficiency("SRA-T0");
            bookEfficiency("SRB-TT");
            bookEfficiency("SRB-TW");
            bookEfficiency("SRB-T0");
            bookEfficiency("SRE");
            bookEfficiency("SRF");

            // Cuts booking section -- do not edit/remove this comment
            bookCut("jet cuts: base");
            bookCut("deltaphi cuts: base");
            bookCut("0 lepton: base");
            bookCut("MET > 250: base");
            bookCut("jet12[0].mass > 120: TT-base");
            bookCut("jet12[1].mass > 120: TT-base");
            bookCut("jet08[0].mass > 60: SRA-TT");
            bookCut("bjets.size() >= 2: SRA-TT");
            bookCut("mT_minb > 200: SRA-TT");
            bookCut("tau veto: SRA-TT");
            bookCut("MET > 400: SRA-TT");
            bookCut("bjets >= 2: SRB-TT");
            bookCut("mT_minb > 200: SRB-TT");
            bookCut("mT_maxb > 200: SRB-TT");
            bookCut("tau veto: SRB-TT");
            bookCut("deltaR_b > 1.2: SRB-TT");
            bookCut("MET > 250: SRB-TT");
            bookCut("jet12[0].mass > 120: TW-base");
            bookCut("jet12[1].mass [60,120]: TW-base");
            bookCut("jet08[0].mass > 60: SRA-TW");
            bookCut("bjets.size() >= 2: SRA-TW");
            bookCut("mT_minb > 200: SRA-TW");
            bookCut("tau veto: SRA-TW");
            bookCut("MET > 450: SRA-TW");
            bookCut("bjets >= 2: SRB-TW");
            bookCut("mT_minb > 200: SRB-TW");
            bookCut("mT_maxb > 200: SRB-TW");
            bookCut("tau veto: SRB-TW");
            bookCut("deltaR_b > 1.2: SRB-TW");
            bookCut("MET > 250: SRB-TW");
            bookCut("jet12[0].mass > 120: T0-base");
            bookCut("jet12[1].mass < 60: T0-base");
            bookCut("jet08[0].mass > 60: SRA-T0");
            bookCut("bjets.size() >= 2: SRA-T0");
            bookCut("mT_minb > 200: SRA-T0");
            bookCut("tau veto: SRA-T0");
            bookCut("MET > 400: SRA-T0");
            bookCut("bjets >= 2: SRB-T0");
            bookCut("mT_minb > 200: SRB-T0");
            bookCut("mT_maxb > 200: SRB-T0");
            bookCut("tau veto: SRB-T0");
            bookCut("deltaR_b > 1.2: SRB-T0");
            bookCut("MET > 500: SRB-T0");
            bookCut("bjets >= 2: SREF base");
            bookCut("jet12[0].mass > 140: SRE");
            bookCut("jet12[1].mass > 60: SRE");
            bookCut("mT_minb > 200: SRE");
            bookCut("tau veto: SRE");
            bookCut("deltaR_b > 1.5: SRE");
            bookCut("MET > 300: SRE");
            bookCut("MET/sqrtHT > 14: SRE");
            bookCut("jet08[0].mass > 120: SRF");
            bookCut("jet08[1].mass > 60: SRF");
            bookCut("mT_minb > 175: SRF");
            bookCut("MET > 250: SRF");
            bookCut("HT > 1100: SRF");
            bookCut("MET/sqrtHT > 15: SRF");
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
            const Particles& jets08 = applyProjection<NearIsoParticle>(event, "Jets08").particlesByPt();
            const Particles& jets12 = applyProjection<NearIsoParticle>(event, "Jets12").particlesByPt();

            const Particles& taus = applyProjection<NearIsoParticle>(event, "TAUJets").particlesByPt();      
            const Particles& leps = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt();      
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();
            // Projection application section -- do not edit/remove this comment
            /// @todo apply projections
            // const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);

            int njet= jets.size();
            double ht = 0;
            for (int i=0; i <njet; i++){
                ht += jets[i].pT();
            }
            double sqrtht= sqrt(ht);

            bool basecut0 = false;
            bool delphicut = false;
            if (njet > 3 and bjets.size() > 0){
                if (jets[0].pT() > 80 and jets[1].pT() > 80 and jets[2].pT() > 40 and jets[3].pT() > 40 and bjets[0].pT() >= jets[3].pT()){
                    basecut0 = true;
                }
                if (abs(deltaPhi(jets[0].momentum(),met)) > 0.4 and abs(deltaPhi(jets[1].momentum(),met)) > 0.4){
                    delphicut = true;
                }
            }

            // deltaphi(jet,MET_track) cut not implemented yet
            if (!cut(basecut0,"jet cuts: base")) vetoEvent;
            if (!cut(delphicut,"deltaphi cuts: base")) vetoEvent;
            if (!cut(leps.size(),CUT_EQ,0,"0 lepton: base")) vetoEvent;
            if (!cut(MET,CUT_GT,250,"MET > 250: base")) vetoEvent;

            
            // SRC and SRD too difficult to implement
            if (jets12.size() > 1){
                //SRA/B
                if (cut(jets12[0].momentum().mass(),CUT_GT,120,"jet12[0].mass > 120: TT-base")){
                    if (cut(jets12[1].momentum().mass(),CUT_GT,120,"jet12[1].mass > 120: TT-base")){
                        if (cut(jets08[0].momentum().mass(),CUT_GT,60,"jet08[0].mass > 60: SRA-TT")){
                            if (cut(bjets.size(),CUT_GE,2,"bjets.size() >= 2: SRA-TT")){
                                double delphibmin = 10000;
                                double mT_minb = 10000;
                                for(int i=0; i< bjets.size(); i++){
                                double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                                double dm = get_mT(bjets[i].momentum(), met);
                                if(delphib < delphibmin) {
                                    delphibmin = delphib;
                                    mT_minb = dm;
                                }
                                 }
                                if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRA-TT")){
                                    if (cut(taus.size(),CUT_EQ,0,"tau veto: SRA-TT")){
                                        if(cut(MET,CUT_GT,400, "MET > 400: SRA-TT")){
                                            pass("SRA-TT");
                                        }
                                    }
                                }
                            }

                        }

                        if (cut(bjets.size(),CUT_GE,2,"bjets >= 2: SRB-TT")){
                            double delphibmin = 10000;
                            double mT_minb = 10000;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib < delphibmin) {
                                delphibmin = delphib;
                                mT_minb = dm;
                            }
                             }
                            if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRB-TT")){
                                double delphibmax = 0;
                            double mT_maxb = 0;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib > delphibmax) {
                                delphibmax = delphib;
                                mT_maxb = dm;
                            }
                             }
                            if(cut(mT_maxb, CUT_GT , 200, "mT_maxb > 200: SRB-TT")){
                                if(cut(taus.size(),CUT_EQ,0,"tau veto: SRB-TT")){
                                    if(cut(deltaR(bjets[0],bjets[1]),CUT_GT,1.2,"deltaR_b > 1.2: SRB-TT")){
                                        if(cut(MET,CUT_GT,250,"MET > 250: SRB-TT")){
                                            pass("SRB-TT");
                                        }
                                    }
                                }
                            }
                            }
                        }
                    }
                }

                if (cut(jets12[0].momentum().mass(),CUT_GT,120,"jet12[0].mass > 120: TW-base")){
                    if (cut(jets12[1].momentum().mass(),CUT_IN,make_pair(60,120),"jet12[1].mass [60,120]: TW-base")){
                        if (cut(jets08[0].momentum().mass(),CUT_GT,60,"jet08[0].mass > 60: SRA-TW")){
                            if (cut(bjets.size(),CUT_GE,2,"bjets.size() >= 2: SRA-TW")){
                                double delphibmin = 10000;
                                double mT_minb = 10000;
                                for(int i=0; i< bjets.size(); i++){
                                double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                                double dm = get_mT(bjets[i].momentum(), met);
                                if(delphib < delphibmin) {
                                    delphibmin = delphib;
                                    mT_minb = dm;
                                }
                                 }
                                if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRA-TW")){
                                    if (cut(taus.size(),CUT_EQ,0,"tau veto: SRA-TW")){
                                        if(cut(MET,CUT_GT,450, "MET > 450: SRA-TW")){
                                            pass("SRA-TW");
                                        }
                                    }
                                }
                            }

                        }
                        if (cut(bjets.size(),CUT_GE,2,"bjets >= 2: SRB-TW")){
                            double delphibmin = 10000;
                            double mT_minb = 10000;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib < delphibmin) {
                                delphibmin = delphib;
                                mT_minb = dm;
                            }
                             }
                            if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRB-TW")){
                                double delphibmax = 0;
                            double mT_maxb = 0;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib > delphibmax) {
                                delphibmax = delphib;
                                mT_maxb = dm;
                            }
                             }
                            if(cut(mT_maxb, CUT_GT , 200, "mT_maxb > 200: SRB-TW")){
                                if(cut(taus.size(),CUT_EQ,0,"tau veto: SRB-TW")){
                                    if(cut(deltaR(bjets[0],bjets[1]),CUT_GT,1.2,"deltaR_b > 1.2: SRB-TW")){
                                        if(cut(MET,CUT_GT,250,"MET > 250: SRB-TW")){
                                            pass("SRB-TW");
                                        }
                                    }
                                }
                            }
                            }
                        }
                    }
                }
                if (cut(jets12[0].momentum().mass(),CUT_GT,120,"jet12[0].mass > 120: T0-base")){
                    if (cut(jets12[1].momentum().mass(),CUT_LT,60,"jet12[1].mass < 60: T0-base")){
                        if (cut(jets08[0].momentum().mass(),CUT_GT,60,"jet08[0].mass > 60: SRA-T0")){
                            if (cut(bjets.size(),CUT_GE,2,"bjets.size() >= 2: SRA-T0")){
                                double delphibmin = 10000;
                                double mT_minb = 10000;
                                for(int i=0; i< bjets.size(); i++){
                                double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                                double dm = get_mT(bjets[i].momentum(), met);
                                if(delphib < delphibmin) {
                                    delphibmin = delphib;
                                    mT_minb = dm;
                                }
                                 }
                                if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRA-T0")){
                                    if (cut(taus.size(),CUT_EQ,0,"tau veto: SRA-T0")){
                                        if(cut(MET,CUT_GT,400, "MET > 400: SRA-T0")){
                                            pass("SRA-T0");
                                        }
                                    }
                                }
                            }

                        }

                        if (cut(bjets.size(),CUT_GE,2,"bjets >= 2: SRB-T0")){
                            double delphibmin = 10000;
                            double mT_minb = 10000;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib < delphibmin) {
                                delphibmin = delphib;
                                mT_minb = dm;
                            }
                             }
                            if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRB-T0")){
                                double delphibmax = 0;
                            double mT_maxb = 0;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib > delphibmax) {
                                delphibmax = delphib;
                                mT_maxb = dm;
                            }
                             }
                            if(cut(mT_maxb, CUT_GT , 200, "mT_maxb > 200: SRB-T0")){
                                if(cut(taus.size(),CUT_EQ,0,"tau veto: SRB-T0")){
                                    if(cut(deltaR(bjets[0],bjets[1]),CUT_GT,1.2,"deltaR_b > 1.2: SRB-T0")){
                                        if(cut(MET,CUT_GT,500,"MET > 500: SRB-T0")){
                                            pass("SRB-T0");
                                        }
                                    }
                                }
                            }
                            }
                        }
                    }
                }

                //SR E/F
                if (cut(bjets.size(),CUT_GE,2,"bjets >= 2: SREF base")){

                    if (cut(jets12[0].momentum().mass(),CUT_GT,140,"jet12[0].mass > 140: SRE")){
                        if (cut(jets12[1].momentum().mass(),CUT_GT,60,"jet12[1].mass > 60: SRE")){
                            double delphibmin = 10000;
                            double mT_minb = 10000;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)) ;
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib < delphibmin) {
                                delphibmin = delphib;
                                mT_minb = dm;
                            }
                             }
                            if (cut(mT_minb, CUT_GT , 200, "mT_minb > 200: SRE")){
                                if(cut(taus.size(),CUT_EQ,0,"tau veto: SRE")){
                                    if(cut(deltaR(bjets[0],bjets[1]),CUT_GT,1.5,"deltaR_b > 1.5: SRE")){
                                        if(cut(MET,CUT_GT,300,"MET > 300: SRE")){
                                            if (cut(MET/sqrtht,CUT_GT,14,"MET/sqrtHT > 14: SRE")){
                                                pass("SRE");

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(jets08.size()>1){
                        if (cut(jets08[0].momentum().mass(),CUT_GT,120,"jet08[0].mass > 120: SRF")){
                            if (cut(jets08[1].momentum().mass(),CUT_GT,60,"jet08[1].mass > 60: SRF")){
                            double delphibmin = 10000;
                            double mT_minb = 10000;
                            for(int i=0; i< bjets.size(); i++){
                            double delphib = abs(deltaPhi(jets[0].momentum(),met)); 
                            double dm = get_mT(bjets[i].momentum(), met);
                            if(delphib < delphibmin) {
                                delphibmin = delphib;
                                mT_minb = dm;
                            }
                             }
                            if (cut(mT_minb, CUT_GT , 175, "mT_minb > 175: SRF")){
                                if(cut(MET,CUT_GT,250,"MET > 250: SRF")){
                                    if(cut(ht,CUT_GT,1100,"HT > 1100: SRF")){
                                        if (cut(MET/sqrtht,CUT_GT,15,"MET/sqrtHT > 15: SRF")){
                                            pass("SRF");
                                        }
                                    }
                                }
                            }
                            }
                        }
                    }
                }
            }


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
    AtomPlugin(ATLAS_CONF_2016_077)
}
