///
/// @file  ATLAS_CONF_2016_066.cc
/// @brief Implementation of ATLAS_CONF_2016_066 analysis
/// @author Seng Pei Liew
/// @date created 12/11/2016
/// @date last revision 12/11/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

    class ATLAS_CONF_2016_066 : public Analysis {
    public:

        ATLAS_CONF_2016_066()
            : Analysis("ATLAS_CONF_2016_066") {
        //  setNeedsCrossSection　not working?
    //        setNeedsCrossSection(true);
        }

        /// @name Analysis methods
        //@{

        /// Book histograms and initialise projections before the run
        void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

            // Projection booking section -- do not edit/remove this comment
            useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );
            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            // fix isolation and filter later

            IsoPhoton gam( Range(PT > 25.) & Range(abseta, 0.0, 2.37) - Range(abseta, 1.37, 1.56) );
          //  gam.addIso(CALO_ISO_ET,  0.4,  0.0, 6, 0.0, CALO_ALL);                        
          //  gam.addIso(TRACK_ISO_PT, 0.2,  0.0,  2.6, 0.0, CALO_ALL, 1.0);
            gam.addConeIso(CALO_ISO_ET,  0.4,  0.0, 6, 0.01, CALO_ALL);                        
            gam.addConeIso(TRACK_ISO_PT, 0.2,  0.0,  2.6, 0.01, CALO_ALL, 1.0);
            gam.setSmearingParams( getPhotonSim( "Photon_Smear_2013_ATLAS" ) );
            gam.setEfficiencyParams( getPhotonEff( "Photon_Ident_Tight_2011_ATLAS" ) );

            IsoElectron ele( Range(PT > 10. & abseta < 2.47) );
            ele.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
            ele.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );
            IsoMuon mu( Range(PT > 10. & abseta < 2.5) );
            mu.addConeIso(TRACK_ISO_PT, 0.2,  0.15,  0.0, 0.01, CALO_ALL);
            mu.addConeIso(CALO_ISO_PT,  0.3,  0.15,  0.0, 0.01, CALO_ALL);            
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );
            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            
            FastJets jets(fsbase, 
                    hadRange & Range( PT > 30 & abseta < 2.5 ), 
                    muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );
            // Overlap removal
            NearIsoParticle jets_clean(jets);
            jets_clean.addFilter(ele, 0.2);
            
            addProjection(gam,"Photon");
      //      MergedFinalState lep_base(ele, mu);
            NearIsoParticle electrons(ele);
     //       
            NearIsoParticle muons(mu);
            addProjection(electrons,     "Electrons");
            addProjection(muons,     "Muons");

            MergedFinalState lep_base(ele, mu);
            NearIsoParticle leptons(lep_base);
         //   leptons.addFilter(jets_clean, 0.4);
            addProjection(leptons,     "Leptons");

            addProjection(jets_clean, "Jets");
         //   addProjection(bjets,       "BJets");
            MergedFinalState met_seed(jets_clean, leptons);
            MissingMomentum met( fsbase, met_seed );
            met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            // Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

            // Efficiency booking section -- do not edit/remove this comment
            bookEfficiency("SRL");
            bookEfficiency("SRH");

            bookCut("Nphoton > 0: base");

            bookCut("photon[0].pt() > 145: SRL");
            bookCut("Nleps = 0: SRL");
            bookCut("Njet > 4: SRL");
            bookCut("delta(j,met) > 0.4: SRL");
            bookCut("delta(photon,met) > 0.4: SRL");
            bookCut("MET > 200: SRL");
            bookCut("meff > 2000: SRL");
            bookCut("RT4 < 0.9 : SRL");


            bookCut("photon[0].pt() > 400: SRH");
            bookCut("Nleps = 0: SRH");
            bookCut("Njet > 2: SRH");
            bookCut("delta(j,met) > 0.4: SRH");
            bookCut("delta(photon,met) > 0.4: SRH");
            bookCut("MET > 400: SRH");
            bookCut("meff > 2000: SRH");
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

            // Projection application section -- do not edit/remove this comment
            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
            const Particles& eles = applyProjection<NearIsoParticle>(event, "Electrons").particlesByPt();      
            const Particles& mus = applyProjection<NearIsoParticle>(event, "Muons").particlesByPt();      
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();      
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(); 
            const Particles& photon = applyProjection<IsoPhoton>(event, "Photon").particlesByPt();
           // const Particles& Ljets = bjproj.getUntaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();

            // Analysis body section -- do not edit/remove this comment
            if (!cut(photon.size(),CUT_GT, 0,"Nphoton > 0: base")) vetoEvent;

         //   if (photon.size()==0) vetoEvent;
         //   if (photon[0].pT() <= 145) vetoEvent;
         //   if (leps.size() > 0 ) vetoEvent;
         //   if (jets.size() < 3) vetoEvent;
         //   if (MET <= 200) vetoEvent;


            double dPhimin= 1000;
            for(int i=0; i < 2; i++){
                double dPhi= deltaPhi(jets[i].momentum(), met);
                if (dPhi < dPhimin) dPhimin = dPhi;
                double dPhip= deltaPhi(photon[0].momentum(), met);
                if (dPhip < dPhimin) dPhimin = dPhip;
            }

         //   if (dPhimin <= 0.4) vetoEvent;

            double HT= photon[0].pT();
            for(int i=0; i < jets.size(); i++){
                HT += jets[i].pT();
            }

            double meff= HT+MET;

        //    if (meff <= 2000) vetoEvent;


            /// @todo fill efficiencies

            // SRL 
            if (cut(photon[0].pT(),CUT_GT,145,"photon[0].pt() > 145: SRL")){
                if (cut(leps.size(),CUT_EQ,0,"Nleps = 0: SRL")){
                    if (cut(jets.size(),CUT_GT,4,"Njet > 4: SRL")){
                            double dPhiminj= 1000;
                                for(int i=0; i < 2; i++){
                                    double dPhi= deltaPhi(jets[i].momentum(), met);
                                    if (dPhi < dPhiminj) dPhiminj = dPhi;
                            }
                            double dPhiminp= 1000;
                            double dPhip= deltaPhi(photon[0].momentum(), met);
                                if (dPhip < dPhiminp) dPhiminp = dPhip;                        
                        if (cut(dPhiminj,CUT_GT,0.4,"delta(j,met) > 0.4: SRL")){
                            if (cut(dPhiminp,CUT_GT,0.4,"delta(photon,met) > 0.4: SRL")){
                                if (cut(MET,CUT_GT,200,"MET > 200: SRL")){
                                    if (cut(meff,CUT_GT,2000,"meff > 2000: SRL")){
                                        double RT4 = 0;
                                           double ptj4 = 0;
                                           double ptjall = 0;
                                           for(int i=0; i < 4; i++){
                                           ptj4 += jets[i].pT();
                                       }
                                           for(int i=0; i < jets.size(); i++){
                                           ptjall += jets[i].pT();
                                       }
                                           RT4 = ptj4/ptjall;
                                        if (cut(RT4,CUT_LT,0.9,"RT4 < 0.9 : SRL")){
                                            pass ("SRL");
                                        }
                                    }
                                }
                            }
                        } 
                    }
                }
            }
            


          if (cut(photon[0].pT(),CUT_GT,400,"photon[0].pt() > 400: SRH")){
                if (cut(leps.size(),CUT_EQ,0,"Nleps = 0: SRH")){
                    if (cut(jets.size(),CUT_GT,4,"Njet > 4: SRH")){
                            double dPhiminj= 1000;
                                for(int i=0; i < 2; i++){
                                    double dPhi= deltaPhi(jets[i].momentum(), met);
                                    if (dPhi < dPhiminj) dPhiminj = dPhi;
                            }
                            double dPhiminp= 1000;
                            double dPhip= deltaPhi(photon[0].momentum(), met);
                                if (dPhip < dPhiminp) dPhiminp = dPhip;                        
                        if (cut(dPhiminj,CUT_GT,0.4,"delta(j,met) > 0.4: SRH")){
                            if (cut(dPhiminp,CUT_GT,0.4,"delta(photon,met) > 0.4: SRH")){
                                if (cut(MET,CUT_GT,400,"MET > 400: SRH")){
                                    if (cut(meff,CUT_GT,2000,"meff > 2000: SRH")){
                                        
                                            pass("SRH");
                                        
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
    AtomPlugin(ATLAS_CONF_2016_066)
}
