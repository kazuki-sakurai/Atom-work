///
/// @file  ATLAS_CONF_2016_038.cc
/// @brief Implementation of ATLAS_CONF_2016_038 analysis
/// @author Seng Pei Liew
/// @date created 12/10/2016
/// @date last revision 12/10/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

    class ATLAS_CONF_2016_038 : public Analysis {
    public:

        ATLAS_CONF_2016_038()
            : Analysis("ATLAS_CONF_2016_038") {
        //  setNeedsCrossSection　not working?
        //    setNeedsCrossSection(true);
        }

        /// @name Analysis methods
        //@{

        /// Book histograms and initialise projections before the run
        void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );
            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );
            // fix isolation and filter later
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
            HeavyFlavorJets bjets(jets_clean);
            bjets.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
            bjets.setCurrentWorkingPoint( 0.77 );
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
            addProjection(bjets,       "BJets");
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
            bookEfficiency("SRL");
            bookEfficiency("SRH");
            bookEfficiency("SRE");
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
           // const Particles& Ljets = bjproj.getUntaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();
          //  double mZ = 91.2;
            // Analysis body section -- do not edit/remove this comment
            int Nlep = leps.size();
            int Ne = eles.size();
            int Nm = mus.size();
            if (Nlep < 2) vetoEvent;
            if (leps[0].pT() <= 40) vetoEvent;
            if (bjets.size() == 0) vetoEvent;
            if (MET <= 100) vetoEvent;

            bool SFOS = false;
            double minve = -1;
            double minvm = -1;
            double ptee =0;
            double ptmm =0;
            if (Ne > 1) {
               minve= (eles[0].momentum()+eles[1].momentum()).mass();
               ptee= (eles[0].momentum()+eles[1].momentum()).pT();
            }
            if (Nm > 1) {
               minvm= (mus[0].momentum()+mus[1].momentum()).mass();
               ptmm= (mus[0].momentum()+mus[1].momentum()).pT();
            }
            if (minve < 106.2 and minve > 76.2) SFOS = true;
            if (minvm < 106.2 and minvm > 76.2) SFOS = true;

            if (!SFOS) vetoEvent;

            //SR

            if (jets.size() > 5 and jets[0].pT() > 60 and MET > 100) pass("SRL");
            if (jets.size() > 4 and jets[0].pT() > 100 and MET > 100) {
                if (ptee > 200 or ptmm > 200) pass("SRH");
            }
            if (jets.size() > 4 and jets[0].pT() > 80 and MET > 100) {
                if (ptee > 100 or ptmm > 100) pass("SRE");
            }


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
    AtomPlugin(ATLAS_CONF_2016_038)
}
