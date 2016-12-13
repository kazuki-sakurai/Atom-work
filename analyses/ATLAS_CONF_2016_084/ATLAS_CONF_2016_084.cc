///
/// @file  ATLAS_CONF_2016_084.cc
/// @brief Implementation of ATLAS_CONF_2016_084 analysis
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

    class ATLAS_CONF_2016_084 : public Analysis {
    public:

        ATLAS_CONF_2016_084()
            : Analysis("ATLAS_CONF_2016_084") {
            setNeedsCrossSection(true);
        }

        /// @name Analysis methods
        //@{

        /// Book histograms and initialise projections before the run
        void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2016" );

            // Projection booking section -- do not edit/remove this comment
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
            addProjection(leps,     "Leptons");

            addProjection(jets_clean, "Jets");
            addProjection(bjets,       "BJets");
            MergedFinalState met_seed(jets_clean, leptons);
            MissingMomentum met( fsbase, met_seed );
            met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            // Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            // bookHisto1D(1,1,1, "Meff");

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

            // Projection application section -- do not edit/remove this comment
            const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
            const Particles& eles = applyProjection<NearIsoParticle>(event, "Electrons").particlesByPt();      
            const Particles& mus = applyProjection<NearIsoParticle>(event, "Muons").particlesByPt();      
            const Particles& leps = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt();      
            const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets = bjproj.getTaggedJets(); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();
            // Analysis body section -- do not edit/remove this comment
            /// @todo apply cuts
            if (jets.size() < 4) vetoEvent;
            

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
    AtomPlugin(ATLAS_CONF_2016_084)
}
