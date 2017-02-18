///
/// @file  ATLAS_1308_1841.cc
/// @brief Implementation of ATLAS_1308_1841 analysis
/// @author Kazuki
/// @date created 04/25/2015
/// @date last revision 04/25/2015
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1308_1841 : public Analysis {
	public:

		ATLAS_1308_1841()
			: Analysis("ATLAS_1308_1841") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

		/// Book histograms and initialise projections before the run
		void init() {

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2013" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT, 20., 8000.) & Range(ETA, -4.5, 4.5), 
                            muDetRange, FastJets::ANTIKT, 0.4 );
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            IsoElectron base_ele( Range(PT, 10., 8000.) & Range(ETA, -2.47, 2.47) );
            base_ele.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            base_ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Medium_2012_ATLAS" ) );

            IsoMuon base_mu(Range(PT, 10., 8000.) & Range(ETA, -2.5, 2.5));
            base_mu.addIso(TRACK_ISO_PT, 0.01,  1.0,  0.0, 0.01, CALO_ALL);
            base_mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            base_mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            // Overlap removal

            MergedFinalState base_leptons(base_ele, base_mu);
            NearIsoParticle base_lep_clean(base_leptons);
            base_lep_clean.addFilter(jets, 0.4);
            addProjection(base_lep_clean, "Leptons");

            // SmearingParams& metsmear = metSim("Smear_MissingET_ATLAS");
            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets, base_lep_clean);
            MissingMomentum met( fsbase, met_seed );
            //MissingMomentum met( fsbase );
            //met.setSmearingParams( FinalState::SELECTED, &dp.metEff( "Jet_PlaceHolder" ) );
            addProjection(met, "MissingEt");

            NearIsoParticle jets50_eta2(jets, Range(PT, 50., 8000.) & Range(ETA, -2., 2.));
            addProjection(jets50_eta2, "Jets50_eta2");

            NearIsoParticle jets80_eta2(jets, Range(PT, 80., 8000.) & Range(ETA, -2., 2.));
            addProjection(jets80_eta2, "Jets80_eta2");

            NearIsoParticle jets40_eta28(jets, Range(PT, 40., 8000.) & Range(ETA, -2.8, 2.8));
            addProjection(jets40_eta28, "Jets40_eta28");

            NearIsoParticle jets50_eta28(jets, Range(PT, 50., 8000.) & Range(ETA, -2.8, 2.8));
            addProjection(jets50_eta28, "Jets50_eta28");

            HeavyFlavorJets bjets(jets, Range(PT, 40., 8000.) & Range(ETA, -2.5, 2.5));
            addProjection(bjets, "BJets");

            NearIsoParticle megajet_seeds(jets, Range(ETA, -2.8, 2.8));
            FastJets megajets(megajet_seeds, Range(PT, 100., 8000.) & Range(ETA, -1.5, 1.5), FastJets::ANTIKT, 1.0 );
            addProjection(megajets, "MegaJets");

            /// @book the efficiencies
            bookEfficiency("FS_8j50_b0");
            bookEfficiency("FS_8j50_b1");
            bookEfficiency("FS_8j50_b2");
            bookEfficiency("FS_9j50_b0");
            bookEfficiency("FS_9j50_b1");
            bookEfficiency("FS_9j50_b2");
            bookEfficiency("FS_10j50");

            bookEfficiency("FS_7j80_b0");
            bookEfficiency("FS_7j80_b1");
            bookEfficiency("FS_7j80_b2");
            bookEfficiency("FS_8j80_b0");
            bookEfficiency("FS_8j80_b1");
            bookEfficiency("FS_8j80_b2");

            bookEfficiency("MS_8j50_340");
            bookEfficiency("MS_8j50_420");
            bookEfficiency("MS_9j50_340");
            bookEfficiency("MS_9j50_420");
            bookEfficiency("MS_10j50_340");
            bookEfficiency("MS_10j50_420");

            /// @book the cuts
	        bookCut("Lepton Veto");
            bookCut("Nj50 = 8: FS_8j50");
            bookCut("MET/sqrt(HT) > 4: FS_8j50");
            bookCut("Nb = 0: FS_8j50");
            bookCut("Nb = 1: FS_8j50");
            bookCut("Nb >= 2: FS_8j50");
            bookCut("Nj50 = 9: FS_9j50");
            bookCut("MET/sqrt(HT) > 4: FS_9j50");
            bookCut("Nb = 0: FS_9j50");
            bookCut("Nb = 1: FS_9j50");
            bookCut("Nb >= 2: FS_9j50");
            bookCut("Nj50 >= 10: FS_10j50");
            bookCut("MET/sqrt(HT) > 4: FS_10j50");
            bookCut("Nj80 = 7: FS_7j80");
            bookCut("MET/sqrt(HT) > 4: FS_7j80");
            bookCut("Nb = 0: FS_7j80");
            bookCut("Nb = 1: FS_7j80");
            bookCut("Nb >= 2: FS_7j80");
            bookCut("Nj80 >= 8: FS_8j80");
            bookCut("MET/sqrt(HT) > 4: FS_8j80");
            bookCut("Nb = 0: FS_8j80");
            bookCut("Nb = 1: FS_8j80");
            bookCut("Nb >= 2: FS_8j80");
            bookCut("Nj50 >= 8: MS_8j50");
            bookCut("MET/sqrt(HT) > 4: MS_8j50");
            bookCut("MjS > 340: MS_8j50");
            bookCut("MjS > 420: MS_8j50");
            bookCut("Nj50 >= 9: MS_9j50");
            bookCut("MET/sqrt(HT) > 4: MS_9j50");
            bookCut("MjS > 340: MS_9j50");
            bookCut("MjS > 420: MS_9j50");
            bookCut("Nj50 >= 10: MS_10j50");
            bookCut("MET/sqrt(HT) > 4: MS_10j50");
            bookCut("MjS > 340: MS_10j50");
            bookCut("MjS > 420: MS_10j50");

		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

            const Particles& jets50_eta2   = applyProjection<NearIsoParticle>(event, "Jets50_eta2").particlesByPt(&event);
            const Particles& jets80_eta2   = applyProjection<NearIsoParticle>(event, "Jets80_eta2").particlesByPt(&event);
            const Particles& jets40_eta28  = applyProjection<NearIsoParticle>(event, "Jets40_eta28").particlesByPt(&event);            
            const Particles& jets50_eta28  = applyProjection<NearIsoParticle>(event, "Jets50_eta28").particlesByPt(&event);
            const Particles& megajets     = applyProjection<FastJets>(event, "MegaJets").particlesByPt(&event);
            const Particles& leps          = applyProjection<NearIsoParticle>(event, "Leptons").particlesByPt(&event);
            const HeavyFlavorJets& bjproj  = applyProjection<HeavyFlavorJets>(event, "BJets");
            const Particles& bjets         = bjproj.getTaggedJets(&event); 
            const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet    = applyProjection<MissingMomentum>(event, "MissingEt");

            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();

            //=============================//
            //         preselecton         //
            //=============================//

            if(!cut(leps.size(), CUT_EQ, 0, "Lepton Veto")) vetoEvent;

            //=============================//
            //    Variable calculation     //
            //=============================//

            int Nb = bjets.size();

            double HT = 0.;
            for(int i=0; i<jets40_eta28.size(); i++){
                HT += jets40_eta28[i].pT();
            }

            double MjS = 0.;
            for(int i=0; i<megajets.size(); i++){
                MjS += megajets[i].mass();
            }

            //=================================//
            //    Multijet + flavour stream    //
            //=================================//

            int Nj50_eta2 = jets50_eta2.size();

            if(cut(Nj50_eta2, CUT_EQ, 8, "Nj50 = 8: FS_8j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: FS_8j50") ){
                    if(cut(Nb, CUT_EQ, 0, "Nb = 0: FS_8j50")) pass("FS_8j50_b0");
                    if(cut(Nb, CUT_EQ, 1, "Nb = 1: FS_8j50")) pass("FS_8j50_b1");
                    if(cut(Nb, CUT_GE, 2, "Nb >= 2: FS_8j50")) pass("FS_8j50_b2");
                }
            }

            if(cut(Nj50_eta2, CUT_EQ, 9, "Nj50 = 9: FS_9j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: FS_9j50") ){
                    if(cut(Nb, CUT_EQ, 0, "Nb = 0: FS_9j50")) pass("FS_9j50_b0");
                    if(cut(Nb, CUT_EQ, 1, "Nb = 1: FS_9j50")) pass("FS_9j50_b1");
                    if(cut(Nb, CUT_GE, 2, "Nb >= 2: FS_9j50")) pass("FS_9j50_b2");
                }
            }

            if(cut(Nj50_eta2, CUT_GE, 10, "Nj50 >= 10: FS_10j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: FS_10j50") ){
                    pass("FS_10j50");
                }
            }

            int Nj80_eta2 = jets80_eta2.size();

            if(cut(Nj80_eta2, CUT_EQ, 7, "Nj80 = 7: FS_7j80") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: FS_7j80") ){
                    if(cut(Nb, CUT_EQ, 0, "Nb = 0: FS_7j80")) pass("FS_7j80_b0");
                    if(cut(Nb, CUT_EQ, 1, "Nb = 1: FS_7j80")) pass("FS_7j80_b1");
                    if(cut(Nb, CUT_GE, 2, "Nb >= 2: FS_7j80")) pass("FS_7j80_b2");
                }
            }

            if(cut(Nj80_eta2, CUT_GE, 8, "Nj80 >= 8: FS_8j80") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: FS_8j80") ){
                    if(cut(Nb, CUT_EQ, 0, "Nb = 0: FS_8j80")) pass("FS_8j80_b0");
                    if(cut(Nb, CUT_EQ, 1, "Nb = 1: FS_8j80")) pass("FS_8j80_b1");
                    if(cut(Nb, CUT_GE, 2, "Nb >= 2: FS_8j80")) pass("FS_8j80_b2");
                }
            }

            //=================================//
            //       Multijet + MjS stream     //
            //=================================//

            int Nj50_eta28 = jets50_eta28.size();

            if(cut(Nj50_eta28, CUT_GE, 8, "Nj50 >= 8: MS_8j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: MS_8j50") ){
                    if(cut( MjS, CUT_GT, 340., "MjS > 340: MS_8j50")) pass("MS_8j50_340");
                    if(cut( MjS, CUT_GT, 420., "MjS > 420: MS_8j50")) pass("MS_8j50_420");
                }
            }

            if(cut(Nj50_eta28, CUT_GE, 9, "Nj50 >= 9: MS_9j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: MS_9j50") ){
                    if(cut( MjS, CUT_GT, 340., "MjS > 340: MS_9j50")) pass("MS_9j50_340");
                    if(cut( MjS, CUT_GT, 420., "MjS > 420: MS_9j50")) pass("MS_9j50_420");
                }
            }

            if(cut(Nj50_eta28, CUT_GE, 10, "Nj50 >= 10: MS_10j50") ){
                if(cut(MET/sqrt(HT), CUT_GT, 4.0, "MET/sqrt(HT) > 4: MS_10j50") ){
                    if(cut( MjS, CUT_GT, 340., "MjS > 340: MS_10j50")) pass("MS_10j50_340");
                    if(cut( MjS, CUT_GT, 420., "MjS > 420: MS_10j50")) pass("MS_10j50_420");
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
	AtomPlugin(ATLAS_1308_1841)
}
