///
/// @file  ATLAS_1410_7238.cc
/// @brief Implementation of ATLAS_1410_7238 analysis
/// @author Kazuki
/// @date created 05/17/2016
/// @date last revision 05/17/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class ATLAS_1410_7238 : public Analysis {
	public:

		ATLAS_1410_7238()
			: Analysis("ATLAS_1410_7238") {
			setNeedsCrossSection(true);
		}

		/// @name Analysis methods
		//@{

        double get_mT(FourMomentum p1, FourMomentum p2){
            double mTsq = 2. * ( p1.pT()*p2.pT() - p1.px()*p2.px() - p1.py()*p2.py() );
            double mT = sqrt(mTsq);
            return mT;
        }


        double n75_80, n80_85;

		/// Book histograms and initialise projections before the run
		void init() {

            n75_80 = 0; n80_85 = 0;

            useDetector( "ATLAS_CMS_all" ); // "ATLAS2014" );

            FinalState fsbase( getRange( "Full_Range_ATLAS" ) );

            IsoElectron ele( Range(PT > 25.) & (Range(abseta < 2.47) - Range(abseta, 1.37, 1.52)) );
            //                         cone  frac  abs  inner 
            ele.addIso(CALO_ISO_ET,  0.3,  0.14,  0.0, 0.01, CALO_ALL);
            ele.addIso(TRACK_ISO_PT, 0.3,  0.13,  0.0, 0.01, CALO_ALL);            
            ele.setVariableThreshold(0.0);
            ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
            ele.setEfficiencyParams( getElectronEff( "Electron_Ident_Tight_2012_ATLAS" ) );

            IsoMuon mu(Range(PT > 25.) & Range(abseta < 2.4));
            mu.addIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.01, CALO_ALL);
            mu.addIso(CALO_ISO_ET, 0.3,  0.14,  0.0, 0.01, CALO_ALL);            
            mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
            mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

            // IsoMuon mu_0(Range(PT > 25.) & Range(abseta < 2.4));
            // mu_0.addIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.01, CALO_ALL);
            // mu_0.addIso(CALO_ISO_ET, 0.3,  0.14,  0.0, 0.01, CALO_ALL);            
            // mu_0.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );
            // addProjection(mu_0, "Muons_uns");

            Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
            Range hadRange   = getRange( "HCal_Range_ATLAS" );
            FastJets jets(fsbase, 
                            hadRange & Range(PT > 25.) & Range(abseta < 2.8), 
                            muDetRange, FastJets::ANTIKT, 0.4 );            
            //jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
            jets.setSmearingParams( getJetSim( "Jet_Smear_LCW_JES_8tev_ATLAS" ) );            
            jets.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

            // FastJets jets_0(fsbase, 
            //                 hadRange & Range(PT > 25.) & Range(abseta < 2.8), 
            //                 muDetRange, FastJets::ANTIKT, 0.4 );
            // addProjection(jets_0, "Jets_unsmear");

            addProjection(ele, "Electrons");
            addProjection(mu, "Muons");
            addProjection(jets, "Jets");

            MergedFinalState leptons(ele, mu);
            addProjection(leptons, "Leptons");

            // FastSimParameterization metsim = createFastSimParam(metsmear);
            MergedFinalState met_seed(jets, leptons);            
            MissingMomentum met( fsbase, met_seed );            
            //MissingMomentum met( fsbase );                        
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Const_PlaceHolder" ) );
            //met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Interpolation_PlaceHolder" ) );            
            addProjection(met, "MissingEt");

			// Projection booking section -- do not edit/remove this comment
            /// @todo define projections (see examples and manual)
			// FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
			// addProjection(jets, "Jets");

			// Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
            //bookHisto1D(1,1,1, "Mjj");
            bookHisto1D("Mjj", 45, 25, 250, "Mjj", "Mjj [GeV]", "Arbitrary Unit");            
            bookHisto1D("R_pT", 30, 0, 2, "R_pT", "pT/pTtrue", "Arbitrary Unit");            

            bookCut("Nlep > 0");
            bookCut("MET > 30");
            bookCut("mT > 40");
            bookCut("Njet == 2");
            bookCut("dphi_j_met > 0.8");
            bookCut("jet requirement");
            bookCut("dEta_jj < 1.5");
            bookCut("dRjj > 0.7 if pTjj < 250");
            bookCut("25 < mjj < 250");
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

            const Particles& jets = applyProjection<FastJets>(event, "Jets").particlesByPt();
            //const Particles& jets_uns = applyProjection<FastJets>(event, "Jets_unsmear").particlesByPt();            
            const Particles& mus = applyProjection<IsoMuon>(event, "Muons").particlesByPt();
            //const Particles& mus_uns = applyProjection<IsoMuon>(event, "Muons_uns").particlesByPt();            
            const Particles& eles = applyProjection<IsoElectron>(event, "Electrons").particlesByPt();
            const Particles& leps = applyProjection<MergedFinalState>(event, "Leptons").particlesByPt();
            //const HeavyFlavorJets& bjproj = applyProjection<HeavyFlavorJets>(event, "BJets");
            //const Particles& bjets = bjproj.getTaggedJets(&event); 
            //const Particles& untagged_jets = bjproj.getUntaggedJets(&event); 
            const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
            const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero
            double MET = met.pT();


            //if( jets.size() > 0 ) cout << jets[0].pT() << endl;

            //vetoEvent;
            //cout << jets.size() <<"  "<< mus.size() <<"  "<< eles.size() <<"  "<< leps.size() <<"  "<< met.pT() << endl;

			// Projection application section -- do not edit/remove this comment
            /// @todo apply projections
			// const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);

			// Analysis body section -- do not edit/remove this comment
            /// @todo apply cuts

           // // cout << "here 1" << endl;
           //  if( mus.size() > 0 && mus_uns.size() > 0 ){
           //      cout << endl;                
           //      for(int i=0; i<mus.size(); i++){
           //          for(int j=0; j<mus_uns.size(); j++){
           //              ///cout << "here 2" << endl;

           //              double dR = deltaR( mus[i].momentum(), mus_uns[j].momentum() );
           //              //cout << "here 3" << endl;

           //              double R_pT = mus[i].pT() / mus_uns[j].pT();
           //              //cout << "here 4" << endl;

           //              if( dR < 1. ){
           //                  //cout << "here 5" << endl;

           //                  cout << dR <<"  "<< R_pT << endl;
           //              }
           //          }
           //      }
           //  }
           //  vetoEvent;

            // if( jets.size() > 0 && jets_uns.size() > 0 ){
            //     cout << endl;                
            //     for(int i=0; i<jets.size(); i++){
            //         for(int j=0; j<jets_uns.size(); j++){
            //             double dR = deltaR( jets[i].momentum(), jets_uns[j].momentum() );
            //             double R_pT = jets[i].pT() / jets_uns[j].pT();
            //             if( dR < 1. ){
            //                 cout << dR <<"  "<< R_pT << endl;
            //             }
            //         }
            //     }
            // }
            // vetoEvent;

            if( cut( leps.size(), CUT_GT, 0, "Nlep > 0" ) ){                       
                if( cut( MET, CUT_GT, 30., "MET > 30" ) ){
                    double mT = get_mT( leps[0].momentum(), met );
                    if( cut( mT, CUT_GT, 40., "mT > 40" ) ){
                        if( cut( jets.size(), CUT_EQ, 2, "Njet == 2") ){
                            double dphi_j_met = deltaPhi(jets[0].momentum(), met);
                            if( cut( dphi_j_met, CUT_GT, 0.8, "dphi_j_met > 0.8") ){
                                bool jet_pass = false;
                                if( jets[0].pT() > 30. && jets[0].abseta() < 2.0 && jets[1].abseta() < 2.0 ) 
                                    jet_pass = true;
                                    if( cut( jet_pass, "jet requirement") ){
                                    double mjj = (jets[0].momentum() + jets[1].momentum()).mass();
                                    double dEta_jj = fabs(jets[0].eta() - jets[1].eta());
                                    double pT_jj = (jets[0].momentum() + jets[1].momentum()).pT();
                                    double dR_jj = deltaR( jets[0].momentum(), jets[1].momentum() );
                                    if( cut( dEta_jj, CUT_LT, 1.5, "dEta_jj < 1.5") ){
                                        jet_pass = false;
                                        if( dR_jj > 0.7 && pT_jj < 250. ) jet_pass = true;
                                        if( pT_jj > 250. ) jet_pass = true;
                                        if( cut( jet_pass, "dRjj > 0.7 if pTjj < 250") ){
                                            if( cut( mjj, CUT_IN, make_pair(25., 250.), "25 < mjj < 250") ){
                                                fillPlot("Mjj", mjj);
                                                if( 75. < mjj && mjj < 80.) n75_80++; 
                                                if( 80. < mjj && mjj < 85.) n80_85++;                                                 
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }


			/// @todo fill histograms 

			// End analyze section -- do not edit/remove this comment
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			// Histogram normalization section -- do not edit/remove this comment
			/// @todo normalize the histograms
			scale("Mjj");

            cout << "75-80: " << n75_80 << endl;
            cout << "80-85: " << n80_85 << endl;

			// End finalize section -- do not edit/remove this comment
		}

		//@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(ATLAS_1410_7238)
}
