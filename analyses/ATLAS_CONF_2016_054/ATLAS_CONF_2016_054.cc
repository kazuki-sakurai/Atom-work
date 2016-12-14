///
/// @file  ATLAS_CONF_2016_054.cc
/// @brief Implementation of ATLAS_CONF_2016_054 analysis
/// @author Seng Pei Liew
/// @date created 09/08/2016
/// @date last revision 09/08/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"
//#include "Rivet/Projections/Sphericity.hh"

using namespace std;

namespace Atom {

	class ATLAS_CONF_2016_054 : public Analysis {
	public:

		ATLAS_CONF_2016_054()
			: Analysis("ATLAS_CONF_2016_054") {
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

      Range muDetRange = getRange( "Muon_Range_Detector_ATLAS" );
      Range hadRange   = getRange( "HCal_Range_ATLAS" );

      FastJets jets_pre(fsbase, 
              hadRange & Range( PT > 20 & abseta < 4.5 ), 
              muDetRange, FastJets::ANTIKT, 0.4 );
      jets_pre.setSmearingParams( getJetSim( "Jet_Smear_Topo_ATLAS" ) );
      jets_pre.setEfficiencyParams( getJetEff( "Jet_Ident_PlaceHolder" ) );

      HeavyFlavorJets bjets_pre(jets_pre, Range(PT > 30. & abseta < 2.5));
      bjets_pre.setTaggingParams( getBJetTag("BJet_Ident_MV1_ATLAS") );
      bjets_pre.setCurrentWorkingPoint( 0.77 );
      

      IsoMuon mu_base(Range(PT > 6. & abseta < 2.5));
      mu_base.addConeIso(TRACK_ISO_PT, 0.2,  0.05,  0.0, 0.01, CALO_ALL);            
      mu_base.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu_base.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      IsoElectron ele_base( Range(PT > 7. & abseta < 2.5) );
      ele_base.addConeIso(TRACK_ISO_PT, 0.2,  0.05,  0.0, 0.01, CALO_ALL);      
      ele_base.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );

      IsoMuon mu(Range(PT > 6. & abseta < 2.5));
      mu.addConeIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.01, CALO_ALL,1.);
      mu.addConeIso(CALO_ISO_ET, 0.2,  0.1,  0.0, 0.01, CALO_ALL);
      mu.setSmearingParams  ( getMuonSim( "Muon_Smear_ID-MS_ATLAS" ) );
      mu.setEfficiencyParams( getMuonEff( "Muon_Ident_CB-ST_ATLAS" ) );

      IsoElectron ele( Range(PT > 7. & abseta < 2.5) );
      ele.addConeIso(TRACK_ISO_PT, 0.3,  0.15,  0.0, 0.01, CALO_ALL,1.);
      ele.addConeIso(CALO_ISO_ET, 0.2,  0.1,  0.0, 0.01, CALO_ALL);
      ele.setSmearingParams  ( getElectronSim( "Electron_Smear_run1_ATLAS" ) );
			
      // Overlap removal

      NearIsoParticle ele_pre(ele_base);
      ele_pre.addFilter(jets_pre,0.2);

      NearIsoParticle ele_sig(ele);
      ele_sig.addFilter(jets_pre,0.2);

      NearIsoParticle bjets(bjets_pre);
      bjets.addFilter(mu_base,0.2);
      addProjection(bjets, "BJets");

      NearIsoParticle mu_pre(mu_base);
      mu_pre.addFilter(jets_pre,0.2);
      mu_pre.addFilter(ele_base,0.01);

      MergedFinalState leps_pre(ele_base,mu_base);
      
      NearIsoParticle mu_sig(mu);
      mu_sig.addFilter(jets_pre,0.2);
      mu_sig.addFilter(ele_base,0.01);

      MergedFinalState leps_sig(ele_sig,mu_sig);
      addProjection(leps_sig,"Signal_Leptons");

      NearIsoParticle jets_sig(jets_pre,Range(PT > 30. & abseta < 2.8));
      addProjection(jets_sig,"Jets");

      MergedFinalState met_seed(jets_pre,leps_pre);
      MissingMomentum met( fsbase, met_seed );
      met.setSmearingParams( getMETSim( "MissingET_Smear_ETOnly_Grid_PlaceHolder" ) );
      addProjection(met, "MissingEt");

      // implement aplanarity

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

			// Cuts booking section -- do not edit/remove this comment
      /// @todo book the cuts
      // bookCut("CutNJets");
      // bookCut("CutPTJ1","description goes here");
      // bookCut("CutEtaJet","this is a control region cut", true);
      //================================// 
			bookEfficiency("GG2J");
			bookEfficiency("GG6J-bulk");
			bookEfficiency("GG6J-high-mass");
			bookEfficiency("GG4J-high-x");
      bookEfficiency("GG4J-low-x-b-veto");
			bookEfficiency("GG4J-low-x");
			bookEfficiency("SS4J-x=1/2");
			bookEfficiency("SS5J-x=1/2");
			bookEfficiency("SS4J-low-x");
			bookEfficiency("SS5J-high-x");
     	//================================// 

      bookCut("Nlep = 1 : base");
      bookCut("lepton pt cut: SS base");

      bookCut("lepton pt cut : GG2J");
      bookCut("Njet >= 2 : GG2J");
      bookCut("jets[0].pt() > 200 : GG2J");
      bookCut("jets[1].pt() > 30 : GG2J");
      bookCut("mT > 100 : GG2J");
      bookCut("MET > 460 : GG2J");
      bookCut("MET/minc_eff > 0.35: GG2J");

      bookCut("lepton pt cut : GG6J-bulk");
      bookCut("Njet >= 6 : GG6J-bulk");
      bookCut("jets[0].pt() > 125 : GG6J-bulk");
      bookCut("jet2/3 pt > 30 : GG6J-bulk");
      bookCut("jets[3].pt() > 30 : GG6J-bulk");
      bookCut("jet5/6 pt > 30 : GG6J-bulk");
      bookCut("mT > 225 : GG6J-bulk");
      bookCut("MET > 250 : GG6J-bulk");
      bookCut("meff_inc > 1000: GG6J-bulk");
      bookCut("MET/meff_inc > 0.2: GG6J-bulk");
      bookCut("jet aplanarity > 0.04: GG6J-bulk");

      bookCut("lepton pt cut : GG6J-high-mass");
      bookCut("Njet >= 6 : GG6J-high-mass");
      bookCut("jets[0].pt() > 125 : GG6J-high-mass");
      bookCut("jet2/3 pt > 30 : GG6J-high-mass");
      bookCut("jets[3].pt() > 30 : GG6J-high-mass");
      bookCut("jet5/6 pt > 30 : GG6J-high-mass");
      bookCut("mT > 225 : GG6J-high-mass");
      bookCut("MET > 250 : GG6J-high-mass");
      bookCut("meff_inc > 1000: GG6J-high-mass");
      bookCut("MET/meff_inc > 0.1: GG6J-high-mass");
      bookCut("jet aplanarity > 0.04: GG6J-high-mass");
      bookCut("GG6J-high-mass");

      bookCut("lepton pt cut : GG4j-low-x");
      bookCut("Njet >= 4 : GG4j-low-x");
      bookCut("jets[0].pt() > 100 : GG4j-low-x");
      bookCut("jet2/3 pt > 100 : GG4j-low-x");
      bookCut("jets[3].pt() > 100 : GG4j-low-x");
      bookCut("mT > 125 : GG4j-low-x");
      bookCut("MET > 250 : GG4j-low-x");
      bookCut("meff_inc > 2000: GG4j-low-x");
      bookCut("jet aplanarity > 0.06: GG4j-low-x");

      bookCut("lepton pt cut : GG4j-high-x");
      bookCut("Njet >= 4 : GG4j-high-x");
      bookCut("jets[0].pt() > 400 : GG4j-high-x");
      bookCut("jet2/3 pt > 100 : GG4j-high-x");
      bookCut("jets[3].pt() [30,100] : GG4j-high-x");
      bookCut("mT > 475 : GG4j-high-x");
      bookCut("MET > 250 : GG4j-high-x");
      bookCut("meff_inc > 1600: GG4j-high-x");
      bookCut("MET/meff_inc > 0.1: GG4j-high-x");

      bookCut("Njet >= 4 : SS4J-x=1/2");
      bookCut("jet1/2 pt > 50 : SS4J-x=1/2");
      bookCut("jet3/4 pt > 50 : SS4J-x=1/2");
      bookCut("bjet veto: SS4J-x=1/2");
      bookCut("mT > 175 : SS4J-x=1/2");
      bookCut("MET > 300 : SS4J-x=1/2");
      bookCut("meff_inc > 1200: SS4J-x=1/2");
      bookCut("lepton aplanarity > 0.08: SS4J-x=1/2");

      bookCut("Njet >= 4 : SS4J-low-x");
      bookCut("jet1/2 pt > 250 : SS4J-low-x");
      bookCut("jet3/4 pt > 30 : SS4J-low-x");
      bookCut("bjet veto: SS4J-low-x");
      bookCut("mT [150,400] : SS4J-low-x");
      bookCut("MET > 250 : SS4J-low-x");
      bookCut("lepton aplanarity > 0.03: SS4J-low-x");

      bookCut("Njet >= 5 : SS5J-x=1/2");
      bookCut("jet1/2 pt > 50 : SS5J-x=1/2");
      bookCut("jet3/4 pt > 50 : SS5J-x=1/2");
      bookCut("jets[4].pT() > 50 : SS5J-x=1/2");
      bookCut("bjet veto: SS5J-x=1/2");
      bookCut("mT > 175 : SS5J-x=1/2");
      bookCut("MET > 300 : SS5J-x=1/2");
      bookCut("MET/meff_inc > 0.1: SS5J-x=1/2");

      bookCut("Njet >= 5 : SS5J-high-x");
      bookCut("jet1/2 pt > 30 : SS5J-high-x");
      bookCut("jet3/4 pt > 30 : SS5J-high-x");
      bookCut("jets[4].pT() > 30 : SS5J-high-x");
      bookCut("bjet veto: SS5J-high-x");
      bookCut("mT > 400 : SS5J-high-x");
      bookCut("MET > 300 : SS5J-high-x");
      bookCut("lepton aplanarity > 0.03: SS5J-high-x");
			// End init section -- do not edit/remove this comment
      bookHisto1D("GG2JMET", 6, 250, 670, "GG2JMET", "MET", "Arbitrary");
      bookHisto1D("GG6J-bulkMT", 6, 65, 545, "GG6J-bulkMT", "MT", "Arbitrary");
      bookHisto1D("GG6J-high-massaplan", 6, 0, 0.12, "GG6J-high-massaplan", "APLAN", "Arbitrary");
      bookHisto1D("GG4J-high-xMT", 4, 175, 575, "GG4J-high-xMT", "MT", "Arbitrary");
      bookHisto1D("GG4J-low-xMT", 4, 25, 425, "GG4J-low-xMT", "MT", "Arbitrary");
      bookHisto1D("GG4J-low-x-b-vetoaplan", 6, 0, 0.09, "GG4J-low-x-b-vetoaplan", "APLAN", "Arbitrary");
      bookHisto1D("SS4J-x=1/2aplan", 4, 0, 0.16, "SS4J-x=1/2aplan", "APLAN", "Arbitrary");
      bookHisto1D("SS5J-x=1/2MET", 6, 200, 500, "SS5J-x=1/2MET", "MET", "Arbitrary");
      bookHisto1D("SS4J-low-xmt", 7, 50, 400, "SS4J-low-xmt", "MT", "Arbitrary");
      bookHisto1D("SS5J-high-xaplan", 3, 0, 0.09, "SS5J-high-xaplan", "APLAN", "Arbitrary");


		}


		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

			// temporary
 	    const Particles& jets = applyProjection<NearIsoParticle>(event, "Jets").particlesByPt();
      // const Particles& large_jets = applyProjection<FastJets>(event, "Large_Jets").particlesByPt();      
      const Particles& leps = applyProjection<MergedFinalState>(event, "Signal_Leptons").particlesByPt();
      // const Particles& leps_cand = applyProjection<MergedFinalState>(event, "Candidate_Leptons").particlesByPt();            
      const HeavyFlavorJets& bjproj = applyProjection<NearIsoParticle>(event, "BJets");
      const Particles& bjets = bjproj.getTaggedJets(); 
      const Particles& untagged_jets = bjproj.getUntaggedJets(); 
      const MissingMomentum& pmet = applyProjection<MissingMomentum>(event, "MissingEt");
      const FourMomentum met = pmet.missingEt(); // met is four-momentum but pz and E is set zero

      double MET = met.pT(); 
      int Njet = jets.size(); 
      int Nb = bjets.size();
      int Nlep = leps.size();
      // int Nlep_cand = leps_cand.size();      
      // int N_J = large_jets.size();
      // temporary

      if (!cut(Nlep,CUT_EQ,1,"Nlep = 1 : base")) vetoEvent; 
      		
      double mT = -100;
     	mT = get_mT(leps[0].momentum(), met);

     	double meff_inc = MET+leps[0].pT();
     	for(int i = 0; i < Njet; i++){
				meff_inc += jets[i].pT();
      } 
		  //redifine aplanarity once bug is fixed
			double aplan_jet = 1;
			double aplan_lep = 1;

		  //GG2J
      if(cut(leps[0].pT(),CUT_LE,35,"lepton pt cut : GG2J")){
        if(cut(Njet,CUT_GT,1,"Njet >= 2 : GG2J")){
          if(cut(jets[0].pT(),CUT_GT,200,"jets[0].pt() > 200 : GG2J")){
            if(cut(jets[1].pT(),CUT_GT,30,"jets[1].pt() > 30 : GG2J")){
              if(cut(mT,CUT_GT,100,"mT > 100 : GG2J")){
                if(cut(MET,CUT_GT,460,"MET > 460 : GG2J")){
                  if(cut(MET/meff_inc,CUT_GT,0.35,"MET/meff_inc > 0.35: GG2J")){
                    pass("GG2J");
                  }
                }
              }
            }
          }
        }
      }
     

      //GG6J
      if(cut(leps[0].pT(),CUT_GT,35,"lepton pt cut : GG6J-bulk")){
        if(cut(Njet,CUT_GT,5,"Njet >= 6 : GG6J-bulk")){
          if(cut(jets[0].pT(),CUT_GT,125,"jets[0].pt() > 125 : GG6J-bulk")){
            bool jet23pt = false;
            if(jets[1].pT()>30 and jets[2].pT()>30) jet23pt = true;
            if(cut(jet23pt,"jet2/3 pt > 30 : GG6J-bulk")){
              if(cut(jets[3].pT(),CUT_GT,30,"jets[3].pt() > 30 : GG6J-bulk")){
                bool jet56pt = false;
                if(jets[4].pT()>30 and jets[5].pT()>30) jet56pt = true;
                if(cut(jet56pt,"jet5/6 pt > 30 : GG6J-bulk")){
                  if(cut(mT,CUT_GT,225,"mT > 225 : GG6J-bulk")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : GG6J-bulk")){
                      if(cut(meff_inc,CUT_GT,1000,"meff_inc > 1000: GG6J-bulk")){
                      if(cut(MET/meff_inc,CUT_GT,0.2,"MET/meff_inc > 0.2: GG6J-bulk")){
                        if(cut(aplan_jet,CUT_GT,0.04,"jet aplanarity > 0.04: GG6J-bulk")){
                        pass("GG6J-bulk");
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

  if(cut(leps[0].pT(),CUT_GT,35,"lepton pt cut : GG6J-high-mass")){
        if(cut(Njet,CUT_GT,5,"Njet >= 6 : GG6J-high-mass")){
          if(cut(jets[0].pT(),CUT_GT,125,"jets[0].pt() > 125 : GG6J-high-mass")){
            bool jet23pt = false;
            if(jets[1].pT()>30 and jets[2].pT()>30) jet23pt = true;
            if(cut(jet23pt,"jet2/3 pt > 30 : GG6J-high-mass")){
              if(cut(jets[3].pT(),CUT_GT,30,"jets[3].pt() > 30 : GG6J-high-mass")){
                bool jet56pt = false;
                if(jets[4].pT()>30 and jets[5].pT()>30) jet56pt = true;
                if(cut(jet56pt,"jet5/6 pt > 30 : GG6J-high-mass")){
                  if(cut(mT,CUT_GT,225,"mT > 225 : GG6J-high-mass")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : GG6J-high-mass")){
                      if(cut(meff_inc,CUT_GT,1000,"meff_inc > 1000: GG6J-high-mass")){
                      if(cut(MET/meff_inc,CUT_GT,0.1,"MET/meff_inc > 0.1: GG6J-high-mass")){
                        if(cut(aplan_jet,CUT_GT,0.04,"jet aplanarity > 0.04: GG6J-high-mass")){
                        pass("GG6J-high-mass");
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
    	      //GG4j

  if(cut(leps[0].pT(),CUT_GT,6,"lepton pt cut : GG4j-low-x")){
        if(cut(Njet,CUT_GT,3,"Njet >= 4 : GG4j-low-x")){
          if(cut(jets[0].pT(),CUT_GT,100,"jets[0].pt() > 100 : GG4j-low-x")){
            bool jet23pt = false;
            if(jets[1].pT()>100 and jets[2].pT()>100) jet23pt = true;
            if(cut(jet23pt,"jet2/3 pt > 100 : GG4j-low-x")){
              if(cut(jets[3].pT(),CUT_GT,100,"jets[3].pt() > 100 : GG4j-low-x")){           
                  if(cut(mT,CUT_GT,125,"mT > 125 : GG4j-low-x")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : GG4j-low-x")){
                      if(cut(meff_inc,CUT_GT,2000,"meff_inc > 2000: GG4j-low-x")){
                        if(cut(aplan_jet,CUT_GT,0.06,"jet aplanarity > 0.06: GG4j-low-x")){
                        pass("GG4j-low-x");                   
                   
                  }
                }
              }
            }
          }
        }
      }
    }
  }

    if(cut(leps[0].pT(),CUT_GT,6,"lepton pt cut : GG4j-low-x-b-veto")){
        if(cut(Njet,CUT_GT,3,"Njet >= 4 : GG4j-low-x-b-veto")){
          if(cut(jets[0].pT(),CUT_GT,100,"jets[0].pt() > 100 : GG4j-low-x-b-veto")){
            bool jet23pt = false;
            if(jets[1].pT()>100 and jets[2].pT()>100) jet23pt = true;
            if(cut(jet23pt,"jet2/3 pt > 100 : GG4j-low-x-b-veto")){
              if(cut(jets[3].pT(),CUT_GT,100,"jets[3].pt() > 100 : GG4j-low-x-b-veto")){   
                if(cut(Nb,CUT_EQ,0,"bjet veto: GG4J-low-x-b-veto")){        
                  if(cut(mT,CUT_GT,125,"mT > 125 : GG4j-low-x-b-veto")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : GG4j-low-x-b-veto")){
                      if(cut(meff_inc,CUT_GT,2000,"meff_inc > 2000: GG4j-low-x-b-veto")){
                        if(cut(aplan_jet,CUT_GT,0.06,"jet aplanarity > 0.06: GG4j-low-x-b-veto")){
                        pass("GG4j-low-x-b-veto");                   
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

  if(cut(leps[0].pT(),CUT_GT,35,"lepton pt cut : GG4j-high-x")){
        if(cut(Njet,CUT_GT,3,"Njet >= 4 : GG4j-high-x")){
          if(cut(jets[0].pT(),CUT_GT,400,"jets[0].pt() > 400 : GG4j-high-x")){
            bool jet23pt = false;
            if(jets[1].pT()>30 and jets[2].pT()>30) jet23pt = true;
            if(cut(jet23pt,"jet2/3 pt > 100 : GG4j-high-x")){
              if(cut(jets[3].pT(),CUT_IN,make_pair(30.,100.),"jets[3].pt() [30,100] : GG4j-high-x")){           
                  if(cut(mT,CUT_GT,475,"mT > 475 : GG4j-high-x")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : GG4j-high-x")){
                      if(cut(meff_inc,CUT_GT,1600,"meff_inc > 1600: GG4j-high-x")){
                        if(cut(MET/meff_inc,CUT_GT,0.3,"MET/meff_inc > 0.1: GG4j-high-x")){
                        pass("GG4j-high-x");                   
                   
                  }
                }
              }
            }
          }
        }
      }
    }
  }
      

if (cut(leps[0].pT(),CUT_GT,35,"lepton pt cut : SS base")){

      //SS4j

  if(cut(Njet,CUT_GT,3,"Njet >= 4 : SS4J-x=1/2")){
          
            bool jet12pt = false;
            if(jets[0].pT()>50 and jets[1].pT()>50) jet12pt = true;
            if(cut(jet12pt,"jet1/2 pt > 50 : SS4J-x=1/2")){
              bool jet34pt = false;
            if(jets[2].pT()>50 and jets[3].pT()>50) jet34pt = true;
              if(cut(jet34pt,"jet3/4 pt > 50 : SS4J-x=1/2")){
                if (cut(Nb,CUT_EQ,0,"bjet veto: SS4J-x=1/2")){
                  if(cut(mT,CUT_GT,175,"mT > 175 : SS4J-x=1/2")){
                    if(cut(MET,CUT_GT,300,"MET > 300 : SS4J-x=1/2")){
                      if(cut(meff_inc,CUT_GT,1200,"meff_inc > 1200: SS4J-x=1/2")){
                        if(cut(aplan_lep,CUT_GT,0.08,"lepton aplanarity > 0.08: SS4J-x=1/2")){
                        pass("SS4J-x=1/2");                   
                   
                  
                }
              }
            }
          }
        }
      }
    }
}
  if(cut(Njet,CUT_GT,3,"Njet >= 4 : SS4J-low-x")){
          
            bool jet12pt = false;
            if(jets[0].pT()>250 and jets[1].pT()>250) jet12pt = true;
            if(cut(jet12pt,"jet1/2 pt > 250 : SS4J-low-x")){
              bool jet34pt = false;
            if(jets[2].pT()>30 and jets[3].pT()>30) jet34pt = true;
              if(cut(jet34pt,"jet3/4 pt > 30 : SS4J-low-x")){
                if (cut(Nb,CUT_EQ,0,"bjet veto: SS4J-low-x")){
                  if(cut(mT,CUT_IN,make_pair(150,400),"mT [150,400] : SS4J-low-x")){
                    if(cut(MET,CUT_GT,250,"MET > 250 : SS4J-low-x")){
                        if(cut(aplan_lep,CUT_GT,0.03,"lepton aplanarity > 0.03: SS4J-low-x")){
                        pass("SS4J-low-x");                   
                   
              }
            }
          }
        }
      }
    }
}

  
      //SS5j
  if(cut(Njet,CUT_GT,4,"Njet >= 5 : SS5J-x=1/2")){
          
            bool jet12pt = false;
            if(jets[0].pT()>50 and jets[1].pT()>50) jet12pt = true;
            if(cut(jet12pt,"jet1/2 pt > 50 : SS5J-x=1/2")){
              bool jet34pt = false;
            if(jets[2].pT()>50 and jets[3].pT()>50) jet34pt = true;
              if(cut(jet34pt,"jet3/4 pt > 50 : SS5J-x=1/2")){
                if(cut(jets[4].pT(),CUT_GT,50,"jets[4].pT() > 50 : SS5J-x=1/2")){
                if (cut(Nb,CUT_EQ,0,"bjet veto: SS5J-x=1/2")){
                  if(cut(mT,CUT_GT,175,"mT > 175 : SS5J-x=1/2")){
                    if(cut(MET,CUT_GT,300,"MET > 300 : SS5J-x=1/2")){
                        if(cut(MET/meff_inc,CUT_GT,0.2,"MET/meff_inc > 0.1: SS5J-x=1/2")){
                        pass("SS5J-x=1/2");                   
                }
              }
            }
          }
        }
      }
    }
}

    if(cut(Njet,CUT_GT,4,"Njet >= 5 : SS5J-high-x")){
          
            bool jet12pt = false;
            if(jets[0].pT()>30 and jets[1].pT()>30) jet12pt = true;
            if(cut(jet12pt,"jet1/2 pt > 30 : SS5J-high-x")){
              bool jet34pt = false;
            if(jets[2].pT()>30 and jets[3].pT()>30) jet34pt = true;
              if(cut(jet34pt,"jet3/4 pt > 30 : SS5J-high-x")){
                if(cut(jets[4].pT(),CUT_GT,30,"jets[4].pT() > 30 : SS5J-high-x")){
                if (cut(Nb,CUT_EQ,0,"bjet veto: SS5J-high-x")){
                  if(cut(mT,CUT_GT,400,"mT > 400 : SS5J-high-x")){
                    if(cut(MET,CUT_GT,300,"MET > 300 : SS5J-high-x")){
                        if(cut(aplan_lep,CUT_GT,0.03,"lepton aplanarity > 0.03: SS5J-high-x")){
                        pass("SS5J-high-x");                   
                }
              }
            }
          }
        }
      }
    }
}

}     				
    	
			// Projection application section -- do not edit/remove this comment
      /// @todo apply projections
			// const Particles& jets = applyProjection< FastJets >(event, "Jets").particlesByPt(&event);

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
	AtomPlugin(ATLAS_CONF_2016_054)
}
