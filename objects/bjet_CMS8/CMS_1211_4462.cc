///
/// @file  CMS_1211_4462.cc
/// @brief Implementation of CMS_1211_4462 analysis
/// @author Pi
/// @date created 07/21/2016
/// @date last revision 07/21/2016
///

// -*- C++ -*-
#include "Atom/Analysis.hh"
// use the stock set or include only the projections you use from "include/Atom/Projections/" directory if compiling takes too long
#include "Atom/Tools/CommonProjections.hh"

using namespace std;

namespace Atom {

	class CMS_1211_4462 : public Analysis {
	public:

	  CMS_1211_4462()
	    : Analysis("CMS_1211_4462") {
	    setNeedsCrossSection(true);
	  }
	  
	  /// @name Analysis methods
	  //@{
	  
	  /// Book histograms and initialise projections before the run
	  void init() {
	    
            useDetector( "ATLAS_CMS_all" ); // "CMS2012" );
	    
	    FinalState base( getRange( "Full_Range_CMS" ) );
	    
	    Range jet_range = Range(PT, 30., 7000.) & 
	      Range(ETA, -2.4, 2.4);
	    
	    FastJets jets(base, 
			  getRange( "HCal_Range_CMS" ) & jet_range, 
			  getRange( "Muon_Range_Detector_CMS" ), 
			  FastJets::ANTIKT, 0.5 );

            jets.useMuons(true);
            jets.setSmearingParams( getJetSim( "Jet_Smear_Topo_CMS" ) );


            HeavyFlavorJets bjets(jets, *getBJetEff("BJet_Ident_CSVM_CMS"));

	    // Projection booking section -- do not edit/remove this comment
            /// @todo define projections (see examples and manual)
	    // FastJets jets(fsbase, hadRange & jet_range, muDetRange, FastJets::ANTIKT, 0.6 );
	    // addProjection(jets, "Jets");
	    addProjection(jets, "Jets");
	    addProjection(bjets, "BJets");	    

	    // Histogram booking section -- do not edit/remove this comment
            /// @todo book histograms
            // triplets of numbers correspond to HepData notation d??-x??-y??
	    _h_eff = bookProfile1D(1, 1, 1, "Eff_PT");
	    
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

	    bookEfficiency("Efficiency");
	    
	    bookCut("NumJets");
	    bookCut("LeadJ_pT");
	    bookCut("LeadJ_etaRange");
	    
	    // End init section -- do not edit/remove this comment
	  }
	  

		/// Perform the per-event analysis
		/// param[in]   event    the event to be analyzed
		void analyze(const Event& event) {

		  const Particles& jets = 
		    applyProjection< FastJets >(event, "Jets").particlesByPt();

		  const HeavyFlavorJets& bjets=applyProjection< HeavyFlavorJets >(event, "BJets");

		  const Particles& bjets_tagged = 
		    bjets.getTaggedJets();

		  const Particles& bjets_untagged = 
		    bjets.getUntaggedJets();

		  cout<<endl;
		  cout<<"New Event"<<endl;
		  for(int i=0; i<jets.size(); i++){
		    cout<<"jets: "<<i<<endl;
		    cout<<jets[i].tag().to_string()<<endl;
		  }

		  // apply cuts		  
		  if (!cut(jets.size(), CUT_GT, 0, "NumJets")) {
		    vetoEvent;
		  }


		  const FourMomentum j0(jets[0].momentum());
		  
		  if (!cut(j0.pT(), CUT_GT, 60.0, "LeadJ_pT")) {
		    vetoEvent;
		  }
		  std::pair<double, double> etarange = 
		    std::make_pair(-2.4, 2.4);
		  
		  if (cut(fabs(j0.eta()), CUT_OUT, etarange, "LeadJ_etaRange")) {
		    vetoEvent;
		  }
		  
		  pass("Efficiency");

		  
		  // Projection application section -- do not edit/remove this comment

		  //bjets_tagged always empty!
		  for(int i=0; i<bjets_tagged.size(); ++i){
		    cout<<"not empty"<<endl;
		    fillProfile1D(_h_eff, bjets_tagged[i].pT(),1.0);
		  }
		  
		  for(int i=0; i<bjets_untagged.size(); ++i){
		    fillProfile1D(_h_eff, bjets_untagged[i].pT(), 0.0);
		  }

		  // End analyze section -- do not edit/remove this comment
		}


		/// Normalise histograms etc., after the run
		void finalize() {
			// Histogram normalization section -- do not edit/remove this comment
			/// @todo normalize the histograms
			// scale("Mjj");

			// End finalize section -- do not edit/remove this comment
		}

	  
	  Profile1DPtr _h_eff;
	  //@}

	};

	// This global object acts as a hook for the plugin system
	AtomPlugin(CMS_1211_4462)
}
