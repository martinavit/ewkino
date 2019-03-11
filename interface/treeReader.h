#ifndef treeReader_h
#define treeReader_h

//include ROOT classes
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

//include other parts of code
#include "Reweighter_old.h"
#include "Sample.h"
#include "HistInfo.h"

class treeReader {
    public :
        //Declare leaf types
        static const unsigned nL_max = 20;
        static const unsigned nJets_max = 20;
        static const unsigned gen_nL_max = 20;
        static const unsigned gen_nPh_max = 10;
   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   ULong64_t       _eventNb;
   UChar_t         _nVertex;
   Bool_t          _passMETFilters;
   Bool_t          _Flag_goodVertices;
   Bool_t          _Flag_HBHENoiseFilter;
   Bool_t          _Flag_HBHENoiseIsoFilter;
   Bool_t          _Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          _Flag_BadPFMuonFilter;
   Bool_t          _Flag_BadChargedCandidateFilter;
   Bool_t          _Flag_eeBadScFilter;
   Bool_t          _updated_ecalBadCalibFilter;
   Bool_t          _passTrigger_1l;
   Bool_t          _HLT_IsoMu24;
   Int_t           _HLT_IsoMu24_prescale;
   Bool_t          _HLT_IsoMu27;
   Int_t           _HLT_IsoMu27_prescale;
   Bool_t          _HLT_Ele32_WPTight_Gsf;
   Int_t           _HLT_Ele32_WPTight_Gsf_prescale;
   UInt_t          _nL;
   UInt_t          _nMu;
   UInt_t          _nEle;
   UInt_t          _nLight;
   UInt_t          _nTau;
   Double_t        _pvX;
   Double_t        _pvY;
   Double_t        _pvZ;
   Double_t        _pvXErr;
   Double_t        _pvYErr;
   Double_t        _pvZErr;
   UChar_t         _nMu;
   UChar_t         _nEle;
   UChar_t         _nLight;
   UChar_t         _nTau;
   UInt_t          _nVFit_os;
   UInt_t          _nVFit;
   UInt_t          _nGoodLeading;
   UInt_t          _nGoodDisplaced;
   Double_t        _vertices_os[4][12];   //[_nVFit_os]
   Double_t        _lDisplaced_os[4][24];   //[_nVFit_os]
   Double_t        _vertices[12][12];   //[_nVFit]
   Double_t        _lDisplaced[12][24];   //[_nVFit]
   UInt_t          _lHasTrigger[20];   //[_nL]
   Double_t        _lPt[20];   //[_nL]
   Double_t        _lEta[20];   //[_nL]
   Double_t        _lEtaSC[4];   //[_nLight]
   Double_t        _lPhi[20];   //[_nL]
   Double_t        _lE[20];   //[_nL]
   UInt_t          _lFlavor[20];   //[_nL]
   Int_t           _lCharge[20];   //[_nL]
   Double_t        _dxy[20];   //[_nL]
   Double_t        _dz[20];   //[_nL]
   Double_t        _3dIP[20];   //[_nL]
   Double_t        _3dIPSig[20];   //[_nL]
   Double_t        _2dIP[20];   //[_nL]
   Double_t        _2dIPSig[20];   //[_nL]
   Bool_t          _lElectronPassEmu[4];   //[_nLight]
   Bool_t          _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[20];   //[_nL]
   Bool_t          _lElectronPassConvVeto[4];   //[_nLight]
   Bool_t          _lElectronChargeConst[4];   //[_nLight]
   UInt_t          _lElectronMissingHits[4];   //[_nLight]
   Bool_t          _lPOGVeto[20];   //[_nL]
   Bool_t          _lPOGLoose[20];   //[_nL]
   Bool_t          _lPOGMedium[20];   //[_nL]
   Bool_t          _lPOGTight[20];   //[_nL]
   Bool_t          _lGlobalMuon[4];   //[_nMu]
   Bool_t          _lTrackerMuon[4];   //[_nMu]
   Double_t        _lInnerTrackValidFraction[4];   //[_nMu]
   Double_t        _lGlobalTrackNormalizeChi2[4];   //[_nMu]
   Double_t        _lCQChi2Position[4];   //[_nMu]
   Double_t        _lCQTrackKink[4];   //[_nMu]
   UInt_t          _lNumberOfMatchedStation[4];   //[_nMu]
   UInt_t          _lNumberOfValidPixelHits[4];   //[_nMu]
   UInt_t          _lTrackerLayersWithMeasurement[4];   //[_nMu]
   Int_t           _lSimType[4];   //[_nMu]
   Int_t           _lSimExtType[4];   //[_nMu]
   Int_t           _lSimFlavour[4];   //[_nMu]
   Int_t           _muDTStationsWithValidHits[4];   //[_nMu]
   Int_t           _muCSCStationsWithValidHits[4];   //[_nMu]
   Int_t           _muRPCStationsWithValidHits[4];   //[_nMu]
   Int_t           _muMuonStationsWithValidHits[4];   //[_nMu]
   Int_t           _lMuRPCTimenDof[4];   //[_nMu]
   Int_t           _lMuTimenDof[4];   //[_nMu]
   Double_t        _lMuRPCTime[4];   //[_nMu]
   Double_t        _lMuRPCTimeErr[4];   //[_nMu]
   Double_t        _lMuTime[4];   //[_nMu]
   Double_t        _lMuTimeErr[4];   //[_nMu]
   UInt_t          _muNumberInnerHits[4];   //[_nMu]
   Bool_t          _lEleIsEB[4];   //[_nLight]
   Bool_t          _lEleIsEE[4];   //[_nLight]
   Double_t        _lEleSuperClusterOverP[4];   //[_nLight]
   Double_t        _lEleEcalEnergy[4];   //[_nLight]
   Double_t        _lElefull5x5SigmaIetaIeta[4];   //[_nLight]
   Double_t        _lEleDEtaInSeed[4];   //[_nLight]
   Double_t        _lEleDeltaPhiSuperClusterTrackAtVtx[4];   //[_nLight]
   Double_t        _lElehadronicOverEm[4];   //[_nLight]
   Double_t        _lEleInvMinusPInv[4];   //[_nLight]
   Double_t        _puCorr[4];   //[_nLight]
   Double_t        _absIso03[20];   //[_nL]
   Double_t        _absIso04[4];   //[_nMu]
   Double_t        _sumNeutralHadronEt04[4];   //[_nMu]
   Double_t        _sumChargedHadronPt04[4];   //[_nMu]
   Double_t        _sumPhotonEt04[4];   //[_nMu]
   Double_t        _sumNeutralHadronEt03[4];   //[_nLight]
   Double_t        _sumChargedHadronPt03[4];   //[_nLight]
   Double_t        _sumPhotonEt03[4];   //[_nLight]
   Double_t        _trackIso[4];   //[_nLight]
   Double_t        _ecalIso[4];   //[_nLight]
   Double_t        _hcalIso[4];   //[_nLight]
   Double_t        _ecalPFClusterIso[4];   //[_nLight]
   Double_t        _hcalPFClusterIso[4];   //[_nLight]
   Bool_t          _tauMuonVeto[20];   //[_nL]
   Double_t        _relIso[4];   //[_nLight]
   Double_t        _relIso0p4[4];   //[_nLight]
   Double_t        _relIso0p4MuDeltaBeta[4];   //[_nMu]
   Double_t        _ptRel[4];   //[_nLight]
   Double_t        _ptRatio[4];   //[_nLight]
   Double_t        _closestJetCsvV2[4];   //[_nLight]
   Double_t        _closestJetDeepCsv_b[4];   //[_nLight]
   Double_t        _closestJEC[4];   //[_nLight]
   Double_t        _closest_lepAwareJetE[4];   //[_nLight]
   Double_t        _closest_lepAwareJetPx[4];   //[_nLight]
   Double_t        _closest_lepAwareJetPy[4];   //[_nLight]
   Double_t        _closest_lepAwareJetPz[4];   //[_nLight]
   Double_t        _closest_l1JetE[4];   //[_nLight]
   Double_t        _closest_l1JetPx[4];   //[_nLight]
   Double_t        _closest_l1JetPy[4];   //[_nLight]
   Double_t        _closest_l1JetPz[4];   //[_nLight]
   Double_t        _closest_lJetE[4];   //[_nLight]
   Double_t        _closest_lJetPx[4];   //[_nLight]
   Double_t        _closest_lJetPy[4];   //[_nLight]
   Double_t        _closest_lJetPz[4];   //[_nLight]
   Double_t        _closestJetDeepCsv_bb[4];   //[_nLight]
   UInt_t          _selectedTrackMult[4];   //[_nLight]
   Double_t        _lMuonSegComp[4];   //[_nMu]
   Double_t        _lMuonTrackPt[4];   //[_nMu]
   Double_t        _lMuonTrackPtErr[4];   //[_nMu]
   UInt_t          _nJets;
   Double_t        _jetPt[20];   //[_nJets]
   Double_t        _jetPt_JECDown[20];   //[_nJets]
   Double_t        _jetPt_JECUp[20];   //[_nJets]
   Double_t        _jetSmearedPt[20];   //[_nJets]
   Double_t        _jetSmearedPt_JECDown[20];   //[_nJets]
   Double_t        _jetSmearedPt_JECUp[20];   //[_nJets]
   Double_t        _jetSmearedPt_JERDown[20];   //[_nJets]
   Double_t        _jetSmearedPt_JERUp[20];   //[_nJets]
   Double_t        _jetPt_Uncorrected[20];   //[_nJets]
   Double_t        _jetPt_L1[20];   //[_nJets]
   Double_t        _jetPt_L2[20];   //[_nJets]
   Double_t        _jetPt_L3[20];   //[_nJets]
   Double_t        _jetEta[20];   //[_nJets]
   Double_t        _jetPhi[20];   //[_nJets]
   Double_t        _jetE[20];   //[_nJets]
   Double_t        _jetCsvV2[20];   //[_nJets]
   Double_t        _jetDeepCsv_udsg[20];   //[_nJets]
   Double_t        _jetDeepCsv_b[20];   //[_nJets]
   Double_t        _jetDeepCsv_c[20];   //[_nJets]
   Double_t        _jetDeepCsv_bb[20];   //[_nJets]
   UInt_t          _jetHadronFlavor[20];   //[_nJets]
   Bool_t          _jetIsLoose[20];   //[_nJets]
   Bool_t          _jetIsTight[20];   //[_nJets]
   Bool_t          _jetIsTightLepVeto[20];   //[_nJets]
   Double_t        _jetNeutralHadronFraction[20];   //[_nJets]
   Double_t        _jetChargedHadronFraction[20];   //[_nJets]
   Double_t        _jetNeutralEmFraction[20];   //[_nJets]
   Double_t        _jetChargedEmFraction[20];   //[_nJets]
   Double_t        _jetHFHadronFraction[20];   //[_nJets]
   Double_t        _jetHFEmFraction[20];   //[_nJets]
   Double_t        _met;
   Double_t        _metRaw;
   Double_t        _metJECDown;
   Double_t        _metJECUp;
   Double_t        _metUnclDown;
   Double_t        _metUnclUp;
   Double_t        _metPhi;
   Double_t        _metRawPhi;
   Double_t        _metPhiJECDown;
   Double_t        _metPhiJECUp;
   Double_t        _metPhiUnclDown;
   Double_t        _metPhiUnclUp;
   Double_t        _metSignificance;


        //Constructor
        treeReader(TTree *tree = nullptr);

        //set up tree for reading and writing
        void initTree(TTree *tree, const bool isData = false);
        void setOutputTree(TTree*, const bool isData = false);

        //skim tree
        void skimTree(const std::string&, std::string outputDirectory = "", const bool isData = false);
        void combinePD(std::vector<std::string>& datasets, const bool is2017, std::string outputDirectory = "");

        //set up tree for analysis
        void readSamples(const std::string& list, const std::string& directory); //read sample list from file
        void readSamples2016(const std::string&, const std::string&);
        void readSamples2017(const std::string&, const std::string&);

        void initSample();                              //event weights will be set according to is2016() ( or equally is2017() ) flag
        void initSample(const Sample&);  

        //functions to analyze tree
        void GetEntry(long unsigned entry);
        void GetEntry(const Sample&, long unsigned entry);
        void Analyze();
        void Analyze(const std::string&, long unsigned, long unsigned);
        void Analyze(const Sample&, long unsigned, long unsigned);
        void Analyze(const std::string&);
        void Analyze(const Sample&);
        void setup();
        void splitJobs();
        void Loop(const std::string& sample, const double xSection);

        //new functions for parallel plotting
        void plot(const std::string&);
        void splitPlots();

        //functions for event selection
        void orderByPt(std::vector<unsigned>&, const double*, const unsigned) const;
        unsigned dilFlavorComb(const std::vector<unsigned>&) const;
        double coneCorr(const unsigned) const;
        void applyConeCorrection();
        bool lepIsLoose(const unsigned) const;
        bool lepIsGood(const unsigned) const;
        bool lepIsTight(const unsigned) const;
        bool lepFromMEExtConversion(const unsigned) const;
        bool eleIsClean(const unsigned) const;
        double closestJetDeepCSV(const unsigned) const;

        unsigned selectLep(std::vector<unsigned>&) const;
        unsigned selectLepConeCorr(std::vector<unsigned>&);
        unsigned tightLepCount(const std::vector<unsigned>&, const unsigned) const;

        bool passPtCuts(const std::vector<unsigned>&) const;
        bool jetIsClean(const unsigned) const;
        bool jetIsGood(const unsigned, const double ptCut = 25., const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;
        unsigned nJets(const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;                                   //without jet pt ordering
        unsigned nJets(std::vector<unsigned>& jetInd, const unsigned unc = 0, const bool clean = true, const bool allowForward = false) const;    //with jet pt ordering
        double deepCSV(const unsigned) const;
        bool bTagged(const unsigned ind, const unsigned wp = 1, const bool deepCSV = true) const;
        unsigned nBJets(const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, const unsigned wp = 1) const;
        unsigned nBJets(std::vector<unsigned>& bJetInd, const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, const unsigned wp = 1) const;

        //baseline selection for leptonMva training
        bool lepPassBaseline(const unsigned) const;

        //trigger decitions
        bool passSingleLeptonTriggers() const;
        bool passDileptonTriggers() const;
        bool passTrileptonTriggers() const;
        bool passTriggerCocktail() const;
        bool passMETFilters() const;

        //overlap removal between samples
        bool photonOverlap(const bool mcNonprompt = true) const;                                            //sample overlap due to photons
        bool photonOverlap(const Sample&, const bool mcNonprompt = true) const;
        bool htOverlap() const;                                                                             //sample overlap due to HT binning
        bool htOverlap(const Sample&) const;

        //check if leptons are prompt in MC
        bool promptLeptons() const;

        //compute b-tagging efficiency
        void computeBTagEff(const std::string& analysis, const bool clean, const bool deepCSV, const bool is2016);

        //event weights
        //pileup reweighting
        double puWeight(const unsigned unc = 0) const;

        //b-tag reweighting
        double bTagWeight_cut_singleJet(const unsigned jetIndex, const unsigned unc = 0) const;
        double bTagWeight_reshaping_singleJet(const unsigned jetIndex, const unsigned unc = 0) const;
        double bTagWeight_base(const unsigned jetFlavor, const unsigned unc, double (treeReader::*jetWeight)(const unsigned, const unsigned) const ) const;
        double bTagWeight_cut( const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight_reshaping( const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight(const unsigned jetFlavor, const unsigned unc = 0) const;
        double bTagWeight(const std::vector<unsigned>& jetInd, const unsigned jetFlavor, const unsigned unc = 0) const; //more efficient version if jets were already selected 
        double bTagWeight_udsg(const unsigned unc = 0) const;
        double bTagWeight_c(const unsigned unc = 0) const;
        double bTagWeight_b(const unsigned unc = 0) const;
        double bTagWeight(const unsigned unc = 0) const;

        //lepton reweighting
        double leptonWeight(const std::string& unc = "") const;

        //fake-rate
        double fakeRateWeight(const unsigned unc = 0);
        double sfWeight();
        double jetPrefiringWeight(const unsigned unc = 0) const;

    private:
        TTree* fChain;                                                          //current Tree
        std::shared_ptr<TFile> sampleFile;                                      //current sample
        std::vector<Sample> samples;                                            //combined list of samples
        std::vector<Sample> samples2016;                                        //2016 data and MC samples
        std::vector<Sample> samples2017;                                        //2017 data and MC samples
        Sample currentSample;                                                   //reference to current sample, needed to check what era sample belongs to
        std::vector<HistInfo> histInfo;                                         //histogram info
        int currentSampleIndex = -1;                                                 //current index in list
        //bool isData = false;
        double scale = 0;
        double weight = 1;                                                      //weight of given event
        unsigned long nEntries = 0;
        const double lumi2017 = 41.53;                                          //in units of 1/fb
        const double lumi2016 = 35.867;                 
        std::shared_ptr<Reweighter> reweighter;                                 //instance of reweighter class used for reweighting functions

        //check whether sample is 2017 or not
        bool is2017() const { return currentSample.is2017(); }
        bool is2016() const { return currentSample.is2016(); }                  //if sample is not 2017 it is automatically 2016
        bool isData() const { return currentSample.isData(); }
        bool isMC() const { return currentSample.isMC(); } 
        bool isSMSignal() const{ return currentSample.isSMSignal(); }
        bool isNewPhysicsSignal() const{ return currentSample.isNewPhysicsSignal(); }

        //check lepton flavors 
        bool isElectron(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 0); }
        bool isMuon(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 1); }
        bool isTau(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 2); }

        //era-specific event selection functions
        bool lepIsLooseBase(const unsigned) const;
        bool lepIsLoose2016(const unsigned) const;
        bool lepIsLoose2017(const unsigned) const;

        bool lepIsGoodBase(const unsigned) const;
        bool lepIsGood2016(const unsigned) const;
        bool lepIsGood2017(const unsigned) const;

        bool lepIsTightBase(const unsigned) const;
        bool lepIsTight2016(const unsigned) const;
        bool lepIsTight2017(const unsigned) const;

        bool eleIsCleanBase(const unsigned, bool (treeReader::*looseMuon)(const unsigned) const) const;
        bool eleIsClean2016(const unsigned) const;
        bool eleIsClean2017(const unsigned) const;

        bool jetIsCleanBase(const unsigned, bool (treeReader::*leptonIsFO)(const unsigned) const) const;

        bool bTaggedDeepCSVBase(const unsigned, const unsigned wp, const double cuts[3]) const;
        bool bTaggedDeepCSV2016(const unsigned, const unsigned wp = 1) const;
        bool bTaggedDeepCSV2017(const unsigned, const unsigned wp = 1) const;
        bool bTaggedDeepCSV(const unsigned, const unsigned wp = 1) const;

        bool bTaggedCSVv2Base(const unsigned, const unsigned wp, const double cuts[3]) const;
        bool bTaggedCSVv22016(const unsigned, const unsigned wp = 1) const;
        bool bTaggedCSVv22017(const unsigned, const unsigned wp = 1) const;
        bool bTaggedCSVv2(const unsigned, const unsigned wp = 1) const;


        //lepton selection for different parts of tZq and ttV analysis (to compute bTag efficiencies for everyone)
        bool lepIsGoodtZq(const unsigned) const;

        bool lepIsGoodttZ3l2016(const unsigned) const;
        bool lepIsGoodttZ3l2017(const unsigned) const;
        bool lepIsGoodttZ3l(const unsigned) const;

        bool lepIsGoodttZ4l2016(const unsigned) const;
        bool lepIsGoodttZ4l2017(const unsigned) const;
        bool lepIsGoodttZ4l(const unsigned) const;

        bool lepIsGoodttW2016(const unsigned) const;
        bool lepIsGoodttW2017(const unsigned) const;
        bool lepIsGoodttW(const unsigned) const;

        bool lepIsGoodMultiAnalysis(const std::string&, const unsigned) const;

        //initialize SF weights
        void initializeWeights();
        
        //some safety-checks for errors 
        void checkSampleEraConsistency() const;  //make sure a sample is not is2016() AND 2017() 
        void checkEraOrthogonality() const;        //make sure no sample from the wrong era is being used (i.e. no 2016 sample in the list of 2017 samples) 

        //debugging prints
        void printLeptonContent( std::ostream& os = std::cout ) const;
        void printLeptonPairing( std::ostream& os = std::cout ) const;

        //general function to read a list of samples
        void readSamples(const std::string&, const std::string&, std::vector<Sample>&);

        //list of branches
        TBranch        *b__runNb;   //!
   TBranch        *b__lumiBlock;   //!
   TBranch        *b__eventNb;   //!
   TBranch        *b__nVertex;   //!
   TBranch        *b__passMETFilters;   //!
   TBranch        *b__Flag_goodVertices;   //!
   TBranch        *b__Flag_HBHENoiseFilter;   //!
   TBranch        *b__Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b__Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b__Flag_BadPFMuonFilter;   //!
   TBranch        *b__Flag_BadChargedCandidateFilter;   //!
   TBranch        *b__Flag_eeBadScFilter;   //!
   TBranch        *b__updated_ecalBadCalibFilter;   //!
   TBranch        *b__passTrigger_1l;   //!
   TBranch        *b__HLT_IsoMu24;   //!
   TBranch        *b__HLT_IsoMu24_prescale;   //!
   TBranch        *b__HLT_IsoMu27;   //!
   TBranch        *b__HLT_IsoMu27_prescale;   //!
   TBranch        *b__HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b__HLT_Ele32_WPTight_Gsf_prescale;   //!
   TBranch        *b__nL;   //!
   TBranch        *b__nMu;   //!
   TBranch        *b__nEle;   //!
   TBranch        *b__nLight;   //!
   TBranch        *b__nTau;   //!
   TBranch        *b__pvX;   //!
   TBranch        *b__pvY;   //!
   TBranch        *b__pvZ;   //!
   TBranch        *b__pvXErr;   //!
   TBranch        *b__pvYErr;   //!
   TBranch        *b__pvZErr;   //!
   TBranch        *b__nMu;   //!
   TBranch        *b__nEle;   //!
   TBranch        *b__nLight;   //!
   TBranch        *b__nTau;   //!
   TBranch        *b__nVFit_os;   //!
   TBranch        *b__nVFit;   //!
   TBranch        *b__nGoodLeading;   //!
   TBranch        *b__nGoodDisplaced;   //!
   TBranch        *b__vertices_os;   //!
   TBranch        *b__lDisplaced_os;   //!
   TBranch        *b__vertices;   //!
   TBranch        *b__lDisplaced;   //!
   TBranch        *b__lHasTrigger;   //!
   TBranch        *b__lPt;   //!
   TBranch        *b__lEta;   //!
   TBranch        *b__lEtaSC;   //!
   TBranch        *b__lPhi;   //!
   TBranch        *b__lE;   //!
   TBranch        *b__lFlavor;   //!
   TBranch        *b__lCharge;   //!
   TBranch        *b__dxy;   //!
   TBranch        *b__dz;   //!
   TBranch        *b__3dIP;   //!
   TBranch        *b__3dIPSig;   //!
   TBranch        *b__2dIP;   //!
   TBranch        *b__2dIPSig;   //!
   TBranch        *b__lElectronPassEmu;   //!
   TBranch        *b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto;   //!
   TBranch        *b__lElectronPassConvVeto;   //!
   TBranch        *b__lElectronChargeConst;   //!
   TBranch        *b__lElectronMissingHits;   //!
   TBranch        *b__lPOGVeto;   //!
   TBranch        *b__lPOGLoose;   //!
   TBranch        *b__lPOGMedium;   //!
   TBranch        *b__lPOGTight;   //!
   TBranch        *b__lGlobalMuon;   //!
   TBranch        *b__lTrackerMuon;   //!
   TBranch        *b__lInnerTrackValidFraction;   //!
   TBranch        *b__lGlobalTrackNormalizeChi2;   //!
   TBranch        *b__lCQChi2Position;   //!
   TBranch        *b__lCQTrackKink;   //!
   TBranch        *b__lNumberOfMatchedStation;   //!
   TBranch        *b__lNumberOfValidPixelHits;   //!
   TBranch        *b__lTrackerLayersWithMeasurement;   //!
   TBranch        *b__lSimType;   //!
   TBranch        *b__lSimExtType;   //!
   TBranch        *b__lSimFlavour;   //!
   TBranch        *b__muDTStationsWithValidHits;   //!
   TBranch        *b__muCSCStationsWithValidHits;   //!
   TBranch        *b__muRPCStationsWithValidHits;   //!
   TBranch        *b__muMuonStationsWithValidHits;   //!
   TBranch        *b__lMuRPCTimenDof;   //!
   TBranch        *b__lMuTimenDof;   //!
   TBranch        *b__lMuRPCTime;   //!
   TBranch        *b__lMuRPCTimeErr;   //!
   TBranch        *b__lMuTime;   //!
   TBranch        *b__lMuTimeErr;   //!
   TBranch        *b__muNumberInnerHits;   //!
   TBranch        *b__lEleIsEB;   //!
   TBranch        *b__lEleIsEE;   //!
   TBranch        *b__lEleSuperClusterOverP;   //!
   TBranch        *b__lEleEcalEnergy;   //!
   TBranch        *b__lElefull5x5SigmaIetaIeta;   //!
   TBranch        *b__lEleDEtaInSeed;   //!
   TBranch        *b__lEleDeltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b__lElehadronicOverEm;   //!
   TBranch        *b__lEleInvMinusPInv;   //!
   TBranch        *b__puCorr;   //!
   TBranch        *b__absIso03;   //!
   TBranch        *b__absIso04;   //!
   TBranch        *b__sumNeutralHadronEt04;   //!
   TBranch        *b__sumChargedHadronPt04;   //!
   TBranch        *b__sumPhotonEt04;   //!
   TBranch        *b__sumNeutralHadronEt03;   //!
   TBranch        *b__sumChargedHadronPt03;   //!
   TBranch        *b__sumPhotonEt03;   //!
   TBranch        *b__trackIso;   //!
   TBranch        *b__ecalIso;   //!
   TBranch        *b__hcalIso;   //!
   TBranch        *b__ecalPFClusterIso;   //!
   TBranch        *b__hcalPFClusterIso;   //!
   TBranch        *b__relIso;   //!
   TBranch        *b__relIso0p4;   //!
   TBranch        *b__relIso0p4MuDeltaBeta;   //!
   TBranch        *b__ptRel;   //!
   TBranch        *b__ptRatio;   //!
   TBranch        *b__closestJetCsvV2;   //!
   TBranch        *b__closestJetDeepCsv_b;   //!
   TBranch        *b__closestJEC;   //!
   TBranch        *b__closest_lepAwareJetE;   //!
   TBranch        *b__closest_lepAwareJetPx;   //!
   TBranch        *b__closest_lepAwareJetPy;   //!
   TBranch        *b__closest_lepAwareJetPz;   //!
   TBranch        *b__closest_l1JetE;   //!
   TBranch        *b__closest_l1JetPx;   //!
   TBranch        *b__closest_l1JetPy;   //!
   TBranch        *b__closest_l1JetPz;   //!
   TBranch        *b__closest_lJetE;   //!
   TBranch        *b__closest_lJetPx;   //!
   TBranch        *b__closest_lJetPy;   //!
   TBranch        *b__closest_lJetPz;   //!
   TBranch        *b__closestJetDeepCsv_bb;   //!
   TBranch        *b__selectedTrackMult;   //!
   TBranch        *b__lMuonSegComp;   //!
   TBranch        *b__lMuonTrackPt;   //!
   TBranch        *b__lMuonTrackPtErr;   //!
   TBranch        *b__nJets;   //!
   TBranch        *b__jetPt;   //!
   TBranch        *b__jetPt_JECDown;   //!
   TBranch        *b__jetPt_JECUp;   //!
   TBranch        *b__jetSmearedPt;   //!
   TBranch        *b__jetSmearedPt_JECDown;   //!
   TBranch        *b__jetSmearedPt_JECUp;   //!
   TBranch        *b__jetSmearedPt_JERDown;   //!
   TBranch        *b__jetSmearedPt_JERUp;   //!
   TBranch        *b__jetPt_Uncorrected;   //!
   TBranch        *b__jetPt_L1;   //!
   TBranch        *b__jetPt_L2;   //!
   TBranch        *b__jetPt_L3;   //!
   TBranch        *b__jetEta;   //!
   TBranch        *b__jetPhi;   //!
   TBranch        *b__jetE;   //!
   TBranch        *b__jetCsvV2;   //!
   TBranch        *b__jetDeepCsv_udsg;   //!
   TBranch        *b__jetDeepCsv_b;   //!
   TBranch        *b__jetDeepCsv_c;   //!
   TBranch        *b__jetDeepCsv_bb;   //!
   TBranch        *b__jetHadronFlavor;   //!
   TBranch        *b__jetIsLoose;   //!
   TBranch        *b__jetIsTight;   //!
   TBranch        *b__jetIsTightLepVeto;   //!
   TBranch        *b__jetNeutralHadronFraction;   //!
   TBranch        *b__jetChargedHadronFraction;   //!
   TBranch        *b__jetNeutralEmFraction;   //!
   TBranch        *b__jetChargedEmFraction;   //!
   TBranch        *b__jetHFHadronFraction;   //!
   TBranch        *b__jetHFEmFraction;   //!
   TBranch        *b__met;   //!
   TBranch        *b__metRaw;   //!
   TBranch        *b__metJECDown;   //!
   TBranch        *b__metJECUp;   //!
   TBranch        *b__metUnclDown;   //!
   TBranch        *b__metUnclUp;   //!
   TBranch        *b__metPhi;   //!
   TBranch        *b__metRawPhi;   //!
   TBranch        *b__metPhiJECDown;   //!
   TBranch        *b__metPhiJECUp;   //!
   TBranch        *b__metPhiUnclDown;   //!
   TBranch        *b__metPhiUnclUp;   //!
   TBranch        *b__metSignificance;   //!
};
#endif
