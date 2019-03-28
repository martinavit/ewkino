#include <iostream>
#include <fstream>

#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"

treeReader::treeReader(TTree *tree) : fChain(nullptr) 
{
    if (tree != nullptr){
        initTree(tree);
    }
}

void treeReader::readSamples(const std::string& list, const std::string& directory, std::vector<Sample>& sampleVector){

    //clean current sample list 
    sampleVector.clear();

    //read list of samples from file
    sampleVector = readSampleList(list, directory);

    //print sample information
    for(auto& sample : sampleVector){
        std::cout << sample << std::endl;
    }
}

void treeReader::readSamples(const std::string& list, const std::string& directory){
    readSamples(list, directory, this->samples);
}

void treeReader::readSamples2016(const std::string& list, const std::string& directory){
    std::cout << "########################################" << std::endl;
    std::cout << "         2016 samples                   " << std::endl;
    std::cout << "########################################" << std::endl;

    readSamples(list, directory, this->samples2016);

    //add the 2016 samples to the total sample list 
    this->samples.insert(samples.end(), samples2016.begin(), samples2016.end() );

    //check for errors
    checkSampleEraConsistency();

}

void treeReader::readSamples2017(const std::string& list, const std::string& directory){
    std::cout << "########################################" << std::endl;
    std::cout << "         2017 samples                   " << std::endl;
    std::cout << "########################################" << std::endl;

    readSamples(list, directory, this->samples2017);

    //add the 2017 samples to the total sample list
    this->samples.insert(samples.end(), samples2017.begin(), samples2017.end() );

    //check for errors 
    checkSampleEraConsistency();
}

void treeReader::initSample(const Sample& samp){ 

    //update current sample
    currentSample = samp;
    sampleFile = samp.getFile();
    sampleFile->cd("blackJackAndHookers");
    fChain = (TTree*) sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree");
    initTree(fChain, samp.isData());
    nEntries = fChain->GetEntries();
    if(!samp.isData()){

        //read sum of simulated event weights
        TH1D* hCounter = new TH1D("hCounter", "Events counter", 1, 0, 1);
        hCounter->Read("hCounter"); 
        double sumSimulatedEventWeights = hCounter->GetBinContent(1);
        delete hCounter;

        //event weights set with lumi depending on sample's era 
        double dataLumi;
        if( is2016() ){
            dataLumi = lumi2016;
        } else {
            dataLumi = lumi2017;
        } 
        scale = samp.getXSec()*dataLumi*1000/sumSimulatedEventWeights;       //xSec*lumi divided by total sum of simulated event weights
    }
}

void treeReader::initSample(){ //initialize the next sample in the list 
    initSample(samples[++currentSampleIndex]);
}


void treeReader::GetEntry(const Sample& samp, long unsigned entry)
{
    if (!fChain) return;
    fChain->GetEntry(entry);
    //Set up correct weights
    if(!samp.isData() ) weight = _weight*scale; //MC
    else weight = 1;                            //data
}

void treeReader::GetEntry(long unsigned entry){    //currently initialized sample when running serial
    GetEntry(samples[currentSampleIndex], entry);
}

void treeReader::initTree(TTree *tree, const bool isData, const bool Is2016)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
   fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
   fChain->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);
   
   
   fChain->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
   fChain->SetBranchAddress("_Flag_goodVertices", &_Flag_goodVertices, &b__Flag_goodVertices);
   fChain->SetBranchAddress("_Flag_HBHENoiseFilter", &_Flag_HBHENoiseFilter, &b__Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("_Flag_HBHENoiseIsoFilter", &_Flag_HBHENoiseIsoFilter, &b__Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("_Flag_EcalDeadCellTriggerPrimitiveFilter", &_Flag_EcalDeadCellTriggerPrimitiveFilter, &b__Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("_Flag_BadPFMuonFilter", &_Flag_BadPFMuonFilter, &b__Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("_Flag_BadChargedCandidateFilter", &_Flag_BadChargedCandidateFilter, &b__Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("_updated_ecalBadCalibFilter", &_updated_ecalBadCalibFilter, &b__updated_ecalBadCalibFilter);
   fChain->SetBranchAddress("_passTrigger_1l", &_passTrigger_1l, &b__passTrigger_1l);
   if (!Is2016)fChain->SetBranchAddress("_HLT_IsoMu27", &_HLT_IsoMu27, &b__HLT_IsoMu27);
   if (!Is2016)fChain->SetBranchAddress("_HLT_IsoMu27_prescale", &_HLT_IsoMu27_prescale, &b__HLT_IsoMu27_prescale);
   if (Is2016) fChain->SetBranchAddress("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, &b__HLT_Ele27_WPTight_Gsf); 
   fChain->SetBranchAddress("_HLT_IsoMu24", &_HLT_IsoMu24, &b__HLT_IsoMu24);
   fChain->SetBranchAddress("_HLT_IsoMu24_prescale", &_HLT_IsoMu24_prescale, &b__HLT_IsoMu24_prescale);
   if (!Is2016)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf", &_HLT_Ele32_WPTight_Gsf, &b__HLT_Ele32_WPTight_Gsf);
   if (!Is2016)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf_prescale", &_HLT_Ele32_WPTight_Gsf_prescale, &b__HLT_Ele32_WPTight_Gsf_prescale);
   if (!Is2016)fChain->SetBranchAddress("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, &b__HLT_Ele35_WPTight_Gsf);
   if (!Is2016)fChain->SetBranchAddress("_HLT_Ele35_WPTight_Gsf_prescale", &_HLT_Ele35_WPTight_Gsf_prescale, &b__HLT_Ele35_WPTight_Gsf_prescale);
   fChain->SetBranchAddress("_nL", &_nL, &b__nL);
   fChain->SetBranchAddress("_nMu", &_nMu, &b__nMu);
   fChain->SetBranchAddress("_nEle", &_nEle, &b__nEle);
   fChain->SetBranchAddress("_nLight", &_nLight, &b__nLight);
   fChain->SetBranchAddress("_nTau", &_nTau, &b__nTau);
   fChain->SetBranchAddress("_pvX", &_pvX, &b__pvX);
   fChain->SetBranchAddress("_pvY", &_pvY, &b__pvY);
   fChain->SetBranchAddress("_pvZ", &_pvZ, &b__pvZ);
   fChain->SetBranchAddress("_pvXErr", &_pvXErr, &b__pvXErr);
   fChain->SetBranchAddress("_pvYErr", &_pvYErr, &b__pvYErr);
   fChain->SetBranchAddress("_pvZErr", &_pvZErr, &b__pvZErr);
   fChain->SetBranchAddress("_nVFit_os", &_nVFit_os, &b__nVFit_os);
   fChain->SetBranchAddress("_nVFit", &_nVFit, &b__nVFit);
   fChain->SetBranchAddress("_nGoodLeading", &_nGoodLeading, &b__nGoodLeading);
   fChain->SetBranchAddress("_nGoodDisplaced", &_nGoodDisplaced, &b__nGoodDisplaced);
   fChain->SetBranchAddress("_vertices_os", _vertices_os, &b__vertices_os);
   fChain->SetBranchAddress("_lDisplaced_os", _lDisplaced_os, &b__lDisplaced_os);
   fChain->SetBranchAddress("_vertices", _vertices, &b__vertices);
   fChain->SetBranchAddress("_lDisplaced", _lDisplaced, &b__lDisplaced);
   fChain->SetBranchAddress("_lHasTrigger", _lHasTrigger, &b__lHasTrigger);
   fChain->SetBranchAddress("_lPt", _lPt, &b__lPt);
   fChain->SetBranchAddress("_lEta", _lEta, &b__lEta);
   fChain->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
   fChain->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
   fChain->SetBranchAddress("_lE", _lE, &b__lE);
   fChain->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
   fChain->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
   fChain->SetBranchAddress("_dxy", _dxy, &b__dxy);
   fChain->SetBranchAddress("_dz", _dz, &b__dz);
   fChain->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
   fChain->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
   fChain->SetBranchAddress("_2dIP", _2dIP, &b__2dIP);
   fChain->SetBranchAddress("_2dIPSig", _2dIPSig, &b__2dIPSig);
   fChain->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
   fChain->SetBranchAddress("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, &b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto);
   fChain->SetBranchAddress("_lElectronPassConvVeto", _lElectronPassConvVeto, &b__lElectronPassConvVeto);
   fChain->SetBranchAddress("_lElectronChargeConst", _lElectronChargeConst, &b__lElectronChargeConst);
   fChain->SetBranchAddress("_lElectronMissingHits", _lElectronMissingHits, &b__lElectronMissingHits);
   fChain->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
   fChain->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
   fChain->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
   fChain->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
   fChain->SetBranchAddress("_lGlobalMuon", _lGlobalMuon, &b__lGlobalMuon);
   fChain->SetBranchAddress("_lTrackerMuon", _lTrackerMuon, &b__lTrackerMuon);
   fChain->SetBranchAddress("_lInnerTrackValidFraction", _lInnerTrackValidFraction, &b__lInnerTrackValidFraction);
   fChain->SetBranchAddress("_lGlobalTrackNormalizeChi2", _lGlobalTrackNormalizeChi2, &b__lGlobalTrackNormalizeChi2);
   fChain->SetBranchAddress("_lCQChi2Position", _lCQChi2Position, &b__lCQChi2Position);
   fChain->SetBranchAddress("_lCQTrackKink", _lCQTrackKink, &b__lCQTrackKink);
   fChain->SetBranchAddress("_lNumberOfMatchedStation", _lNumberOfMatchedStation, &b__lNumberOfMatchedStation);
   fChain->SetBranchAddress("_lNumberOfValidPixelHits", _lNumberOfValidPixelHits, &b__lNumberOfValidPixelHits);
   fChain->SetBranchAddress("_lTrackerLayersWithMeasurement", _lTrackerLayersWithMeasurement, &b__lTrackerLayersWithMeasurement);
   fChain->SetBranchAddress("_lSimType", _lSimType, &b__lSimType);
   fChain->SetBranchAddress("_lSimExtType", _lSimExtType, &b__lSimExtType);
   fChain->SetBranchAddress("_lSimFlavour", _lSimFlavour, &b__lSimFlavour);
   fChain->SetBranchAddress("_muDTStationsWithValidHits", _muDTStationsWithValidHits, &b__muDTStationsWithValidHits);
   fChain->SetBranchAddress("_muCSCStationsWithValidHits", _muCSCStationsWithValidHits, &b__muCSCStationsWithValidHits);
   fChain->SetBranchAddress("_muRPCStationsWithValidHits", _muRPCStationsWithValidHits, &b__muRPCStationsWithValidHits);
   fChain->SetBranchAddress("_muMuonStationsWithValidHits", _muMuonStationsWithValidHits, &b__muMuonStationsWithValidHits);
   fChain->SetBranchAddress("_lMuRPCTimenDof", _lMuRPCTimenDof, &b__lMuRPCTimenDof);
   fChain->SetBranchAddress("_lMuTimenDof", _lMuTimenDof, &b__lMuTimenDof);
   fChain->SetBranchAddress("_lMuRPCTime", _lMuRPCTime, &b__lMuRPCTime);
   fChain->SetBranchAddress("_lMuRPCTimeErr", _lMuRPCTimeErr, &b__lMuRPCTimeErr);
   fChain->SetBranchAddress("_lMuTime", _lMuTime, &b__lMuTime);
   fChain->SetBranchAddress("_lMuTimeErr", _lMuTimeErr, &b__lMuTimeErr);
   fChain->SetBranchAddress("_muNumberInnerHits", _muNumberInnerHits, &b__muNumberInnerHits);
   fChain->SetBranchAddress("_lEleIsEB", _lEleIsEB, &b__lEleIsEB);
   fChain->SetBranchAddress("_lEleIsEE", _lEleIsEE, &b__lEleIsEE);
   fChain->SetBranchAddress("_lEleSuperClusterOverP", _lEleSuperClusterOverP, &b__lEleSuperClusterOverP);
   fChain->SetBranchAddress("_lEleEcalEnergy", _lEleEcalEnergy, &b__lEleEcalEnergy);
   fChain->SetBranchAddress("_lElefull5x5SigmaIetaIeta", _lElefull5x5SigmaIetaIeta, &b__lElefull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("_lEleDEtaInSeed", _lEleDEtaInSeed, &b__lEleDEtaInSeed);
   fChain->SetBranchAddress("_lEleDeltaPhiSuperClusterTrackAtVtx", _lEleDeltaPhiSuperClusterTrackAtVtx, &b__lEleDeltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("_lElehadronicOverEm", _lElehadronicOverEm, &b__lElehadronicOverEm);
   fChain->SetBranchAddress("_lEleInvMinusPInv", _lEleInvMinusPInv, &b__lEleInvMinusPInv);
   fChain->SetBranchAddress("_puCorr", _puCorr, &b__puCorr);
   fChain->SetBranchAddress("_absIso03", _absIso03, &b__absIso03);
   fChain->SetBranchAddress("_absIso04", _absIso04, &b__absIso04);
   fChain->SetBranchAddress("_sumNeutralHadronEt04", _sumNeutralHadronEt04, &b__sumNeutralHadronEt04);
   fChain->SetBranchAddress("_sumChargedHadronPt04", _sumChargedHadronPt04, &b__sumChargedHadronPt04);
   fChain->SetBranchAddress("_sumPhotonEt04", _sumPhotonEt04, &b__sumPhotonEt04);
   fChain->SetBranchAddress("_sumNeutralHadronEt03", _sumNeutralHadronEt03, &b__sumNeutralHadronEt03);
   fChain->SetBranchAddress("_sumChargedHadronPt03", _sumChargedHadronPt03, &b__sumChargedHadronPt03);
   fChain->SetBranchAddress("_sumPhotonEt03", _sumPhotonEt03, &b__sumPhotonEt03);
   fChain->SetBranchAddress("_trackIso", _trackIso, &b__trackIso);
   fChain->SetBranchAddress("_ecalIso", _ecalIso, &b__ecalIso);
   fChain->SetBranchAddress("_hcalIso", _hcalIso, &b__hcalIso);
   fChain->SetBranchAddress("_ecalPFClusterIso", _ecalPFClusterIso, &b__ecalPFClusterIso);
   fChain->SetBranchAddress("_hcalPFClusterIso", _hcalPFClusterIso, &b__hcalPFClusterIso);
   fChain->SetBranchAddress("_relIso", _relIso, &b__relIso);
   fChain->SetBranchAddress("_relIso0p4", _relIso0p4, &b__relIso0p4);
   fChain->SetBranchAddress("_relIso0p4MuDeltaBeta", _relIso0p4MuDeltaBeta, &b__relIso0p4MuDeltaBeta);
   fChain->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
   fChain->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
   fChain->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
   fChain->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
   fChain->SetBranchAddress("_closestJEC", _closestJEC, &b__closestJEC);
   fChain->SetBranchAddress("_closest_lepAwareJetE", _closest_lepAwareJetE, &b__closest_lepAwareJetE);
   fChain->SetBranchAddress("_closest_lepAwareJetPx", _closest_lepAwareJetPx, &b__closest_lepAwareJetPx);
   fChain->SetBranchAddress("_closest_lepAwareJetPy", _closest_lepAwareJetPy, &b__closest_lepAwareJetPy);
   fChain->SetBranchAddress("_closest_lepAwareJetPz", _closest_lepAwareJetPz, &b__closest_lepAwareJetPz);
   fChain->SetBranchAddress("_closest_l1JetE", _closest_l1JetE, &b__closest_l1JetE);
   fChain->SetBranchAddress("_closest_l1JetPx", _closest_l1JetPx, &b__closest_l1JetPx);
   fChain->SetBranchAddress("_closest_l1JetPy", _closest_l1JetPy, &b__closest_l1JetPy);
   fChain->SetBranchAddress("_closest_l1JetPz", _closest_l1JetPz, &b__closest_l1JetPz);
   fChain->SetBranchAddress("_closest_lJetE", _closest_lJetE, &b__closest_lJetE);
   fChain->SetBranchAddress("_closest_lJetPx", _closest_lJetPx, &b__closest_lJetPx);
   fChain->SetBranchAddress("_closest_lJetPy", _closest_lJetPy, &b__closest_lJetPy);
   fChain->SetBranchAddress("_closest_lJetPz", _closest_lJetPz, &b__closest_lJetPz);
   fChain->SetBranchAddress("_closestJetDeepCsv_bb", _closestJetDeepCsv_bb, &b__closestJetDeepCsv_bb);
   fChain->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
   fChain->SetBranchAddress("_lMuonSegComp", _lMuonSegComp, &b__lMuonSegComp);
   fChain->SetBranchAddress("_lMuonTrackPt", _lMuonTrackPt, &b__lMuonTrackPt);
   fChain->SetBranchAddress("_lMuonTrackPtErr", _lMuonTrackPtErr, &b__lMuonTrackPtErr);   
   fChain->SetBranchAddress("_lPtCorr", _lPtCorr, &b__lPtCorr);
   fChain->SetBranchAddress("_lPtScaleUp", _lPtScaleUp, &b__lPtScaleUp);
   fChain->SetBranchAddress("_lPtScaleDown", _lPtScaleDown, &b__lPtScaleDown);
   fChain->SetBranchAddress("_lPtResUp", _lPtResUp, &b__lPtResUp);
   fChain->SetBranchAddress("_lPtResDown", _lPtResDown, &b__lPtResDown);
   fChain->SetBranchAddress("_lECorr", _lECorr, &b__lECorr);
   fChain->SetBranchAddress("_lEScaleUp", _lEScaleUp, &b__lEScaleUp);
   fChain->SetBranchAddress("_lEScaleDown", _lEScaleDown, &b__lEScaleDown);
   fChain->SetBranchAddress("_lEResUp", _lEResUp, &b__lEResUp);
   fChain->SetBranchAddress("_lEResDown", _lEResDown, &b__lEResDown);
   fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
   fChain->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
   fChain->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
   fChain->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
   fChain->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
   fChain->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
   fChain->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
   fChain->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
   fChain->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
   fChain->SetBranchAddress("_jetPt_Uncorrected", _jetPt_Uncorrected, &b__jetPt_Uncorrected);
   fChain->SetBranchAddress("_jetPt_L1", _jetPt_L1, &b__jetPt_L1);
   fChain->SetBranchAddress("_jetPt_L2", _jetPt_L2, &b__jetPt_L2);
   fChain->SetBranchAddress("_jetPt_L3", _jetPt_L3, &b__jetPt_L3);
   fChain->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
   fChain->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
   fChain->SetBranchAddress("_jetE", _jetE, &b__jetE);
   fChain->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
   fChain->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
   fChain->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
   fChain->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
   fChain->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
   fChain->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
   fChain->SetBranchAddress("_jetIsLoose", _jetIsLoose, &b__jetIsLoose);
   fChain->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
   fChain->SetBranchAddress("_jetIsTightLepVeto", _jetIsTightLepVeto, &b__jetIsTightLepVeto);
   fChain->SetBranchAddress("_jetNeutralHadronFraction", _jetNeutralHadronFraction, &b__jetNeutralHadronFraction);
   fChain->SetBranchAddress("_jetChargedHadronFraction", _jetChargedHadronFraction, &b__jetChargedHadronFraction);
   fChain->SetBranchAddress("_jetNeutralEmFraction", _jetNeutralEmFraction, &b__jetNeutralEmFraction);
   fChain->SetBranchAddress("_jetChargedEmFraction", _jetChargedEmFraction, &b__jetChargedEmFraction);
   fChain->SetBranchAddress("_jetHFHadronFraction", _jetHFHadronFraction, &b__jetHFHadronFraction);
   fChain->SetBranchAddress("_jetHFEmFraction", _jetHFEmFraction, &b__jetHFEmFraction);
   fChain->SetBranchAddress("_met", &_met, &b__met);
   fChain->SetBranchAddress("_metRaw", &_metRaw, &b__metRaw);
   fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
   fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
   fChain->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
   fChain->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
   fChain->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
   fChain->SetBranchAddress("_metRawPhi", &_metRawPhi, &b__metRawPhi);
   fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
   fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
   fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
   fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
   fChain->SetBranchAddress("_metSignificance", &_metSignificance, &b__metSignificance);

    if(!isData){
        fChain->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
        fChain->SetBranchAddress("_weight", &_weight, &b__weight);
        fChain->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
        fChain->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
        fChain->SetBranchAddress("_nLheTau", &_nLheTau, &b__nLheTau);
        fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
        fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
        fChain->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
        fChain->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
        fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
        fChain->SetBranchAddress("_gen_pdgID", _gen_pdgID, &b__gen_pdgID);
        fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
        fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
        fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
        fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
        fChain->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
        fChain->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
        fChain->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
        fChain->SetBranchAddress("_gen_vertex_x", _gen_vertex_x, &b__gen_vertex_x);
        fChain->SetBranchAddress("_gen_vertex_y", _gen_vertex_y, &b__gen_vertex_y);
        fChain->SetBranchAddress("_gen_vertex_z", _gen_vertex_z, &b__gen_vertex_z);
        fChain->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
        fChain->SetBranchAddress("_gen_lMinDeltaR", _gen_lMinDeltaR, &b__gen_lMinDeltaR);
        fChain->SetBranchAddress("_gen_lPassParentage", _gen_lPassParentage, &b__gen_lPassParentage);
        
         fChain->SetBranchAddress("_lGenIndex", _lGenIndex, &b__lGenIndex);
        fChain->SetBranchAddress("_lMatchType", _lMatchType, &b__lMatchType);
        fChain->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
        fChain->SetBranchAddress("_lIsPromptFinalState", _lIsPromptFinalState, &b__lIsPromptFinalState);
         fChain->SetBranchAddress("_lIsPromptDecayed", _lIsPromptDecayed, &b__lIsPromptDecayed);
         fChain->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
        fChain->SetBranchAddress("_lMomPdgId", _lMomPdgId, &b__lMomPdgId);
        fChain->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
        fChain->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
         fChain->SetBranchAddress("_lProvenanceConversion", _lProvenanceConversion, &b__lProvenanceConversion);
        fChain->SetBranchAddress("_lMatchPt", _lMatchPt, &b__lMatchPt);
         fChain->SetBranchAddress("_lMatchEta", _lMatchEta, &b__lMatchEta);
         fChain->SetBranchAddress("_lMatchPhi", _lMatchPhi, &b__lMatchPhi);
        fChain->SetBranchAddress("_lMatchVertexX", _lMatchVertexX, &b__lMatchVertexX);
        fChain->SetBranchAddress("_lMatchVertexY", _lMatchVertexY, &b__lMatchVertexY);
        fChain->SetBranchAddress("_lMatchVertexZ", _lMatchVertexZ, &b__lMatchVertexZ);
        
        
        
    }
}

void treeReader::setOutputTree(TTree* outputTree, const bool isData, const bool Is2016){
    outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
    outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
    outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
    outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/b");
    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");
    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    outputTree->Branch("_metSignificance",              &_metSignificance,              "_metSignificance/D");
    outputTree->Branch("_passTrigger_1l", &_passTrigger_1l, "_passTrigger_1l/O");
    outputTree->Branch("_HLT_IsoMu24", &_HLT_IsoMu24, "_HLT_IsoMu24/O");
    if (Is2016) utputTree->Branch("_HLT_IsoTkMu24", &_HLT_IsoTkMu24, "_HLT_IsoTkMu24/O");
    if (Is2016) outputTree->Branch("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, "_HLT_Ele27_WPTight_Gsf/O");
        
    if (!Is2016)outputTree->Branch("_HLT_IsoMu27", &_HLT_IsoMu27, "_HLT_IsoMu27/O");
    if (!Is2016)outputTree->Branch("_HLT_Ele32_WPTight_Gsf", &_HLT_Ele32_WPTight_Gsf, "_HLT_Ele32_WPTight_Gsf/O");
    if (!Is2016)outputTree->Branch("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, "_HLT_Ele35_WPTight_Gsf/O");
    if (!Is2016)outputTree->Branch("_HLT_Ele32_WPTight_Gsf_L1DoubleEG", &_HLT_Ele32_WPTight_Gsf_L1DoubleEG, "_HLT_Ele32_WPTight_Gsf_L1DoubleEG/O");

    outputTree->Branch("_passMETFilters", &_passMETFilters, "_passMETFilters/O");
    //TEMPORARY FOR CHECK, CAN BE REMOVED LATER
    outputTree->Branch("_Flag_BadPFMuonFilter", &_Flag_BadPFMuonFilter, "_Flag_BadPFMuonFilter/O");
    outputTree->Branch("_Flag_BadChargedCandidateFilter", &_Flag_BadChargedCandidateFilter, "_Flag_BadChargedCandidateFilter/O");
    //////////////////////////////////////////
    outputTree->Branch("_nL",                           &_nL,                           "_nL/i");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/i");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/i");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/i");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/i");
    outputTree->Branch("_pvX",                          &_pvX,                          "_pvX/D");       // displaced specific [TODO: it seems these should be moved to multilep.cc or so]
    outputTree->Branch("_pvY",                          &_pvY,                          "_pvY/D");       // "
    outputTree->Branch("_pvZ",                          &_pvZ,                          "_pvZ/D");       // "
    outputTree->Branch("_pvXErr",                       &_pvXErr,                       "_pvXErr/D");    // "
    outputTree->Branch("_pvYErr",                       &_pvYErr,                       "_pvYErr/D");    // "
    outputTree->Branch("_pvZErr",                       &_pvZErr,                       "_pvZErr/D");    // "
    outputTree->Branch("_nVFit",                        &_nVFit,                        "_nVFit/i");                   // displaced specific
   outputTree->Branch("_nVFit_os",                        &_nVFit_os,                        "_nVFit_os/i");                   // displaced specific
    outputTree->Branch("_nGoodLeading",                 &_nGoodLeading,                 "_nGoodLeading/i");            // "
    outputTree->Branch("_nGoodDisplaced",               &_nGoodDisplaced,               "_nGoodDisplaced/i");          // "
    outputTree->Branch("_vertices_os",                     &_vertices_os,                     "_vertices_os[_nVFit_os][12]/D");    // displaced specific
    outputTree->Branch("_lDisplaced_os",                   &_lDisplaced_os,                   "_lDisplaced_os[_nVFit_os][24]/D");  // "
    outputTree->Branch("_vertices",                     &_vertices,                     "_vertices[_nVFit][12]/D");    // displaced specific
    outputTree->Branch("_lDisplaced",                   &_lDisplaced,                   "_lDisplaced[_nVFit][24]/D");  // "
    outputTree->Branch("_lHasTrigger",                  &_lHasTrigger,                  "_lHasTrigger[_nL]/i");        // "
    outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
    outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
    outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
    outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
    outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
    outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
    outputTree->Branch("_2dIP",                         &_2dIP,                         "_2dIP[_nL]/D");          // displaced specific
    outputTree->Branch("_2dIPSig",                      &_2dIPSig,                      "_2dIPSig[_nL]/D");       // "
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", &_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, "_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[_nL]/O"); // displaced specific (still needed ???)
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");       // DELETED (same as the original lElectronMissingHits)
    outputTree->Branch("_lGlobalMuon",                  &_lGlobalMuon,                  "_lGlobalMuon[_nMu]/O");                    // displaced specific
    outputTree->Branch("_lTrackerMuon",                 &_lTrackerMuon,                 "_lTrackerMuon[_nMu]/O");                   // "
    outputTree->Branch("_lInnerTrackValidFraction",     &_lInnerTrackValidFraction,     "_lInnerTrackValidFraction[_nMu]/D");       // "
    outputTree->Branch("_lGlobalTrackNormalizeChi2",    &_lGlobalTrackNormalizeChi2,    "_lGlobalTrackNormalizeChi2[_nMu]/D");      // "
    outputTree->Branch("_lCQChi2Position",              &_lCQChi2Position,              "_lCQChi2Position[_nMu]/D");                // "
    outputTree->Branch("_lCQTrackKink",                 &_lCQTrackKink,                 "_lCQTrackKink[_nMu]/D");                   // "
    outputTree->Branch("_lNumberOfMatchedStation",      &_lNumberOfMatchedStation,      "_lNumberOfMatchedStation[_nMu]/i");        // "
    outputTree->Branch("_lNumberOfValidPixelHits",      &_lNumberOfValidPixelHits,      "_lNumberOfValidPixelHits[_nMu]/i");        // "
    outputTree->Branch("_lTrackerLayersWithMeasurement",&_lTrackerLayersWithMeasurement,"_lTrackerLayersWithMeasurement[_nMu]/i");  // "
    outputTree->Branch("_lSimType",                     &_lSimType,                     "_lSimType[_nMu]/I");                       // "
    outputTree->Branch("_lSimExtType",                  &_lSimExtType,                  "_lSimExtType[_nMu]/I");                    // "
    outputTree->Branch("_lSimFlavour",                  &_lSimFlavour,                  "_lSimFlavour[_nMu]/I");                    // "
    outputTree->Branch("_muDTStationsWithValidHits",    &_muDTStationsWithValidHits,    "_muDTStationsWithValidHits[_nMu]/I");      // "
    outputTree->Branch("_muCSCStationsWithValidHits",   &_muCSCStationsWithValidHits,   "_muCSCStationsWithValidHits[_nMu]/I");     // "
    outputTree->Branch("_muRPCStationsWithValidHits",   &_muRPCStationsWithValidHits,   "_muRPCStationsWithValidHits[_nMu]/I");     // "
    outputTree->Branch("_muMuonStationsWithValidHits",  &_muMuonStationsWithValidHits,  "_muMuonStationsWithValidHits[_nMu]/I");    // "
    outputTree->Branch("_lMuRPCTimenDof",               &_lMuRPCTimenDof,               "_lMuRPCTimenDof[_nMu]/I");                 // "
    outputTree->Branch("_lMuTimenDof",                  &_lMuTimenDof,                  "_lMuTimenDof[_nMu]/I");                    // "
    outputTree->Branch("_lMuRPCTime",                   &_lMuRPCTime,                   "_lMuRPCTime[_nMu]/D");                     // "
    outputTree->Branch("_lMuRPCTimeErr",                &_lMuRPCTimeErr,                "_lMuRPCTimeErr[_nMu]/D");                  // "
    outputTree->Branch("_lMuTime",                      &_lMuTime,                      "_lMuTime[_nMu]/D");                        // "
    outputTree->Branch("_lMuTimeErr",                   &_lMuTimeErr,                   "_lMuTimeErr[_nMu]/D");                     // "
    outputTree->Branch("_muNumberInnerHits",            &_muNumberInnerHits,            "_muNumberInnerHits[_nMu]/i");              // "
    outputTree->Branch("_lEleIsEB",                     &_lEleIsEB ,                    "_lEleIsEB[_nLight]/O");                    // "
    outputTree->Branch("_lEleIsEE",                     &_lEleIsEE ,                    "_lEleIsEE[_nLight]/O");                    // "
    outputTree->Branch("_lEleSuperClusterOverP",        &_lEleSuperClusterOverP ,       "_lEleSuperClusterOverP[_nLight]/D");       // "
    outputTree->Branch("_lEleEcalEnergy",               &_lEleEcalEnergy ,              "_lEleEcalEnergy[_nLight]/D");              // "
    outputTree->Branch("_lElefull5x5SigmaIetaIeta",     &_lElefull5x5SigmaIetaIeta ,    "_lElefull5x5SigmaIetaIeta[_nLight]/D");    // "
    outputTree->Branch("_lEleDEtaInSeed",               &_lEleDEtaInSeed ,              "_lEleDEtaInSeed[_nLight]/D");              // "
    outputTree->Branch("_lEleDeltaPhiSuperClusterTrackAtVtx", &_lEleDeltaPhiSuperClusterTrackAtVtx , "_lEleDeltaPhiSuperClusterTrackAtVtx[_nLight]/D"); // "
    outputTree->Branch("_lElehadronicOverEm",           &_lElehadronicOverEm ,          "_lElehadronicOverEm[_nLight]/D");          // "
    outputTree->Branch("_lEleInvMinusPInv",             &_lEleInvMinusPInv ,            "_lEleInvMinusPInv[_nLight]/D");            // "
    outputTree->Branch("_puCorr",                       &_puCorr,                       "_puCorr[_nLight]/D");                      // "
    outputTree->Branch("_absIso03",                     &_absIso03,                     "_absIso03[_nL]/D");                        // "
    outputTree->Branch("_absIso04",                     &_absIso04,                     "_absIso04[_nMu]/D");                       // "
    outputTree->Branch("_sumNeutralHadronEt04",         &_sumNeutralHadronEt04,         "_sumNeutralHadronEt04[_nMu]/D");           // "
    outputTree->Branch("_sumChargedHadronPt04",         &_sumChargedHadronPt04,         "_sumChargedHadronPt04[_nMu]/D");           // "
    outputTree->Branch("_sumPhotonEt04",                &_sumPhotonEt04,                "_sumPhotonEt04[_nMu]/D");                  // "
    outputTree->Branch("_sumNeutralHadronEt03",         &_sumNeutralHadronEt03,         "_sumNeutralHadronEt03[_nLight]/D");        // "
    outputTree->Branch("_sumChargedHadronPt03",         &_sumChargedHadronPt03,         "_sumChargedHadronPt03[_nLight]/D");        // "
    outputTree->Branch("_sumPhotonEt03",                &_sumPhotonEt03,                "_sumPhotonEt03[_nLight]/D");               // "
    outputTree->Branch("_trackIso",                     &_trackIso ,                    "_trackIso[_nLight]/D");                    // "
    outputTree->Branch("_ecalIso",                      &_ecalIso ,                     "_ecalIso[_nLight]/D");                     // "
    outputTree->Branch("_hcalIso",                      &_hcalIso ,                     "_hcalIso[_nLight]/D");                     // "
    outputTree->Branch("_ecalPFClusterIso",             &_ecalPFClusterIso ,            "_ecalPFClusterIso[_nLight]/D");            // displaced specific
    outputTree->Branch("_hcalPFClusterIso",             &_hcalPFClusterIso ,            "_hcalPFClusterIso[_nLight]/D");            // "
    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_relIso0p4",                    &_relIso0p4,                    "_relIso0p4[_nLight]/D");
    outputTree->Branch("_relIso0p4MuDeltaBeta",         &_relIso0p4MuDeltaBeta,         "_relIso0p4MuDeltaBeta[_nMu]/D");  
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJEC",          &_closestJEC,          "_closestJEC[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetE",          &_closest_lepAwareJetE,          "_closest_lepAwareJetE[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPx",          &_closest_lepAwareJetPx,          "_closest_lepAwareJetPx[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPy",          &_closest_lepAwareJetPy,          "_closest_lepAwareJetPy[_nLight]/D");
     outputTree->Branch("_closest_lepAwareJetPz",          &_closest_lepAwareJetPz,          "_closest_lepAwareJetPz[_nLight]/D");
    outputTree->Branch("_closest_l1JetE",          &_closest_l1JetE,          "_closest_l1JetE[_nLight]/D");
    outputTree->Branch("_closest_l1JetPx",          &_closest_l1JetPx,          "_closest_l1JetPx[_nLight]/D");
    outputTree->Branch("_closest_l1JetPy",          &_closest_l1JetPy,          "_closest_l1JetPy[_nLight]/D");
    outputTree->Branch("_closest_l1JetPz",          &_closest_l1JetPz,          "_closest_l1JetPz[_nLight]/D");
    outputTree->Branch("_closest_lJetE",          &_closest_lJetE,          "_closest_lJetE[_nLight]/D");
    outputTree->Branch("_closest_lJetPx",          &_closest_lJetPx,          "_closest_lJetPx[_nLight]/D");
    outputTree->Branch("_closest_lJetPy",          &_closest_lJetPy,          "_closest_lJetPy[_nLight]/D");
    outputTree->Branch("_closest_lJetPz",          &_closest_lJetPz,          "_closest_lJetPz[_nLight]/D");  
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/i");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetPt_Uncorrected",         &_jetPt_Uncorrected,        "_jetPt_Uncorrected[_nJets]/D");
    outputTree->Branch("_jetPt_L1",                  &_jetPt_L1,                 "_jetPt_L1[_nJets]/D");
    outputTree->Branch("_jetPt_L2",                  &_jetPt_L2,                 "_jetPt_L2[_nJets]/D");
    outputTree->Branch("_jetPt_L3",                  &_jetPt_L3,                 "_jetPt_L3[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    outputTree->Branch("_jetIsLoose",                &_jetIsLoose,               "_jetIsLoose[_nJets]/O");
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");
    outputTree->Branch("_jetNeutralHadronFraction",  &_jetNeutralHadronFraction, "_jetNeutralHadronFraction[_nJets]/D");
    outputTree->Branch("_jetChargedHadronFraction",  &_jetChargedHadronFraction, "_jetChargedHadronFraction[_nJets]/D");
    outputTree->Branch("_jetNeutralEmFraction",      &_jetNeutralEmFraction,     "_jetNeutralEmFraction[_nJets]/D");
    outputTree->Branch("_jetChargedEmFraction",      &_jetChargedEmFraction,     "_jetChargedEmFraction[_nJets]/D");
    outputTree->Branch("_jetHFHadronFraction",       &_jetHFHadronFraction,      "_jetHFHadronFraction[_nJets]/D");
    outputTree->Branch("_jetHFEmFraction",           &_jetHFEmFraction,          "_jetHFEmFraction[_nJets]/D");

    if(!isData){
        outputTree->Branch("_nTrueInt",                  &_nTrueInt,                  "_nTrueInt/F");
        outputTree->Branch("_nLheWeights",               &_nLheWeights,               "_nLheWeights/i");
        outputTree->Branch("_lheWeight",                 &_lheWeight,                 "_lheWeight[_nLheWeights]/D");
        outputTree->Branch("_weight",                    &_weight,                    "_weight/D");
        outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
        outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/b");
        outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
        outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
        outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
        outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
        outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/i");
        outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
        outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
        outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");
        outputTree->Branch("_lGenIndex",                  &_lGenIndex,                    "_lGenIndex[_nL]/i");
      outputTree->Branch("_lMatchType",                 &_lMatchType,                   "_lMatchType[_nL]/i");
      outputTree->Branch("_lIsPrompt",                  &_lIsPrompt,                    "_lIsPrompt[_nL]/O");
      outputTree->Branch("_lIsPromptFinalState",        &_lIsPromptFinalState,          "_lIsPromptFinalState[_nL]/O");
      outputTree->Branch("_lIsPromptDecayed",           &_lIsPromptDecayed,             "_lIsPromptDecayed[_nL]/O");
      outputTree->Branch("_lMatchPdgId",                &_lMatchPdgId,                  "_lMatchPdgId[_nL]/I");
      outputTree->Branch("_lMomPdgId",                  &_lMomPdgId,                    "_lMomPdgId[_nL]/I");
      outputTree->Branch("_lProvenance",                &_lProvenance,                  "_lProvenance[_nL]/i");
      outputTree->Branch("_lProvenanceCompressed",      &_lProvenanceCompressed,        "_lProvenanceCompressed[_nL]/i");
      outputTree->Branch("_lProvenanceConversion",      &_lProvenanceConversion,        "_lProvenanceConversion[_nL]/i");
      outputTree->Branch("_lMatchPt",                   &_lMatchPt,                     "_lMatchPt[_nL]/D");
      outputTree->Branch("_lMatchEta",                  &_lMatchEta,                    "_lMatchEta[_nL]/D");
      outputTree->Branch("_lMatchPhi",                  &_lMatchPhi,                    "_lMatchPhi[_nL]/D");
      outputTree->Branch("_lMatchVertexX",              &_lMatchVertexX,                "_lMatchVertexX[_nL]/D");
      outputTree->Branch("_lMatchVertexY",              &_lMatchVertexY,                "_lMatchVertexY[_nL]/D");
      outputTree->Branch("_lMatchVertexZ",              &_lMatchVertexZ,                "_lMatchVertexZ[_nL]/D");
    }
}
