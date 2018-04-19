//include ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <algorithm>

//include other parts of the code
#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"


void treeReader::skimTree(const std::string& fileName, std::string outputDirectory, const bool isData){//std::string outputFileName){
    //Read tree	
    std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>( (const TString& ) fileName,"read");	
    sampleFile->cd("blackJackAndHookers");
    //Determine hcounter and lheCounter for MC cross section scaling and uncertainties
    TH1D* hCounter;
    TH1D* lheCounter;
    TH1D* nTrueInt;
    TH1D* nVertices;
    if(!isData){
        hCounter = (TH1D*) sampleFile->Get("blackJackAndHookers/hCounter");
        lheCounter = (TH1D*) sampleFile->Get("blackJackAndHookers/lheCounter");
        nTrueInt = (TH1D*) sampleFile->Get("blackJackAndHookers/nTrueInteractions");
    }
    nVertices = (TH1D*) sampleFile->Get("blackJackAndHookers/nVertices");
    //Get Tree
    TTree* sampleTree = (TTree*) (sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree"));
    initTree(sampleTree, isData);
    outputDirectory = (outputDirectory == "") ? "~/Work/ntuples_temp/" : outputDirectory;
    std::string outputFileName = fileName;
    /*
    auto it = outputFileName.find_last_of("/");
    outputFileName.erase(0, it + 1);
    */
    outputFileName.erase(std::remove(outputFileName.begin(), outputFileName.end(), '/'), outputFileName.end());
    outputFileName.insert(0, outputDirectory);
    auto it = outputFileName.find(".root");
    outputFileName.erase(it, outputFileName.size());
    outputFileName.append("_trilepSkim.root");
    std::cout << "output file : " << outputFileName << std::endl;
    
    std::shared_ptr<TFile> outputFile = std::make_shared<TFile>((const TString&) outputFileName ,"RECREATE");
    outputFile->mkdir("blackJackAndHookers");
    outputFile->cd("blackJackAndHookers"); 
    TTree* outputTree = new TTree("blackJackAndHookersTree","blackJackAndHookersTree");
    setOutputTree(outputTree, isData);

    //TEMPORARY MVA READER: REMOVE THIS LATER
    float pt, eta, trackMultClosestJet, miniIsoCharged, miniIsoNeutral, ptRel, ptRatio, relIso, deepCsvClosestJet, sip3d, dxy, dz;
    std::shared_ptr<TMVA::Reader> leptonMvaReader[2][2];
    for(unsigned era = 0; era < 2; ++era){
        for(unsigned flavor = 0; flavor < 2; ++flavor){
            leptonMvaReader[era][flavor] = std::make_shared<TMVA::Reader>( "!Color:!Silent");
            leptonMvaReader[era][flavor]->AddVariable( "pt", &pt );
            leptonMvaReader[era][flavor]->AddVariable( "eta", &eta );
            leptonMvaReader[era][flavor]->AddVariable( "trackMultClosestJet", &trackMultClosestJet );
            leptonMvaReader[era][flavor]->AddVariable( "miniIsoCharged", &miniIsoCharged );
            leptonMvaReader[era][flavor]->AddVariable( "miniIsoNeutral", &miniIsoNeutral);
            leptonMvaReader[era][flavor]->AddVariable( "pTRel", &ptRel);
            leptonMvaReader[era][flavor]->AddVariable( "ptRatio", &ptRatio);
            leptonMvaReader[era][flavor]->AddVariable( "relIso", &relIso);
            leptonMvaReader[era][flavor]->AddVariable( "deepCsvClosestJet", &deepCsvClosestJet);
            leptonMvaReader[era][flavor]->AddVariable( "sip3d", &sip3d);
            leptonMvaReader[era][flavor]->AddVariable( "dxy", &dxy);
            leptonMvaReader[era][flavor]->AddVariable( "dz", &dz);
            if(flavor == 0){
                if(era == 0){
                    leptonMvaReader[era][flavor]->AddVariable("electronMva", &electronMva);
                } else {
                    leptonMvaReader[era][flavor]->AddVariable("electronMvaFall17NoIso", &electronMvaFall17NoIso);
                }    
            } else {
                leptonMvaReader[era][flavor]->AddVariable("segmentCompatibility", &LepGood_segmentCompatibility);
            }
            if(flavor = 0 && era == 0){
                leptonMvaReader[era][flavor]->BookMVA("BDTG method", "weights/leptonMva/el_tZqTTV16_BDTG.weights.xml");
            } else if (flavor == 0 && era == 1){
                leptonMvaReader[era][flavor]->BookMVA("BDTG method", "weights/leptonMva/el_tZqTTV17_BDTG.weights.xml");
            } else if (flavor == 1 && era == 0){
                leptonMvaReader[era][flavor]->BookMVA("BDTG method", "weights/leptonMva/mu_tZqTTV16_BDTG.weights.xml");
            } else {
                leptonMvaReader[era][flavor]->BookMVA("BDTG method", "weights/leptonMva/mu_tZqTTV17_BDTG.weights.xml");
            }
        }
    }
    /////////////////////////////////////////

    double progress = 0; 	//For printing progress bar
    long nEntries = sampleTree->GetEntries();
    for (long it=0; it <nEntries; ++it) {
        if(it%100 == 0 && it != 0){
            progress += (double) (100./ (double) nEntries);
            tools::printProgress(progress);
        } else if(it == nEntries -1){
            progress = 1.;
            tools::printProgress(progress);
        }
        sampleTree->GetEntry(it);

        std::vector<unsigned> ind;
        //Select both TOP-16-020 and FO leptons. Remove the former later.
        unsigned lCount = std::max(selectLep(ind), selectLep_TOP16_020(ind) );
        if(lCount < 3) continue;

        outputTree->Fill();
    }   
    std::cout << std::endl;
    if(!isData){
        hCounter->Write();
        lheCounter->Write();
        nTrueInt->Write();
    }
    nVertices->Write();
    outputTree->Write("",  BIT(2));
    outputFile->Close(); 
}

int main(int argc, char* argv[]){
    treeReader reader;
    bool isData = false;
    if(argc != 0){
        std::vector<std::string> datasets = {"SingleElectron", "SingleMuon", "DoubleEG", "DoubleMuon", "MuonEG"}; 
        for(auto it = datasets.cbegin(); it != datasets.cend(); ++it){
            std::string name(argv[1]);
            auto pos = name.find(*it);
            if(pos < name.size()) isData = true;
        }
    }
    switch(argc){
        case 2:{
                   reader.skimTree(argv[1], "", isData);
                   return 0;
               }
        case 3:{
                   reader.skimTree(argv[1], argv[2], isData);
                   return 0;
               }
        default:{
                    std::cerr << "Error: Wrong number of options given!" << std::endl;
                    return 1;
                }
    }
}



