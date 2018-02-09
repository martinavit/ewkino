//include ROOT classes
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <stdlib.h>

//include other parts of the code
#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"
#include "../interface/ewkinoTools.h"
#include "../interface/HistCollectionSample.h"
#include "../interface/HistCollectionDist.h"
#include "../plotting/plotCode.h"
#include "../plotting/tdrStyle.h"


void treeReader::setup(){
    //Set CMS plotting style
    setTDRStyle();
    gROOT->SetBatch(kTRUE);
    //read samples and cross sections from txt file
    readSamples("sampleLists/samples_dilepCR.txt");
    //info on kinematic distributions to plot
    //name      xlabel    nBins,  min, max
    histInfo = {
        HistInfo("sip3d", "SIP_{3D}", 100, 0, 8),
        HistInfo("dxy", "|d_{xy}| (cm)", 100, 0, 0.05),
        HistInfo("dz", "|d_{z}| (cm)", 100, 0, 0.1),
        HistInfo("miniIso", "miniIso", 100, 0, 0.4),
        HistInfo("leptonMvaSUSY", "SUSY lepton MVA value", 100, -1, 1),
        HistInfo("leptonMvaTTH", "TTH lepton MVA value", 100, -1, 1),
        HistInfo("ptRel", "P_{T}^{rel} (GeV)", 100, 0, 200),
        HistInfo("ptRatio", "P_{T}^{ratio}", 100, 0, 2),
        HistInfo("closestJetCsv", "closest jet CSV", 100, 0, 1),
        HistInfo("chargedTrackMult", "closest jet track multiplicity", 20, 0, 20),
        HistInfo("electronMvaGP", "electron GP MVA value", 100, -1, 1),
        HistInfo("muonSegComp", "muon segment compatibility", 100, 0, 1),
        HistInfo("met", "E_{T}^{miss} (GeV)", 100, 0, 300),
        HistInfo("mll", "M_{ll} (GeV)", 200, 12, 200),
        HistInfo("leadPt", "P_{T}^{leading} (GeV)", 100, 25, 200),
        HistInfo("trailPt", "P_{T}^{trailing} (GeV)", 100, 15, 150),
        HistInfo("nVertex", "number of vertices", 100, 0, 100),
        HistInfo("nJets", "number of jets", 10, 0, 10),
        HistInfo("nBJets_CSVv2", "number of b-jets (CSVv2)", 8, 0, 8),
        HistInfo("nBJets_DeepCSV", "number of b-jets (Deep CSV)", 8, 0, 8)
    };
}

void treeReader::Analyze(const std::string& sampName, const long unsigned begin, const long unsigned end){
    auto samIt = std::find_if(samples.cbegin(), samples.cend(), [sampName](const Sample& s) { return s.getFileName() == sampName; } );
    Analyze(*samIt, begin, end);
}

void treeReader::Analyze(const Sample& samp, const long unsigned begin, const long unsigned end){
    //set up histogram collection for particular sample
    HistCollectionSample histCollection(histInfo, samp, { {"all2017", "RunB", "RunC", "RunD", "RunE", "RunF"}, {"inclusive", "ee", "em", "mm"}, {"nJetsInclusive", "1pt40Jet"}, {"noPuW", "PuW"} });

    //read pu weights for every period
    TFile* puFile = TFile::Open("weights/puWeights2017.root");
    const std::string eras[6] = {"Inclusive", "B", "C", "D", "E", "F"};
    TH1D* puWeights[6];
    for(unsigned e = 0; e < 6; ++e){
        puWeights[e] = (TH1D*) puFile->Get( (const TString&) "puw_Run" + eras[e]);
    }

    const unsigned nDist = histInfo.size();
    const unsigned nRuns = histCollection.categoryRange(0);
    const unsigned nFlav = histCollection.categoryRange(1);
    const unsigned nJetCat = histCollection.categoryRange(2);
    const unsigned nPuRew = histCollection.categoryRange(3);

    initSample(samp, 1);  //use 2017 lumi
    for(long unsigned it = begin; it < end; ++it){
        GetEntry(samp, it);
        //vector containing good lepton indices
        std::vector<unsigned> ind;
        //select leptons
        const unsigned lCount = selectLep(ind);
        if(lCount != 2) continue;
        //require pt cuts (25, 20) to be passed
        if(!passPtCuts(ind)) continue;
        //make lorentzvectors for leptons
        TLorentzVector lepV[lCount];
        for(unsigned l = 0; l < lCount; ++l) lepV[l].SetPtEtaPhiE(_lPt[ind[l]], _lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
        //Cut of Mll at 12 GeV
        if((lepV[0] + lepV[1]).M() < 12) continue;
        //reject SS events
        if(_lCharge[ind[0]] == _lCharge[ind[1]]) continue;
        //determine flavor compositions
        unsigned flav = dilFlavorComb(ind) + 1; //reserve 0 for inclusive
        //order jets and check the number of jets
        std::vector<unsigned> jetInd;
        unsigned jetCount = nJets(jetInd);

        //Extra category selection: for DY select onZ, for ttbar select 1 b-jet 2-jets
        if(flav == 1 || flav == 3){ //OSSF
            if(fabs((lepV[0] + lepV[1]).M() - 91) > 10) continue;
            if(nBJets(0,  true, false, 0) != 0) continue;
        } else if(flav == 2){
            if(_met < 50) continue;
            if(jetCount < 2 || nBJets() < 1) continue;
        }
        //determine run perios
        unsigned run;
        run = ewk::runPeriod(_runNb) + 1 - 1; //reserve 0 for inclusive // -1 because we start at run B 
        //max pu to extract
        float max = (run == 0 || run > 3) ? 80 : 80;
        //float max = 100;
        //loop over leading leptons
        for(unsigned l = 0; l < 2; ++l){
            double fill[12] = {_3dIPSig[ind[l]], fabs(_dxy[ind[l]]), fabs(_dz[ind[l]]), _miniIso[ind[l]], _leptonMvaSUSY[ind[l]], _leptonMvaTTH[ind[l]], _ptRel[ind[l]], _ptRatio[ind[l]], _closestJetCsvV2[ind[l]], (double) _selectedTrackMult[ind[l]], (_lFlavor[ind[l]] == 0) ? _lElectronMva[ind[l]] : 0,  (_lFlavor[ind[l]] == 1) ? _lMuonSegComp[ind[l]] : 0};
            //fill histograms
            for(unsigned j = 0; j < nJetCat; ++j){
                if(j == 1 && ( (jetCount == 0) ? false :_jetPt[jetInd[0]] <= 40 ) ) continue;
                for(unsigned pu = 0; pu < nPuRew; ++pu){
                    for(unsigned dist = 0; dist < 12; ++dist){
                        if(_lFlavor[ind[l]] == 0 && dist == 11) continue;  //do not plot muonSegComp for electrons
                        if(_lFlavor[ind[l]] == 1 && dist == 10) continue;  //do not plot electronMva for muons
                        for(unsigned r = 0; r < nRuns; ++r){
                            if(!samp.isData() || r == run || r == 0){
                                double puw = 1.;
                                if(!samp.isData() && pu == 1){
                                    puw = puWeights[run]->GetBinContent(puWeights[run]->FindBin( std::min(_nTrueInt, max) ) );
                                }
                                histCollection.access(dist, {r, flav, j, pu})->Fill(std::min(fill[dist], histInfo[dist].maxBinCenter()), weight*puw); 
                                histCollection.access(dist, {r, 0,    j, pu})->Fill(std::min(fill[dist], histInfo[dist].maxBinCenter()), weight*puw);

                            }
                        }
                    }
                }
            }
        }
        double fill[nDist - 12] = {_met, (lepV[0] + lepV[1]).M(), _lPt[ind[0]], _lPt[ind[1]], (double) _nVertex, (double) jetCount, (double) nBJets(0, false), (double) nBJets()}; //replace 0 by _met for correct trees
        for(unsigned j = 0; j < nJetCat; ++j){
            if(j == 1 && ( (jetCount == 0) ? false :_jetPt[jetInd[0]] <= 40 ) ) continue;
            for(unsigned pu = 0; pu < nPuRew; ++pu){
                for(unsigned dist = 12; dist < nDist; ++dist){
                    for(unsigned r = 0; r < nRuns; ++r){
                        if(!samp.isData() || r == run || r == 0){
                            double puw = 1.;
                            if(!samp.isData() && pu == 1){
                                puw = puWeights[run]->GetBinContent(puWeights[run]->FindBin( std::min(_nTrueInt, max) ) );
                            }
                            histCollection.access(dist, {r, flav, j, pu})->Fill(std::min(fill[dist - 12], histInfo[dist].maxBinCenter()), weight*puw);
                            histCollection.access(dist, {r, 0,    j, pu})->Fill(std::min(fill[dist - 12], histInfo[dist].maxBinCenter()), weight*puw);
                        }
                    }
                }
            } 
        }
    }
    //store histograms for later merging
    histCollection.store("tempHists/", begin, end);
}

void treeReader::splitJobs(){
    for(unsigned sam = 0; sam < samples.size(); ++sam){
        initSample(1);
        for(long unsigned it = 0; it < nEntries; it+=1000000){
            long unsigned begin = it;
            long unsigned end = std::min(nEntries, it + 1000000);
            //make temporary job script 
            std::ofstream script("runTuples.sh");
            tools::initScript(script);
            script << "./dilepCR " << samples[currentSample].getFileName() << " " << std::to_string(begin) << " " << std::to_string(end);
            script.close();
            //submit job
            tools::submitScript("runTuples.sh", "01:00:00");
         }
    }
}

void treeReader::plot(const std::string& distName){
    //loop over all distributions and find the one to plot
    for(size_t d = 0; d < histInfo.size(); ++d){
        if(histInfo[d].name() == distName){
            //read collection for this distribution from files
            HistCollectionDist col("inputList.txt", histInfo[d], samples, { {"all2017", "RunB", "RunC", "RunD", "RunE", "RunF"}, {"inclusive", "ee", "em", "mm"}, {"nJetsInclusive", "1pt40Jet"}, {"noPuW", "PuW"} });
            for(unsigned c = 0; c < col.categorySize(); ++c){
                std::string header;
                std::string categoryName = col.categoryName(c);

                if(categoryName.find("_mm_") != std::string::npos) header += "#mu#mu : ";
                else if(categoryName.find("_em_") != std::string::npos) header += "e#mu : ";
                else if(categoryName.find("_ee_") != std::string::npos) header += "ee : ";
                
                if(categoryName.find("RunB") != std::string::npos) header += "2017 Run B";
                else if(categoryName.find("RunC") != std::string::npos) header += "2017 Run C";
                else if(categoryName.find("RunD") != std::string::npos) header += "2017 Run D";
                else if(categoryName.find("RunE") != std::string::npos) header += "2017 Run E";
                else if(categoryName.find("RunF") != std::string::npos) header += "2017 Run F";
                else( header += "42 fb^{-1}");
                header += " (13 TeV)";
                if(d == 10 && categoryName.find("_mm_") != std::string::npos) continue; //do not plot electronMva for mm category
                if(d == 11 && categoryName.find("_ee_") != std::string::npos) continue; //do not plot muonSegComp for ee category
                col.getPlot(c).draw("plots/ewkino/dilepCR", "ewkinoDilep", true, true, header);
                //col.getPlot(c).draw("ewkinoDilep", true, true, "");
            }
        }
    }
}

void treeReader::readPlots(){
    //histCollection.read("tempHists");
}

void treeReader::splitPlots(){
    //histCollection = HistCollection("tempHists", histInfo, samples, { {"all2017", "RunB", "RunC", "RunD", "RunE", "RunF"}, {"inclusive", "ee", "em", "mm"}, {"nJetsInclusive", "1pt40Jet"}, {"noPuW", "PuW"} });
    std::system("touch inputList.txt");
    std::system("for f in tempHists/*; do echo $f >> inputList.txt; done");
    for(auto& h: histInfo){
        std::ofstream script("printPlots.sh");
        tools::initScript(script);
        script << "./dilepCR plot " << h.name();
        script.close();
        tools::submitScript("printPlots.sh", "00:30:00");
    }
}

int main(int argc, char* argv[]){
    treeReader reader;
    reader.setup();
    std::vector<std::string> argvStr;
    for(unsigned i = 0; i < argc; ++i){
        argvStr.push_back(std::string(argv[i]));
    }
    //submit single job if sample and range given
    if(argc == 4){
        long unsigned begin = std::stoul(argvStr[2]);
        long unsigned end = std::stoul(argvStr[3]);
        //sample, first entry, last entry:
        reader.Analyze(argvStr[1], begin, end);
    }
    else if(argc > 1 && argvStr[1] == "plot"){
        reader.readPlots();
        if(argc > 2){
            std::string dist(argv[2]);
            reader.plot(dist);
        } else{
            reader.splitPlots();           
        }
    }
    else if(argc > 1 && argvStr[1] == "run"){
        //Analyze all, or split jobs
        reader.splitJobs();
    } else{
        std::cerr << "no valid input given" << std::endl;
    }
    return 0;
}
