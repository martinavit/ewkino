#include "../interface/HistCollection.h"

//include c++ library classes
#include <iostream>
#include <set>

HistCollectionSample::HistCollectionSample(const std::vector<HistInfo>& infoList, std::shared_ptr< Category > categorization, std::shared_ptr<Sample> sam): cat(categorization), sample(sam){
    for(auto infoIt = infoList.cbegin(); infoIt != infoList.cend(); ++infoIt){
        collection.push_back(std::vector< std::shared_ptr<TH1D> >() );
        size_t counter = 0;
        for(auto catIt = cat->getCat().cbegin(); catIt != cat->getCat().cend(); ++catIt){
            collection[counter].push_back(infoIt->makeHist(*catIt + sample->getFileName() ) );
            ++counter;
        }
    }
}

HistCollectionSample::HistCollectionSample(const std::vector<HistInfo>& infoList, const std::vector < std::vector < std::string > >& categorization, std::shared_ptr<Sample> sam):
    HistCollectionSample(infoList, std::make_shared<Category>(categorization), sam) {}

std::shared_ptr<TH1D> HistCollectionSample::access(size_t infoIndex, const std::vector<size_t>& catIndices) const{
    size_t catIndex = cat->getIndex(catIndices);
    return collection[infoIndex][catIndex];
}

void HistCollectionSample::setNegZero(){
    for(auto dIt = collection.begin(); dIt != collection.cend(); ++dIt){
        for(auto cIt = dIt->cbegin(); cIt != dIt->cend(); ++cIt){
            for(unsigned b = 1; b < (*cIt)->GetNbinsX() + 1; ++b){
                if((*cIt)->GetBinContent(b) < 0.) (*cIt)->SetBinContent(b, 0.);
            }
        }
    }
}

HistCollectionSample& HistCollectionSample::operator+=(const HistCollectionSample& rhs){
    if(collection.size() != rhs.collection.size() || cat->getCat().size() != rhs.cat->getCat().size()){
        std::cerr << "HistCollection of incompatible dimensions can not be added: returning left hand side!" << std::endl;
    } else{
        for(size_t dist = 0; dist < collection.size(); ++dist){
            for(size_t c = 0; c < collection[dist].size(); ++c){
                collection[dist][c]->Add(rhs.collection[dist][c].get());
            }
        }
    }
    return *this;
};

HistCollectionSample operator+(const HistCollectionSample& lhs, const HistCollectionSample& rhs){
    HistCollectionSample ret(lhs);
    ret += rhs;
    return ret;
}


HistCollection::HistCollection(const std::vector<HistInfo>& infoList, std::shared_ptr< Category > categorization, const std::vector<Sample>& samList){ 
    for(auto samIt = samList.cbegin(); samIt != samList.cend(); ++samIt){
        fullCollection.push_back(HistCollectionSample(infoList, categorization, std::make_shared<Sample>(Sample(*samIt) ) ) );
    }
}

HistCollection::HistCollection(const std::vector<HistInfo>& infoList, const std::vector < std::vector < std::string > >& categorization, const std::vector<Sample>& samList):
    HistCollection(infoList, std::make_shared<Category>(Category(categorization)), samList) {}

std::shared_ptr<TH1D> HistCollection::access(size_t samIndex, size_t infoIndex, const std::vector<size_t>& catIndices) const{
    return fullCollection[samIndex].access(infoIndex, catIndices);
}

void HistCollection::setNegZero(){
    for(HistCollectionSample& samCol : fullCollection){
        samCol.setNegZero();
    }
}

void HistCollection::mergeProcesses(){
    setNegZero();
    HistCollection tempCol;
    std::set<std::string> usedProcesses;
    for(auto it = fullCollection.cbegin(); it != fullCollection.cend(); ++it){
        if(usedProcesses.find(it->sample->getProc()) != usedProcesses.end()){
            HistCollectionSample tempSam = *it;
            for(auto jt = it + 1; jt != fullCollection.cend(); ++jt){
                if(it->sample->getProc() == jt->sample->getProc()){
                    tempSam += *jt;
                }
            }
            tempCol.fullCollection.push_back(tempSam);
        }
    }
    *this = tempCol;
}
