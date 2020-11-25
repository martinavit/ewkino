#!/bin/bash

#include bash functiosn to set up CMSSW
source setCMSSW.sh

cwd=$(pwd)                                          #current working directory needed to locate code 
echo $pwd
skimSample(){                                           #function to skim one sample
    name="${1%/*}"                                      #remove everything before the last "/" in the path to the sample

    if [[ $1 = *"_realistic_v10"* ]]; then
        name="${name}_realistic_v10"
    elif [[ $1 = *"_realistic_v14"* ]]; then
        name="${name}_realistic_v14"
    elif [[ $1 = *"_realistic_v11"* ]]; then
        name="${name}_realistic_v11"
    fi

    if [[ $1 = *"_RECOPF"* ]]; then
        name="${name}_RECOPF"
    fi

    if [[ $1 = *"Fall17"* ]] || [[ $1 = *"Run2017"* ]] || [[ $1 = *"2017"* ]]; then
        #name="${name}_Fall17"
        name="${name}_2017"
    #if [[ $1 = *"18Mini"* ]]; then
        #name="${name}_Fall17"
        #name="${name}_2018"    
    else 
        name="${name}_2016"
    fi
    echo "$name"
    echo "this is it the pwd $pwd"

    outputDir=~/Work/ntuples_temp_${name}
    echo "~/Work/ntuples_temp_${name}"
    if [ ! -d "$outputDir" ]; then                      #make output directory if it doesn't exist 
        mkdir -p $outputDir
    fi
    submit=~/skimJob.sh
    makeSubmit $submit $2                               #make temporary submission script

    count=0                                             #file counter
    subdir=$(ls $1 | sort -r | tail -1)
    files=$1/$subdir/*/*.root
    #files=${1}/*/*/*root
    for f in $files
        do if (( $count % 50 == 0)); then
            submitJob $submit "12:00:00"
            makeSubmit $submit $2
        fi
        #filename=${f##*/}                               
        filename=${f///}
        filename=${filename%.*}
        echo "${cwd}/../skimTree $f $outputDir/ > ${outputDir}/${filename}_log.txt 2> ${outputDir}/${filename}_err.txt" >> $submit
        echo "${cwd}/../skimTree $f $outputDir/ > ${outputDir}/${filename}_log.txt 2> ${outputDir}/${filename}_err.txt"
        count=$((count+1))
    done
    submitJob $submit "12:00:00"
    $submit                                          #remove temporary submit file
}

#baseFolder=/pnfs/iihe/cms/store/user/mvit/heavyNeutrino/2018_fromTom
baseFolder=/pnfs/iihe/cms/store/user/tomc/heavyNeutrino

cd $baseFolder
folderTestmu=HeavyNeutrino_trilepton_M-*mu_Dirac_massiveAndCKM_LO*/*17*displaced_signals_v4
folderTeste=HeavyNeutrino_trilepton_M-*e_Dirac_massiveAndCKM_LO*/*17*displaced_signals_v4
#folderTest=SingleMuon/*displaced_2018_v2
folderTomMC=*/*displaced_2018_v1
foldersData=*/*2016_legacy9March
foldersMC=2016_94Mc9March2/*/*2016_94Mc9March2
foldersData17=*/*2017_rereco9March
foldersMC17=*/*2016_fromTom
foldersLeptonMva16=*CUETP8M1*/*leptonMvaTrainingList-v5
foldersLeptonMva17=TTTo*CP5*/*leptonMvaTrainingList-v5
foldersFR_dataEG_2017=*/*2017_FR_EG22MarchFR
foldersFR_SingleM_2017=*/*2017_FRSingleM22MarchFR
foldersFR_DoubleM_2017=*/*2017_FRDoubleM22MarchFR
foldersFR_mc_2017=*/*2017_FRmc22MarchFR
#for d in $foldersMC $foldersMC17 $foldersData $foldersData17                        #skim all samples 
for d in $folderTestmu $folderTeste
    do skimSample $d $baseFolder
done
