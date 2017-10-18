#!/bin/bash

#include bash functiosn to set up CMSSW
source setCMSSW.sh

skimSample(){                                           #function to skim one sample
    name="${1%/*}"                                      #remove everything before the last "/" in the path to the sample
    echo "$name"
    outputDir=~/Work/ntuples_temp_${name}
    if [ ! -d "$outputDir" ]; then                      #make output directory if it doesn't exist 
        mkdir $outputDir
    fi
    submit=~/skimJob.sh
    makeSubmit $submit $2                               #make temporary submission script
    
    count=0                                             #file counter
    files=${1}/*/*/*root
    for f in $files
        do if (( $count % 50 == 0))
            then qsub $submit -l walltime=04:00:00;
            makeSubmit $submit $2
        fi
        #filename=${f##*/}                               
        filename=${f///}
        filename=${filename%.*}
        echo "~/Work/AnalysisCode/ewkino/skimTree $f $outputDir/ > ${outputDir}/${filename}_log.txt 2> ${outputDir}/${filename}_err.txt" >> $submit
        count=$((count+1))
    done
    qsub $submit -l walltime=04:00:00;
    rm $submit                                          #remove temporary submit file
}

baseFolder=/pnfs/iihe/cms/store/user/wverbeke/heavyNeutrino
cd $baseFolder
foldersMC=*/*ewkinoMCList                               #add suffix for newer versions
foldersData=*/*ewkinoDataList
for d in $foldersMC $foldersData                        #skim all samples
    do skimSample $d $baseFolder
done
