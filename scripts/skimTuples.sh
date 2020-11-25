#!/bin/bash

#include bash functiosn to set up CMSSW
source setCMSSW.sh
echo "at the beginning"
echo "before cwd"
echo "$(pwd)"
cwd=$(pwd)                                          #current working directory needed to locate code 
echo "cwd"
echo "${cwd}"
#skimSample(){                                           #function to skim one sample
    #name="${1%/*}"                                      #remove everything before the last "/" in the path to the sample

    #if [[ $1 = *"_realistic_v10"* ]]; then
        #name="${name}_realistic_v10"
    #elif [[ $1 = *"_realistic_v14"* ]]; then
        #name="${name}_realistic_v14"
    #elif [[ $1 = *"_realistic_v11"* ]]; then
        name="${name}_realistic_v11"
    #fi

    #if [[ $1 = *"_RECOPF"* ]]; then
        #name="${name}_RECOPF"
    #fi

    #if [[ $1 = *"Fall17"* ]] || [[ $1 = *"Run2017"* ]] || [[ $1 = *"2017"* ]]; then
        #name="${name}_Fall17"
        #name="${name}_2017"
    #if [[ $1 = *"18Mini"* ]]; then
        #name="${name}_Fall17"
        #name="${name}_2018"    
    #else 
        #name="${name}_2016"
    #fi
    #echo "name -->"
    #echo "$name"
    #echo "this is it the cwd" 
    #echo "${cwd}" 

    #outputDir=~/Work/ntuples_temp_${name}
    #echo "~/Work/ntuples_temp_${name}"
    #if [ ! -d "$outputDir" ]; then                      #make output directory if it doesn't exist 
        mkdir -p $outputDir
    #fi
    #echo "before submit = skimjob.sh"
    #submit=~/skimJob.sh
    #makeSubmit $submit $2                               #make temporary submission script

    #count=0                                             #file counter
    #subdir=$(ls $1 | sort -r | tail -1)
    #files=$1/$subdir/*/*.root
    #files=${1}/*/*/*root
    #or f in $files
        #do if (( $count % 50 == 0)); then
            #submitJob $submit "12:00:00"
            #makeSubmit $submit $2
        #fi
        #filename=${f##*/}                               
        #filename=${f///}
        #filename=${filename%.*}
        #echo "${cwd}/../skimTree $f $outputDir/ > ${outputDir}/${filename}_log.txt 2> ${outputDir}/${filename}_err.txt" >> $submit
        #count=$((count+1))
    #done
    #echo "before submitJob $submit"
    #submitJob $submit "12:00:00"
    #rm $submit                                          #remove temporary submit file
#}

