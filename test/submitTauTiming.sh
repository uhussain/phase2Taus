#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat runTauTiming.py > SUBtiming.py
cat submit.py >> SUBtiming.py


for dir in TTBar-RelVal-0  TTBar-RelVal-140  TTBar-RelVal-200 ZTT-RelVal-0  ZTT-RelVal-140   ZTT-RelVal-200  QCDMuEnriched_realistic
do 
    echo " "
    echo "====================" $dir "========================="

    rm -r /data/ojalvo/$1-$dir-SUBtiming/
    mkdir /data/ojalvo/$1-$dir-SUBtiming/
    
#make dag dir
    mkdir -p /data/ojalvo/$1-$dir-SUBtiming/dags
    mkdir -p /data/ojalvo/$1-$dir-SUBtiming/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/$1-$dir-SUBtiming/
    
    farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=$dir.txt \
	--submit-dir=/data/ojalvo/$1-$dir-SUBtiming/submit \
	--output-dag-file=/data/ojalvo/$1-$dir-SUBtiming/dags/dag \
	$1-$dir  \
	$CMSSW_BASE  \
	$CMSSW_BASE/src/RecoTauTag/phase2Taus/test/SUBtiming.py 

done

rm SUBtiming.py