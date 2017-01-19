```
SCRAM_ARCH=slc6_amd64_gcc600
mkdir phase2-tutorial-taus
cd phase2-tutorial-taus
cmsrel CMSSW_8_1_0_pre12
cd CMSSW_8_1_0_pre12/src
cmsenv
git cms-init
git cms-merge-topic 16450 #this should work, but please pay attention in case there is an issue here
git cms-addpkg RecoTauTag/RecoTau
cd RecoTauTag
git clone --recursive -b exerciseNov16 https://github.com/isobelojalvo/phase2Taus.git
cd ../
scram b -j 5
cd RecoTauTag/phase2Taus/test
cmsRun make_ntuple_from_miniaod.py


#OR

cmsRun runTauEfficiency.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root
root -l -b -q plotEfficiency.C


#.... wait for the build
```
