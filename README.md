```
#!/bin/bash
scramv1 project CMSSW CMSSW_8_1_0_pre16
cd CMSSW_8_1_0_pre16/src
eval `scramv1 runtime -sh`
git cms-init
git  cms-merge-topic 16450
git cms-addpkg RecoTauTag/RecoTau
cd RecoTauTag
git clone git@github.com:isobelojalvo/phase2Taus.git
cd ../
scram b -j 5
cd RecoTauTag/phase2Taus/test

#.... wait for the build
```
