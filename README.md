```
#!/bin/bash
mkdir phase2-tutorial
cmsrel CMSSW_8_1_0_pre12
cd CMSSW_8_1_0_pre12/src
git cms-init
git cms-addpkg RecoTauTag/RecoTau
cd RecoTauTag
git clone git@github.com:isobelojalvo/phase2Taus.git
scram b -j 5
cd RecoTauTag/phase2Taus/test
cmsRun make_ntuple_from_miniaod.py


#OR

cmsRun runTauEfficiency.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root
root -l -b -q MiniAOD_effi_80x_DYtoLL.root


#.... wait for the build
```
