Checkout Instructions
phase2-tutorial/CMSSW_8_1_0_pre12/src/RecoTauTag/phase2Taus
mkdir phase2-tutorial
cmsrel CMSSW_8_1_0_pre12
cd CMSSW_8_1_0_pre12/src
git cms-init
git cms-addpkg RecoTauTag/RecoTau
git clone 
scram b -j 5

#.... wait for the build
