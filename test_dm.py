from ROOT import TH1F
from DataFormats.FWLite import Events, Handle
# events = Events('root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DarkMatter_MonoZToLL_V_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CE49A2DE-0DB7-E611-9B9F-141877344D39.root')

# events = Events('root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DarkMatter_MonoZToLL_V_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/AACE521C-32B7-E611-BDE3-0025905B860E.root ')

# events = Events('root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DarkMatter_MonoZToLL_V_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0889AC17-27B9-E611-AD12-0CC47A4D7644.root')

events = Events('root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DarkMatter_MonoZToLL_V_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0013148D-9DAC-E611-9D9D-20CF3027A5CD.root')

handleGen = Handle('std::vector<reco::GenParticle>')
labelGen = 'prunedGenParticles'

z_pt = TH1F('z_pt', '', 100, 0., 1000.)

for i_ev, ev in enumerate(events):
    if i_ev % 10 == 1:
        print i_ev
    ev.getByLabel(labelGen, handleGen)
    gens = handleGen.product()
    for gen in gens:
        if gen.pdgId() == 23:
            if gen.daughter(0).pdgId() != 23:
                z_pt.Fill(gen.pt())

z_pt.Draw()


