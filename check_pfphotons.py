from DataFormats.FWLite import Events, Handle
from CMGTools.H2TauTau.proto.analyzers.HTTGenAnalyzer import HTTGenAnalyzer
p4sum = HTTGenAnalyzer.p4sum

if __name__ == '__main__':

    signal = True

    f_name = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/AToZhToLLTauTau_M-220_13TeV_madgraph_4f_LO/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/60000/8A590B2F-DDD0-E611-83A6-002590D9D8C2.root'

    events = Events(f_name)

    handleGen = Handle('std::vector<reco::GenParticle>')
    labelGen = 'prunedGenParticles'

    n_ev = 0
    for n_ev, ev in enumerate(events):
        if n_ev > 1000: 
            break
        
        ev.getByLabel(labelGen, handleGen)
        gens = handleGen.product()

        abosons = [p for p in gens if p.pdgId() == 36]
        for aboson in abosons:
            print 'm of A in gen list', aboson.mass()

        leptons_prompt = [p for p in gens if abs(p.pdgId()) in [11, 12, 13, 14] and p.fromHardProcessFinalState()]
        taus_prompt = [p for p in gens if p.statusFlags().isDirectHardProcessTauDecayProduct()]
        m_a = p4sum(leptons_prompt + taus_prompt)
        print 'm from leptons and taus', m_a.mass()

        if abosons and abs(abosons[0].mass() - m_a.mass()) > 1:
            import pdb; pdb.set_trace()

