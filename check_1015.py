import math
from DataFormats.FWLite import Events, Handle

def deltaPhi(p1, p2):
    dphi = p1 - p2
    while dphi > math.pi:
        dphi -= 2*math.pi
    while dphi < -math.pi:
        dphi += 2*math.pi
    return dphi

def deltaR(t1, t2):
    deta = abs(t1.eta() - t2.eta())
    dphi = deltaPhi(t1.phi(), t2.phi())
    return math.sqrt(deta**2 + dphi**2)

class GenTau(object):
    def __init__(self, p4):
        self.p4 = p4

    def __getattr__(self,name):
        return getattr(self.p4, name)


events = Events('root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_1_7/RelValTenTau_15_500/MINIAODSIM/PU25ns_101X_upgrade2018_realistic_HEmiss_v1-v1/10000/C070CA2D-2B80-E811-A723-0CC47A7C3610.root')
# events = Events('root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_10_1_7/RelValZpTT_1500_13/MINIAODSIM/PU25ns_101X_upgrade2018_realistic_HEmiss_v1-v1/10000/E66A71A9-E87F-E811-813C-0CC47A4C8ECA.root')

# events = Events('root://cms-xrd-global.cern.ch///store/mc/RunIISummer18MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/FlatPU0to70_101X_upgrade2018_realistic_v7-v2/70000/FA946AC7-E370-E811-A369-0CC47A4C8F08.root')

h_vertex = Handle('vector<reco::Vertex> ')
l_vertex = 'offlineSlimmedPrimaryVertices'

h_tau = Handle('vector<pat::Tau>')
l_tau = 'slimmedTaus'

h_bs = Handle('reco::BeamSpot')
l_bs = 'offlineBeamSpot'

h_gen = Handle('vector<reco::GenParticle>')
l_gen = 'prunedGenParticles'

for i_ev, ev in enumerate(events):
    print '\nAnalysing event #', i_ev
    ev.getByLabel(l_vertex, h_vertex)
    ev.getByLabel(l_tau, h_tau)
    ev.getByLabel(l_bs, h_bs)
    ev.getByLabel(l_gen, h_gen)

    gens = h_gen.product()
    print 'GEN'
    for gen in gens:
        if gen.pdgId() != 2212:
            print '{:.2f} {:.2f}'.format(gen.vx(), gen.vy())


    vertices = h_vertex.product()
    print 'RECO VERTICES'
    for v in vertices:
        print '{:.2f} {:.2f} {:.2f} {}'.format(v.x(), v.y(), v.ndof(), v.nTracks())
    
    bs = h_bs.product()
    print 'Offline BS x, y {:.2f} {:.2f}'.format(bs.x(0.), bs.y(0.))

    taus = h_tau.product()
    for tau in taus:
        print 'tau dxy, dz', tau.leadChargedHadrCand().get().dxy(), tau.leadChargedHadrCand().get().dz() 

