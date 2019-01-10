import math
import ROOT
import numpy
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

class GenJet(object):
    def __init__(self, p4):
        self.p4 = p4
        self.dm = -99

    def __getattr__(self,name):
        return getattr(self.p4, name)

class Var(object):
    def __init__(self, name, type, len=1, len_var=None, default=-999):
        self.name = name
        self.type = type
        self.len = len
        self.len_var = len_var
        self.storage = None
        self.default = default

    def reset(self):
        for i in xrange(self.len):
            self.storage[i] = self.default

    def fill(self, val, i_val=0):
        if i_val is not None:
            self.storage[i_val] = val
        else:
            for i in self.len:
                self.storage[i] = val[i]

    def __str__(self):
        return 'Var: name={}, type={}, val={:.2f}'.format(self.name, self.type, self.storage[0])

def daughters(p):
    ds = []
    for i_d in xrange(p.numberOfDaughters()):
        ds.append(p.daughter(i_d))
    return ds

def getNChargedAndNeutral(ds, n_charged=0, n_pizero=0):
    for d in ds:
        if d.charge():
            n_charged += 1
        elif d.pdgId() == 111:
            n_pizero += 1
        else:
            drs = daughters(d)
            if drs:
                n_charged, n_pizero = getNChargedAndNeutral(drs, n_charged, n_pizero)
            else:
                if d.pdgId() not in [310, 130, 22]:
                    import pdb; pdb.set_trace()

    return n_charged, n_pizero

def findGenTaus(gens):
    vis_had_taus = []
    for p in gens:
        if abs(p.pdgId()) == 15:
            # Check if not final tau or leptonic tau decay
            ds = daughters(p)
            if any(abs(d.pdgId()) in [11, 13, 15] for d in ds):
                continue
            ds = [d for d in ds if abs(d.pdgId()) not in [12, 14, 16]] # keep visible products
            if ds:
                vis_had_taus.append(GenTau(ds[-1].p4()))
                vis_tau = vis_had_taus[-1]
                for i_d in xrange(len(ds) - 1):
                    vis_tau.p4 += ds[i_d].p4()
                
                n_charged, n_pizero = getNChargedAndNeutral(ds)

                if n_charged == 1:
                    vis_tau.dm = n_pizero
                elif n_charged == 3:
                    vis_tau.dm = 10 + n_pizero
                else:
                    print 'Weird gen DM:', n_charged, n_pizero
                    vis_tau.dm = 100
                    # import pdb; pdb.set_trace()

                vis_tau.charged = [d for d in ds if d.charge()]
                vis_tau.pizeros = [daughters(d) for d in ds if d.pdgId()==111]
    return vis_had_taus

if __name__ == '__main__':

    signal = True

    f_name = 'root://eoscms.cern.ch//eos/cms/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/50000/D40255DF-7EE1-E711-A0E9-4C79BA181477.root' if signal else 'root://eoscms.cern.ch//eos/cms/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/50000/D40255DF-7EE1-E711-A0E9-4C79BA181477.root'

    events = Events(f_name)
    handle = Handle('std::vector<pat::Tau>')
    label = ('slimmedTaus', '', 'PAT')

    handlePU =  Handle('std::vector<PileupSummaryInfo>')
    labelPU = 'slimmedAddPileupInfo'

    handleGen = Handle('std::vector<reco::GenParticle>')
    labelGen = 'prunedGenParticles'

    handleGenJets = Handle('vector<reco::GenJet>')
    labelGenJets = 'slimmedGenJets'

    n_gen_objs = 0
    n_matched_taus = 0

    out_file = ROOT.TFile('per_tau.root', 'recreate')
    tau_tree = ROOT.TTree('per_tau', 'per_tau')

    all_vars = [
        Var('tau_dm', int),
        Var('tau_pt', float),
        Var('tau_eta', float),
        Var('tau_phi', float),
        Var('tau_mass', float),
        Var('tau_z', float),
        Var('tau_gendm', int),
        Var('tau_genpt', float),
        Var('tau_geneta', float),
        Var('tau_genphi', float),
        Var('n_true_interactions', int),
        Var('n_iso_ch', int, default=0),
        Var('tau_ch_dz', float, 99, 'n_iso_ch'),
        Var('tau_ch_pt', float, 99, 'n_iso_ch')
    ]

    all_var_dict = {var.name: var for var in all_vars}
    avd = all_var_dict

    for var in all_vars:
        var.storage = numpy.zeros(var.len, dtype=var.type)
        addendum = '[{}]'.format(var.len_var) if var.len > 1 else ''
        tau_tree.Branch(var.name, var.storage, var.name + addendum +
                        '/'+('I' if var.type == int else 'D'))

    n_ev = 0
    for n_ev, ev in enumerate(events):
        if n_ev > 1000: 
            break
        ev.getByLabel(label, handle)
        taus = handle.product()

        ev.getByLabel(labelGen, handleGen)
        gens = handleGen.product()

        ev.getByLabel(labelGenJets, handleGenJets)
        genJets = handleGenJets.product()

        ev.getByLabel(labelPU, handlePU)
        for puInfo in handlePU.product():
            if puInfo.getBunchCrossing() == 0:
                nTrue = puInfo.getTrueNumInteractions()

        # Find visible generated taus
        vis_had_taus = findGenTaus(gens)

        if not signal:
            vis_had_taus = [GenJet(genJet) for genJet in genJets if not any(deltaR(genJet, vis_had_tau) < 0.5 for vis_had_tau in vis_had_taus)]

        # Find reconstructed taus matched to gen taus
        for gen_obj in vis_had_taus:
            if gen_obj.pt() > 25. and abs(gen_obj.eta()) < 2.3:
                n_gen_objs += 1
                avd['tau_gendm'].fill(gen_obj.dm)
                avd['tau_genpt'].fill(gen_obj.pt())
                avd['tau_geneta'].fill(gen_obj.eta())
                avd['tau_genphi'].fill(gen_obj.phi())

                avd['n_true_interactions'].fill(nTrue)

                matched_tau = None
                for tau in taus:
                    if deltaR(tau, gen_obj) < 0.2:
                        matched_tau = tau
                        break
                if matched_tau:
                    n_matched_taus += 1
                    # print matched_tau.tauID('chargedIsoPtSum')
                    avd['n_iso_ch'].fill(matched_tau.isolationChargedHadrCands().size())
                    for n_iso, iso_ch in enumerate(matched_tau.isolationChargedHadrCands()):
                        avd['tau_ch_dz'].fill(iso_ch.dz(matched_tau.vertex()), n_iso)
                        avd['tau_ch_pt'].fill(iso_ch.pt(), n_iso)
                    avd['tau_dm'].fill(matched_tau.decayMode())
                    avd['tau_pt'].fill(matched_tau.pt())
                    avd['tau_eta'].fill(matched_tau.eta())
                    avd['tau_phi'].fill(matched_tau.phi())
                    avd['tau_z'].fill(matched_tau.vertex().z())
            
                tau_tree.Fill()
                for var in avd.values():
                    var.reset()
        

    print 'Found generated, matched:', n_gen_objs, n_matched_taus
    out_file.Write()
    out_file.Close()
