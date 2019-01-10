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


events = Events('no_eg_filters.root')
# events = Events('ori.root')
handle = Handle('std::vector<pat::Tau>')
label = 'slimmedTaus'

handleGen = Handle('std::vector<reco::GenParticle>')
labelGen = 'prunedGenParticles'

handlePF = Handle('std::vector<pat::PackedCandidate>')
labelPF = 'packedPFCandidates'

countGenMatched1 = 0
countGenUnMatched1 = 0

# These lists mean "missing in ..."

rereco_but_not_ori = [(6342L, 1), (6393L, 0), (7623L, 0), (7773L, 1), (7842L, 1)]

ori_but_not_rereco = [(6331L, 1), (6348L, 0), (6398L, 0), (6537L, 1), (6544L, 0), (7040L, 1), (7049L, 0), (7110L, 1), (7145L, 0), (7164L, 0), (7168L, 0), (7194L, 0), (7249L, 0), (7325L, 0), (7340L, 0), (7334L, 0), (7357L, 0), (7364L, 1), (7365L, 0), (7383L, 1), (7425L, 0), (7424L, 1), (7623L, 1), (7752L, 1), (7795L, 0), (7804L, 0), (7802L, 1), (7863L, 1), (7894L, 1), (7950L, 1), (7938L, 0), (7952L, 1), (7986L, 0)]

ev_list = [] # Pairs of ev, i(gen tau) that aren't matched

for ev in events:
    ev.getByLabel(label, handle)
    taus = handle.product()
    taus = [tau for tau in taus]# if tau.pt() > 18.]
    
    ev.getByLabel(labelGen, handleGen)
    gens = handleGen.product()

    ev.getByLabel(labelPF, handlePF)
    pfcands = handlePF.product()

    iev = ev.eventAuxiliary().id().event() # lumi()  run()

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

                    
    
    for i_gen_tau, gen_tau in enumerate(vis_had_taus):
        if gen_tau.pt() > 25. and abs(gen_tau.eta()) < 2.3:
            matched_tau1 = None
            best_tau1 = None
            for tau1 in taus:
                if deltaR(tau1, gen_tau) < 0.3:
                    best_tau1 = tau1
                    if tau1.pt() > 0.8 * gen_tau.pt() and tau.tauID('decayModeFinding'): # and tau.tauID('againstElectronVLooseMVA6'):
                        matched_tau1 = tau1
            if matched_tau1:
                countGenMatched1 += 1
            else:
                countGenUnMatched1 += 1
                ev_list.append((iev, i_gen_tau))
                
            if (iev, i_gen_tau) in ori_but_not_rereco + rereco_but_not_ori:

                print '####   Event, i_tau ', iev, i_gen_tau
                print '\nDid not match to vis tau', '{:.2f} {:.2f} {:.2f} {}'.format(gen_tau.pt(), gen_tau.eta(), gen_tau.phi(), gen_tau.dm)
                for d in gen_tau.charged:
                    print ' charged {:.2f} {:.2f} {:.2f} {}'.format(d.pt(), d.eta(), d.phi(), d.pdgId())
                for pi0 in gen_tau.pizeros:
                    print ' pizero'
                    for d in pi0:
                        print '  gamma {:.2f} {:.2f} {:.2f} {}'.format(d.pt(), d.eta(), d.phi(), d.pdgId())
                if best_tau1:
                    print ' ## best tau {:.2f} {:.2f} {:.2f} {:.2f} {}'.format(best_tau1.pt(), best_tau1.eta(), best_tau1.phi(), best_tau1.mass(), best_tau1.decayMode())
                # Can investigate why here
                print '  ## PF candidates in cone'
                for pf in pfcands:
                    if deltaR(gen_tau, pf) < 0.3:
                        print '    {:.2f} {:.2f} {:.2f} {}'.format(pf.pt(), pf.eta(), pf.phi(), pf.pdgId())
                    # print '  ## Combo taus in cone'
                    # for c_tau in c_taus:
                    #     if deltaR(gen_tau, c_tau) < 0.3 and c_tau.pt() > 20.:
                    #         print '    {:.2f} {:.2f} {:.2f} {:.2f} {}'.format(c_tau.pt(), c_tau.eta(), c_tau.phi(), c_tau.mass(), c_tau.decayMode())
                    # import pdb; pdb.set_trace()


print 'Matched (re-reco to gen)', countGenMatched1
print 'Unmatched (re-reco to gen)', countGenUnMatched1

print ev_list
