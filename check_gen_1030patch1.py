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


events = Events('root://cms-xrd-global.cern.ch////store/relval/CMSSW_10_3_0/RelValZTT_13/MINIAODSIM/PUpmx25ns_103X_upgrade2018_realistic_v7-v1/10000/1D5B2A91-8EE2-E84C-9A2C-70C59F816C51.root')

# files = [
# 'root://cms-xrd-global.cern.ch////store/relval/CMSSW_10_3_0/RelValTenTau_15_500/MINIAODSIM/PU25ns_103X_upgrade2018_realistic_v7-v1/10000/14B5CC4F-1B90-174E-9989-415E30E47EF5.root',
# 'root://cms-xrd-global.cern.ch////store/relval/CMSSW_10_3_0/RelValTenTau_15_500/MINIAODSIM/PU25ns_103X_upgrade2018_realistic_v7-v1/10000/7EED20DF-C4C4-F94E-BABF-52870CCD9972.root'
# ]
# events = Events(files)
# handle = Handle('std::vector<pat::Tau>')
# label = ('patTaus', '', 'TAURECO')
handle = Handle('std::vector<pat::Tau>')
label = ('slimmedTaus', '', 'RECO')

handlePF = Handle('std::vector<pat::PackedCandidate>')
labelPF = 'packedPFCandidates'

handleLost = Handle('std::vector<pat::PackedCandidate>')
labelLost = 'lostTracks'

# handleComb = Handle('std::vector<reco::PFTau')
# labelComb = 'combinatoricRecoTaus'

# handleCH = Handle('reco::PFJetChargedHadronAssociation')
# labelCH = 'ak4PFJetsRecoTauChargedHadrons'

handleJets = Handle('std::vector<pat::Jet>')
labelJets = 'slimmedJets'


handleGen = Handle('std::vector<reco::GenParticle>')
labelGen = 'prunedGenParticles'

countMatched = 0
countUnMatched = 0

countMatched2 = 0
countUnMatched2 = 0

countGenMatched1 = 0
countGenUnMatched1 = 0

countGenMatched2 = 0
countGenUnMatched2 = 0

count_eg_in_cone = 0
count_g_in_cone = 0

count_electrons = 0
count_electrons_eg = 0
count_muons = 0
count_neutral = 0
count_unknown = 0
count_missing = 0


for ev in events:
    ev.getByLabel(label, handle)
    taus = handle.product()
    taus = [tau for tau in taus]# if tau.pt() > 18.]

    # ev.getByLabel(labelComb, handleComb)
    # c_taus = handleComb.product()

    ev.getByLabel(labelJets, handleJets)
    jets = handleJets.product()

    # ev.getByLabel(labelCH, handleCH)
    # tau_chs = handleCH.product()

    ev.getByLabel(labelPF, handlePF)
    pfcands = handlePF.product()

    ev.getByLabel(labelLost, handleLost)
    pflost = handleLost.product()

    ev.getByLabel(labelGen, handleGen)
    gens = handleGen.product()

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

                    

    chs = pfcands
    # chs = [p for p in pfcands if abs(p.pdgId()) in [11, 13, 211]]
    # chs = [p for p in pfcands if p.charge()!=0]

    # print 're-reco:', [tau.pt() for tau in taus]
    # print 'PAT ori:', [tau.pt() for tau in taus2], '\n'
    # print 'Re-miniAOD reco:'
    # for tau in taus:
    #   print 'pT {:.2f} eta {:.2f} phi {:.2f} dm {}'.format(tau.pt(), tau.eta(), tau.phi(), tau.decayMode())
    #   import pdb; pdb.set_trace()
    # print 'PAT original:'
    # for tau in taus2:
    #   print 'pT {:.2f} eta {:.2f} phi {:.2f} dm {}'.format(tau.pt(), tau.eta(), tau.phi(), tau.decayMode())

    
    # for tau2 in taus2:
    #     if not (tau2.genJet() and tau2.genJet().pt() > 20.):
    #         continue
    #     for tau1 in taus:
    #         if deltaR(tau1, tau2) < 0.3:
    #             countMatched += 1
    #             break
    #     else:
    #         countUnMatched += 1
    
    for tau1 in taus:    
        if not (tau1.genJet() and tau1.genJet().pt() > 20.):
            continue
        # for tau2 in taus2:    
        #     if deltaR(tau1, tau2) < 0.3:
        #         countMatched2 += 1
        #         if tau1.leadChargedHadrCand().pdgId() == 22:
        #             print 'Photon-induced tau with pT', tau1.pt(), tau1.decayMode()
        #             print 'Matched to tau with pT', tau2.pt(), tau2.decayMode()
        #             print 'Gen tau pT', tau1.genJet().pt()
        #         break
        # else:
        #     countUnMatched2 += 1
        #     # import pdb; pdb.set_trace()
        #     print 'Found a photon-induced tau?', tau1.leadChargedHadrCand().pdgId()
        #     print tau1.decayMode(), tau1.pt(), tau1.genJet().pt()

    for gen_tau in vis_had_taus:
        if gen_tau.pt() > 25. and abs(gen_tau.eta()) < 2.3:
            if gen_tau.dm not in [0, 1, 2, 3, 10, 11]:
                continue
            matched_tau1 = None
            matched_tau2 = None
            best_tau1 = None
            for tau1 in taus:
                if deltaR(tau1, gen_tau) < 0.3:
                    best_tau1 = tau1
                    if tau1.pt() > 0.8 * gen_tau.pt():
                        matched_tau1 = tau1
            # for tau2 in taus2:
            #     if deltaR(tau2, gen_tau) < 0.3 and tau2.pt() > 0.8 * gen_tau.pt():
            #         matched_tau2 = tau2

            if matched_tau1:
                countGenMatched1 += 1
                for ch in matched_tau1.signalChargedHadrCands():
                    pdg_id = abs(ch.pdgId())
                    if pdg_id == 11:
                        count_electrons += 1
                        if ch.isGoodEgamma():
                            count_electrons_eg += 1
                    elif pdg_id == 13:
                        count_muons += 1
                    elif pdg_id == 130:
                        count_neutral += 1
                    elif pdg_id != 211:
                        count_unknown += 1
                if matched_tau1.decayMode() in [0, 1, 2]:
                    if len(matched_tau1.signalChargedHadrCands()) != 1:
                        count_missing += 1
                elif matched_tau1.decayMode() in [10, 11]:
                    if len(matched_tau1.signalChargedHadrCands()) != 3:
                        count_missing += 1
            else:
                countGenUnMatched1 += 1
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
                        print '    {:.2f} {:.2f} {:.2f} {} {}'.format(pf.pt(), pf.eta(), pf.phi(), pf.pdgId(), pf.isGoodEgamma())
                # print '  ## Combo taus in cone'
                # for c_tau in c_taus:
                #     if deltaR(gen_tau, c_tau) < 0.3 and c_tau.pt() > 20.:
                #         print '    {:.2f} {:.2f} {:.2f} {:.2f} {}'.format(c_tau.pt(), c_tau.eta(), c_tau.phi(), c_tau.mass(), c_tau.decayMode())
                eg_in_cone = [pf for pf in pfcands if deltaR(gen_tau, pf)<0.3 and pf.isGoodEgamma()] 
                if eg_in_cone:
                    count_eg_in_cone += 1
                if eg_in_cone and any(eg.pdgId() == 22 for eg in eg_in_cone):
                    count_g_in_cone += 1

            if matched_tau2:
                countGenMatched2 += 1
            else:
                countGenUnMatched2 += 1


print 'Matched (re-reco to reco)', countMatched
print 'Unmatched (re-reco to reco)', countUnMatched

print 'Matched (reco to re-reco)', countMatched2
print 'Unmatched (reco to re-reco)', countUnMatched2

print 'Matched (re-reco to gen)', countGenMatched1
print 'Unmatched (re-reco to gen)', countGenUnMatched1

print 'Of which with good e/g', count_eg_in_cone
print 'Of which with good g', count_g_in_cone

print 'Matched (reco to gen)', countGenMatched2
print 'Unmatched (reco to gen)', countGenUnMatched2

print 'N electrons', count_electrons
print 'Of which with good e/g', count_electrons_eg
print 'N muons', count_muons
print 'N neutral', count_neutral
print 'N unknown', count_unknown
print 'N missing', count_missing
