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

events = Events('outputFULL.root')
handle = Handle('std::vector<pat::Tau>')
label = ('patTaus', '', 'TAURECO')
handle2 = Handle('std::vector<pat::Tau>')
label2 = ('slimmedTaus', '', 'PAT')

handlePF = Handle('std::vector<pat::PackedCandidate>')
labelPF = 'packedPFCandidates'

handleLost = Handle('std::vector<pat::PackedCandidate>')
labelLost = 'lostTracks'

handleComb = Handle('std::vector<reco::PFBaseTau')
labelComb = 'combinatoricRecoTaus'

handleCH = Handle('reco::PFJetChargedHadronAssociation')
labelCH = 'ak4PFJetsRecoTauChargedHadrons'

handleJets = Handle('std::vector<pat::Jet>')
labelJets = 'slimmedJets'


countMatched = 0
countUnMatched = 0

for ev in events:
    ev.getByLabel(label, handle)
    ev.getByLabel(label2, handle2)
    taus = handle.product()
    taus2 = handle2.product()
    taus = [tau for tau in taus if tau.pt() > 18.]
    
    ev.getByLabel(labelComb, handleComb)
    c_taus = handleComb.product()

    ev.getByLabel(labelJets, handleJets)
    jets = handleJets.product()

    ev.getByLabel(labelCH, handleCH)
    tau_chs = handleCH.product()

    ev.getByLabel(labelPF, handlePF)
    pfcands = handlePF.product()

    ev.getByLabel(labelLost, handleLost)
    pflost = handleLost.product()
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

    
    for tau2 in taus2:
        for tau1 in taus:
            if deltaR(tau1, tau2) < 0.3:
                if False:
                    print
                    print 'Comparing miniAOD tau (first) and reco tau (second)'
                    print '{:.2f} {:.2f} {}'.format(tau1.pt(), tau2.pt(), 'pt')
                    print '{:.2f} {:.2f} {}'.format(tau1.eta(), tau2.eta(), 'eta')
                    print '{:.2f} {:.2f} {}'.format(tau1.phi(), tau2.phi(), 'phi')
                    print '{:.2f} {:.2f} {}'.format(tau1.vertex().z(), tau2.vertex().z(), 'z')
                    print '{} {} {}'.format(tau1.decayMode(), tau2.decayMode(), 'decayMode')
                    print '{:.3f} {:.3f} {}'.format(tau1.tauID('chargedIsoPtSum'), tau2.tauID('chargedIsoPtSum'), 'chargedIsoPtSum')
                    
                    print '  REGULAR PF Collection'
                    print '   Printing pdgId, pt(), normalized chi2, dz, dxy for all charged PF cands'
                    for ch in chs:
                        # if deltaR(ch, tau1) < 0.5 and abs(tau1.vertex().z() - ch.vertex().z())<0.2:
                        # #     print ch.pdgId(), ch.pt(), ch.pseudoTrack().normalizedChi2() if ch.hasTrackDetails() else -1.
                        # if deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [11,13,211] and ch.hasTrackDetails():
                        #     print '    ', ch.pdgId(), '{:.2f}'.format(ch.pt()), ch.pseudoTrack().normalizedChi2(), '{:.2f}'.format(ch.pseudoTrack().dz(tau1.vertex())), '{:.2f}'.format(ch.pseudoTrack().dxy(tau1.vertex()))

                        if deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [11,13,211]:
                            print '    ', ch.pdgId(), '{:.2f}'.format(ch.pt()),  '{:.2f}'.format(ch.dz(tau1.vertex())), '{:.2f}'.format(ch.dxy(tau1.vertex())), ch.trackHighPurity(), ch.hcalFraction(), ch.rawCaloFraction(), ch.numberOfPixelHits(), ch.pixelLayersWithMeasurement()
                    
                    print 'Sum CH', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(tau1.vertex().z() - ch.vertex().z())<0.2 and abs(ch.pdgId()) in [11, 13, 211]) - sum(p.pt() for p in tau1.signalChargedHadrCands())
                    print 'Sum NH', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [130])
                    print 'Sum ph', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [22])

                    print '   LOST tracks collection - printing everything in deltaR(tau, track)<0.5'
                    for ch_lost in pflost:
                        if deltaR(ch_lost, tau1) < 0.5: # and abs(tau1.vertex().z() - ch_lost.vertex().z())<0.2:
                            print ch_lost.pdgId(), ch_lost.pt()

                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('footprintCorrection'), tau2.tauID('footprintCorrection'), 'footprintCorrection')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('photonPtSumOutsideSignalCone'), tau2.tauID('photonPtSumOutsideSignalCone'), 'photonPtSumOutsideSignalCone')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('neutralIsoPtSum'), tau2.tauID('neutralIsoPtSum'), 'neutralIsoPtSum')
                    
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits'), tau2.tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits'), 'byCombinedIsolationDeltaBetaCorrRaw3Hits')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw'), tau2.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw'), 'byIsolationMVArun2v1DBdR03oldDMwLTraw')
                    print
                countMatched += 1
                break
        else:
            print '\n    No match found for reco tau', '{:.2f} {:.2f} {:.2f} {:.2f}'.format(tau2.pt(),tau2.eta(), tau2.phi(), tau2.mass()), tau2.decayMode()
            print '    Matching jet', '{:.2f} {:.2f} {:.2f}'.format(tau2.p4Jet().pt(), tau2.p4Jet().eta(), tau2.p4Jet().phi())
            lch = tau2.leadChargedHadrCand()
            # import pdb; pdb.set_trace()
            print '    Lead CH {:.2f} {:.2f} {:.2f}'.format(lch.pt(),lch.eta(), lch.phi()), lch.pdgId()
            countUnMatched += 1
            printCHInfo = True
            if printCHInfo:
                print '    REGULAR'
                for ch in chs:
                    if deltaR(ch, tau2) < 0.5 and (ch.charge()==0 or abs(tau2.vertex().z() - ch.vertex().z())<0.2):
                        print ch.pdgId(), ch.pt()
                        # if abs(ch.pdgId())==11:
                            # import pdb; pdb.set_trace()
                    
                print '    LOST'
                for ch_lost in pflost:
                    if deltaR(ch_lost, tau2) < 0.5 and abs(tau2.vertex().z() - ch_lost.vertex().z())<0.2:
                        print '    ', ch_lost.pdgId(), ch_lost.pt()

                print '    Tau reco CH'
                for i_t, tau_ch in enumerate(tau_chs):
                    if deltaR(tau_ch.first, tau2) < 0.5:
                        print '    Matching jet', '{:.2f} {:.2f} {:.2f}'.format(tau_ch.first.pt(), tau_ch.first.eta(), tau_ch.first.phi())
                        for jet in jets:
                            if deltaR(jet, tau2) < 0.5:
                                print jet.eta(), jet.phi(), jet.pt()
                                break
                        # import pdb; pdb.set_trace()
                        for the_ch in tau_chs.value(i_t):
                            print the_ch.pdgId(), the_ch.pt(), the_ch.eta(), the_ch.phi()
                        break
                else:
                    for jet in jets:
                        if deltaR(jet, tau2) < 0.5:
                            print jet.eta(), jet.phi(), jet.pt()
                            break
                    else:
                        print 'No CHS jet found - hence probably a seeding issue!'

            print '    Combo taus'
            for c_tau in c_taus:
                if deltaR(c_tau, tau2) < 0.5:
                    print '{:.2f} {:.2f} {:.2f} {:.2f}'.format(c_tau.pt(),c_tau.eta(), c_tau.phi(), c_tau.mass()), c_tau.decayMode()
                    lch = c_tau.leadPFChargedHadrCand()
                    print '    Lead CH {:.2f} {:.2f} {:.2f}'.format(lch.pt(),lch.eta(), lch.phi()), lch.pdgId()
                    # print 'Jet', '{:.2f} {:.2f} {:.2f}'.format(c_tau.jetRef().pt(), c_tau.jetRef().eta(), c_tau.jetRef().phi())


print 'Matched', countMatched
print 'Unmatched', countUnMatched


