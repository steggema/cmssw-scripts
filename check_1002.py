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

events = Events('miniAOD_TauReco_ak4PFJets_ggH.root')
handle = Handle('std::vector<pat::Tau>')
label = ('selectedPatTaus', '', 'TAURECO')
handle2 = Handle('std::vector<pat::Tau>')
label2 = ('slimmedTaus', '', 'RECO')

handlePF = Handle('std::vector<pat::PackedCandidate>')
labelPF = 'packedPFCandidates'

handleLost = Handle('std::vector<pat::PackedCandidate>')
labelLost = 'lostTracks'

handleJets = Handle('std::vector<pat::Jet>')
labelJets = 'slimmedJets'

handleCombo = Handle('vector<reco::PFTau>')
labelCombo = 'combinatoricRecoTaus'

handleCH = Handle('reco::PFJetChargedHadronAssociation')
labelCH = 'ak4PFJetsRecoTauChargedHadrons'

countMatched = 0
countDifferentReco = 0
countDifferentMVA = 0
countUnMatched = 0

def n_photons_tau(tau):
    n_ph = 0
    for ph in tau.signalGammaCands():
        if ph.pt() > 0.5:
            n_ph += 1
    for ph in tau.isolationGammaCands():
        if ph.pt() > 0.5:
            n_ph += 1
    return n_ph

def e_over_h(tau):
    total = tau.ecalEnergy() + tau.hcalEnergy()
    return tau.ecalEnergy()/total if total > 0. else -1.

def pt_weighted_dx(tau, mode=0, var=0, decaymode=-1):
    sum_pt = 0.;
    sum_dx_pt = 0.;
    signalrad = max(0.05, min(0.1, 3./max(1., tau.pt())))
    cands = getGammas(tau, mode < 2)
    for cand in cands:
        if cand.pt() < 0.5:
          continue;
        dr = deltaR(cand, tau)
        deta = abs(cand.eta() - tau.eta())
        dphi = abs(deltaPhi(cand.phi(), tau.phi()))
        pt = cand.pt()
        flag = isInside(pt, deta, dphi)
        if decaymode != 10:
            if mode == 2 or (mode == 0 and dr < signalrad) or (mode == 1 and dr > signalrad):
                sum_pt += pt
                if var == 0:
                  sum_dx_pt += pt * dr
                elif var == 1:
                  sum_dx_pt += pt * deta
                elif var == 2:
                  sum_dx_pt += pt * dphi
        else:
            if (mode == 2 and flag == False) or (mode == 1 and flag == True) or mode==0:
                sum_pt += pt
                if var == 0:
                  sum_dx_pt += pt * dr
                elif var == 1:
                  sum_dx_pt += pt * deta
                elif var == 2:
                  sum_dx_pt += pt * dphi
    if sum_pt > 0.:
        return sum_dx_pt/sum_pt;
    return 0.

def tau_pt_weighted_dr_iso(tau):
    return pt_weighted_dx(tau, 2, 0, tau.decayMode())

def tau_pt_weighted_dphi_strip(tau):
    dm = tau.decayMode()
    return pt_weighted_dx(tau, 2 if dm == 10 else 1, 2, dm)

def tau_pt_weighted_deta_strip(tau):
    dm = tau.decayMode()
    return pt_weighted_dx(tau, 2 if dm == 10 else 1, 1, dm)

def tau_pt_weighted_dr_signal(tau):
    return pt_weighted_dx(tau, 0, 0, tau.decayMode())

def getGammas(tau, signal=True):
    if signal:
        return tau.signalGammaCands()
    return tau.isolationGammaCands()


def isInside(photon_pt, deta, dphi):
    stripEtaAssociationDistance_0p95_p0 = 0.197077
    stripEtaAssociationDistance_0p95_p1 = 0.658701
    stripPhiAssociationDistance_0p95_p0 = 0.352476
    stripPhiAssociationDistance_0p95_p1 = 0.707716
    if photon_pt == 0.:
        return False
    if dphi < 0.3 and dphi < max(0.05, stripPhiAssociationDistance_0p95_p0*pow(photon_pt, -stripPhiAssociationDistance_0p95_p1)) and  deta<0.15 and deta<max(0.05, stripEtaAssociationDistance_0p95_p0*pow(photon_pt, -stripEtaAssociationDistance_0p95_p1)):
        return True
    return False

for ev in events:
    ev.getByLabel(label, handle)
    ev.getByLabel(label2, handle2)
    taus = handle.product()
    taus2 = handle2.product()
    taus = [tau for tau in taus if tau.pt() > 18.]
    
    ev.getByLabel(labelJets, handleJets)
    jets = handleJets.product()

    ev.getByLabel(labelPF, handlePF)
    pfcands = handlePF.product()

    ev.getByLabel(labelLost, handleLost)
    pflost = handleLost.product()
    chs = pfcands

    ev.getByLabel(labelCombo, handleCombo)
    combos = handleCombo.product()
    
    ev.getByLabel(labelCH, handleCH)
    tau_chs = handleCH.product()


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
                if abs(tau1.pt() - tau2.pt()) > 0.02 or tau1.decayMode() != tau2.decayMode() or abs(tau1.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw') -tau2.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw')) > 0.02:
                    if abs(tau1.pt() - tau2.pt()) > 0.02 or tau1.decayMode() != tau2.decayMode():
                        countDifferentReco += 1
                    else:
                        countDifferentMVA += 1
                    # print
                    print '\nComparing miniAOD tau (first) and reco tau (second)'
                    print '{:.2f} {:.2f} {}'.format(tau1.pt(), tau2.pt(), 'pt')
                    print '{:.2f} {:.2f} {}'.format(tau1.eta(), tau2.eta(), 'eta')
                    print '{:.2f} {:.2f} {}'.format(tau1.phi(), tau2.phi(), 'phi')
                    print '{:.3f} {:.3f} {}'.format(tau1.mass(), tau2.mass(), 'mass')
                    print '{} {} {}'.format(tau1.decayMode(), tau2.decayMode(), 'decayMode')
                    print '{} {} {}'.format(n_photons_tau(tau1), n_photons_tau(tau2), 'n(photons)')

                    print '{:.2f} {:.2f} {}'.format(tau1.vertex().z(), tau2.vertex().z(), 'z')

                    print '{:.3f} {:.3f} {}'.format(tau1.tauID('chargedIsoPtSum'), tau2.tauID('chargedIsoPtSum'), 'chargedIsoPtSum')
                    print '{:.3f} {:.3f} {}'.format(tau1.tauID('puCorrPtSum'), tau2.tauID('puCorrPtSum'), 'puCorrPtSum')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('footprintCorrection'), tau2.tauID('footprintCorrection'), 'footprintCorrection')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('photonPtSumOutsideSignalCone'), tau2.tauID('photonPtSumOutsideSignalCone'), 'photonPtSumOutsideSignalCone')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('neutralIsoPtSum'), tau2.tauID('neutralIsoPtSum'), 'neutralIsoPtSum')
                    
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits'), tau2.tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits'), 'byCombinedIsolationDeltaBetaCorrRaw3Hits')
                    print '{:.2f} {:.2f} {}'.format(tau1.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw'), tau2.tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw'), 'byIsolationMVArun2v1DBdR03oldDMwLTraw')
                    
                    print '{:.4f} {:.4f} {}'.format(tau1.flightLength().r(), tau2.flightLength().r(), 'flightLength')
                    print '{:.2f} {:.2f} {}'.format(tau1.flightLengthSig(), tau2.flightLengthSig(), 'flightLengthSig')
                    print '{:.4f} {:.4f} {}'.format(tau1.dxy(), tau2.dxy(), 'dxy')
                    print '{:.4f} {:.4f} {}'.format(tau1.dxy_error(), tau2.dxy_error(), 'dxy_error')
                    print '{:.2f} {:.2f} {}'.format(tau1.dxy_Sig(), tau2.dxy_Sig(), 'dxy_sig')
                    print '{:.4f} {:.4f} {}'.format(tau1.ip3d(), tau2.ip3d(), 'ip3d')
                    print '{:.4f} {:.4f} {}'.format(tau1.ip3d_error(), tau2.ip3d_error(), 'ip3d_error')
                    print '{:.2f} {:.2f} {}'.format(tau1.ip3d_Sig(), tau2.ip3d_Sig(), 'ip3d_sig')
                    print '{:.2f} {:.2f} {}'.format(tau1.leadingTrackNormChi2(), tau2.leadingTrackNormChi2(), 'leadingTrackNormChi2')

                    print '{:.2f} {:.2f} {}'.format(e_over_h(tau1), e_over_h(tau2), 'E_over_H')
                    print '{:.2f} {:.2f} {}'.format(tau_pt_weighted_dr_iso(tau1), tau_pt_weighted_dr_iso(tau2), 'ptWeightedDRIsolation')
                    print '{:.2f} {:.2f} {}'.format(tau_pt_weighted_dr_signal(tau1), tau_pt_weighted_dr_signal(tau2), 'ptWeightedDRSignal')
                    print '{:.2f} {:.2f} {}'.format(tau_pt_weighted_deta_strip(tau1), tau_pt_weighted_deta_strip(tau2), 'ptWeightedDetaStrip')
                    print '{:.2f} {:.2f} {}'.format(tau_pt_weighted_dphi_strip(tau1), tau_pt_weighted_dphi_strip(tau2), 'ptWeightedDphiStrip')
                    

                    

                    # import pdb; pdb.set_trace()
                    
                    # print '  MiniAOD tau charged hadrons'
                    # for ch in tau1.signalChargedHadrCands():
                    #     print '    {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.charge(), ch.pdgId(), '(pt eta charge pdg)')
                    # print '  RECO tau charged hadrons'
                    # for ch in tau2.signalChargedHadrCands():
                    #     print '    {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.charge(), ch.pdgId(), '(pt eta charge pdg)')
                    # print '  MiniAOD tau photons'
                    # for ch in tau1.signalGammaCands():
                    #     print '    {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.charge(), '(pt eta charge)')
                    # print '  RECO tau photons'
                    # for ch in tau2.signalGammaCands():
                    #     print '    {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.charge(), '(pt eta charge)')

                    # print 'Checking all combo taus'
                    # for c_tau in combos:
                    #     if deltaR(c_tau, tau2) < 0.3:
                    #         if abs(c_tau.pt() - tau2.pt()) < 0.02:
                    #             print '### ',
                    #         print '{:.2f} {} {:.3f} {}'.format(c_tau.pt(), c_tau.decayMode(), c_tau.mass(), '(pT, DM, mass)')
                    
                    # if True:
                    #     print '  REGULAR PF Collection'
                    #     print '   Printing pdgId, dR, pt, dz, dxy, highPurity, hcalFraction, rawCaloFraction, n(pixel), n(pixel layers) for all charged PF cands'
                    #     for ch in chs:
                    #         # if deltaR(ch, tau1) < 0.5 and abs(tau1.vertex().z() - ch.vertex().z())<0.2:
                    #         # #     print ch.pdgId(), ch.pt(), ch.pseudoTrack().normalizedChi2() if ch.hasTrackDetails() else -1.
                    #         # if deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [11,13,211] and ch.hasTrackDetails():
                    #         #     print '    ', ch.pdgId(), '{:.2f}'.format(ch.pt()), ch.pseudoTrack().normalizedChi2(), '{:.2f}'.format(ch.pseudoTrack().dz(tau1.vertex())), '{:.2f}'.format(ch.pseudoTrack().dxy(tau1.vertex()))

                    #         if deltaR(ch, tau1) < 0.5 and (abs(ch.pdgId()) in [11,13,211] or (ch.pdgId()==22 and ch.ptTrk()>0.)):
                    #             print '    ', ch.pdgId(), '{:.2f}'.format(deltaR(ch, tau1)), '{:.2f}'.format(ch.pt()),  '{:.2f}'.format(ch.dz(tau1.vertex())), '{:.2f}'.format(ch.dxy(tau1.vertex())), ch.trackHighPurity(), ch.hcalFraction(), ch.rawCaloFraction(), ch.numberOfPixelHits(), ch.pixelLayersWithMeasurement()
                    #             if ch.pdgId()==22 and ch.hasTrackDetails(): import pdb; pdb.set_trace()
                        
                    #     print '   CHECKING RECO TAU CHARGED HADRONS -----'
                    #     print '    Printing PDG id, algo, dr, pt, n(neutral)'
                    #     for i_chs in xrange(tau_chs.size()):
                    #         chs = tau_chs.value(i_chs)
                    #         for i_ch in xrange(chs.size()):
                    #             ch = chs.at(i_ch)
                    #             if deltaR(ch, tau1) < 0.5:
                    #                 print '    ', ch.pdgId(), ch.algo(), '{:.2f}.'.format(deltaR(ch, tau1) < 0.5), '{:.2f}'.format(ch.pt()), ch.getNeutralPFCandidates().size(), ch.getChargedPFCandidate().pdgId() if ch.algo() == 1 else ''

                    #     print 'Sum CH', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(tau1.vertex().z() - ch.vertex().z())<0.2 and abs(ch.pdgId()) in [11, 13, 211]) - sum(p.pt() for p in tau1.signalChargedHadrCands())
                    #     print 'Sum NH', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [130])
                    #     print 'Sum ph', sum(ch.pt() for ch in chs if  deltaR(ch, tau1) < 0.5 and abs(ch.pdgId()) in [22])

                    #     print '   LOST tracks collection - printing everything in deltaR(tau, track)<0.5'
                    #     print '    pdg Id, dR, pT, dz, dxy, high purity, f_hcal, f_rawcalo, n(pixel), n(pixel layers)'
                    #     for ch in pflost:
                    #         if deltaR(ch, tau1) < 0.5: # and abs(tau1.vertex().z() - ch_lost.vertex().z())<0.2:
                    #             print '    ', ch.pdgId(), '{:.2f}'.format(deltaR(ch, tau1)), '{:.2f}'.format(ch.pt()),  '{:.2f}'.format(ch.dz(tau1.vertex())), '{:.2f}'.format(ch.dxy(tau1.vertex())), ch.trackHighPurity(), ch.hcalFraction(), ch.rawCaloFraction(), ch.numberOfPixelHits(), ch.pixelLayersWithMeasurement()

                        
                        # print
                else:
                    countMatched += 1
                break
        else:
            print '\n    No match found for reco tau', '{:.2f} {:.2f} {:.2f} {:.2f}'.format(tau2.pt(),tau2.eta(), tau2.phi(), tau2.mass()), tau2.decayMode()
            # print '    Matching jet', '{:.2f} {:.2f} {:.2f}'.format(tau2.p4Jet().pt(), tau2.p4Jet().eta(), tau2.p4Jet().phi())
            lch = tau2.leadChargedHadrCand()
            # import pdb; pdb.set_trace()
            print '    Lead CH {:.2f} {:.2f} {:.2f}'.format(lch.pt(),lch.eta(), lch.phi()), lch.pdgId()
            countUnMatched += 1

            print '    Checking all combo taus'
            for c_tau in combos:
                if deltaR(c_tau, tau2) < 0.3:
                    if abs(c_tau.pt() - tau2.pt()) < 0.02:
                        print '### ',
                    print '      {:.2f} {} {:.3f} {}'.format(c_tau.pt(), c_tau.decayMode(), c_tau.mass(), '(pT, DM, mass)')

            printCHInfo = True
            if printCHInfo:
                print '    REGULAR'
                for ch in chs:
                    if deltaR(ch, tau2) < 0.5 and (ch.charge()==0 or abs(tau2.vertex().z() - ch.vertex().z())<0.2):
                        print '      {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.pdgId(), '(pt eta PDG)')
                        # if abs(ch.pdgId())==11:
                            # import pdb; pdb.set_trace()
                    
                print '    LOST'
                for ch_lost in pflost:
                    if deltaR(ch_lost, tau2) < 0.5 and abs(tau2.vertex().z() - ch_lost.vertex().z())<0.2:
                        print '      {:.2f} {:.2f} {} {}'.format(ch.pt(), ch.eta(), ch.pdgId(), '(pt eta PDG)')

                for jet in jets:
                    if deltaR(jet, tau2) < 0.5:
                        print '      {:.2f} {:.2f} {:.2f} {}'.format(jet.pt(), jet.eta(), jet.phi(), '(pt eta phi)')
                        break
                else:
                    print 'No CHS jet found - hence probably a seeding issue!'

total = float(countMatched + countDifferentReco + countDifferentMVA + countUnMatched)
print 'Matched', countMatched, round(countMatched/total, 2)
print 'Different reco', countDifferentReco, round(countDifferentReco/total, 2)
print 'Different MVA', countDifferentMVA, round(countDifferentMVA/total, 2)
print 'Unmatched', countUnMatched, round(countUnMatched/total, 2)


