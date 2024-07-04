from DataFormats.FWLite import Events, Handle
# events = Events('/eos/user/s/steggema/pickevents_merged_2022.root')
events = Events('/eos/user/a/amascell/picked_events/pickevents_merged_2022_sideband.root')

h_m = Handle('std::vector<pat::Muon>')
h_d = Handle('std::vector<pat::Muon>')
h_t = Handle('std::vector<reco::Track>')
h_j = Handle('std::vector<pat::Jet>')

for i, ev in enumerate(events):
    print(f'\n  ### EVENT {i} ###\n')
    ev.getByLabel('slimmedMuons', h_m)
    ev.getByLabel('slimmedDisplacedMuons', h_d)
    ev.getByLabel('displacedGlobalMuons', h_t)
    ev.getByLabel('slimmedJets', h_j)
    
    muons = h_m.product()
    displaced = h_d.product()
    tracks = h_t.product()
    jets = h_j.product()

    for muon in muons:
        if muon.pt() < 4.:
            continue
        print('   MUONS')
        print(f'pt: {muon.pt():.2f}, eta: {muon.eta():.2f}, phi: {muon.phi():.2f}, charge: {muon.charge()}')
        print(f'Loose: {muon.isLooseMuon()}, medium: {muon.isMediumMuon()}, normchi2: {muon.bestTrack().normalizedChi2():.2f}, chi2loc: {muon.combinedQuality().chi2LocalPosition:.2f} kink: {muon.combinedQuality().trkKink:.2f}, matched stat: {muon.numberOfMatchedStations()}')
        print(f'vx: {muon.vx():.2f}, vy: {muon.vy():.2f}, vz: {muon.vz():.2f}')
        print(f'db: {muon.dB():.2f}, edb: {muon.edB():.2f}, dz error: {muon.dzError():.2f}')
        print(f'glo: {muon.isGlobalMuon()}, seg_comp: {muon.segmentCompatibility():.2f}, time ndof: {muon.time().nDof:.2f}, time in out: {muon.time().timeAtIpInOut:.2f}')
        if muon.innerTrack().isAvailable():
            print(f'inner val frac: {muon.innerTrack().validFraction():.2f}, n_pixel: {muon.innerTrack().hitPattern().numberOfValidPixelHits():.2f}, layers w/ meas: {muon.innerTrack().hitPattern().trackerLayersWithMeasurement():.2f}, time in out: {muon.time().timeAtIpInOut:.2f}')
        if muon.isGlobalMuon():
            print(f'muon hits in fit: {muon.globalTrack().hitPattern().numberOfValidMuonHits()}')
        print(f'sta: {muon.isStandAloneMuon()}, n_valid_hits: {muon.numberOfValidHits()}, isPF: {muon.isPFMuon()}')
        print(f'isolationR03: {muon.isolationR03().sumPt:.2f}, isolationR05: {muon.isolationR05().sumPt:.2f}')
        print(f'pfisolationR03 ch: {muon.pfIsolationR03().sumChargedParticlePt:.2f}, pfisolationR03 ph: {muon.pfIsolationR03().sumPhotonEt:.2f}, , pfisolationR03 pu: {muon.pfIsolationR03().sumPUPt:.2f}')
        print(f'pfisolationR04 ch: {muon.pfIsolationR04().sumChargedParticlePt:.2f}, pfisolationR04 ph: {muon.pfIsolationR04().sumPhotonEt:.2f}, , pfisolationR04 pu: {muon.pfIsolationR04().sumPUPt:.2f}')
    print('\n')
    for muon in displaced:
        print('   DISPLACED MUONS')
        print(f'pt: {muon.pt():.2f}, eta: {muon.eta():.2f}, phi: {muon.phi():.2f}, charge: {muon.charge()}')
        print(f'Loose: {muon.isLooseMuon():.2f}, medium: {muon.isMediumMuon():.2f}, normchi2: {muon.bestTrack().normalizedChi2():.2f}, chi2loc: {muon.combinedQuality().chi2LocalPosition:.2f} kink: {muon.combinedQuality().trkKink:.2f}, matched stat: {muon.numberOfMatchedStations()}')
        print(f'vx: {muon.vx():.2f}, vy: {muon.vy():.2f}, vz: {muon.vz():.2f}')
        print(f'db: {muon.dB():.2f}, edb: {muon.edB():.2f}, dz error: {muon.dzError():.2f}')
        print(f'glo: {muon.isGlobalMuon()}, tra: {muon.isTrackerMuon()}, seg_comp: {muon.segmentCompatibility():.2f}, time ndof: {muon.time().nDof:.2f}, time in out: {muon.time().timeAtIpInOut:.2f}')
        if muon.innerTrack().isAvailable():
            print(f'inner val frac: {muon.innerTrack().validFraction():.2f}, n_pixel: {muon.innerTrack().hitPattern().numberOfValidPixelHits():.2f}, layers w/ meas: {muon.innerTrack().hitPattern().trackerLayersWithMeasurement():.2f}, time in out: {muon.time().timeAtIpInOut:.2f}, chi2 norm: {muon.innerTrack().normalizedChi2():.2f}')
            print(f'inner charge: {muon.innerTrack().charge()}, inner pt: {muon.innerTrack().pt():.2f}, inner pt error: {muon.innerTrack().ptError():.2f},inner eta: {muon.innerTrack().eta():.2f}, inner phi: {muon.innerTrack().phi():.2f}')
            print(f'inner dxy: {muon.innerTrack().dxy():.2f},inner dxy error: {muon.innerTrack().dxyError():.2f},  inner n hits: {muon.innerTrack().numberOfValidHits():.2f}, inner dz: {muon.innerTrack().dz():.2f}, inner missing: {muon.innerTrack().missingInnerHits()}')
        print(f'sta: {muon.isStandAloneMuon()}, n_valid_hits: {muon.numberOfValidHits()}, isPF: {muon.isPFMuon()}')
        if muon.isStandAloneMuon():
            print(f'outer charge: {muon.outerTrack().charge()}, outer pt: {muon.outerTrack().pt():.2f}, outer eta: {muon.outerTrack().eta():.2f}, outer phi: {muon.outerTrack().phi():.2f}')
        print(f'isolationR03: {muon.isolationR03().sumPt:.2f}, isolationR05: {muon.isolationR05().sumPt:.2f}')
        print(f'pfisolationR03 ch: {muon.pfIsolationR03().sumChargedParticlePt:.2f}, pfisolationR03 ph: {muon.pfIsolationR03().sumPhotonEt:.2f}, , pfisolationR03 pu: {muon.pfIsolationR03().sumPUPt:.2f}')
        print(f'pfisolationR04 ch: {muon.pfIsolationR04().sumChargedParticlePt:.2f}, pfisolationR04 ph: {muon.pfIsolationR04().sumPhotonEt:.2f}, , pfisolationR04 pu: {muon.pfIsolationR04().sumPUPt:.2f}')
    print('\n')
    for jet in jets:
        if jet.pt() > 15.:
            print('    JETS')
            print(f'pt: {jet.pt():.2f}, eta: {jet.eta():.2f}, phi: {jet.phi():.2f}, n dau: {muon.numberOfDaughters()}')