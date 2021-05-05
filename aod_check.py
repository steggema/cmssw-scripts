import ROOT
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR

ROOT.gSystem.Load("libDataFormatsMuonReco")
ROOT.gInterpreter.ProcessLine("#include <DataFormats/MuonReco/interface/MuonSelectors.h>")

events = Events(['file:/eos/user/a/amascell/HNL/pickevents/pickevents_2016Wplus_sigReg.root', 'file:/eos/user/a/amascell/HNL/pickevents/pickevents_2016Wminus_sigReg.root'])

h_gen = Handle('std::vector<reco::GenParticle>')
l_gen = 'genParticles'
h_dsa = Handle('std::vector<reco::Track>')
l_dsa = 'displacedStandAloneMuons'
h_sta = Handle('std::vector<reco::Track>')
l_sta = 'standAloneMuons'
h_glo = Handle('std::vector<reco::Track>')
l_glo = 'globalMuons'
h_muo = Handle('std::vector<reco::Muon>')
l_muo = 'muons'
h_sim = Handle('edm::ValueMap<reco::MuonSimInfo>')
l_sim = "muonSimClassifier" 
    

def finalDaughters(gen, daughters=None):
    if daughters is None:
        daughters = []
    for i in range(gen.numberOfDaughters()):
        daughter = gen.daughter(i)
        if daughter.numberOfDaughters() == 0:
            daughters.append(daughter)
        else:
            finalDaughters(daughter, daughters)

    return daughters


def visibleP4(gen):
    gen.final_ds = finalDaughters(gen)
    return sum((d.p4() for d in gen.final_ds if abs(d.pdgId()) not in [12, 14, 16]), ROOT.math.XYZTLorentzVectorD())


for n_ev, ev in enumerate(events):
    print '\n#### Analysing event #', n_ev
    ev.getByLabel(l_gen, h_gen)
    ev.getByLabel(l_dsa, h_dsa)
    ev.getByLabel(l_sta, h_sta)
    ev.getByLabel(l_glo, h_glo)
    ev.getByLabel(l_muo, h_muo)
    ev.getByLabel(l_sim, h_sim)
    gen = h_gen.product()
    dsa = h_dsa.product()
    glo = h_glo.product()
    muo = h_muo.product()
    sim = h_sim.product()
    gen = [g for g in gen if abs(g.pdgId()) == 13 and g.isLastCopy()]
    # gen = [g for g in gen if g.isLastCopy()]

    for d in dsa:
        if d.pt() > 5. and d.numberOfValidHits() >= 15: 
            print '\nDSA muon: {:.2f} {:.2f} {:.2f}'.format(d.pt(), d.eta(), d.phi())
            for g in gen:
                if deltaR(g, d) < 0.5 and (g.pt() > 1. or g.mass() > 1.):
                    print 'Found matched gen particle: {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {} {}'.format(g.pt(), g.eta(), g.phi(), g.vx(), g.vy(), g.vz(), g.status(), g.pdgId())
            
            # for s in sta:
            #     if deltaR(s, d) < 0.5 and s.pt() > 1.:
            #         print 'Found matched STA muon: {:.2f} {:.2f} {:.2f}'.format(s.pt(), s.eta(), s.phi())

            # for g in glo:
            #     if deltaR(g, d) < 0.5 and g.pt() > 1.:
            #         print 'Found matched GLO muon: {:.2f} {:.2f} {:.2f}'.format(g.pt(), g.eta(), g.phi())

            for i_m, m in enumerate(muo):
                if deltaR(m, d) < 0.5 and m.pt() > 1.:
                    print ' Found matched muon: {:.2f} {:.2f} {:.2f}'.format(m.pt(), m.eta(), m.phi())
                    print '   global? tracker? standalone? type', m.isGlobalMuon(), m.isTrackerMuon(), m.isStandAloneMuon(), m.type()
                    print '   loose? medium? tight? soft cb? in time?', m.passed(1<<0), m.passed(1<<1), m.passed(1<<3), m.passed(1<<13), m.passed(1<<23)

                    print '   chi2, nhits, pt error, seg comp ', m.bestTrack().normalizedChi2(), m.bestTrack().numberOfValidHits(), m.bestTrack().ptError(), ROOT.muon.segmentCompatibility(m)

                    si = sim.get(i_m)
                    print '    Sim classes:', si.primaryClass, si.extendedClass
                    print '    Sim info (1):', si.flavour, si.pdgId, si.g4processType, si.motherPdgId, si.grandMotherPdgId, si.motherFlavour, si.heaviestMotherFlavour
                    print '    Sim info (2):',si.tpId, si.tpEvent, si.tpBX, si.charge

            for i in range(d.recHitsSize()):
                rh = d.recHit(i)
                if rh.isValid():
                    if rh.geographicalId().subdetId() == 3: #RPC
                        print 'RPC: i, Time, Error, clustersize, bunch time', i, rh.time(), rh.timeError(), rh.clusterSize(), rh.BunchX()*25.
                    elif rh.geographicalId().subdetId() == 2: #CSC
                        print 'CSC: time, chi2', rh.time(), rh.chi2()
                    elif rh.geographicalId().subdetId() == 1: #DT
                        phi_segm = rh.phiSegment()
                        z_segm = rh.zSegment()
                        print 'DT: i, t0 phi, t0 z', i, phi_segm.t0() if phi_segm else -1., z_segm.t0() if z_segm else -1.
