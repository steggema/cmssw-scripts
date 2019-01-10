from DataFormats.FWLite import Events, Handle

if __name__ == '__main__':

    signal = True

    f_name = '/eos/cms/store/relval/CMSSW_10_3_0/RelValZTT_13/GEN-SIM-RECO/PUpmx25ns_103X_upgrade2018_realistic_v7-v1/10000/4BBB7042-E996-654D-8A04-6756881BC02B.root'

    events = Events(f_name)

    handlePF = Handle('std::vector<reco::PFCandidate>')
    labelPF = 'particleFlow'

    n_ev = 0
    for n_ev, ev in enumerate(events):
        if n_ev > 1000:
            break

        ev.getByLabel(labelPF, handlePF)
        pfs = handlePF.product()

        for pf in pfs:
            if pf.pdgId() == 22 and abs(pf.eta()) > 3.3:
                print 'eta', '{:.2f}'.format(pf.eta())
                print 'ECAL eta', '{:.2f}'.format(pf.positionAtECALEntrance().eta())
                print 'Vertex z', '{:.2f}'.format(pf.vz())
                print 'ECAL energy', '{:.2f}'.format(pf.ecalEnergy())
                print 'energy', '{:.2f}'.format(pf.energy()), '\n'
                print pf.bestTrack()
