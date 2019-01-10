import numpy as num

from ROOT import TFile, gROOT, TTree

from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR


class GenTau(object):
    def __init__(self, p4):
        self.p4 = p4

    def __getattr__(self,name):
        return getattr(self.p4, name)


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
                    pass
                    # import pdb; pdb.set_trace()

    return n_charged, n_pizero

gROOT.SetBatch(True)

RelVal = 'CMSSW_10_4_0_pre3'

electronH = Handle('vector<reco::GsfElectron>')
vertexH = Handle('std::vector<reco::Vertex>')
handleGen = Handle('std::vector<reco::GenParticle>')
labelGen = 'prunedGenParticles'

ztt_dir = '/eos/cms/store/relval/CMSSW_10_4_0_pre3/RelValZTT_13/GEN-SIM-RECO/PU25ns_103X_mc2017_realistic_v2-v1/20000/'
ztt_files = [
    '1F7B665B-C486-6140-83C4-3FD5B60850DB.root',
    '23BF126F-1763-424E-BF01-BB46ABF71F72.root',
    '26F61E4A-4D61-A548-8F81-3E61A59245F6.root',
    '50B7A671-CA0C-144E-9326-4F95A81308EC.root',
    '5D811999-BBDE-1942-89CB-117EB69F93E2.root',
    'BC9976E6-AB29-6548-8AA9-2F11944B2C5E.root',
]


filelist = [ztt_dir + ztt_file for ztt_file in ztt_files]

if __name__ == '__main__':

    events = Events(filelist)
    print len(filelist), 'files will be analyzed'

    outputname = 'Myroot_tauPOG_' + RelVal + '_tautau.root'

    out_file = TFile(outputname, 'recreate')
    tree = TTree('tree', 'tree')


    gen_tau_pt = num.zeros(1, dtype=float)
    gen_tau_eta = num.zeros(1, dtype=float)
    gen_tau_phi = num.zeros(1, dtype=float)
    gen_tau_dm = num.zeros(1, dtype=float)

    electron_pt = num.zeros(1, dtype=float)
    electron_eta = num.zeros(1, dtype=float)
    electron_phi = num.zeros(1, dtype=float)
    electron_q = num.zeros(1, dtype=int)


    tree.Branch('gen_tau_pt', gen_tau_pt, 'gen_tau_pt/D')
    tree.Branch('gen_tau_eta', gen_tau_eta, 'gen_tau_eta/D')
    tree.Branch('gen_tau_phi', gen_tau_phi, 'gen_tau_phi/D')
    tree.Branch('gen_tau_dm', gen_tau_dm, 'gen_tau_dm/D')

    tree.Branch('electron_pt', electron_pt, 'electron_pt/D')
    tree.Branch('electron_eta', electron_eta, 'electron_eta/D')
    tree.Branch('electron_phi', electron_phi, 'electron_phi/D')
    tree.Branch('electron_q', electron_q, 'electron_q/I')


    evtid = 0

    for event in events:

        evtid += 1

        if evtid % 1000 == 0:
            print 'Event ', evtid, 'processed'
    #        break

        event.getByLabel("gedGsfElectrons", electronH)
        event.getByLabel("offlineSlimmedPrimaryVertices", vertexH)
        event.getByLabel(labelGen, handleGen)
        gens = handleGen.product()

        electrons = electronH.product()
        vertices = vertexH.product()

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
                    vis_tau.pizeros = [daughters(d) for d in ds if d.pdgId() == 111]

        for i_gen_tau, gen_tau in enumerate(vis_had_taus):
            if gen_tau.pt() > 20. and abs(gen_tau.eta()) < 2.3:
                gen_tau_pt[0] = gen_tau.pt()
                gen_tau_eta[0] = gen_tau.eta()
                gen_tau_phi[0] = gen_tau.phi()
                gen_tau_dm[0] = gen_tau.dm

                for electron in electrons:
                    if deltaR(electron, gen_tau) < 0.5:
                        electron_pt[0] = electron.pt()
                        electron_eta[0] = electron.eta()
                        electron_phi[0] = electron.phi()
                        electron_q[0] = electron.charge()
                        break # Take first matching electron
            
                tree.Fill()
                electron_pt[0] = -99.
                electron_eta[0] = -99.
                electron_phi[0] = -99.
                electron_q[0] = -99.



    print evtid, 'events are processed !'

    out_file.Write()
    out_file.Close()
