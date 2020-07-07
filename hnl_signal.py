#!/usr/bin/env python

''' Investigate signal samples.
Authors: Jan Steggemann.
'''

import math
import copy

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, bestMatch, deltaR2
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from Var import Var

ROOT.gROOT.SetBatch(True)

sign = lambda x: math.copysign(1, x)

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


def isGenLepton(lep_cand, pid):
    return (abs(lep_cand.pdgId()) == pid and
            lep_cand.status() == 1 and
            (lep_cand.isPromptFinalState() or lep_cand.isDirectPromptTauDecayProductFinalState()))

def mother(gen):
    mo = gen.mother()
    if mo.pdgId() != gen.pdgId():
        return mo
    else:
        return mother(mo)

if __name__ == '__main__':
    filelist = ['/eos/user/s/steggema/HNL/samples/M3mu.root']

    print len(filelist), "files will be analyzed:", filelist
    events = Events(filelist)
    out_file_name = 'hnl.root'
    print "out_file_name:", out_file_name

    out_file = ROOT.TFile(out_file_name, 'recreate')

    tree = ROOT.TTree('tree', 'tree')

    all_vars = [
        Var('event_no', int),
        Var('event_id', int),
    ]

    lepton_vars = ['pt', 'eta', 'phi', 'd0']
    lepton_colls = ['l0_muon', 'l1_muon', 'l2_muon', 'l0_sta', 'l1_sta', 'l2_sta']
    all_vars += [Var('_'.join([coll, var]), float) for coll in lepton_colls for var in lepton_vars]

    gen_colls = ['l0_gen', 'l1_gen', 'l2_gen'] # l0 = l from W decay, l1 = l from prompt HNL decay, l2 = l from non-prompt HNL decay
    gen_vars = ['pt', 'eta', 'phi', 'id', 'flightlength', 'flightlengthxy']

    all_vars += [Var('_'.join([coll, var]), float) for coll in gen_colls for var in gen_vars]

    all_vars += [
        Var('gen_hnl_pt'),
        Var('gen_neutrino_id')
    ]

    avd = {var.name: var for var in all_vars} # all var dict

    for var in all_vars:
        tree.Branch(var.name, var.storage, var.name +
                        '/' + ('I' if var.type == int else 'D'))

    evtid = 0

    h_muon = Handle('vector<pat::Muon>')
    h_stamuon = Handle('vector<reco::Track>')
    h_gen = Handle('vector<reco::GenParticle>')

    for event in events:
        for var in all_vars:
            var.reset()

        evtid += 1
        eid = event.eventAuxiliary().id().event()
        avd['event_no'].fill(evtid)
        avd['event_id'].fill(eid)


        if evtid % 1000 == 0:
            print 'Event ', evtid, 'processed'

        event.getByLabel("slimmedMuons", h_muon)
        event.getByLabel("displacedStandAloneMuons", h_stamuon)
        event.getByLabel("prunedGenParticles", h_gen)

        muons = h_muon.product()
        stamuons = h_stamuon.product()
        gens = h_gen.product()

        gen_muons = [p for p in gens if isGenLepton(p, 13) or isGenLepton(p, 11)]

        gen_hnl = [p for p in gens if p.pdgId() == 9900012 and p.isLastCopy()]
        if len(gen_hnl) == 1:
            gen_hnl = gen_hnl[0]
        else:
            import pdb; pdb.set_trace()

        hnl_daughters  = [gen_hnl.daughter(i) for i in range(gen_hnl.numberOfDaughters())]
        
        neutrino = [d for d in hnl_daughters if d.pdgId()%2 == 0][0]
        lepton_from_virtualw = [d for d in hnl_daughters if d.pdgId()%2 == 1 and sign(d.pdgId()) != sign(neutrino.pdgId())][0]
        lepton_from_directhnldecay = [d for d in hnl_daughters if d.pdgId()%2 == 1 and sign(d.pdgId()) == sign(neutrino.pdgId())][0]

        avd['gen_hnl_pt'] = gen_hnl.pt()
        avd['gen_neutrino_id'].fill(neutrino.pdgId())
        # print 'Neutrino', neutrino
        # print 'Lepton from virtual W', lepton_from_virtualw
        # print 'Lepton from direct HNL decay', lepton_from_directhnldecay

        gen_w = [p for p in gens if abs(p.pdgId()) == 24 and p.isLastCopy()]
        if len(gen_w) != 1:
            print 'No W found, taking first quark'
            gen_w = [p for p in gens if abs(p.pdgId()) < 6]
        
        w = gen_w[0]
        try:
            prompt_lepton = [p for p in gen_muons if abs(mother(p).pdgId()) in [24, 1, 2, 3, 4, 5]][0]
        except IndexError:
            import pdb; pdb.set_trace()

        for label, gen in [('l0', prompt_lepton), ('l1', lepton_from_directhnldecay), ('l2', lepton_from_virtualw)]:
            avd[label+'_gen_pt'].fill(gen.pt())
            avd[label+'_gen_eta'].fill(gen.eta())
            avd[label+'_gen_phi'].fill(gen.phi())
            avd[label+'_gen_id'].fill(gen.pdgId())

            flight_length = math.sqrt((gen.vx() - w.vx())**2 + (gen.vy() - w.vy())**2 + (gen.vz() - w.vz())**2)
            flight_lengthxy = math.sqrt((gen.vx() - w.vx())**2 + (gen.vy() - w.vy())**2)

            avd[label+'_gen_flightlength'].fill(flight_length)
            avd[label+'_gen_flightlengthxy'].fill(flight_lengthxy)

            for rec_label, coll in [('muon', muons), ('sta', stamuons)]:
                muon, dr = bestMatch(gen, coll)
                if dr < 0.2:
                    avd['_'.join([label, rec_label, 'pt'])].fill(muon.pt())
                    avd['_'.join([label, rec_label, 'eta'])].fill(muon.eta())
                    avd['_'.join([label, rec_label, 'phi'])].fill(muon.phi())
                    avd['_'.join([label, rec_label, 'd0'])].fill(muon.dB() if rec_label == 'muon' else muon.dxy())

        tree.Fill()
    print evtid, 'events are processed !'

    out_file.Write()
    out_file.Close()
