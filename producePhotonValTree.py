#!/usr/bin/env python

''' Produces a flat tree for tau release/data validation.
Authors: Yuta Takahashi, Michal Bluj, Jan Steggemann.
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
            (lep_cand.isPromptFinalState() or lep_cand.isDirectPromptTauDecayProductFinalState()) and
            lep_cand.pt() > 20 and
            abs(lep_cand.eta()) < 2.3)


if __name__ == '__main__':
    runtype = 'electron'
    if runtype == 'tau':
        # ztt_dir = '/eos/cms/store/relval/CMSSW_10_4_0_pre1/RelValZTT_13/GEN-SIM-RECO/PUpmx25ns_103X_upgrade2018_realistic_v7-v1/10000'
        # ztt_files = [
        #     '12A90B88-6C8C-3C4D-82AB-5EDF4912DE5B.root',
        #     '29EEE068-8D53-E540-8F9F-E5207BD9C4B9.root',
        #     '4BD7B325-E91E-3345-B02F-D6ADBAC8AE87.root',
        #     '79A166A9-8B8D-F549-864A-B6E4C723EE1D.root',
        #     '7C162A8E-4B18-BD44-9C18-B63854B5FF90.root',
        #     'AC6EE629-617A-C347-A79F-8E01503AC4B5.root',
        #     'B1AC14B0-6839-7A46-A95F-F54C8EA73052.root',
        #     'C338FB84-2626-D246-84E9-87899D81809D.root',
        #     'E6BEB3A1-8370-4149-B16D-801481807944.root',
        #     'F5B481EF-716D-0C40-BD58-9D9D6BE6842F.root',
        # ]
        ztt_dir = '/eos/cms/store/relval/CMSSW_10_4_0_pre1/RelValTenTau_15_500/GEN-SIM-RECO/PU25ns_103X_upgrade2018_realistic_v7-v1/10000/'
        ztt_files = [
            '36A8D50C-79BD-D14C-9DA7-756A3CF161A2.root',
            '3AEF2D54-047C-5446-A5C4-802DB7402ACC.root',
            '47F7023D-BE6E-B544-AFF8-DFA8DAAF2678.root',
            '486344AE-6FC0-2545-A753-C77C39154BEA.root',
            '59F14666-2784-6F44-9D0B-48DB8D804323.root',
            '9EDE7C05-41EB-1A47-8264-1061E87F2DAF.root',
            'ACD9DC22-08AB-3248-B506-B306D1B60053.root',
            'E7FEBEB5-A1B5-6A45-873E-F73623387E44.root',
            'EF6AD547-A228-CA40-A1F3-38826EAF1239.root',
        ]
        filelist = ['/'.join([ztt_dir, ztt_file]) for ztt_file in ztt_files]
        outputFileName = 'ztt_photons.root'
    elif runtype == 'electron':
        zee_dir = '/eos/cms/store/relval/CMSSW_10_4_0_pre1/RelValZEE_13/GEN-SIM-RECO/PUpmx25ns_103X_upgrade2018_realistic_v7-v1/10000'
        zee_files = [
            '2EEFFE84-90F6-E24B-89D0-91C000B8DB51.root',
            '456D3274-DA0E-C14B-9533-E7A703B3D33D.root',
            '49FA4A7A-FC30-9F48-9C1B-B64C9300CA19.root',
            '71E58EB6-E7BA-CF43-B282-53343B27AA78.root',
            '99361246-41D1-DE4D-A79F-F34D7FB68657.root',
            'BE5BF7E6-F1E4-C14D-9F08-D1CAD10DFB0D.root',
            'C39541D6-7E68-234F-85D7-79069C404234.root',
            'DB65626E-EC8B-714C-8649-CE6C769599B3.root',
            'EAACBF5C-7C70-E543-888B-CDB2AD9CF079.root',
        ]
        filelist = ['/'.join([zee_dir, zee_file]) for zee_file in zee_files]
        tt_dir = '/eos/cms/store/relval/CMSSW_10_4_0_pre1/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_103X_upgrade2018_realistic_v7-v1/10000'
        tt_files = [
            '38344961-08FF-F24C-AA49-50B22EABBDF9.root',
            '78C1BDD8-8DF1-6A40-AF45-A38CFD1C5012.root',
            '88518E15-B797-2243-A2EE-98EED6289597.root',
            '9D9F753B-BE5F-4549-B425-AF6832B4EE34.root',
            'A949E818-7F4B-3841-B9E6-A5D89469B824.root',
            'E3D7D9F5-F152-DE48-8CD0-D8FD60E3659D.root',
            'E4F4B482-6DB2-B340-AC3A-9C65F75D3812.root',
            'ECB460E0-ADD5-5A40-949C-F42C6090C873.root',
            'FC749E84-AF7D-6E47-BF1D-E7B46B962721.root',
        ]
        filelist += ['/'.join([tt_dir, tt_file]) for tt_file in tt_files]

        ten_e_dir = '/eos/cms/store/relval/CMSSW_10_3_0/RelValTenE_0_200/GEN-SIM-RECO/PU25ns_103X_upgrade2018_realistic_v7_resub-v1/10000'

        ten_e_files = [
            '0FF8567C-E108-BE42-960D-1F940A46D6F1.root',
            '1A6FDE4E-9373-0043-9880-99FF3FFE7926.root',
            '2239D7E2-3077-4444-B86B-03BABFEFFC2C.root',
            '273F0B8F-FBF8-4F47-BC12-357F8E107A1D.root',
            '3860091B-8A6B-5441-B929-F7BCDFBA259F.root',
            '4B634622-505A-9940-BAA9-8FE4B969B1F5.root',
            '6207F1C4-92C0-3643-9B06-CDEA4BE4D24E.root',
            'A8A7BC70-CA71-D34C-A0BD-C238FC1C3935.root',
            'D668E414-A18F-014A-AA0F-194C19B22C2B.root',
        ]
        filelist += ['/'.join([ten_e_dir, ten_e_file]) for ten_e_file in ten_e_files]
        outputFileName = 'zee_photons.root'


    print len(filelist), "files will be analyzed:", filelist
    events = Events(filelist)

    
    print "outputFileName:", outputFileName

    out_file = ROOT.TFile(outputFileName, 'recreate')

    tau_tree = ROOT.TTree('per_tau', 'per_tau')

    all_vars = [
        Var('event_id', int),
        Var('ev_no', int),
        Var('tau_dm', int),
        Var('photon_pt', float),
        Var('photon_eta', float),
        Var('photon_phi', float),
        Var('gen_dm', int),
        Var('gen_pt', float),
        Var('gen_eta', float),
        Var('gen_phi', float),
        Var('gen_charged_pt', float),
        Var('gen_neutral_pt', float),
        Var('n_vertices', int),
    ]

    ph_vars = ['sigmaEtaEta', 'sigmaIetaIeta', 'e1x5',  'e5x5', 'r9', 'full5x5_sigmaEtaEta', 'full5x5_sigmaIetaIeta', 'isStandardPhoton', 'hasConversionTracks', 'hasPixelSeed', 'hadronicOverEm', 'hadTowOverEm', 'mipTotEnergy']

    for ph_var in ph_vars:
        all_vars.append(Var(ph_var, float))

    all_var_dict = {var.name: var for var in all_vars}

    for var in all_vars:
        tau_tree.Branch(var.name, var.storage, var.name +
                        '/' + ('I' if var.type == int else 'D'))

    evtid = 0

    NMatchedTaus = 0

    vertexH = Handle('std::vector<reco::Vertex>')
    photonH = Handle('vector<reco::Photon>')
    genParticlesH = Handle('std::vector<reco::GenParticle>')
    pfH = Handle('vector<reco::PFCandidate>')
    tauH = Handle('vector<reco::PFTau>')

    for event in events:
        evtid += 1
        eid = event.eventAuxiliary().id().event()

        if evtid % 1000 == 0:
            print 'Event ', evtid, 'processed'

        event.getByLabel("gedPhotons", photonH)
        event.getByLabel("offlinePrimaryVertices", vertexH)
        event.getByLabel('genParticles', genParticlesH)
        event.getByLabel('particleFlow', pfH)

        photons = photonH.product()
        vertices = vertexH.product()
        genParticles = genParticlesH.product()
        pfs = pfH.product()

        photons = [e for e in photons if any(abs(pf.pdgId()) == 22 and deltaR(e, pf) < 0.005 for pf in pfs)]

        genTaus = [p for p in genParticles if abs(
            p.pdgId()) == 15 and p.isPromptDecayed()]
        genElectrons = [
            p for p in genParticles if isGenLepton(p, 11) and p.pt() > 20.]

        refObjs = []
        if runtype == 'tau':
            for gen_tau in genTaus:
                gen_tau.visP4 = visibleP4(gen_tau)

                gen_dm = tauDecayModes.genDecayModeInt(
                    [d for d in gen_tau.final_ds if abs(d.pdgId()) not in [12, 14, 16]])
                if abs(gen_tau.visP4.eta()) > 2.3:
                    continue
                if gen_tau.visP4.pt() < 20:
                    continue
                if gen_dm == -11 or gen_dm == -13:
                    continue
                # For the 10-tau sample, remove gen taus that have overlap
                if any(deltaR(other_tau, gen_tau) < 0.5 for other_tau in genTaus if other_tau is not gen_tau):
                    continue

                # Only consider gen taus that don't have a matched nice reco tau
                event.getByLabel('hpsPFTauProducer', tauH)
                taus = tauH.product()

                has_matched_tau = False
                for tau in taus:
                    if deltaR(tau, gen_tau) < 0.3 and tau.pt()/gen_tau.visP4.pt() > 0.7 and tau.decayMode() in [0, 1, 2, 3, 10]:
                        has_matched_tau = True
                        break

                if not has_matched_tau:
                    refObjs.append(gen_tau)

        elif runtype == 'electron':
            refObjs = copy.deepcopy(genElectrons)

        for refObj in refObjs:
            for var in all_vars:
                var.reset()
            all_var_dict['ev_no'].fill(evtid)
            all_var_dict['event_id'].fill(eid)
            all_var_dict['n_vertices'].fill(len(vertices))

            if runtype == 'tau':
                gen_dm = tauDecayModes.genDecayModeInt(
                    [d for d in finalDaughters(refObj) if abs(d.pdgId()) not in [12, 14, 16]])
                all_var_dict['gen_dm'].fill(gen_dm)
                all_var_dict['gen_pt'].fill(refObj.visP4.pt())
                all_var_dict['gen_eta'].fill(refObj.visP4.eta())
                all_var_dict['gen_phi'].fill(refObj.visP4.phi())
                charged_p4 = sum((d.p4() for d in refObj.final_ds if d.charge()), ROOT.math.XYZTLorentzVectorD())
                neutral_p4 = sum((d.p4() for d in refObj.final_ds if abs(d.pdgId()) not in [12, 14, 16] and not d.charge()), ROOT.math.XYZTLorentzVectorD())
                all_var_dict['gen_charged_pt'].fill(charged_p4.pt())
                all_var_dict['gen_neutral_pt'].fill(neutral_p4.pt())
            else:
                all_var_dict['gen_dm'].fill(-1)
                all_var_dict['gen_pt'].fill(refObj.pt())
                all_var_dict['gen_eta'].fill(refObj.eta())
                all_var_dict['gen_phi'].fill(refObj.phi())

            photon, dr = bestMatch(refObj, photons)
            if dr < 0.5:
                # Fill reco-tau variables if it exists...
                NMatchedTaus += 1

                all_var_dict['photon_pt'].fill(photon.pt())
                all_var_dict['photon_eta'].fill(photon.eta())
                all_var_dict['photon_phi'].fill(photon.phi())

                for ph_var in ph_vars:
                    all_var_dict[ph_var].fill(getattr(photon, ph_var)())
                tau_tree.Fill()
    print "MATCHED TAUS:", NMatchedTaus
    print evtid, 'events are processed !'

    out_file.Write()
    out_file.Close()
