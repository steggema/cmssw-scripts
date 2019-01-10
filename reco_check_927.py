import ROOT

from DataFormats.FWLite import Events, Handle

events = Events('/eos/cms/store/mc/RunIISummer17DRStdmix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v8-v1/110000/14F6127A-2276-E711-B1C9-D067E5F91E51.root')

handlePF = Handle('std::vector<reco::PFCandidate>')
labelPF = 'particleFlow'

h_raw_to_reco = ROOT.TH1F('raw_to_reco', '', 100, 0., 2.)
h_ratio_vs_eta = ROOT.TH2F('h_ratio_vs_eta', '', 100, 0., 2., 100, -1.4, 1.4)
h_ratio_vs_pt = ROOT.TH2F('h_ratio_vs_pt', '', 100, 0., 2., 100, 0., 50.)
h_raw_to_reco_pt05_to_10 = ROOT.TH1F('h_raw_to_reco_pt05_to_10', '', 100, 0., 2.)
h_ratio_vs_eta_pt05_to_10 = ROOT.TH2F('h_ratio_vs_eta_pt05_to_10', '', 100, 0., 2., 100, -1.4, 1.4)
h_pizero = ROOT.TH1F('pizero', '', 200, 0, 0.5)
h_pizero_ptgr15 = ROOT.TH1F('h_pizero_ptgr15', '', 200, 0, 0.5)

i_ev = 0

for event in events:
	i_ev += 1
	event.getByLabel(labelPF, handlePF)
	pfs = handlePF.product()
	photons = [pf for pf in pfs if pf.pdgId() == 22 and pf.pt() > 0.5 and pf.pt()<1.5 and abs(pf.eta()) < 1.4]

	for pf in photons:
		# print  '{:.2f} {:.2f} {:.2f} {}'.format(pf.pt(), pf.energy(), pf.eta(), pf.phi())
		# import pdb; pdb.set_trace()
		# print '{:.2f}'.format(pf.rawEcalEnergy()/pf.energy())
		ratio = pf.energy()/pf.rawEcalEnergy()
		pt = pf.pt()
		eta = pf.eta()
		# if ratio < 0.2 and pt > 0.5:
		# 	import pdb; pdb.set_trace()
		h_raw_to_reco.Fill(ratio)
		h_ratio_vs_eta.Fill(ratio, eta)
		h_ratio_vs_pt.Fill(ratio, pt)
		if pt > 0.5 and pt < 1.5:
			h_raw_to_reco_pt05_to_10.Fill(ratio)
			h_ratio_vs_eta_pt05_to_10.Fill(ratio, pf.eta())

		for pf2 in photons:
			if pf2.pdgId() == 22 and abs(pf2.eta()) < 1.4:
				if pf2 != pf:
					p4_pi0 = pf.p4() + pf2.p4()
					m_pi0 = p4_pi0.mass()
					pt_pi0 = p4_pi0.pt()
					h_pizero.Fill(m_pi0)
					if pt_pi0 > 1.5:
						h_pizero_ptgr15.Fill(m_pi0)

print 'Processed', i_ev, 'events'
