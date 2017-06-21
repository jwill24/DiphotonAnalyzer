#/usr/bin/python

import ROOT
from DataFormats.FWLite import Events, Handle

filename = 'root:///eos/cms/store/group/phys_pps/diphoton/DoubleEG/crab_pickEvents/170320_152646/0000/pickevents_1.root'
diphoton_vertex = ROOT.TVector3(0.0517226, 0.100587, 3.57555)

events = Events(filename)
vtx_handle = Handle('std::vector<reco::Vertex>')
vtx_label = ('offlinePrimaryVertices')
pho_handle = Handle('std::vector<reco::Photon>')
pho_label = ('photons')
met_handle = Handle('std::vector<reco::PFMET>')
met_label = ('pfMet')
ele_handle = Handle('std::vector<reco::GsfElectron>')
ele_label = ('gedGsfElectrons')
muo_handle = Handle('std::vector<reco::Muon>')
muo_label = ('muons')

for ev in events:
    ev.getByLabel(ele_label, ele_handle)
    electrons = ele_handle.product()

    ev.getByLabel(muo_label, muo_handle)
    muons = muo_handle.product()

    ev.getByLabel(vtx_label, vtx_handle)
    vertices = vtx_handle.product()
    for vtx in vertices:
        vtx_pos = ROOT.TVector3(vtx.x(), vtx.y(), vtx.z())
        if (vtx_pos-diphoton_vertex).Perp()>0.001:
            continue
        print vtx.x(), vtx.y(), vtx.z(), (vtx_pos-diphoton_vertex).Perp()
        print 'number of tracks associated to this vertex:', vtx.tracksSize()
        for i in range(vtx.tracksSize()):
            track = vtx.trackRefAt(i)
            track_kin = ROOT.TVector3(track.px(), track.py(), track.pz())
            for ele in electrons:
                ele_trk = gsfTrack()
                ele_kin = ROOT.TVector3(ele_trk.px(), ele_trk.py(), ele_trk.pz())
                print (ele_kin-track_kin).Perp()

    ev.getByLabel(pho_label, pho_handle)
    photons = pho_handle.product()
    for i in range(len(photons)):
        pho1 = photons[i]
        pho1_kin = ROOT.TLorentzVector(pho1.px(), pho1.py(), pho1.pz(), pho1.mass())
        for j in range(i+1, len(photons)):
            pho2 = photons[j]
            pho2_kin = ROOT.TLorentzVector(pho2.px(), pho2.py(), pho2.pz(), pho2.mass())
            diph_kin = (pho1_kin+pho2_kin)
            print 'diphoton', i, j, 'has mass', diph_kin.M(), 'and pT', diph_kin.Pt()    

    ev.getByLabel(met_label, met_handle)
    met = met_handle.product()[0]
    print 'MEt:', met.et()

