#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TPad.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TGraph2D.h>
#include <TH3D.h>
#include <TTUBE.h>
#include <TStyle.h>

#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "ITSMFTSimulation/Hit.h"

#include <string>
#include <array>
#include <vector>
#include <iostream>
#endif

class trackRef
{
public:
    trackRef() = default;
    trackRef(int id, int mid, int pdg, double vx, double vy, double vz) : mTrackID{id}, mMotherID{mid}, mPDG{pdg}, mVertex{vx, vy, vz} {}
    trackRef &operator=(trackRef &tr)
    {
        this->mTrackID = tr.mTrackID;
        this->mMotherID = tr.mMotherID;
        this->mPDG = tr.mPDG;
        this->mVertex = tr.mVertex;
        this->mHits = tr.mHits;

        return *this;
    }

    int mTrackID;
    int mMotherID;
    int mPDG;
    std::array<double, 3> mVertex;
    std::vector<o2::itsmft::Hit> mHits;
};

void omegaDCA(const std::string basePath = "./",
              const std::string simFileName = "o2sim_Kine.root",
              const std::string itsHitsFileName = "o2sim_HitsITS3.root")
{
    auto simFile = TFile::Open((basePath + simFileName).data(), "r");
    auto hitsFile = TFile::Open((basePath + itsHitsFileName).data(), "r");

    /// Get simulation data
    auto simTree = static_cast<TTree *>(simFile->Get("o2sim"));
    auto hitsTree = static_cast<TTree *>(hitsFile->Get("o2sim"));

    std::vector<o2::MCTrack> *MCTracks = nullptr;
    std::vector<o2::TrackReference> *TrackRefs = nullptr;
    std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

    simTree->SetBranchAddress("MCTrack", &MCTracks);
    simTree->SetBranchAddress("TrackRefs", &TrackRefs);
    hitsTree->SetBranchAddress("ITS3Hit", &ITSHits);

    /// Number of events
    auto nev = simTree->GetEntriesFast();
    std::vector<std::map<int, trackRef>> vecTrackMap;
    vecTrackMap.resize(nev);
    TFile *f = TFile::Open("ITSTracksHits.root", "recreate");
    TNtuple hitdata("hits", "hit ntuple", "ev:pdg:partid:momid:x:y:z:detectorid:vx:vy:vz");

    for (auto iEv{0}; iEv < nev; ++iEv)
    {
        simTree->GetEntry(iEv);
        hitsTree->GetEntry(iEv);

        /// Tracks loop
        for (auto iTrack{0}; iTrack < (int)MCTracks->size(); ++iTrack)
        {
            auto &track = (*MCTracks)[iTrack];
            vecTrackMap[iEv].emplace(std::make_pair(iTrack, trackRef{iTrack, track.getMotherTrackId(), track.GetPdgCode(), track.Vx(), track.Vy(), track.Vz()}));
        }

        for (auto &hit : *ITSHits)
        {
            auto it = vecTrackMap[iEv].find(hit.GetTrackID()); /// if hit -> track must exist
            it->second.mHits.push_back(hit);
            hitdata.Fill(iEv, it->second.mPDG, it->second.mTrackID, it->second.mMotherID, hit.GetX(), hit.GetY(), hit.GetZ(), hit.GetDetectorID(), it->second.mVertex[0], it->second.mVertex[1], it->second.mVertex[2]);
            Printf("\t%d %d %d %d %lf %lf %lf %d %lf %lf %lf", iEv, it->second.mPDG, it->second.mTrackID, it->second.mMotherID, hit.GetX(), hit.GetY(), hit.GetZ(), hit.GetDetectorID(), it->second.mVertex[0], it->second.mVertex[1], it->second.mVertex[2]);
        }
    }
    f->cd();
    hitdata.Write();
}