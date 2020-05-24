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
#include <vector>
#include <iostream>
#endif

struct Omega
{
    int idLambda = -1;
    int idK = -1;
    int idLambdaPi = -1;
    int idLambdaP = -1;
};

void plotDecayChain(const int pdgCode = 3334, // Omega
                    const std::string basePath = "./",
                    const std::string simFileName = "o2sim_Kine.root",
                    const std::string itsHitsFileName = "o2sim_HitsITS.root")
{

    gStyle->SetOptStat(0);
    auto simFile = TFile::Open((basePath + simFileName).data(), "r");
    auto hitsFile = TFile::Open((basePath + itsHitsFileName).data(), "r");

    auto simTree = static_cast<TTree *>(simFile->Get("o2sim"));
    auto hitsTree = static_cast<TTree *>(hitsFile->Get("o2sim"));

    std::vector<o2::MCTrack> *MCTracks = nullptr;
    std::vector<o2::TrackReference> *TrackRefs = nullptr;
    std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

    simTree->SetBranchAddress("MCTrack", &MCTracks);
    simTree->SetBranchAddress("TrackRefs", &TrackRefs);
    hitsTree->SetBranchAddress("ITSHit", &ITSHits);
    TH2D *hitsLayout = new TH2D("hl", "hl", 1500, -50, 50, 1500, -50, 50);
    std::vector<std::map<int, std::vector<o2::itsmft::Hit>>> hitmap; //every element of the vector corresponds to an event, map contains id of the track and corresponding Hit array
    std::vector<std::map<int, Omega>> omegamap;                      //every element of the vector corresponds to an event, map contains id of the omega and corresponding struct (duaghters IDs)

    const size_t nev = simTree->GetEntriesFast();
    std::vector<std::map<int, o2::MCTrack>> intTracksTotal;

    for (size_t iEvent{0}; iEvent < nev; ++iEvent)
    {
        std::cout << "Event: " << iEvent << std::endl;
        simTree->GetEntry(iEvent);
        hitsTree->GetEntry(iEvent);

        std::map<int, Omega> omegaMap; //id of the Omega and corresponding struct (daughters IDs)
        std::map<int, int> lambdaMap;  //id of the mother (Omega) and id of the lambda
        std::map<int, int> kaonMap;
        std::map<int, std::vector<int>> pipMap; //id of the mother (Lambda) and vector of id of pions and protons

        for (size_t iTrack{0}; iTrack < MCTracks->size(); ++iTrack)
        {
            auto &track = (*MCTracks)[iTrack];
            if (track.getMotherTrackId() < 0 && std::abs(track.GetPdgCode()) == pdgCode) // Primary Omega
            {
                omegaMap[iTrack] = Omega();
            }
            else
            {
                if (std::abs(track.GetPdgCode()) == 3122) // lambda
                {
                    lambdaMap[track.getMotherTrackId()] = iTrack;
                }
                if (track.getMotherTrackId() > 0 && std::abs(track.GetPdgCode()) == 321) // kaons
                {
                    kaonMap[track.getMotherTrackId()] = iTrack;
                }
                if (std::abs(track.GetPdgCode()) == 2212 || std::abs(track.GetPdgCode()) == 211) // proton and pions
                {
                    if (pipMap.find(track.getMotherTrackId()) == pipMap.end())
                    {
                        pipMap[track.getMotherTrackId()] = std::vector<int>(1, iTrack);
                    }
                    else
                    {
                        pipMap[track.getMotherTrackId()].push_back(iTrack);
                    }
                }
            }
        }
        for (auto &Omega : omegaMap)
        {
            int omega_id = Omega.first;
            int lambda_id = -1;
            int kaon_id = -1;
            if (lambdaMap.find(omega_id) != lambdaMap.end())
            {
                lambda_id = lambdaMap[omega_id];
                Omega.second.idLambda = lambda_id;
            }
            if (kaonMap.find(omega_id) != kaonMap.end())
            {
                kaon_id = kaonMap[omega_id];
                Omega.second.idK = kaon_id;
            }
            if (pipMap.find(lambda_id) != pipMap.end())
            {
                auto l_daughters = pipMap[lambda_id];
                for (size_t i{0}; i < l_daughters.size(); i++)
                {
                    int id_dau = l_daughters.at(i);
                    if (std::abs(((*MCTracks)[id_dau]).GetPdgCode()) == 211)
                    {
                        Omega.second.idLambdaPi = id_dau;
                    }
                    else if (std::abs(((*MCTracks)[id_dau]).GetPdgCode()) == 2212)
                    {
                        Omega.second.idLambdaP = id_dau;
                    }
                }
            }
            Printf("\tOmega  id = %10d", omega_id);
            Printf("\tkaon   id = %10d  mother = %10d", Omega.second.idK, ((*MCTracks)[Omega.second.idK]).getMotherTrackId());
            Printf("\tlambda id = %10d  mother = %10d", Omega.second.idLambda, ((*MCTracks)[Omega.second.idLambda]).getMotherTrackId());
            Printf("\tpion   id = %10d  mother = %10d", Omega.second.idLambdaPi, ((*MCTracks)[Omega.second.idLambdaPi]).getMotherTrackId());
            Printf("\tproton id = %10d  mother = %10d", Omega.second.idLambdaP, ((*MCTracks)[Omega.second.idLambdaP]).getMotherTrackId());
            Printf("\n");
        }
        omegamap.push_back(omegaMap);
        // hit map per event
        std::map<int, std::vector<o2::itsmft::Hit>> hits_map; // first==track id, second==vector of Hits associated

        for (size_t iHit{0}; iHit < ITSHits->size(); iHit++) // <<< loop over hits array
        {
            const auto &hit = (*ITSHits)[iHit];
            hitsLayout->Fill(hit.GetX(), hit.GetY());
            // get MC info
            int trID = hit.GetTrackID();

            bool is_ok = false;
            for (auto &omega : omegaMap)
            {
                int omega_id = omega.first;
                int omega_k = omega.second.idK;
                int omega_lambda_pi = omega.second.idLambdaPi;
                int omega_lambda_p = omega.second.idLambdaP;
                if (trID == omega_id && omega_id != -1)
                    is_ok = true;
                if (trID == omega_k && omega_k != -1)
                    is_ok = true;
                if (trID == omega_lambda_pi && omega_lambda_pi != -1)
                    is_ok = true;
                if (trID == omega_lambda_p && omega_lambda_p != -1)
                    is_ok = true;
            }
            if (is_ok)
            {
                if (hits_map.find(trID) == hits_map.end())
                    hits_map[trID] = std::vector<o2::itsmft::Hit>(1, hit);
                else
                    hits_map[trID].push_back(hit);
            }
        }

        hitmap.push_back(hits_map);
    }
    printf("\nsize of omegamap = %lu\n", omegamap.size());
    printf("size of Hit = %lu\n", hitmap.size());

    // store data
    std::vector<std::vector<int>> n_omega;
    TFile *f = TFile::Open("OmegaClusters.root", "recreate");
    TNtuple nt("ntc", "cluster ntuple", "ev:omegaid:partid:x:y:z:detectorid");

    for (size_t iev{0}; iev < hitmap.size(); ++iev)
    {
        int iomega = 0;
        std::vector<int> i_omega;
        for (auto &omega : omegamap[iev])
        {
            bool has_hits = false;
            int track_id[4] = {omega.first, omega.second.idK, omega.second.idLambdaPi, omega.second.idLambdaP};
            printf("iev = %zu    iomega: %d \t track id: %d %d %d %d\n", iev, iomega, track_id[0], track_id[1], track_id[2], track_id[3]);
            for (int itr{0}; itr < 4; itr++)
            {
                auto hit_entry = hitmap[iev].find(track_id[itr]);
                if (hit_entry != hitmap[iev].end())
                {
                    std::vector<o2::itsmft::Hit> hit_arr = hit_entry->second;
                    if (hit_arr.size() != 0)
                    {
                        has_hits = true;
                        for (size_t iHit{0}; iHit < hit_arr.size(); iHit++)
                        {
                            auto hit = hit_arr[iHit];
                            int detectID = hit.GetDetectorID();
                            printf("%5zu, %5d, %10d, %5d, %3.2f, %3.2f, %3.2f, %10d, %10d\n", iev, iomega, track_id[0], itr, hit.GetX(), hit.GetY(), hit.GetZ(), detectID, track_id[itr]);
                            nt.Fill(iev, track_id[0], itr, hit.GetX(), hit.GetY(), hit.GetZ(), detectID);
                        }
                    }
                }
            }
            if (has_hits)
            {
                i_omega.push_back(track_id[0]);
                iomega++;
            }
            n_omega.push_back(i_omega);
        }
    }
    f->cd();
    nt.Write();

    // Rendering
    TH2D *himppar_transverse = new TH2D("himppar_transverse", ";imp par rphi;n clusters", 2000, -0.5, 0.5, 10, 0, 10);
    TH2D *himppar_longitudin = new TH2D("himppar_longitudin", ";imp par z;n clusters", 2000, -0.5, 0.5, 10, 0, 10);

    TEllipse *l0 = new TEllipse(0., 0., 1.8);
    l0->SetFillStyle(0);
    l0->SetLineColor(kGray + 1);
    TEllipse *l1 = new TEllipse(0., 0., 2.4);
    l1->SetFillStyle(0);
    l1->SetLineColor(kGray + 1);
    TEllipse *l2 = new TEllipse(0., 0., 3.0);
    l2->SetFillStyle(0);
    l2->SetLineColor(kGray + 1);
    TEllipse *l3 = new TEllipse(0., 0., 7.0);
    l3->SetFillStyle(0);
    l3->SetLineColor(kGray + 1);
    TEllipse *l4 = new TEllipse(0., 0., 19.605);
    l4->SetFillStyle(0);
    l4->SetLineColor(kGray + 1);
    TEllipse *l5 = new TEllipse(0., 0., 24.545);
    l5->SetFillStyle(0);
    l5->SetLineColor(kGray + 1);
    TEllipse *l6 = new TEllipse(0., 0., 34.368);
    l6->SetFillStyle(0);
    l6->SetLineColor(kGray + 1);
    TEllipse *l7 = new TEllipse(0., 0., 39.355);
    l7->SetFillStyle(0);
    l7->SetLineColor(kGray + 1);

    TCanvas *c_ev36_omega = new TCanvas("c_ev3_omega", "c_ev36_omega", 1200, 1200);
    c_ev36_omega->SetTicks();

    TH2D *hframe = new TH2D("hframe", ";x;y", 200, -50, 50, 200, -50, 50);
    c_ev36_omega->cd();
    hframe->DrawCopy();
    hitsLayout->SetMarkerColor(kGray + 1);
    hitsLayout->Draw("same");
    // l0->Draw("same");
    // l1->Draw("same");
    // l2->Draw("same");
    // l3->Draw("same");
    // l4->Draw("same");
    // l5->Draw("same");
    // l6->Draw("same");
    // l7->Draw("same");
    nt.Draw("x:y", "ev==36 && omegaid==18563 && partid==0", "same");
    TGraph *graph_omega = (TGraph *)gPad->GetPrimitive("Graph"); // 2D
    if (graph_omega)
    {
        graph_omega->SetNameTitle(Form("graph_omega_%d_%d", 36, 0), Form("graph_omega_%d_%d", 36, 0));
        graph_omega->SetMarkerColor(kRed);
        graph_omega->SetMarkerStyle(2);
    }
    nt.Draw("x:y", "ev==36 && omegaid==18563 && partid==1", "same");
    TGraph *graph_k = (TGraph *)gPad->GetPrimitive("Graph"); // 2D
    if (graph_k)
    {
        graph_k->SetNameTitle(Form("graph_k_%d_%d", 36, 0), Form("graph_k_%d_%d", 36, 0));
        graph_k->SetMarkerColor(kBlue);
        graph_k->SetMarkerStyle(3);
    }
    nt.Draw("x:y", "ev==36 && omegaid==18563 && partid==2", "same");
    TGraph *graph_piL = (TGraph *)gPad->GetPrimitive("Graph"); // 2D
    if (graph_piL)
    {
        graph_piL->SetNameTitle(Form("graph_pi_%d_%d", 36, 0), Form("graph_pi_%d_%d", 36, 0));
        graph_piL->SetMarkerColor(kGreen);
        graph_piL->SetMarkerStyle(4);
    }
    nt.Draw("x:y", "ev==36 && omegaid==18563 && partid==3", "same");
    TGraph *graph_pL = (TGraph *)gPad->GetPrimitive("Graph"); // 2D
    if (graph_pL)
    {
        graph_pL->SetNameTitle(Form("graph_p_%d_%d", 36, 0), Form("graph_p_%d_%d", 36, 0));
        graph_pL->SetMarkerColor(kBlack);
        graph_pL->SetMarkerStyle(5);
    }
    c_ev36_omega->SaveAs("omega36_2D.png", "r");

    // float x, y, z;
    // float ev, xi, part;
    // nt.SetBranchAddress("ev", &ev);
    // nt.SetBranchAddress("oemgaid", &xi);
    // nt.SetBranchAddress("partid", &part);
    // nt.SetBranchAddress("x", &x);
    // nt.SetBranchAddress("y", &y);
    // nt.SetBranchAddress("z", &z);
    // TH3D *hframe3d = new TH3D("hframe3d", ";z;x;y", 200, -50, 50, 200, -50, 50, 200, -50, 50);
    // TCanvas *c_ev3 = new TCanvas("c_ev36", "c_ev36", 1200, 1200);
    // hframe3d->DrawCopy();
    // c_ev3->SetTicks();

    // TGraph2D *omega_3d = new TGraph2D();
    // omega_3d->SetNameTitle(Form("xi_%d_%d", 36, 0), Form("xi_%d_%d", 36, 0));
    // omega_3d->SetMarkerColor(kRed);
    // omega_3d->SetMarkerStyle(2);
    // TGraph2D *k_3d = new TGraph2D();
    // k_3d->SetNameTitle(Form("pi_%d_%d", 36, 0), Form("pi_%d_%d", 36, 0));
    // k_3d->SetMarkerColor(kBlue);
    // k_3d->SetMarkerStyle(3);
    // TGraph2D *piL_3d = new TGraph2D();
    // piL_3d->SetNameTitle(Form("piL_%d_%d", 36, 0), Form("piL_%d_%d", 36, 0));
    // piL_3d->SetMarkerColor(kGreen + 1);
    // piL_3d->SetMarkerStyle(4);
    // TGraph2D *pL_3d = new TGraph2D();
    // pL_3d->SetNameTitle(Form("pL_%d_%d", 36, 0), Form("pL_%d_%d", 36, 0));
    // pL_3d->SetMarkerColor(kBlack);
    // pL_3d->SetMarkerStyle(5);

    // // 3D plot
    // TPolyLine3D *omega_3d_l = new TPolyLine3D;
    // omega_3d_l->SetLineColor(kRed);
    // TPolyLine3D *k_3d_l = new TPolyLine3D;
    // k_3d_l->SetLineColor(kBlue);
    // TPolyLine3D *piL_3d_l = new TPolyLine3D;
    // piL_3d_l->SetLineColor(kGreen + 1);
    // TPolyLine3D *pL_3d_l = new TPolyLine3D;
    // pL_3d_l->SetLineColor(kBlack);
    // int k_0 = 0;
    // int k_1 = 0;
    // int k_2 = 0;
    // int k_3 = 0;
    // for (int i{0}; i < nt.GetEntries(); i++)
    // {
    //     nt.GetEntry(i);
    //     if (ev == 36)
    //     {
    //         if (part == 0)
    //         {
    //             omega_3d->SetPoint(k_0, z, x, y);
    //             omega_3d_l->SetNextPoint(z, x, y);
    //             k_0++;
    //         }
    //         else if (part == 1)
    //         {
    //             k_3d->SetPoint(k_1, z, x, y);
    //             k_3d_l->SetNextPoint(z, x, y);
    //             k_1++;
    //         }
    //         else if (part == 2)
    //         {
    //             piL_3d->SetPoint(k_2, z, x, y);
    //             piL_3d_l->SetNextPoint(z, x, y);
    //             k_2++;
    //         }
    //         else if (part == 3)
    //         {
    //             pL_3d->SetPoint(k_3, z, x, y);
    //             pL_3d_l->SetNextPoint(z, x, y);
    //             k_3++;
    //         }
    //     }
    // }
    // int n_xi = omega_3d->GetN();
    // int n_pi = k_3d->GetN();
    // int n_piL = piL_3d->GetN();
    // int n_pL = pL_3d->GetN();
    // if (n_xi > 0)
    // {
    //     omega_3d->SetMarkerColor(kRed);
    //     omega_3d->SetMarkerStyle(2);
    //     omega_3d->Draw("psame");
    //     omega_3d_l->Draw("same");
    // }
    // if (n_pi > 0)
    // {
    //     k_3d->SetMarkerColor(kBlue);
    //     k_3d->SetMarkerStyle(3);
    //     k_3d->Draw("psame");
    //     k_3d_l->Draw("same");
    // }
    // if (n_piL > 0)
    // {
    //     piL_3d->SetMarkerColor(kGreen + 1);
    //     piL_3d->SetMarkerStyle(4);
    //     piL_3d->Draw("psame");
    //     piL_3d_l->Draw("same");
    // }
    // if (n_pL > 0)
    // {
    //     pL_3d->SetMarkerColor(kBlack);
    //     pL_3d->SetMarkerStyle(5);
    //     pL_3d->Draw("psame");
    //     pL_3d_l->Draw("same");
    // }
    // TTUBE *its0_3d = new TTUBE("l0", "Layer 0", "Silicon", 1.8, 1.8, 15.f);
    // TTUBE *its1_3d = new TTUBE("l1", "Layer 1", "Silicon", 2.4, 2.4, 15.f);
    // TTUBE *its2_3d = new TTUBE("l2", "Layer 2", "Silicon", 3.0, 3.0, 15.f);
    // TTUBE *its3_3d = new TTUBE("l3", "Layer 3", "Silicon", 7.0, 7.0, 15.f);
    // TTUBE *its4_3d = new TTUBE("l4", "Layer 4", "Silicon", 19.605, 19.605, 14.f);
    // TTUBE *its5_3d = new TTUBE("l5", "Layer 5", "Silicon", 24.545, 24.545, 14.f);
    // TTUBE *its6_3d = new TTUBE("l6", "Layer 6", "Silicon", 34.368, 34.368, 75.f);
    // TTUBE *its7_3d = new TTUBE("l7", "Layer 7", "Silicon", 39.355, 39.355, 75.f);
    // its0_3d->Draw("psame");
    // its1_3d->Draw("psame");
    // its2_3d->Draw("psame");
    // its3_3d->Draw("psame");
    // its4_3d->Draw("psame");
    // its5_3d->Draw("psame");
    // its6_3d->Draw("psame");
    // its7_3d->Draw("psame");
}
