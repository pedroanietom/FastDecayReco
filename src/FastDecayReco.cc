#include "FastDecayReco.h"

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packetv1.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <trackreco/ActsPropagator.h>
#include <Acts/Surfaces/CylinderSurface.hpp>

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTowerv2.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoDefs.h>

#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Acts/Surfaces/CylinderSurface.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <TFile.h>
#include <TH1.h>
#include "TF1.h"
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cmath>
#include <utility>

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;


int FastDecayReco::process_event(PHCompositeNode* topNode)
{

  EventHeader* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  int m_runNumber;
  int m_evtNumber;
  if (evtHeader)
  {
    m_runNumber = evtHeader->get_RunNumber();
    m_evtNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    m_runNumber = m_evtNumber = -1;
  }

  // Global Vertex information
  double Zvtx = -999.;

  for (GlobalVertexMap::Iter iter = m_global_vtxmap->begin(); iter != m_global_vtxmap->end(); ++iter)  {
    GlobalVertex *vtx = iter->second;
    //std::cout << "Vertex ID: " << vtx->get_id() << ", Z: " << vtx->get_z() << std::endl;
    Zvtx = vtx->get_z();
  }
  if(Verbosity() > 0){
    std::cout << "Global vertex Z = " << Zvtx << " " << std::endl;
  }
  if(fabs(Zvtx)>30.){
    std::cout << "No vertex" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Trigger information
  std::vector<int> _triggers;
  bool isMB10 = false;
  bool isMB12 = false;
  bool isPhoton4 = false;
  bool isPhoton5 = false;
  bool isMBPhoton3 = false;
  bool isMBPhoton4 = false;
  bool isMBPhoton5 = false;

  if(gl1Packet)
  {
    auto scaled_vector = gl1Packet->getScaledVector();
    for(int i = 0; i < 40; i++)
    {
      if((scaled_vector & (int)std::pow(2,i)) != 0)
      {
        _triggers.push_back(i);
      }
    }
  }
  else
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (int trigger : _triggers)
  {
    if (trigger == 10)
    {
      isMB10 = true;
      h_trigger->Fill(10);
    }
    else if(trigger == 12)
    {
      isMB12 = true;
      h_trigger->Fill(12);
    }
    else if (trigger == 30)
    {
      isPhoton4 = true;
      h_trigger->Fill(30);
    }
    else if (trigger == 31)
    {
      isPhoton5 = true;
      h_trigger->Fill(31);
    }
    else if (trigger == 36)
    {
      isMBPhoton3 = true;
      h_trigger->Fill(36);
    }
    else if (trigger == 37)
    {
      isMBPhoton4 = true;
      h_trigger->Fill(37);
    }
    else if (trigger == 38)
    {
      isMBPhoton5 = true;
      h_trigger->Fill(38);
    }
  }

  // Loop over tracks and check for close DCA match with all other tracks
  for (auto tr1_it = m_svtxTrackMap->begin(); tr1_it != m_svtxTrackMap->end(); ++tr1_it)
  {
    auto id1 = tr1_it->first;
    auto *tr1 = tr1_it->second;
    if (tr1->get_quality() > _qual_cut)
    {
      continue;
    }
    if (tr1->get_pt() < track_pt_cut)
    {
      continue;
    }

    short int crossing1 = tr1->get_crossing();

    // calculate number silicon tracks
    double this_dca_cut = track_dca_cut;
    TrackSeed* siliconseed = tr1->get_silicon_seed();

    if (!siliconseed)
    {
      this_dca_cut *= 5;
      if (Verbosity() > 2)
      {
        std::cout << "silicon seed not found" << std::endl;
      }
      if (_require_mvtx)
      {
        continue;
      }
    }

    std::vector<unsigned int> nstates1 = getTrackStates(tr1);
    unsigned int track1_mvtx_state_size = nstates1[0];
    unsigned int track1_intt_state_size = nstates1[1];
    // unsigned int track1_tpc_state_size = nstates1[2];
    // unsigned int track1_mms_state_size = nstates1[3];

    unsigned int track1_silicon_cluster_size = std::numeric_limits<unsigned int>::quiet_NaN();
    if (siliconseed)
    {
      track1_silicon_cluster_size = siliconseed->size_cluster_keys();
    }

    std::vector<TrkrDefs::cluskey> ckeys1;
    if (siliconseed)
    {
      ckeys1.insert(ckeys1.end(), siliconseed->begin_cluster_keys(), siliconseed->end_cluster_keys());
    }

    unsigned int track1_mvtx_cluster_size = 0;
    unsigned int track1_intt_cluster_size = 0;
    for (const auto& ckey : ckeys1)
    {
      auto detid = TrkrDefs::getTrkrId(ckey);
      if (detid == TrkrDefs::TrkrId::mvtxId)
      {
        track1_mvtx_cluster_size++;
      }
      else if (detid == TrkrDefs::TrkrId::inttId)
      {
        track1_intt_cluster_size++;
      }
    }

    Acts::Vector3 pos1(tr1->get_x(), tr1->get_y(), tr1->get_z());
    Acts::Vector3 mom1(tr1->get_px(), tr1->get_py(), tr1->get_pz());
    Acts::Vector3 dcaVals1 = calculateDca(tr1, mom1, pos1);
    // first dca cuts
    if (fabs(dcaVals1(0)) < this_dca_cut || fabs(dcaVals1(1)) < this_dca_cut)
    {
      continue;
    }

    // look for close DCA matches with all other such tracks
    for (auto tr2_it = std::next(tr1_it); tr2_it != m_svtxTrackMap->end(); ++tr2_it)
    {
      auto id2 = tr2_it->first;
      auto *tr2 = tr2_it->second;
      if (tr2->get_quality() > _qual_cut)
      {
        continue;
      }
      if (tr2->get_pt() < track_pt_cut)
      {
        continue;
      }

      short int crossing2 = tr2->get_crossing();

      // calculate number silicon tracks
      TrackSeed* siliconseed2 = tr2->get_silicon_seed();

      double this_dca_cut2 = track_dca_cut;

      if (!siliconseed2)
      {
        this_dca_cut2 *= 5;
        if (Verbosity() > 2)
        {
          std::cout << "silicon seed not found" << std::endl;
        }
        if (_require_mvtx)
        {
          continue;
        }
      }

      std::vector<unsigned int> nstates2 = getTrackStates(tr2);
      unsigned int track2_mvtx_state_size = nstates2[0];
      unsigned int track2_intt_state_size = nstates2[1];
      // unsigned int track2_tpc_state_size = nstates2[2];
      // unsigned int track2_mms_state_size = nstates2[3];

      unsigned int track2_silicon_cluster_size = std::numeric_limits<unsigned int>::quiet_NaN();
      if (siliconseed2)
      {
        track2_silicon_cluster_size = siliconseed2->size_cluster_keys();
      }

      std::vector<TrkrDefs::cluskey> ckeys2;
      if (siliconseed2)
      {
        ckeys2.insert(ckeys2.end(), siliconseed2->begin_cluster_keys(), siliconseed2->end_cluster_keys());
      }

      unsigned int track2_mvtx_cluster_size = 0;
      unsigned int track2_intt_cluster_size = 0;
      for (const auto& ckey : ckeys2)
      {
        auto detid = TrkrDefs::getTrkrId(ckey);
        if (detid == TrkrDefs::TrkrId::mvtxId)
        {
          track2_mvtx_cluster_size++;
        }
        else if (detid == TrkrDefs::TrkrId::inttId)
        {
          track2_intt_cluster_size++;
        }
      }

      // dca xy and dca z cut here compare to track dca cut
      Acts::Vector3 pos2(tr2->get_x(), tr2->get_y(), tr2->get_z());
      Acts::Vector3 mom2(tr2->get_px(), tr2->get_py(), tr2->get_pz());
      Acts::Vector3 dcaVals2 = calculateDca(tr2, mom2, pos2);

      if (fabs(dcaVals2(0)) < this_dca_cut2 || fabs(dcaVals2(1)) < this_dca_cut2)
      {
        continue;
      }

      // find DCA of these two tracks
      if (Verbosity() > 3)
      {
        std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
      }

      if (tr1->get_charge() == tr2->get_charge())
      {
        // continue;
      }

      // declare these variables to pass into findPCAtwoTracks and fillHistogram by reference
      double pair_dca;
      Acts::Vector3 pca_rel1;
      Acts::Vector3 pca_rel2;
      double invariantMass;
      double invariantPt;
      float invariantPhi;
      float rapidity;
      float pseudorapidity;
      float decaymassA;
      float decaymassB;

      // Matching calo params
      float track1_cemc_eta;
      float track1_cemc_phi;
      float cemc_eta1;
      float cemc_z1;
      float cemc_phi1;
      float cemc_ecore1;

      float track2_cemc_eta;
      float track2_cemc_phi;
      float cemc_eta2;
      float cemc_z2;
      float cemc_phi2;
      float cemc_ecore2;

      float track1_ihcal_eta;
      float track1_ihcal_phi;
      float ihcal_eta1;
      float ihcal_phi1;
      float ihcal_e3x3_1;

      float track2_ihcal_eta;
      float track2_ihcal_phi;
      float ihcal_eta2;
      float ihcal_phi2;
      float ihcal_e3x3_2;

      float track1_ohcal_eta;
      float track1_ohcal_phi;
      float ohcal_eta1;
      float ohcal_phi1;
      float ohcal_e3x3_1;

      float track2_ohcal_eta;
      float track2_ohcal_phi;
      float ohcal_eta2;
      float ohcal_phi2;
      float ohcal_e3x3_2;

      std::array<std::pair<float, float>, 2> mass_combos = {
        std::make_pair(decaymass1, decaymass2),
        std::make_pair(decaymass2, decaymass1)
      };
      int ncombinations = (decaymass1 == decaymass2) ? 1 : 2; // if decaymass1 == decaymass2 -> (pi pi),  if decaymass1 != decaymass2 -> (pi K, K pi)

      for(int icomb = 0; icomb<ncombinations; icomb++){

        decaymassA = mass_combos[icomb].first;
        decaymassB = mass_combos[icomb].second;

        // Initial calculation of point of closest approach between the two tracks
        // This presently assumes straight line tracks to get a rough answer
        // Should update to use circles instead?
        findPcaTwoTracks(decaymassA, decaymassB, pos1, pos2, mom1, mom2, pca_rel1, pca_rel2, pair_dca);

        // tracks with small relative pca are k short candidates
        if (abs(pair_dca) < pair_dca_cut)
        {
          // Pair pca and dca were calculated with nominal track parameters and are approximate
          // Project tracks to this rough pca
          Eigen::Vector3d projected_pos1;
          Eigen::Vector3d projected_mom1;
          Eigen::Vector3d projected_pos2;
          Eigen::Vector3d projected_mom2;

          bool ret1 = projectTrackToPoint(tr1, pca_rel1, projected_pos1, projected_mom1);
          bool ret2 = projectTrackToPoint(tr2, pca_rel2, projected_pos2, projected_mom2);

          if (!ret1 || !ret2)
          {
            continue;
          }

          // recalculate pca starting with projected position and momentum
          double pair_dca_proj;
          Acts::Vector3 pca_rel1_proj;
          Acts::Vector3 pca_rel2_proj;
          findPcaTwoTracks(decaymassA, decaymassB, projected_pos1, projected_pos2, projected_mom1, projected_mom2, pca_rel1_proj, pca_rel2_proj, pair_dca_proj);

          // if(pair_dca_proj > pair_dca_cut) continue;
          // matching with calos is calculated here
          match_calos(tr1, tr2, track1_cemc_phi, track1_cemc_eta, cemc_phi1, cemc_eta1, cemc_z1, cemc_ecore1, track2_cemc_phi, track2_cemc_eta, cemc_phi2, cemc_eta2, cemc_z2, cemc_ecore2, track1_ihcal_phi, track1_ihcal_eta, ihcal_phi1, ihcal_eta1, ihcal_e3x3_1, track2_ihcal_phi, track2_ihcal_eta, ihcal_phi2, ihcal_eta2, ihcal_e3x3_2, track1_ohcal_phi, track1_ohcal_eta, ohcal_phi1, ohcal_eta1, ohcal_e3x3_1, track2_ohcal_phi, track2_ohcal_eta, ohcal_phi2, ohcal_eta2, ohcal_e3x3_2, Zvtx);

          // invariant mass is calculated in this method
          fillHistogram(projected_mom1, projected_mom2, decaymassA, decaymassB, recomass, invariantMass, invariantPt, invariantPhi, rapidity, pseudorapidity);

          if(crossing1==0 && crossing2==0){ // To do the track-calo matching

            if (_photonconv && invariantMass >= 0.06) continue; // For the Photon Conversion option

            fillNtp(tr1, tr2, decaymassA, decaymassB, dcaVals1, dcaVals2, pca_rel1, pca_rel2, pair_dca, invariantMass, invariantPt, invariantPhi, rapidity, pseudorapidity, projected_pos1, projected_pos2, projected_mom1, projected_mom2, pca_rel1_proj, pca_rel2_proj, pair_dca_proj, track1_silicon_cluster_size, track2_silicon_cluster_size, track1_mvtx_cluster_size, track1_mvtx_state_size, track1_intt_cluster_size, track1_intt_state_size, track2_mvtx_cluster_size, track2_mvtx_state_size, track2_intt_cluster_size, track2_intt_state_size, icomb, track1_cemc_phi, track1_cemc_eta, cemc_phi1, cemc_eta1, cemc_z1, cemc_ecore1, track2_cemc_phi, track2_cemc_eta, cemc_phi2, cemc_eta2, cemc_z2, cemc_ecore2, track1_ihcal_phi, track1_ihcal_eta, ihcal_phi1, ihcal_eta1, ihcal_e3x3_1, track2_ihcal_phi, track2_ihcal_eta, ihcal_phi2, ihcal_eta2, ihcal_e3x3_2, track1_ohcal_phi, track1_ohcal_eta, ohcal_phi1, ohcal_eta1, ohcal_e3x3_1, track2_ohcal_phi, track2_ohcal_eta, ohcal_phi2, ohcal_eta2, ohcal_e3x3_2, Zvtx, m_runNumber, m_evtNumber, isMB10, isMB12, isPhoton4, isPhoton5, isMBPhoton3, isMBPhoton4, isMBPhoton5);

          }

          if (Verbosity() > 1)
          {
            std::cout << " Accepted Track Pair" << std::endl;
            std::cout << " id1 " << id1 << " id2 " << id2 << std::endl;
            std::cout << " crossing1 " << crossing1 << " crossing2 " << crossing2 << std::endl;
            std::cout << " invariant mass: " << invariantMass << std::endl;
            std::cout << " track1 dca_cut: " << this_dca_cut << " track2 dca_cut: " << this_dca_cut2 << std::endl;
            std::cout << " dca3dxy1,dca3dz1,phi1: " << dcaVals1 << std::endl;
            std::cout << " dca3dxy2,dca3dz2,phi2: " << dcaVals2 << std::endl;
            std::cout << "Initial:  pca_rel1: " << pca_rel1 << " pca_rel2: " << pca_rel2 << std::endl;
            std::cout << " Initial: mom1: " << mom1 << " mom2: " << mom2 << std::endl;
            std::cout << "Proj_pca_rel:  proj_pos1: " << projected_pos1 << " proj_pos2: " << projected_pos2 << " proj_mom1: " << projected_mom1 << " proj_mom2: " << projected_mom2 << std::endl;
            std::cout << " Relative PCA = " << abs(pair_dca) << " pca_cut = " << pair_dca_cut << std::endl;
            std::cout << " charge 1: " << tr1->get_charge() << " charge2: " << tr2->get_charge() << std::endl;
            std::cout << "found viable projection" << std::endl;
            std::cout << "Final: pca_rel1_proj: " << pca_rel1_proj << " pca_rel2_proj: " << pca_rel2_proj << " mom1: " << projected_mom1 << " mom2: " << projected_mom2 << std::endl;
            std::cout << "Combination: "  << icomb << " decaymassA = " << decaymassA << " decaymassB " << decaymassB << std::endl
              << std::endl;
          }

          if (m_save_tracks)
          {
            m_output_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_output_trackMap_node_name);
            m_output_trackMap->insertWithKey(tr1, tr1->get_id());
            m_output_trackMap->insertWithKey(tr2, tr2->get_id());
          }

        } // old loop

      }
    }
  }
  return 0;
}

std::vector<unsigned int> FastDecayReco::getTrackStates(SvtxTrack *track)
{
  std::vector<unsigned int> nstates;
  unsigned int nmapsstate = 0;
  unsigned int ninttstate = 0;
  unsigned int ntpcstate = 0;
  unsigned int nmmsstate = 0;

  // the track states from the Acts fit are fitted to fully corrected clusters, and are on the surface
  for (auto state_iter = track->begin_states();
      state_iter != track->end_states();
      ++state_iter)
  {
    SvtxTrackState* tstate = state_iter->second;

    if (tstate->get_pathlength() != 0)  // The first track state is an extrapolation so has no cluster
    {

      auto stateckey = tstate->get_cluskey();
      uint8_t id = TrkrDefs::getTrkrId(stateckey);

      switch (id)
      {
        case TrkrDefs::mvtxId:
          nmapsstate++;
          break;
        case TrkrDefs::inttId:
          ninttstate++;
          break;
        case TrkrDefs::tpcId:
          ntpcstate++;
          break;
        case TrkrDefs::micromegasId:
          nmmsstate++;
          break;
        default:
          //std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
          break;

          /*
             default:
             std::cout << PHWHERE << " unknown key " << stateckey << std::endl;
             gSystem->Exit(1);
             exit(1);
             */

      }
    }
  }
  nstates.push_back(nmapsstate);
  nstates.push_back(ninttstate);
  nstates.push_back(ntpcstate);
  nstates.push_back(nmmsstate);

  return nstates;
}

void FastDecayReco::fillNtp(SvtxTrack* track1, SvtxTrack* track2, float decaymassA, float decaymassB, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float invariantPhi, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2, Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj, unsigned int track1_silicon_cluster_size, unsigned int track2_silicon_cluster_size, unsigned int track1_mvtx_cluster_size,  unsigned int track1_mvtx_state_size, unsigned int track1_intt_cluster_size,  unsigned int track1_intt_state_size, unsigned int track2_mvtx_cluster_size,  unsigned int track2_mvtx_state_size, unsigned int track2_intt_cluster_size,  unsigned int track2_intt_state_size, int icomb, float track1_cemc_phi, float track1_cemc_eta, float cemc_phi1, float cemc_eta1, float cemc_z1, float cemc_ecore1, float track2_cemc_phi, float track2_cemc_eta, float cemc_phi2, float cemc_eta2, float cemc_z2, float cemc_ecore2, float track1_ihcal_phi, float track1_ihcal_eta, float ihcal_phi1, float ihcal_eta1, float ihcal_e3x3_1, float track2_ihcal_phi, float track2_ihcal_eta, float ihcal_phi2, float ihcal_eta2, float ihcal_e3x3_2, float track1_ohcal_phi, float track1_ohcal_eta, float ohcal_phi1, float ohcal_eta1, float ohcal_e3x3_1, float track2_ohcal_phi, float track2_ohcal_eta, float ohcal_phi2, float ohcal_eta2, float ohcal_e3x3_2,  double Zvtx,  int runNumber, int eventNumber, bool isMB10, bool isMB12, bool isPhoton4, bool isPhoton5, bool isMBPhoton3, bool isMBPhoton4, bool isMBPhoton5)
{
  double px1 = track1->get_px();
  double py1 = track1->get_py();
  double pz1 = track1->get_pz();
  auto *tpcSeed1 = track1->get_tpc_seed();
  size_t tpcClusters1 = tpcSeed1->size_cluster_keys();
  double eta1 = asinh(pz1 / sqrt(pow(px1, 2) + pow(py1, 2)));

  double px2 = track2->get_px();
  double py2 = track2->get_py();
  double pz2 = track2->get_pz();
  auto *tpcSeed2 = track2->get_tpc_seed();
  size_t tpcClusters2 = tpcSeed2->size_cluster_keys();
  double eta2 = asinh(pz2 / sqrt(pow(px2, 2) + pow(py2, 2)));

  double trackid1 = track1->get_id();
  double trackid2 = track2->get_id();

  // dEdx calculation
  float calculated_dEdx1 = get_dEdx(track1);
  float calculated_dEdx2 = get_dEdx(track2);

  int PID1 = get_PID(decaymassA, track1->get_charge());
  int PID2 = get_PID(decaymassB, track2->get_charge());

  float qmomentum1 = (track1->get_charge())*sqrt(px1*px1 + py1*py1 + pz1*pz1);
  float qmomentum2 = (track2->get_charge())*sqrt(px2*px2 + py2*py2 + pz2*pz2);

  double expected_dEdx1 = get_dEdx_fitValue(qmomentum1, PID1);
  double expected_dEdx2 = get_dEdx_fitValue(qmomentum2, PID2);

  double matching1 = fabs(calculated_dEdx1 - expected_dEdx1)/expected_dEdx1;
  double matching2 = fabs(calculated_dEdx2 - expected_dEdx2)/expected_dEdx2;

  if(Verbosity()>0){
    std::cout << "calculated_dEdx1: " << calculated_dEdx1 << ", PID1: " << PID1 << ", qmomentum1: " << qmomentum1 << ", expected_dEdx1: " << expected_dEdx1 << ", matching1: " << matching1 << std::endl;
  }

  auto vtxid = track1->get_vertex_id();

  Acts::Vector3 vertex(0, 0, track1->get_z());  // fake primary vertex
  auto *svtxVertex = m_vertexMap->get(vtxid);
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x();
    vertex(1) =  svtxVertex->get_y();
    vertex(2) = svtxVertex->get_z(); 
  }

  Acts::Vector3 pathLength = (pca_rel1 + pca_rel2) * 0.5 - vertex;
  Acts::Vector3 pathLength_proj = (pca_rel1_proj + pca_rel2_proj) * 0.5 - vertex;

  float mag_pathLength = sqrt(pow(pathLength(0), 2) + pow(pathLength(1), 2) + pow(pathLength(2), 2));
  float mag_pathLength_proj = sqrt(pow(pathLength_proj(0), 2) + pow(pathLength_proj(1), 2) + pow(pathLength_proj(2), 2));

  Acts::Vector3 projected_momentum = projected_mom1 + projected_mom2;
  float cos_theta_reco = pathLength_proj.dot(projected_momentum) / (projected_momentum.norm() * pathLength_proj.norm());

  float mb10, mb12, photon4, photon5, mbphoton3, mbphoton4, mbphoton5;

  if(isMB10){mb10 = 1;}else{mb10 = 0;}
  if(isMB12){mb12 = 1;}else{mb12 = 0;}
  if(isPhoton4){photon4 = 1;}else{photon4 = 0;}
  if(isPhoton5){photon5 = 1;}else{photon5 = 0;}
  if(isMBPhoton3){mbphoton3 = 1;}else{mbphoton3 = 0;}
  if(isMBPhoton4){mbphoton4 = 1;}else{mbphoton4 = 0;}
  if(isMBPhoton5){mbphoton5 = 1;}else{mbphoton5 = 0;}

  if (Verbosity() > 0)
  {
    std::cout << " Accepted Track Pair" << std::endl;
    std::cout << " Combination: " << icomb << std::endl;
    std::cout << " decaymassA: " << decaymassA << " decaymassB: " << decaymassB << std::endl;
    std::cout << " qmomentum1: " << qmomentum1 << " qmomentum2: " << qmomentum2 << std::endl;
    std::cout << " PID1: " << PID1 << " PID2: " << PID2 << std::endl;
    std::cout << " calculated_dEdx1: " << calculated_dEdx1 << " calculated_dEdx2: " << calculated_dEdx2 << std::endl;
    std::cout << " expected_dEdx1: " << expected_dEdx1 << " expected_dEdx2: " << expected_dEdx2 << std::endl;
    std::cout << " matching1: " << matching1 << " matching2: " << matching2 << std::endl
      << std::endl;
  }

  float reco_info[] = {(float) track1->get_id(), (float) track1->get_crossing(), track1->get_x(), track1->get_y(), track1->get_z(), track1->get_px(), track1->get_py(), track1->get_pz(), (float) decaymassA, (float) calculated_dEdx1, (float) matching1, (float) dcavals1(0), (float) dcavals1(1), (float) dcavals1(2), (float) pca_rel1(0), (float) pca_rel1(1), (float) pca_rel1(2), (float) eta1, (float) track1->get_charge(), (float) tpcClusters1, (float) track2->get_id(), (float) track2->get_crossing(), track2->get_x(), track2->get_y(), track2->get_z(), track2->get_px(), track2->get_py(), track2->get_pz(), (float) decaymassB, (float) calculated_dEdx2, (float) matching2, (float) dcavals2(0), (float) dcavals2(1), (float) dcavals2(2), (float) pca_rel2(0), (float) pca_rel2(1), (float) pca_rel2(2), (float) eta2, (float) track2->get_charge(), (float) tpcClusters2, (float) vertex(0), (float) vertex(1), (float) vertex(2), (float) pair_dca, (float) invariantMass, (float) invariantPt, invariantPhi, (float) pathLength(0), (float) pathLength(1), (float) pathLength(2), mag_pathLength, rapidity, pseudorapidity, (float) projected_pos1(0), (float) projected_pos1(1), (float) projected_pos1(2), (float) projected_pos2(0), (float) projected_pos2(1), (float) projected_pos2(2), (float) projected_mom1(0), (float) projected_mom1(1), (float) projected_mom1(2), (float) projected_mom2(0), (float) projected_mom2(1), (float) projected_mom2(2), (float) pca_rel1_proj(0), (float) pca_rel1_proj(1), (float) pca_rel1_proj(2), (float) pca_rel2_proj(0), (float) pca_rel2_proj(1), (float) pca_rel2_proj(2), (float) pair_dca_proj, (float) pathLength_proj(0), (float) pathLength_proj(1), (float) pathLength_proj(2), mag_pathLength_proj, track1->get_quality(), track2->get_quality(), cos_theta_reco, (float) track1_silicon_cluster_size, (float) track2_silicon_cluster_size, (float) track1_mvtx_cluster_size, (float) track1_mvtx_state_size, (float) track1_intt_cluster_size,  (float) track1_intt_state_size, (float) track2_mvtx_cluster_size, (float) track2_mvtx_state_size,  (float) track2_intt_cluster_size, (float) track2_intt_state_size, (float) icomb, (float) track1_cemc_phi, (float) track1_cemc_eta, (float) cemc_phi1, (float) cemc_eta1, (float) cemc_z1, (float) cemc_ecore1, (float) track2_cemc_phi, (float) track2_cemc_eta, (float) cemc_phi2, (float) cemc_eta2, (float) cemc_z2, (float) cemc_ecore2, (float) track1_ihcal_phi, (float) track1_ihcal_eta, (float) ihcal_phi1, (float) ihcal_eta1, (float) ihcal_e3x3_1, (float) track2_ihcal_phi, (float) track2_ihcal_eta, (float) ihcal_phi2, (float) ihcal_eta2, (float) ihcal_e3x3_2, (float) track1_ohcal_phi, (float) track1_ohcal_eta, (float) ohcal_phi1, (float) ohcal_eta1, (float) ohcal_e3x3_1, (float) track2_ohcal_phi, (float) track2_ohcal_eta, (float) ohcal_phi2, (float) ohcal_eta2, (float) ohcal_e3x3_2, (float) Zvtx, (float) runNumber, (float) eventNumber, (float) trackid1, (float) trackid2, (float) _runnumber, (float) _segment, (float) mb10, (float) mb12, (float) photon4, (float) photon5, (float) mbphoton3, (float) mbphoton4, (float) mbphoton5};



  ntp_reco_info->Fill(reco_info);
}

void FastDecayReco::match_calos(SvtxTrack* track1, SvtxTrack* track2, float& track1_cemc_phi, float& track1_cemc_eta, float& cemc_phi1, float& cemc_eta1, float& cemc_z1, float& cemc_ecore1, float& track2_cemc_phi, float& track2_cemc_eta, float& cemc_phi2, float& cemc_eta2, float& cemc_z2, float& cemc_ecore2, float& track1_ihcal_phi, float& track1_ihcal_eta, float& ihcal_phi1, float& ihcal_eta1, float& ihcal_e3x3_1, float& track2_ihcal_phi, float& track2_ihcal_eta, float& ihcal_phi2, float& ihcal_eta2, float& ihcal_e3x3_2, float& track1_ohcal_phi, float& track1_ohcal_eta, float& ohcal_phi1, float& ohcal_eta1, float& ohcal_e3x3_1, float& track2_ohcal_phi, float& track2_ohcal_eta, float& ohcal_phi2, float& ohcal_eta2, float& ohcal_e3x3_2, double Zvtx)
{

  std::vector<SvtxTrackState*> states1;
  std::vector<SvtxTrackState*> states2;

  std::vector<double> proj1;
  std::vector<double> proj2;

  for (SvtxTrack::StateIter stateiter = track1->begin_states(); stateiter != track1->end_states(); ++stateiter)
  {
    SvtxTrackState *trackstate = stateiter->second;
    if(trackstate) { proj1.push_back(trackstate->get_pathlength()); }
  }
  for (SvtxTrack::StateIter stateiter = track2->begin_states(); stateiter != track2->end_states(); ++stateiter)
  {
    SvtxTrackState *trackstate = stateiter->second;
    if(trackstate) { proj2.push_back(trackstate->get_pathlength()); }
  }

  double cemc_pathlength1 = 102.9; // from DetailedCalorimeterGeometry, project to inner surface
  double cemc_pathlength2 = 102.9;
  double ihcal_pathlength1 = proj1[proj1.size()-4]; // HCALIN is next to last
  double ihcal_pathlength2 = proj2[proj2.size()-4];
  double ohcal_pathlength1 = proj1[proj1.size()-2]; // HCALOUT is the last
  double ohcal_pathlength2 = proj2[proj2.size()-2];

  if(Verbosity() > 0){

    std::cout << "projections pathlength in track1: " << cemc_pathlength1 << ", projections pathlength in track2: " << cemc_pathlength2 << std::endl;

  }

  // track1
  SvtxTrackState *cemc_trackstate1 = track1->get_state(cemc_pathlength1);
  // track2
  SvtxTrackState *cemc_trackstate2 = track2->get_state(cemc_pathlength2);

  if (!cemc_trackstate1 || !cemc_trackstate2) {
    if (Verbosity() > 0) std::cout << "No CEMC state found → skip this pair\n";
    return;  // exit function early
  }

  // Find the closest EMCal cluster to track1 and track2


  RawCluster* clus1 = MatchClusterCEMC(cemc_trackstate1, m_cemc_clusters, track1_cemc_eta, track1_cemc_phi, Zvtx);
  RawCluster* clus2 = MatchClusterCEMC(cemc_trackstate2, m_cemc_clusters, track2_cemc_eta, track2_cemc_phi, Zvtx);

  double cemc_x1 = clus1->get_x();
  double cemc_y1 = clus1->get_y();
  cemc_z1 = clus1->get_z() - Zvtx; // correct for event vertex position
  double cemc_r1 = sqrt(pow(cemc_x1,2)+pow(cemc_y1,2));
  cemc_eta1 = asinh( cemc_z1 / cemc_r1 );
  cemc_phi1 = atan2( cemc_y1, cemc_x1 );
  cemc_ecore1 = clus1->get_ecore();

  double cemc_x2 = clus2->get_x();
  double cemc_y2 = clus2->get_y();
  cemc_z2 = clus2->get_z() - Zvtx; // correct for event vertex position
  double cemc_r2 = sqrt(pow(cemc_x2,2)+pow(cemc_y2,2));
  cemc_eta2 = asinh( cemc_z2 / cemc_r2 );
  cemc_phi2 = atan2( cemc_y2, cemc_x2 );
  cemc_ecore2 = clus2->get_ecore();

  if(Verbosity() > 0){
    std::cout << "ecore1: " << cemc_ecore1 << ", ecore2: " << cemc_ecore2 << std::endl;
  }

  SvtxTrackState *ihcal_trackstate1 = track1->get_state(ihcal_pathlength1);
  SvtxTrackState *ihcal_trackstate2 = track2->get_state(ihcal_pathlength2);

  if (!ihcal_trackstate1 || !ihcal_trackstate2) {
    if (Verbosity() > 0) std::cout << "No IHCAL state found → skip this pair\n";
    return;
  }

  ihcal_e3x3_1 = Get_CAL_e3x3(ihcal_trackstate1, m_ihcal_towers, track1_ihcal_eta, track1_ihcal_phi, ihcal_eta1, ihcal_phi1, Zvtx, 1);
  ihcal_e3x3_2 = Get_CAL_e3x3(ihcal_trackstate2, m_ihcal_towers, track2_ihcal_eta, track2_ihcal_phi, ihcal_eta2, ihcal_phi2, Zvtx, 1);

  SvtxTrackState *ohcal_trackstate1 = track1->get_state(ohcal_pathlength1);
  SvtxTrackState *ohcal_trackstate2 = track2->get_state(ohcal_pathlength2);

  if (!ohcal_trackstate1 || !ohcal_trackstate2) {
    if (Verbosity() > 0) std::cout << "No OHCAL state found → skip this pair\n";
    return;
  }

  ohcal_e3x3_1 = Get_CAL_e3x3(ohcal_trackstate1, m_ohcal_towers, track1_ohcal_eta, track1_ohcal_phi, ohcal_eta1, ohcal_phi1, Zvtx, 2);
  ohcal_e3x3_2 = Get_CAL_e3x3(ohcal_trackstate2, m_ohcal_towers, track2_ohcal_eta, track2_ohcal_phi, ohcal_eta2, ohcal_phi2, Zvtx, 2);

}

RawCluster* FastDecayReco::MatchClusterCEMC(SvtxTrackState* trackstate, RawClusterContainer* cemc_clusters, float& track_eta, float& track_phi, double Zvtx)
{

  // track coordinates
  double track_x = trackstate->get_x();
  double track_y = trackstate->get_y();
  double track_z = trackstate->get_z() - Zvtx;

  double track_r = sqrt(track_x*track_x+track_y*track_y);
  track_eta = asinh( track_z / track_r );
  track_phi = atan2( track_y, track_x );

  // clusters loop
  double dist = 99999.;
  RawCluster* returnCluster = nullptr;
  RawClusterContainer::Range begin_end = cemc_clusters->getClusters();
  RawClusterContainer::Iterator clusiter;

  for (clusiter = begin_end.first; clusiter != begin_end.second; ++clusiter)
  {
    RawCluster* cluster = clusiter->second;
    if(!cluster) { std::cout << "ERROR: bad cluster pointer = " << cluster << std::endl; continue; }
    else {
      double clus_x = cluster->get_x();
      double clus_y = cluster->get_y();
      double clus_z = cluster->get_z() - Zvtx; // correct for event vertex position
      double clus_r = sqrt(pow(clus_x,2)+pow(clus_y,2));
      double clus_eta = asinh( clus_z / clus_r );
      double clus_phi = atan2( clus_y, clus_x );
      double dPhi = TVector2::Phi_mpi_pi(track_phi - clus_phi);
      double dEta = track_eta-clus_eta;
      double tmpdist = sqrt(pow(dEta,2)+pow(dPhi,2));
      if(tmpdist<dist) { dist = tmpdist; returnCluster = cluster;}
    }
  }
  return returnCluster;
}

double FastDecayReco::Get_CAL_e3x3(SvtxTrackState* trackstate, TowerInfoContainerv4* _towersInfo, float& track_eta, float& track_phi, float& hcal_eta, float& hcal_phi, double Zvtx, int what)
{

  const float eta_mapIH[] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

  const float eta_mapOH[] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

  const float phi_mapIH[] = {0.0936424,0.191817,0.289992,0.388167,0.486341,0.584516,0.682691,0.780866,0.879041,0.977215,1.07539,1.17356,1.27174,1.36991,1.46809,1.56626,1.66444,1.76261,1.86079,1.95896,2.05714,2.15531,2.25349,2.35166,2.44984,2.54801,2.64619,2.74436,2.84254,2.94071,3.03889,3.13706,3.23523,3.33341,3.43158,3.52976,3.62793,3.72611,3.82428,3.92246,4.02063,4.11881,4.21698,4.31516,4.41333,4.51151,4.60968,4.70786,4.80603,4.90421,5.00238,5.10056,5.19873,5.29691,5.39508,5.49325,5.59143,5.6896,5.78778,5.88595,5.98413,6.0823,6.18048,6.27865};

  const float phi_mapOH[] = {0.0731489,0.171324,0.269498,0.367673,0.465848,0.564023,0.662198,0.760372,0.858547,0.956722,1.0549,1.15307,1.25125,1.34942,1.4476,1.54577,1.64395,1.74212,1.84029,1.93847,2.03664,2.13482,2.23299,2.33117,2.42934,2.52752,2.62569,2.72387,2.82204,2.92022,3.01839,3.11657,3.21474,3.31292,3.41109,3.50927,3.60744,3.70562,3.80379,3.90196,4.00014,4.09831,4.19649,4.29466,4.39284,4.49101,4.58919,4.68736,4.78554,4.88371,4.98189,5.08006,5.17824,5.27641,5.37459,5.47276,5.57094,5.66911,5.76729,5.86546,5.96363,6.06181,6.15998,6.25816};

  double track_x = trackstate->get_x();
  double track_y = trackstate->get_y();
  double track_z = trackstate->get_z() - Zvtx;
  double track_r = sqrt(track_x*track_x+track_y*track_y);
  track_eta = asinh( track_z / track_r );
  track_phi = atan2( track_y, track_x );

  double e3x3 = 0.;
  double dist = 9999.;
  //  double dphi = 9999.;
  //  double deta = 9999.;
  int m_channelkey = 0;
  int nchannels = _towersInfo->size();

  for (int channel = 0; channel < nchannels; channel++){
    int channelkey = _towersInfo->encode_key(channel);
    double tower_eta = 9999.;
    double tower_phi = 9999.;
    if(what == 1){
      tower_eta = eta_mapIH[_towersInfo->getTowerEtaBin(channelkey)];
      tower_phi = phi_mapIH[_towersInfo->getTowerPhiBin(channelkey)]-M_PI;
    } else if(what == 2){
      tower_eta = eta_mapOH[_towersInfo->getTowerEtaBin(channelkey)];
      tower_phi = phi_mapOH[_towersInfo->getTowerPhiBin(channelkey)]-M_PI;
    }
    double dPhi = TVector2::Phi_mpi_pi(track_phi - tower_phi);
    double dEta = track_eta-tower_eta;
    double tmpdist = sqrt(pow(dEta,2)+pow(dPhi,2));
    if(tmpdist<dist) { dist = tmpdist; m_channelkey = channelkey; hcal_eta = tower_eta; hcal_phi = tower_phi; }
  }

  //int m_channel = _towersInfo->decode_key(m_channelkey);
  unsigned int ieta = _towersInfo->getTowerEtaBin(m_channelkey);
  unsigned int jphi = _towersInfo->getTowerPhiBin(m_channelkey);

  if( Verbosity() > 0){
    std::cout << "dist: " << dist << std::endl;
  }

  TowerInfo* thetower = _towersInfo->get_tower_at_key(m_channelkey);
  if(!thetower) { return e3x3; }

  unsigned int maxbinphi = 63; //if(what==0) maxbinphi = 255;
  unsigned int maxbineta = 23; //if(what==0) maxbineta = 93;

  for(unsigned int i=0; i<=2; i++) {
    for(unsigned int j=0; j<=2; j++) {
      unsigned int itmp = ieta-1+i;
      unsigned int jtmp = 0;
      if(jphi==0 && j==0) { jtmp = maxbinphi; }      // wrap around
      else if(jphi==maxbinphi && j==2) { jtmp = 0; } // wrap around
      else { jtmp = jphi-1+j; }
      if(itmp<=maxbineta) {
        unsigned int tmpchannelkey = TowerInfoDefs::encode_hcal(itmp, jtmp);
        TowerInfo *tmptower = _towersInfo->get_tower_at_key(tmpchannelkey);
        if(tmptower) { e3x3 += tmptower->get_energy(); }
      }
    }
  }
  return e3x3;

}

void FastDecayReco::fillHistogram(Eigen::Vector3d mom1, Eigen::Vector3d mom2, float& decaymassA, float& decaymassB, TH1* massreco, double& invariantMass, double& invariantPt, float& invariantPhi, float& rapidity, float& pseudorapidity)
{

  double E1 = sqrt(pow(mom1(0), 2) + pow(mom1(1), 2) + pow(mom1(2), 2) + pow(decaymassA, 2));
  double E2 = sqrt(pow(mom2(0), 2) + pow(mom2(1), 2) + pow(mom2(2), 2) + pow(decaymassB, 2));

  TLorentzVector v1(mom1(0), mom1(1), mom1(2), E1);
  TLorentzVector v2(mom2(0), mom2(1), mom2(2), E2);

  TLorentzVector tsum;
  tsum = v1 + v2;

  rapidity = tsum.Rapidity();
  pseudorapidity = tsum.Eta();
  invariantMass = tsum.M();
  invariantPt = tsum.Pt();
  invariantPhi = tsum.Phi();

  if (Verbosity() > 2)
  {
    std::cout << "px1: " << mom1(0) << " py1: " << mom1(1) << " pz1: " << mom1(2) << " E1: " << E1 << std::endl;
    std::cout << "px2: " << mom2(0) << " py2: " << mom2(1) << " pz2: " << mom2(2) << " E2: " << E2 << std::endl;
    std::cout << "tsum: " << tsum(0) << " " << tsum(1) << " " << tsum(2) << " " << tsum(3) << std::endl;
    std::cout << "invariant mass: " << invariantMass << " invariant Pt: " << invariantPt << " invariantPhi: " << invariantPhi << std::endl;
  }

  if (invariantPt > invariant_pt_cut)
  {
    massreco->Fill(invariantMass);
  }
}

bool FastDecayReco::projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom)
{
  bool ret = true;

  /// create perigee surface
  ActsPropagator actsPropagator(_tGeometry);
  auto perigee = actsPropagator.makeVertexSurface(PCA);  // PCA is in cm here
  auto params = actsPropagator.makeTrackParams(track, m_vertexMap);
  if (!params.ok())
  {
    return false;
  }
  auto result = actsPropagator.propagateTrack(params.value(), perigee);

  if (result.ok())
  {
    auto projectionPos = result.value().second.position(_tGeometry->geometry().getGeoContext());
    const auto momentum = result.value().second.momentum();
    pos(0) = projectionPos.x() / Acts::UnitConstants::cm;
    pos(1) = projectionPos.y() / Acts::UnitConstants::cm;
    pos(2) = projectionPos.z() / Acts::UnitConstants::cm;

    if (Verbosity() > 2)
    {
      std::cout << "                 Input PCA " << PCA << "  projection out " << pos << std::endl;
    }

    mom(0) = momentum.x();
    mom(1) = momentum.y();
    mom(2) = momentum.z();
  }
  else
  {
    pos(0) = track->get_x();
    pos(1) = track->get_y();
    pos(2) = track->get_z();

    mom(0) = track->get_px();
    mom(1) = track->get_py();
    mom(2) = track->get_pz();

    if(Verbosity() > 0)
    {
      std::cout << result.error() << std::endl;
      std::cout << result.error().message() << std::endl;
      std::cout << " Failed projection of track with: " << std::endl;
      std::cout << " x,y,z = " << track->get_x() << "  " << track->get_y() << "  " << track->get_z() << std::endl;
      std::cout << " px,py,pz = " << track->get_px() << "  " << track->get_py() << "  " << track->get_pz() << std::endl;
      std::cout << " to point (x,y,z) = " << PCA(0) / Acts::UnitConstants::cm << "  " << PCA(1) / Acts::UnitConstants::cm << "  " << PCA(2) / Acts::UnitConstants::cm << std::endl;
    }

    //    ret = false;
  }

  return ret;
}

bool FastDecayReco::projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom)
{
  // Make a cylinder surface at the radius and project the track to that
  bool ret = true;
  const double eta = 2.0;
  const double theta = 2. * atan(exp(-eta));
  const double halfZ = Radius / tan(theta) * Acts::UnitConstants::cm;
  Radius *= Acts::UnitConstants::cm;

  /// Make a cylindrical surface at (0,0,0) aligned along the z axis
  auto transform = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> cylSurf =
    Acts::Surface::makeShared<Acts::CylinderSurface>(transform,
        Radius,
        halfZ);
  ActsPropagator actsPropagator(_tGeometry);
  auto params = actsPropagator.makeTrackParams(track, m_vertexMap);
  if (!params.ok())
  {
    return false;
  }

  auto result = actsPropagator.propagateTrack(params.value(), cylSurf);
  if (result.ok())
  {
    auto projectionPos = result.value().second.position(_tGeometry->geometry().getGeoContext());
    const auto momentum = result.value().second.momentum();
    pos(0) = projectionPos.x() / Acts::UnitConstants::cm;
    pos(1) = projectionPos.y() / Acts::UnitConstants::cm;
    pos(2) = projectionPos.z() / Acts::UnitConstants::cm;

    mom(0) = momentum.x();
    mom(1) = momentum.y();
    mom(2) = momentum.z();
  }
  else
  {
    ret = false;
  }

  return ret;
}

Acts::Vector3 FastDecayReco::getVertex(SvtxTrack* track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
  Acts::Vector3 vertex = Acts::Vector3::Zero();
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
    vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
    vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
  }

  return vertex;
}

void FastDecayReco::findPcaTwoTracks(float& decaymassA, float& decaymassB, const Acts::Vector3& pos1, const Acts::Vector3& pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca) const
{
  TLorentzVector v1;
  TLorentzVector v2;

  double px1 = mom1(0);
  double py1 = mom1(1);
  double pz1 = mom1(2);
  double px2 = mom2(0);
  double py2 = mom2(1);
  double pz2 = mom2(2);

  Float_t E1 = sqrt(pow(px1, 2) + pow(py1, 2) + pow(pz1, 2) + pow(decaymassA, 2));
  Float_t E2 = sqrt(pow(px2, 2) + pow(py2, 2) + pow(pz2, 2) + pow(decaymassB, 2));

  v1.SetPxPyPzE(px1, py1, pz1, E1);
  v2.SetPxPyPzE(px2, py2, pz2, E2);

  // calculate lorentz vector
  const Eigen::Vector3d& a1 = pos1;
  const Eigen::Vector3d& a2 = pos2;

  Eigen::Vector3d b1(v1.Px(), v1.Py(), v1.Pz());
  Eigen::Vector3d b2(v2.Px(), v2.Py(), v2.Pz());

  // The shortest distance between two skew lines described by
  //  a1 + c * b1
  //  a2 + d * b2
  // where a1, a2, are vectors representing points on the lines, b1, b2 are direction vectors, and c and d are scalars
  // dca = (b1 x b2) .(a2-a1) / |b1 x b2|

  // bcrossb/mag_bcrossb is a unit vector perpendicular to both direction vectors b1 and b2
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2 - a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  dca = 999;
  if (mag_bcrossb != 0)
  {
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  }
  else
  {
    return;  // same track, skip combination
  }

  // get the points at which the normal to the lines intersect the lines, where the lines are perpendicular
  double X = b1.dot(b2) - (b1.dot(b1) * b2.dot(b2) / b2.dot(b1));
  double Y = (a2.dot(b2) - a1.dot(b2)) - ((a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1));
  double c = Y / X;

  double F = b1.dot(b1) / b2.dot(b1);
  double G = -(a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = (c * F) + G;

  // then the points of closest approach are:
  pca1 = a1 + c * b1;
  pca2 = a2 + d * b2;

  return;
}

FastDecayReco::FastDecayReco(const std::string& name)
  : SubsysReco(name)
{
}

Acts::Vector3 FastDecayReco::calculateDca(SvtxTrack* track, const Acts::Vector3& momentum, Acts::Vector3 position)
{
  // For the purposes of this module, we set default values to prevent this track from being rejected if the dca calc fails
  Acts::Vector3 r = momentum.cross(Acts::Vector3(0., 0., 1.));
  float phi = atan2(r(1), r(0));
  Acts::Vector3 outVals(track_dca_cut*1.1, track_dca_cut*1.1, phi);
  auto vtxid = track->get_vertex_id();
  if (!m_vertexMap)
  {
    //std::cout << "Could not find m_vertexmap " << std::endl;
    return outVals;
  }
  auto *svtxVertex = m_vertexMap->get(vtxid);
  if (!svtxVertex)
  {
    //std::cout << "Could not find vtxid in m_vertexMap " << vtxid << std::endl;
    return outVals;
  }
  Acts::Vector3 vertex(svtxVertex->get_x(), svtxVertex->get_y(), svtxVertex->get_z());
  position -= vertex;

  Acts::RotationMatrix3 rot;
  rot(0, 0) = std::cos(phi);
  rot(0, 1) = -std::sin(phi);
  rot(0, 2) = 0;
  rot(1, 0) = std::sin(phi);
  rot(1, 1) = std::cos(phi);
  rot(1, 2) = 0;
  rot(2, 0) = 0;
  rot(2, 1) = 0;
  rot(2, 2) = 1;

  Acts::Vector3 pos_R = rot * position;
  double dca3dxy = pos_R(0);
  double dca3dz = pos_R(2);

  outVals(0) = abs(dca3dxy);
  outVals(1) = abs(dca3dz);
  outVals(2) = phi;

  if (Verbosity() > 4)
  {
    std::cout << " pre-position: " << position << std::endl;
    std::cout << " vertex: " << vertex << std::endl;
    std::cout << " vertex subtracted-position: " << position << std::endl;
  }

  return outVals;
}

/*
   float FastDecayReco::get_dEdx(SvtxTrack* track)
   {

   if(!m_cluster_map || !m_geom_container)
   {
   std::cout << "Can't continue in KFParticle_Tools::get_dEdx, returning -1" << std::endl;
   return -1.0;
   }

   TrackSeed *tpcseed = track->get_tpc_seed();
   std::vector<TrkrDefs::cluskey> clusterKeys;
   clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys());
   std::vector<float> dedxlist;
   for (unsigned long cluster_key : clusterKeys)
   {
   unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
   if(TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::TrkrId::tpcId)
   {
   continue;
   }
   TrkrCluster* cluster = m_cluster_map->findCluster(cluster_key);

   float adc = cluster->getAdc();
   PHG4TpcCylinderGeom* GeoLayer_local = m_geom_container->GetLayerCellGeom(layer_local);
   float thick = GeoLayer_local->get_thickness();

   float r = GeoLayer_local->get_radius();
   float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
   float beta = atan(tpcseed->get_slope());

   float alphacorr = cos(alpha);
   if(alphacorr<0||alphacorr>4)
   {
   alphacorr=4;
   }

   float betacorr = cos(beta);
   if(betacorr<0||betacorr>4)
   {
   betacorr=4;
   }

   adc/=thick;
   adc*=alphacorr;
   adc*=betacorr;
   dedxlist.push_back(adc);
   sort(dedxlist.begin(), dedxlist.end());
   }

   int trunc_min = 0;
   int trunc_max = (int)dedxlist.size()*0.7;
   float sumdedx = 0;
   int ndedx = 0;
   for(int j = trunc_min; j<=trunc_max;j++)
   {
   sumdedx+=dedxlist.at(j);
   ndedx++;
   }

   sumdedx/=ndedx;
   return sumdedx;
   }
   */

float FastDecayReco::get_dEdx(SvtxTrack* track)
{

  if(!m_cluster_map || !m_geom_container)
  {
    std::cout << "Can't continue in KFParticle_Tools::get_dEdx, returning -1" << std::endl;
    return -1.0;
  }

  TrackSeed *tpcseed = track->get_tpc_seed();

  float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
  layerThicknesses[0] = m_geom_container->GetLayerCellGeom(7)->get_thickness();
  layerThicknesses[1] = m_geom_container->GetLayerCellGeom(8)->get_thickness();
  layerThicknesses[2] = m_geom_container->GetLayerCellGeom(27)->get_thickness();
  layerThicknesses[3] = m_geom_container->GetLayerCellGeom(50)->get_thickness();

  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys());
  std::vector<float> dedxlist;

  for (unsigned long cluster_key : clusterKeys)
  {
    auto detid = TrkrDefs::getTrkrId(cluster_key);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    TrkrCluster* cluster = m_cluster_map->findCluster(cluster_key);
    float adc = cluster->getAdc();
    float thick = 0;
    if (layer_local < 23)
    {
      if (layer_local % 2 == 0)
      {
        thick = layerThicknesses[1];
      }
      else
      {
        thick = layerThicknesses[0];
      }
    }
    else if (layer_local < 39)
    {
      thick = layerThicknesses[2];
    }
    else
    {
      thick = layerThicknesses[3];
    }
    auto cglob = _tGeometry->getGlobalPosition(cluster_key, cluster);
    float r = std::sqrt(cglob(0) * cglob(0) + cglob(1) * cglob(1));
    float alpha = (r * r) / (2 * r * std::abs(1.0 / tpcseed->get_qOverR()));
    float beta = std::atan(tpcseed->get_slope());
    float alphacorr = std::cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = std::cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  if (dedxlist.empty())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;

}

int FastDecayReco::get_PID(float& mass, int charge)
{
  int PID = 0;
  if (std::abs(mass - m_pion) < mass_tolerance || std::abs(mass - m_electron) < mass_tolerance)
  {
    PID = (charge > 0) ? 211 : -211;
  }
  else if (std::abs(mass - m_kaon) < mass_tolerance)
  {
    PID = (charge > 0) ? 321 : -321;
  }
  else if (std::abs(mass - m_proton) < mass_tolerance)
  {
    PID = (charge > 0) ? 2212 : -2212;
  }
  else
  {
    std::cerr << "Unknown particle with mass = " << mass << ", charge = " << charge << std::endl;
  }

  return PID;

}

double FastDecayReco::get_dEdx_fitValue(float& qmomentum, int& PID)
{

  static bool initialized = false;

  if (!initialized)
  {
    std::string dedx_fitparams = CDBInterface::instance()->getUrl("TPC_DEDX_FITPARAM");
    TFile *filefit = TFile::Open(dedx_fitparams.c_str());

    if (!filefit->IsOpen())
    {
      std::cerr << "Error opening filefit!" << std::endl;
      return -1;
    }

    filefit->GetObject("f_piband", f_pion_plus);
    filefit->GetObject("f_Kband", f_kaon_plus);
    filefit->GetObject("f_pband", f_proton_plus);
    filefit->GetObject("f_piminus_band", f_pion_minus);
    filefit->GetObject("f_Kminus_band", f_kaon_minus);
    filefit->GetObject("f_pbar_band", f_proton_minus);

    pidMap[ 211]  = f_pion_plus;
    pidMap[ 321]  = f_kaon_plus;
    pidMap[2212]  = f_proton_plus;
    pidMap[-211]  = f_pion_minus;
    pidMap[-321]  = f_kaon_minus;
    pidMap[-2212] = f_proton_minus;

    initialized = true;
  }

  // Safety check
  if (pidMap.find(PID) == pidMap.end() || !pidMap[PID])
  {
    std::cerr << "No dEdx fit function for PID = " << PID << std::endl;
    return -1;
  }

  if( Verbosity() >  0)
  {
    std::cout << "qmomentum: " << qmomentum << " PID " << PID << " EvalFit: " << pidMap[PID]->Eval(qmomentum) << std::endl;
  }

  return pidMap[PID]->Eval(qmomentum);
}

int FastDecayReco::InitRun(PHCompositeNode* topNode)
{
  std::cout << "L=1" << std::endl; 
  const char* cfilepath = filepath.c_str();
  std::cout << "L=2" << std::endl;
  fout = new TFile(cfilepath, "recreate");
  std::cout << "L=3" << std::endl;
  ntp_reco_info = new TNtuple("ntp_reco_info", "decay_pairs", "id1:crossing1:x1:y1:z1:px1:py1:pz1:decaymassA:dEdx1:matching1:dca3dxy1:dca3dz1:phi1:pca_rel1_x:pca_rel1_y:pca_rel1_z:eta1:charge1:tpcClusters_1:id2:crossing2:x2:y2:z2:px2:py2:pz2:decaymassB:dEdx2:matching2:dca3dxy2:dca3dz2:phi2:pca_rel2_x:pca_rel2_y:pca_rel2_z:eta2:charge2:tpcClusters_2:vertex_x:vertex_y:vertex_z:pair_dca:invariant_mass:invariant_pt:invariantPhi:pathlength_x:pathlength_y:pathlength_z:pathlength:rapidity:pseudorapidity:projected_pos1_x:projected_pos1_y:projected_pos1_z:projected_pos2_x:projected_pos2_y:projected_pos2_z:projected_mom1_x:projected_mom1_y:projected_mom1_z:projected_mom2_x:projected_mom2_y:projected_mom2_z:projected_pca_rel1_x:projected_pca_rel1_y:projected_pca_rel1_z:projected_pca_rel2_x:projected_pca_rel2_y:projected_pca_rel2_z:projected_pair_dca:projected_pathlength_x:projected_pathlength_y:projected_pathlength_z:projected_pathlength:quality1:quality2:cosThetaReco:track1_silicon_clusters:track2_silicon_clusters:track1_mvtx_clusters:track1_mvtx_states:track1_intt_clusters:track1_intt_states:track2_mvtx_clusters:track2_mvtx_states:track2_intt_clusters:track2_intt_states:icomb:proj_cemc_phi1:proj_cemc_eta1:cemc_phi1:cemc_eta1:cemc_z1:cemc_ecore1:proj_cemc_phi2:proj_cemc_eta2:cemc_phi2:cemc_eta2:cemc_z2:cemc_ecore2:proj_ihcal_phi1:proj_ihcal_eta1:ihcal_phi1:ihcal_eta1:ihcal_e3x3_1:proj_ihcal_phi2:proj_ihcal_eta2:ihcal_phi2:ihcal_eta2:ihcal_e3x3_2:proj_ohcal_phi1:proj_ohcal_eta1:ohcal_phi1:ohcal_eta1:ohcal_e3x3_1:proj_ohcal_phi2:proj_ohcal_eta2:ohcal_phi2:ohcal_eta2:ohcal_e3x3_2:Zvtx:runNumber:eventNumber:trackid1:trackid2:run:segment:mb10:mb12:photon4:photon5:mbphoton3:mbphoton4:mbphoton5");
  std::cout << "L=4" << std::endl;
  getNodes(topNode);

  std::cout << "L=5" << std::endl;
  recomass = new TH1D("recomass", "recomass", 1000, 0.0, 1);  // root histogram arguments: name,title,bins,minvalx,maxvalx

  h_trigger = new TH1D("h_trigger","h_trigger", 41, -0.5, 40.5);

  std::cout << "L=6" << std::endl;
  //Add new track map to save selected tracks
  if (m_save_tracks)
  {
    PHNodeIterator nodeIter(topNode);

    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
      std::cout << "DST node added" << std::endl;
    }

    m_output_trackMap = new SvtxTrackMap_v2();
    PHIODataNode<PHObject> *outputTrackNode = new PHIODataNode<PHObject>(m_output_trackMap, m_output_trackMap_node_name, "PHObject");
    dstNode->addNode(outputTrackNode);
    if (Verbosity() > 1) { std::cout << m_output_trackMap_node_name << " node added" << std::endl; }
  }

  return 0;
}

int FastDecayReco::End(PHCompositeNode* /**topNode*/)
{
  fout->cd();
  ntp_reco_info->Write();
  recomass->Write();
  h_trigger->Write();
  /*
     for (size_t i = 0; i < m_debug_canvases.size(); ++i)
     {
     m_debug_canvases[i]->Write();  // Writes canvas with name "event_debug_i"
     }
     */
  fout->Close();

  return 0;
}

int FastDecayReco::getNodes(PHCompositeNode* topNode)
{
  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_svtxTrackMap)
  {
    std::cout << PHWHERE << "No SvtxTrackMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No vertexMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_global_vtxmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if(!m_global_vtxmap) {
    std::cerr << PHWHERE << " ERROR: Can not find GlobalVertexMap node." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_cemc_clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if(!m_cemc_clusters) {
    std::cerr << PHWHERE << " ERROR: Can not find CLUSTERINFO_CEMC node." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(Verbosity() > 0){
    std::cout << " Number of CEMC clusters: " << m_cemc_clusters->size() << std::endl;
  }

  m_ihcal_towers = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_HCALIN");
  if(!m_ihcal_towers) {
    std::cerr << PHWHERE << " ERROR: Can not find TOWERINFO_CALIB_HCALIN node." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(Verbosity() > 0){
    std::cout << " Number of HCALIN towers: " << m_ihcal_towers->size() << std::endl;
  }

  m_ohcal_towers = findNode::getClass<TowerInfoContainerv4>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if(!m_ohcal_towers) {
    std::cerr << PHWHERE << " ERROR: Can not find TOWERINFO_CALIB_HCALIN node." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(Verbosity() > 0){
    std::cout << " Number of HCALOUT towers: " << m_ohcal_towers->size() << std::endl;
  }

  m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  m_geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  gl1Packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if(!gl1Packet) {
    std::cerr << PHWHERE << " ERROR: Can not find GL1RAWHIT node." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
