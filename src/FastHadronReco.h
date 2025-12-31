#ifndef KSHORTRECONSTRUCTION_H
#define KSHORTRECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <globalvertex/SvtxVertexMap.h>
#include <ffarawobjects/Gl1Packetv1.h>

#include <trackbase_historic/SvtxTrackMap.h>

#include <Acts/Definitions/Algebra.hpp>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <Eigen/Dense>

class TFile;
class TH1;
class TNtuple;
class TF1;
class TCanvas;

class ActsGeometry;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class RawCluster;
class RawClusterContainer;
class TowerInfoContainerv4;
class GlobalVertexMap;
class TrkrClusterContainer;
class PHG4TpcCylinderGeomContainer;
class Gl1Packet;

class FastHadronReco : public SubsysReco
{
 public:
  FastHadronReco(const std::string& name = "FastHadronReco");
  virtual ~FastHadronReco() = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

  void setPtCut(double ptcut) { invariant_pt_cut = ptcut; }
  void setTrackPtCut(double ptcut) { track_pt_cut = ptcut; }
  void setTrackQualityCut(double cut) { _qual_cut = cut; }
  void setPairDCACut(double cut) { pair_dca_cut = cut; }
  void setTrackDCACut(double cut) { track_dca_cut = cut; }
  void setRequireMVTX(bool set) { _require_mvtx = set; }
  void setDecayMass(float decayMassSet) { decaymass = decayMassSet; }  //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511)
  void setDecayMass1(Float_t decayMassSet1) { decaymass1 = decayMassSet1; }
  void setDecayMass2(Float_t decayMassSet2) { decaymass2 = decayMassSet2; }
  void setRunNumber(int runnumber){ _runnumber = runnumber;}
  void setSegment(int segment){ _segment = segment;}
  void setPhotonConv(bool photonconv){ _photonconv = photonconv;}
  void set_output_file(const std::string& outputfile) { filepath = outputfile; }
  void save_tracks(bool save = true) { m_save_tracks = save; }

 private:

  void fillNtp(SvtxTrack* track1, SvtxTrack* track2, float decaymassA, float decaymassB, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float invariantPhi, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2, Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj,unsigned int track1_silicon_cluster_size, unsigned int track2_silicon_cluster_size, unsigned int track1_mvtx_cluster_size, unsigned int track1_mvtx_state_size, unsigned int track1_intt_cluster_size, unsigned int track1_intt_state_size, unsigned int track2_mvtx_cluster_size, unsigned int track2_mvtx_state_size, unsigned int track2_intt_cluster_size, unsigned int track2_intt_state_size, int icomb, float track1_cemc_phi, float track1_cemc_eta, float cemc_phi1, float cemc_eta1, float cemc_z1, float cemc_ecore1, float track2_cemc_phi, float track2_cemc_eta, float cemc_phi2, float cemc_eta2, float cemc_z2, float cemc_ecore2, float track1_ihcal_phi, float track1_ihcal_eta, float ihcal_phi1, float ihcal_eta1, float ihcal_e3x3_1, float track2_ihcal_phi, float track2_ihcal_eta, float ihcal_phi2, float ihcal_eta2, float ihcal_e3x3_2, float track1_ohcal_phi, float track1_ohcal_eta, float ohcal_phi1, float ohcal_eta1, float ohcal_e3x3_1, float track2_ohcal_phi, float track2_ohcal_eta, float ohcal_phi2, float ohcal_eta2, float ohcal_e3x3_2, double Zvtx, int runNumber, int eventNumber, bool isMB10, bool isMB12, bool isPhoton4, bool isPhoton5, bool isMBPhoton3, bool isMBPhoton4, bool isMBPhoton5);

  void fillHistogram(Eigen::Vector3d mom1, Eigen::Vector3d mom2, float& decaymassA, float& decaymassB, TH1* massreco, double& invariantMass, double& invariantPt, float& invariantPhi, float& rapidity, float& pseudorapidity);

  void match_calos(SvtxTrack* track1, SvtxTrack* track2, float& track1_cemc_phi, float& track1_cemc_eta, float& cemc_phi1, float& cemc_eta1, float& cemc_z1, float& cemc_ecore1, float& track2_cemc_phi, float& track2_cemc_eta, float& cemc_phi2, float& cemc_eta2, float& cemc_z2, float& cemc_ecore2, float& track1_ihcal_phi, float& track1_ihcal_eta, float& ihcal_phi1, float& ihcal_eta1, float& ihcal_e3x3_1, float& track2_ihcal_phi, float& track2_ihcal_eta, float& ihcal_phi2, float& ihcal_eta2, float& ihcal_e3x3_2, float& track1_ohcal_phi, float& track1_ohcal_eta, float& ohcal_phi1, float& ohcal_eta1, float& ohcal_e3x3_1, float& track2_ohcal_phi, float& track2_ohcal_eta, float& ohcal_phi2, float& ohcal_eta2, float& ohcal_e3x3_2, double Zvtx);

  RawCluster* MatchClusterCEMC(SvtxTrackState* trackstate, RawClusterContainer* cemc_clusters, float& track_phi, float& track_eta, double Zvtx);

  double Get_CAL_e3x3(SvtxTrackState* trackstate, TowerInfoContainerv4* m_ohcal_towers, float& track_eta, float& track_phi, float& hcal_eta, float& hcal_phi, double Zvtx, int what);

  // void findPcaTwoTracks(SvtxTrack *track1, SvtxTrack *track2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca);
  void findPcaTwoTracks(float& decaymassA, float& decaymassB, const Acts::Vector3& pos1, const Acts::Vector3& pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca) const;

  float get_dEdx(SvtxTrack* track);
  //void init_dEdx_fits();
  int get_PID(float& mass, int charge);
  double get_dEdx_fitValue(float& qmomentum, int& PID);

  int getNodes(PHCompositeNode* topNode);

  Acts::Vector3 calculateDca(SvtxTrack* track, const Acts::Vector3& momentum, Acts::Vector3 position);

  bool projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom);
  bool projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom);

  Acts::Vector3 getVertex(SvtxTrack* track);
  static std::vector<unsigned int> getTrackStates(SvtxTrack *track);
  
  TNtuple* ntp_reco_info {nullptr};
  ActsGeometry* _tGeometry {nullptr};
  SvtxTrackMap* m_svtxTrackMap {nullptr};
  SvtxVertexMap* m_vertexMap {nullptr};
  SvtxTrackMap* m_dst_trackmap {nullptr};
  Gl1Packet* gl1Packet {nullptr};
  TrkrClusterContainer* m_cluster_map {nullptr};
  PHG4TpcCylinderGeomContainer* m_geom_container {nullptr};
  GlobalVertexMap* m_global_vtxmap {nullptr};
  RawClusterContainer* m_cemc_clusters {nullptr};
  TowerInfoContainerv4* m_ihcal_towers {nullptr};
  TowerInfoContainerv4* m_ohcal_towers {nullptr};
  TowerInfoContainerv4* _geomIH {nullptr};
  TowerInfoContainerv4* _geomOH {nullptr};

  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};

  std::map<int, TF1*> pidMap;
  
  std::string filepath {""};
  std::string m_trk_map_node_name;
  float decaymass {0.13957};  // pion decay mass
  float decaymass1 {0.13957}; // pions as default
  float decaymass2 {0.13957}; // pions as default
  bool _require_mvtx {true};
  bool _photonconv {false};
  double _qual_cut {1000.0};
  double pair_dca_cut {0.05};  // kshort relative cut 500 microns
  double track_dca_cut {0.01};
  double invariant_pt_cut {0.1};
  double track_pt_cut {0.2};
  int _projection {0};
  int _runnumber {53877};
  int _segment {0};

 // std::vector<TCanvas*> m_debug_canvases;

  double m_electron {0.000511};
  double m_pion   {0.13957};
  double m_kaon   {0.49367};
  double m_proton {0.93827};
  double mass_tolerance {0.05};

  TFile* fout {nullptr};
  TH1* recomass {nullptr};
  TH1* h_trigger {nullptr};

  bool m_save_tracks {false};
  SvtxTrackMap *m_output_trackMap {nullptr};
  std::string m_output_trackMap_node_name {"FastHadronReco_SvtxTrackMap"};
};

#endif  // KSHORTRECONSTRUCTION_H
