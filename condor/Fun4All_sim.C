#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <trackreco/PHTpcDeltaZCorrection.h>

#include <tpccalib/PHTpcResiduals.h>

#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TpcSiliconQA.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

#include <ffamodules/CDBInterface.h>

#include <globalvertex/GlobalVertexReco.h>
#include <trackreco/PHActsTrkFitter.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <globalvertex/GlobalVertexReco.h>
#include <calovalid/CaloValid.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/DeadHotMapLoader.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/TowerInfoDeadHotMask.h>
#include <fastdecayreco/FastDecayReco.h>
#include <Calo_Calib.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <phool/recoConsts.h>
#include <string.h>

#include <fastdecayreco/FastDecayReco.h>

#include <MbdReco.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libzdcinfo.so)

R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)

void Fun4All_sim(int ct = 0, string inputFile = "/sphenix/user/pnietomar/simulation/generations/outputs/mb_MBDC1/DST_mb_0.root")
{
  char ListFile[99];
  char outputFile[99];
  sprintf(outputFile,"outputs/MB_sim_%d.root", ct);

  int runnumber = 0;

  auto *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "2025p009");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
    
  Fun4AllServer *se = Fun4AllServer::instance();

  TrackingInit();

  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  se->registerInputManager(in);

  in->AddFile(inputFile);

  Global_Reco();

  //MbdReco* mbdreco = new MbdReco();
  //se->registerSubsystem(mbdreco);

  std::string faststring = outputFile;
  auto fastreco = new FastDecayReco();
  fastreco->Verbosity(0);
  fastreco->setMC(true);
  fastreco->setTrackQualityCut(1000);
  fastreco->setPtCut(0.25);
  fastreco->setPairDCACut(5.0);
  fastreco->setRequireMVTX(true);
  fastreco->setTrackDCACut(0.0);  // requires fabs(dca) > this
  fastreco->setTrackCaloMatching(true); // if this is set as true then crossing == 0
  fastreco->setPhotonConv(false); // if this is set as true then set DecayMass1 and DecayMass2 as the electron mass
  fastreco->setDecayMass1(0.13957);    //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511) (kaons = 0.49367)
  fastreco->setDecayMass2(0.13957);  //(muons decaymass = 0.1057) (pions = 0.13957) (electron = 0.000511) (kaons = 0.49367)
  fastreco->setRunNumber(runnumber);
  fastreco->setSegment(ct);
  fastreco->set_output_file(faststring);
  se->registerSubsystem(fastreco);

  se->run(-1);
  se->End();
  se->PrintTimer();
  gSystem->Exit(0);

}

