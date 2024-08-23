// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   AgingLaserTask.h
/// \author Andreas Molander <andreas.molander@cern.ch>, Sandor Lokos <sandor.lokos@cern.ch>, Edmundo Garcia-Solis <edmundo.garcia@cern.ch>
///

#ifndef QC_MODULE_FT0_VISUALIZATIONTASK_H
#define QC_MODULE_FT0_VISUALIZATIONTASK_H

#include "QualityControl/TaskInterface.h"

#include <CommonConstants/LHCConstants.h>
#include <FT0Base/Constants.h>

#include <TH1I.h>
#include <TH2I.h>
#include <TH2Poly.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <vector>

using namespace o2::quality_control::core;

namespace o2::quality_control_modules::ft0
{

class VisualizationTask final : public TaskInterface
{
 public:
  VisualizationTask() = default;
  ~VisualizationTask() override;

  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(const Activity& activity) override;
  void startOfCycle() override;
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(const Activity& activity) override;
  void reset() override;

 private:
  /// Check if the laser was triggered for this BC
  /// \param bc      BC to check
  /// \param bcDelay Expected BC delay from trigger to signal
  /// \return        True if the laser was triggered
  bool bcIsTrigger(int bc, int bcDelay) const;

  /// Check if a detector signal is expected for this BC
  /// \param bc BC to check
  /// \return   True if a detector signal is expected
  bool bcIsDetector(int bc) const;

  /// Check if the first reference signal is expected for this BC
  /// \param bc      BC to check
  /// \param refChId Rerernce channel where signal is seen
  /// \return        True if the first reference signal is expected for this BC
  bool bcIsPeak1(int bc, int refChId) const;

  /// Check if the second reference signal is expected for this BC
  /// \param bc      BC to check
  /// \param refChId Rerernce channel where signal is seen
  /// \return        True if the second reference signal is expected for this BC
  bool bcIsPeak2(int bc, int refChId) const;

  void buildFT0Frames();

  // Constants
  
  constexpr static std::size_t sNCHANNELS_PM = o2::ft0::Constants::sNCHANNELS_PM; ///< Max number of FT0 channels
  constexpr static std::size_t sMaxBC = o2::constants::lhc::LHCMaxBunches;        ///< Max number of BCs
  
    float XChannelPosition[208] = {
    59.75,86.25,86.25,59.75,74.25,47.75,74.25,47.75,13.25,13.25,-13.25,-13.25,
    -47.75,-47.75,-74.25,-74.25,-59.75,-86.25,-86.25,-59.75,-74.25,-47.75,-74.25,
    -47.75,-13.25,-13.25,13.25,13.25,47.75,47.75,74.25,74.25,120.75,147.25,147.25,
    120.75,135.25,108.75,135.25,108.75,135.25,108.75,135.25,108.75,74.25,74.25,
    47.75,47.75,13.25,13.25,-13.25,-13.25,-47.75,-47.75,-74.25,-74.25,-108.75,-108.75,
    -135.25,-135.25,-108.75,-135.25,-108.75,-135.25,-120.75,-147.25,-147.25,-120.75,-135.25,-108.75,
    -135.25,-108.75,-135.25,-108.75,-135.25,-108.75,-74.25,-74.25,-47.75,-47.75,-13.25,-13.25,
    13.25,13.25,47.75,47.75,74.25,74.25,108.75,108.75,135.25,135.25,108.75,135.25,
    108.75,135.25,103.2,76.9,103.1,76.8,103.2,76.8,103.2,76.8,43.2,43.2,
    16.8,16.8,-16.8,-16.8,-43.2,-43.2,-76.8,-76.8,-103.2,-103.2,-76.8,-103.1,
    -76.9,-103.2,-103.2,-76.9,-103.1,-76.8,-103.2,-76.8,-103.2,-76.8,-43.2,-43.2,
    -16.8,-16.8,16.8,16.8,43.2,43.2,76.8,76.8,103.2,103.2,76.8,103.1,76.9,103.2,163,137,163,137,163,137,
    163.0,137.0,103.4,102.9,77.1,76.6,43.3,43.2,16.9,16.7,-16.7,-16.9,-43.2,-43.3,-76.6,-77.1,-102.9,
    -103.4,-137,-163,-137,-163,-137,-163,-137,-163,-163,-137,-163,-137,-163,-137,-163,-137,-103.4,
    -102.9,-77.1,-76.6,-43.3,-43.2,-16.9,-16.7,16.7,16.9,43.2,43.3,76.6,77.1,102.9,103.4,137,163,137,
    163.0,137.0,163.0,137.0,163.0};
  
  float YChannelPosition[208] = {
    -13.25,-13.25,13.25,13.25,47.75,47.75,74.25,74.25,60.75,87.25,87.25,60.75,
    74.25,47.75,74.25,47.75,13.25,13.25,-13.25,-13.25,-47.75,-47.75,-74.25,
    -74.25,-60.75,-87.25,-87.25,-60.75,-74.25,-47.75,-74.25,-47.75,-13.25,-13.25,
    13.25,13.25,47.75,47.75,74.25,74.25,108.75,108.75,135.25,135.25,108.75,135.25,
    108.75,135.25,121.75,148.25,148.25,121.75,135.25,108.75,135.25,108.75,135.25,
    108.75,135.25,108.75,74.25,74.25,47.75,47.75,13.25,13.25,-13.25,-13.25,-47.75,
    -47.75,-74.25,-74.25,-108.8,-108.8,-135.3,-135.3,-108.8,-135.3,-108.8,-135.3,
    -121.8,-148.3,-148.3,-121.8,-135.3,-108.8,-135.3,-108.8,-135.3,-108.8,-135.3,
    -108.8,-74.25,-74.25,-47.75,-47.75,17.8,17.8,44.2,44.2,78.7,79,105,105.3,78.8,
    105.1,78.9,105.2,105.2,78.9,105.1,78.8,105.3,79,105,78.7,44.2,44.2,17.8,17.8,
    -17.8,-17.8,-44.2,-44.2,-78.7,-79,-105,-105.3,-78.8,-105.1,-78.9,-105.2,-105.2,
    -78.9,-105.1,-78.8,-105.3,-79,-105,-78.7,-44.2,-44.2,-17.8,-17.8,18.7,18.9,45.2,
    45.3,78.6,79.1,104.9,105.4,138,164,138,164,139,165,139,165,165,139,165,139,164,
    138,164,138,105.4,104.9,79.1,78.6,45.3,45.2,18.9,18.7,-18.7,-18.9,-45.2,-45.3,-78.6,
    -79.1,-104.9,-105.4,-138,-164,-138,-164,-139,-165,-139,-165,-165,-139,-165,-139,-164,
    -138,-164,-138,-105.4,-104.9,-79.1,-78.6,-45.3,-45.2,-18.9,-18.7};
  
  // Task parameters
  std::vector<uint8_t> mDetectorChIDs;            ///< Enabled detector channels
  std::vector<uint8_t> mReferenceChIDs;           ///< Enabled reference channels
  int mDetectorAmpCut;                            ///< Amplitude cut for detector channels
  int mReferenceAmpCut;                           ///< Amplitude cut for reference channels
  std::vector<int> mLaserTriggerBCs;              ///< BCs in which the laser is triggered
  int mDetectorBCdelay;                           ///< BC delay for detector channels (same for all)
  std::map<uint8_t, int> mReferencePeak1BCdelays; ///< BC delays for reference channel peak 1 (per channel)
  std::map<uint8_t, int> mReferencePeak2BCdelays; ///< BC delays for reference channel peak 2 (per channel)
  std::vector<uint8_t> mvisualizationOption;      ///< Visualization option parameters {1,0,0} Amplitude {0,1,0} Occupancy, {0,0,1} DeadChannelMap
  
  bool mDebug = false; ///< Enable more histograms in debug mode

  // Histograms
  std::unique_ptr<TH1I> visnormalize;
  std::unique_ptr<TH1I> vischannelMap;
  std::unique_ptr<TH1F> vischannelAmplitude;
  std::unique_ptr<TH1F> vischannelCDFTime;
  std::unique_ptr<TH2Poly> mFT0AFrame;      ///< Frame for FT0A visualization
  std::unique_ptr<TH2Poly> mFT0CFrame;      ///< Frame for FT0C visualization
  
 
};

} // namespace o2::quality_control_modules::ft0

#endif // QC_MODULE_FT0_VISUALIZATIONTASK_H
