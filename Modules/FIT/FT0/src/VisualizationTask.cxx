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
/// \file   VisualizationTask.cxx
/// \author Edmundo Garcia-Solis <edmundo.garcia@cern.ch>
///

#include "FT0/VisualizationTask.h"

#include "Common/Utils.h"
#include "FITCommon/HelperCommon.h"
#include "QualityControl/QcInfoLogger.h"

#include <DataFormatsFT0/ChannelData.h>
#include <DataFormatsFT0/Digit.h>
#include <Framework/InputRecordWalker.h>
#include <Framework/DataRefUtils.h>

#include <TH1I.h>
#include <TH2I.h>
#include <TH2Poly.h>
#include <cstdint>
#include <memory>
#include <vector>

namespace o2::quality_control_modules::ft0
{

VisualizationTask::~VisualizationTask()
{
}

void VisualizationTask::initialize(o2::framework::InitContext&)
{
  // Read task parameters
  
  // Enabled detector channels
  const std::string detectorChannels = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "detectorChannelIDs", "");
  if (detectorChannels.size()) {
    mDetectorChIDs = fit::helper::parseParameters<uint8_t>(detectorChannels, ",");
  } else {
    // Not specified, enable all
    for (uint8_t chId = 0; chId < sNCHANNELS_PM; chId++) {
      mDetectorChIDs.push_back(chId);
    }
  }

  // Enabled reference channels
  const std::string referenceChannels = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "referenceChannelIDs", "");
  if (referenceChannels.size()) {
    mReferenceChIDs = fit::helper::parseParameters<uint8_t>(referenceChannels, ",");
  } else {
    // Not specified, enable all
    // TODO: return with fatal if not specified, to avoid hard coded numbers?
    for (uint8_t chId = 208; chId < 211; chId++) {
      mReferenceChIDs.push_back(chId);
    }
  }

  mDetectorAmpCut = o2::quality_control_modules::common::getFromConfig<int>(mCustomParameters, "detectorAmpCut", 0);
  mReferenceAmpCut = o2::quality_control_modules::common::getFromConfig<int>(mCustomParameters, "referenceAmpCut", 100);

  // BCs
  // Laser trigger BCs
  const std::string laserTriggerBCs = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "laserTriggerBCs", "");
  if (laserTriggerBCs.size()) {
    const auto vecParams = fit::helper::parseParameters<int>(laserTriggerBCs, ",");
    for (const int bc : vecParams) {
      mLaserTriggerBCs.push_back(bc);
    }
  }
  if (mLaserTriggerBCs.size() == 0) {
    //LOG(fatal) << "No laser trigger BCs specified in QC config!";
  }

  // BC delay for detector channels
  mDetectorBCdelay = o2::quality_control_modules::common::getFromConfig<int>(mCustomParameters, "detectorBCdelay", -1);
  if (mDetectorBCdelay < 0) {
    //LOG(fatal) << "No detector BC delay specified in QC config!";
  }

  // BC delay for reference channels peak1
  const std::string referencePeak1BCdelays = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "referencePeak1BCdelays", "");
  if (referencePeak1BCdelays.size()) {
    const auto vecParams = fit::helper::parseParameters<int>(referencePeak1BCdelays, ",");
    if (vecParams.size() != mReferenceChIDs.size()) {
      //LOG(fatal) << "Number of reference channels and reference peak 1 BC delays do not match!";
    }
    for (int i = 0; i < mReferenceChIDs.size(); i++) {
      mReferencePeak1BCdelays.insert({ mReferenceChIDs.at(i), vecParams.at(i) });
    }
  } else {
    //LOG(fatal) << "No reference peak 1 BC delays specified in QC config!";
  }

  // BC delay for reference channels peak2
  const std::string referencePeak2BCdelays = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "referencePeak2BCdelays", "");
  if (referencePeak2BCdelays.size()) {
    const auto vecParams = fit::helper::parseParameters<int>(referencePeak2BCdelays, ",");
    if (vecParams.size() != mReferenceChIDs.size()) {
      //LOG(fatal) << "Number of reference channels and reference peak 2 BC delays do not match!";
    }
    for (int i = 0; i < mReferenceChIDs.size(); i++) {
      mReferencePeak2BCdelays.insert({ mReferenceChIDs.at(i), vecParams.at(i) });
    }
  } else {
    //LOG(fatal) << "No reference peak 2 BC delays specified in QC config!";
  }

  //Visualization option
  const std::string visualizationOption = o2::quality_control_modules::common::getFromConfig<std::string>(mCustomParameters, "visualizationOption", "");
  if (visualizationOption.size()) mvisualizationOption = fit::helper::parseParameters<uint8_t>(visualizationOption, ",");
    else {
      // Not specified, enable amplitude
      mvisualizationOption[0] = 1;
      mvisualizationOption[1] = 0;
      mvisualizationOption[2] = 0;
    }
  
  mDebug = o2::quality_control_modules::common::getFromConfig<bool>(mCustomParameters, "debug", false);
  if (mDebug) {
    LOG(warning) << "Running in debug mode!";
  }

  // Initialize histograms
  visnormalize= std::make_unique<TH1I>("visnormalize","Normalize",5,0,5);
  vischannelMap = std::make_unique<TH1I>("vischannelMap","Channel Map",208,0,208);
  vischannelAmplitude = std::make_unique<TH1F>("vischannelAmplitude","Channel Amplitude",208,0.,208.);
  vischannelCDFTime = std::make_unique<TH1F>("vischannelCDFTime","Channel CDF Time",208,0.,208.);
  mFT0AFrame= std::make_unique<TH2Poly>("mFT0AFrame","FT0-A",-20.0,20.0,-20.0,20.0);
  mFT0CFrame= std::make_unique<TH2Poly>("mFT0CFrame","FT0-C",-20.0,20.0,-20.0,20.0);
 
  mFT0AFrame->SetStats(0);
  mFT0CFrame->SetStats(0);
 
  buildFT0Frames();
  getObjectsManager()->startPublishing(visnormalize.get());
  getObjectsManager()->startPublishing(vischannelMap.get());
  getObjectsManager()->startPublishing(vischannelAmplitude.get());
   getObjectsManager()->startPublishing(vischannelCDFTime.get());
  getObjectsManager()->startPublishing(mFT0AFrame.get());
  getObjectsManager()->setDefaultDrawOptions(mFT0AFrame.get(), "TEXT COLZ L");
  getObjectsManager()->startPublishing(mFT0CFrame.get());
  getObjectsManager()->setDefaultDrawOptions(mFT0CFrame.get(), "TEXT COLZ L");

}

void VisualizationTask::startOfActivity(const Activity& activity)
{
  reset();
}

void VisualizationTask::startOfCycle()
{ 
}

void VisualizationTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  auto channels = ctx.inputs().get<gsl::span<o2::ft0::ChannelData>>("channels");
  auto digits = ctx.inputs().get<gsl::span<o2::ft0::Digit>>("digits");
  buildFT0Frames();
  // Loop over digits
  for (const auto& digit : digits) {
    const int bc = digit.getIntRecord().bc;
    const auto& digitChannelData = digit.getBunchChannelData(channels);
    
    // Conditions wether to fill BC histograms for this BC. For debug use only
    // 'AmpCut' means there was at least one channel with chAmp > mReferenceAmpCut
    // 'ADCX' means there was at least one channel data with ADCX
    // 'Detector' means there was at least one detector channel
    // 'Reference' means there was at least one reference channel
    bool bcHasAmpCut = false;
    bool bcHasAmpCutADC0 = false;
    bool bcHasAmpCutADC1 = false;
    bool bcHasDetectorCh = false;
    bool bcHasDetectorChAmpCut = false;
    bool bcHasDetectorChAmpCutADC0 = false;
    bool bcHasDetectorChAmpCutADC1 = false;
    bool bcHasReferenceCh = false;
    bool bcHasReferenceChAmpCut = false;
    bool bcHasReferenceChAmpCutADC0 = false;
    bool bcHasReferenceChAmpCutADC1 = false;
    
    
    // Loop over channels
    for (const auto& chData : digitChannelData) {
      const int chId = chData.ChId;
      const int chAmp = chData.QTCAmpl;
      const int chTime = chData.CFDTime;
      const bool isRef = std::find(mReferenceChIDs.begin(), mReferenceChIDs.end(), chId) != mReferenceChIDs.end(); // TODO: optimize
      const bool isDet = !isRef;
      const bool isADC0 = !chData.getFlag(o2::ft0::ChannelData::kNumberADC);
      const bool isADC1 = !isADC0;
      const bool isAmpCutOk = chAmp > mReferenceAmpCut;
      const bool isDetAmpCutOk = chAmp > mDetectorAmpCut; // TODO: use this
      const bool isRefAmpCutOk = chAmp > mReferenceAmpCut;
      
      // Use var = condition || var, so that var is never set to back to false if once true
      bcHasAmpCut = isAmpCutOk || bcHasAmpCut;
      bcHasAmpCutADC0 = (isAmpCutOk && isADC0) || bcHasAmpCutADC0;
      bcHasAmpCutADC1 = (isAmpCutOk && isADC1) || bcHasAmpCutADC1;
      bcHasDetectorCh = isDet || bcHasDetectorCh;
      bcHasDetectorChAmpCut = (isDet && isAmpCutOk) || bcHasDetectorChAmpCut;
      bcHasDetectorChAmpCutADC0 = (isDet && isAmpCutOk && isADC0) || bcHasDetectorChAmpCutADC0;
      bcHasDetectorChAmpCutADC1 = (isDet && isAmpCutOk && isADC1) || bcHasDetectorChAmpCutADC1;
      bcHasReferenceCh = isRef || bcHasReferenceCh;
      bcHasReferenceChAmpCut = (isRef && isAmpCutOk) || bcHasReferenceChAmpCut;
      bcHasReferenceChAmpCutADC0 = (isRef && isAmpCutOk && isADC0) || bcHasReferenceChAmpCutADC0;
      bcHasReferenceChAmpCutADC1 = (isRef && isAmpCutOk && isADC1) || bcHasReferenceChAmpCutADC1;


      visnormalize->Fill(1,1);

      if(mvisualizationOption[0]==1){
	//channel map option
	for(int ii=0; ii<209; ii++)vischannelMap->Fill(ii,ii);
	for(int ii=0; ii<96; ii++)
	  mFT0AFrame->SetBinContent(mFT0AFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),(vischannelMap->GetBinContent(vischannelMap->FindBin(ii+1)))/visnormalize->GetBinContent(2));
	for(int ii=96; ii<209; ii++)
	  mFT0CFrame->SetBinContent(mFT0CFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),(vischannelMap->GetBinContent(vischannelMap->FindBin(ii+1)))/visnormalize->GetBinContent(2));
      }
      if(mvisualizationOption[0]==2){
	// channel amplitud
	vischannelAmplitude->Fill(chId,chAmp);
	for(int ii=0; ii<96; ii++)
	  mFT0AFrame->SetBinContent(mFT0AFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),vischannelAmplitude->GetBinContent(ii+1)/visnormalize->GetBinContent(2));
	for(int ii=96; ii<209; ii++)
	  mFT0CFrame->SetBinContent(mFT0CFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),vischannelAmplitude->GetBinContent(ii+1)/visnormalize->GetBinContent(2));
      }
      if(mvisualizationOption[0]==3){
	//channel CDF time
	vischannelCDFTime->Fill(chId,chTime);
	for(int ii=0; ii<96; ii++)
	  mFT0AFrame->SetBinContent(mFT0AFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),vischannelCDFTime->GetBinContent(ii+1)/visnormalize->GetBinContent(2));
	for(int ii=96; ii<209; ii++)
	  mFT0CFrame->SetBinContent(mFT0CFrame->FindBin(XChannelPosition[ii]/10.0, YChannelPosition[ii]/10.0),vischannelCDFTime->GetBinContent(ii+1)/visnormalize->GetBinContent(2));
      }
    }
  } // digit loop
  
 
} // monitorData



  
  //Channel map check
  // monitorData

void VisualizationTask::endOfCycle()
{
  ILOG(Debug, Devel) << "endOfCycle" << ENDM;
}

void VisualizationTask::endOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "endOfActivity" << ENDM;
}

void VisualizationTask::reset()
{
  // Reset histograms
  ILOG(Debug, Devel) << "Resetting the histograms" << ENDM;
  
  mFT0AFrame->TH2Poly::Reset("content");
  mFT0CFrame->TH2Poly::Reset("content");
  visnormalize->Reset("content");
  vischannelMap->Reset("content");
  vischannelAmplitude->Reset("content");
  vischannelCDFTime->Reset("content");
}

bool VisualizationTask::bcIsTrigger(int bc, int bcDelay) const
{
  for (const int bcTrg : mLaserTriggerBCs) {
    if (bc == bcTrg + bcDelay) {
      return true;
    }
  }
  return false;
}

bool VisualizationTask::bcIsDetector(int bc) const
{
  return bcIsTrigger(bc, mDetectorBCdelay);
}

bool VisualizationTask::bcIsPeak1(int bc, int refChId) const
{
  return bcIsTrigger(bc, mReferencePeak1BCdelays.at(refChId));
}

bool VisualizationTask::bcIsPeak2(int bc, int refChId) const
{
  return bcIsTrigger(bc, mReferencePeak2BCdelays.at(refChId));
}

void VisualizationTask::buildFT0Frames(){
  //Construct mFT0AFrame
  float X=-14.46;
  float Y=-14.56;
  float sizeCell=2.66;
  float sizeAlum=0.6;
  float sizeShft=0.9;
  //row1
  mFT0AFrame->AddBin(X,                      Y,  X+sizeCell,              Y+sizeCell);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y,  X+sizeCell*2.0,          Y+sizeCell*1.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*1.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*1.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y-sizeShft,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*1.0-sizeShft);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y-sizeShft,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*1.0-sizeShft);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*1.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*1.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*1.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*1.0);  //bin10
  //row2
  mFT0AFrame->AddBin(X,                      Y+sizeCell*1.0,  X+sizeCell,              Y+sizeCell*2.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*1.0,  X+sizeCell*2.0,          Y+sizeCell*2.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*1.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*2.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*1.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*2.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y-sizeShft+sizeCell*1.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*2.0-sizeShft);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y-sizeShft+sizeCell*1.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*2.0-sizeShft);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*1.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*2.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*1.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*2.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*1.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*2.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*1.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*2.0);  //bin10  
  //row3
  mFT0AFrame->AddBin(X,                      Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell,              Y+sizeCell*3.0+sizeAlum*1.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*2.0,          Y+sizeCell*3.0+sizeAlum*1.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*3.0+sizeAlum*1.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*3.0+sizeAlum*1.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y-sizeShft+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*3.0-sizeShft+sizeAlum*1.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y-sizeShft+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*3.0-sizeShft+sizeAlum*1.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*3.0+sizeAlum*1.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*3.0+sizeAlum*1.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*3.0+sizeAlum*1.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*2.0+sizeAlum*1.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*3.0+sizeAlum*1.0);  //bin10  
  //row4
  mFT0AFrame->AddBin(X,                      Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell,              Y+sizeCell*4.0+sizeAlum*1.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*2.0,          Y+sizeCell*4.0+sizeAlum*1.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*4.0+sizeAlum*1.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*4.0+sizeAlum*1.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y-sizeShft+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*4.0-sizeShft+sizeAlum*1.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y-sizeShft+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*4.0-sizeShft+sizeAlum*1.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*4.0+sizeAlum*1.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*4.0+sizeAlum*1.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*4.0+sizeAlum*1.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*3.0+sizeAlum*1.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*4.0+sizeAlum*1.0);  //bin10 
  //row5
  mFT0AFrame->AddBin(X-sizeShft,                      Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell-sizeShft,              Y+sizeCell*5.0+sizeAlum*2.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0-sizeShft,         Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*2.0-sizeShft,          Y+sizeCell*5.0+sizeAlum*2.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum-sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*3.0+sizeAlum-sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum-sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*4.0+sizeAlum-sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0+sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*7.0+sizeAlum*3.0+sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0+sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*8.0+sizeAlum*3.0+sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0+sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*9.0+sizeAlum*4.0+sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0+sizeShft,Y+sizeCell*4.0+sizeAlum*2.0,  X+sizeCell*10.0+sizeAlum*4.0+sizeShft, Y+sizeCell*5.0+sizeAlum*2.0);  //bin10 
  //row6
  mFT0AFrame->AddBin(X-sizeShft,                      Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell-sizeShft,              Y+sizeCell*6.0+sizeAlum*2.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0-sizeShft,         Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*2.0-sizeShft,          Y+sizeCell*6.0+sizeAlum*2.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum-sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*3.0+sizeAlum-sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum-sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*4.0+sizeAlum-sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0+sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*7.0+sizeAlum*3.0+sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0+sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*8.0+sizeAlum*3.0+sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0+sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*9.0+sizeAlum*4.0+sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0+sizeShft,Y+sizeCell*5.0+sizeAlum*2.0,  X+sizeCell*10.0+sizeAlum*4.0+sizeShft, Y+sizeCell*6.0+sizeAlum*2.0);  //bin10 
  //row7
  mFT0AFrame->AddBin(X,                      Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell,              Y+sizeCell*7.0+sizeAlum*3.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*2.0,          Y+sizeCell*7.0+sizeAlum*3.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*7.0+sizeAlum*3.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*7.0+sizeAlum*3.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y+sizeShft+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*7.0+sizeShft+sizeAlum*3.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y+sizeShft+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*7.0+sizeShft+sizeAlum*3.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*7.0+sizeAlum*3.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*7.0+sizeAlum*3.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*7.0+sizeAlum*3.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*6.0+sizeAlum*3.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*7.0+sizeAlum*3.0);  //bin10
  //row8
  mFT0AFrame->AddBin(X,                      Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell,              Y+sizeCell*8.0+sizeAlum*3.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*2.0,          Y+sizeCell*8.0+sizeAlum*3.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*8.0+sizeAlum*3.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*8.0+sizeAlum*3.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y+sizeShft+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*8.0+sizeShft+sizeAlum*3.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y+sizeShft+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*8.0+sizeShft+sizeAlum*3.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*8.0+sizeAlum*3.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*8.0+sizeAlum*3.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*8.0+sizeAlum*3.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*7.0+sizeAlum*3.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*8.0+sizeAlum*3.0);  //bin10
  //row9
  mFT0AFrame->AddBin(X,                      Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell,              Y+sizeCell*9.0+sizeAlum*4.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*2.0,          Y+sizeCell*9.0+sizeAlum*4.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*9.0+sizeAlum*4.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*9.0+sizeAlum*4.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y+sizeShft+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*9.0+sizeShft+sizeAlum*4.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y+sizeShft+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*9.0+sizeShft+sizeAlum*4.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*9.0+sizeAlum*4.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*9.0+sizeAlum*4.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*9.0+sizeAlum*4.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*8.0+sizeAlum*4.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*9.0+sizeAlum*4.0);  //bin10
  //row10
  mFT0AFrame->AddBin(X,                      Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell,              Y+sizeCell*10.0+sizeAlum*4.0);     //bin1
  mFT0AFrame->AddBin(X+sizeCell*1.0,         Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*2.0,          Y+sizeCell*10.0+sizeAlum*4.0); //bin2
  mFT0AFrame->AddBin(X+sizeCell*2.0+sizeAlum,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*3.0+sizeAlum, Y+sizeCell*10.0+sizeAlum*4.0);  //bin3
  mFT0AFrame->AddBin(X+sizeCell*3.0+sizeAlum,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*4.0+sizeAlum, Y+sizeCell*10.0+sizeAlum*4.0);  //bin4
  mFT0AFrame->AddBin(X+sizeCell*4.0+sizeAlum*2.0,Y+sizeShft+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*5.0+sizeAlum*2.0, Y+sizeCell*10.0+sizeShft+sizeAlum*4.0);  //bin5-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*5.0+sizeAlum*2.0,Y+sizeShft+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*6.0+sizeAlum*2.0, Y+sizeCell*10.0+sizeShft+sizeAlum*4.0);  //bin6-shifted-down
  mFT0AFrame->AddBin(X+sizeCell*6.0+sizeAlum*3.0,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*7.0+sizeAlum*3.0, Y+sizeCell*10.0+sizeAlum*4.0);  //bin7
  mFT0AFrame->AddBin(X+sizeCell*7.0+sizeAlum*3.0,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*8.0+sizeAlum*3.0, Y+sizeCell*10.0+sizeAlum*4.0);  //bin8
  mFT0AFrame->AddBin(X+sizeCell*8.0+sizeAlum*4.0,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*9.0+sizeAlum*4.0, Y+sizeCell*10.0+sizeAlum*4.0);  //bin9
  mFT0AFrame->AddBin(X+sizeCell*9.0+sizeAlum*4.0,Y+sizeCell*9.0+sizeAlum*4.0,  X+sizeCell*10.0+sizeAlum*4.0, Y+sizeCell*10.0+sizeAlum*4.0);  //bin10

  //Construct mFT0CFrame
  X=-17.46;
  Y=-17.46;
  //row1
  mFT0CFrame->AddBin(X+sizeCell*2+sizeAlum, Y,X+sizeCell*3.0+sizeAlum,Y+sizeCell);     //bin1
  mFT0CFrame->AddBin(X+sizeCell*3+sizeAlum, Y,X+sizeCell*4.0+sizeAlum,Y+sizeCell);     //bin2
  mFT0CFrame->AddBin(X+sizeCell*4+2*sizeAlum, Y,X+sizeCell*5.0+2*sizeAlum,Y+sizeCell); //bin3
  mFT0CFrame->AddBin(X+sizeCell*5+2*sizeAlum, Y,X+sizeCell*6.0+2*sizeAlum,Y+sizeCell); //bin4
  mFT0CFrame->AddBin(X+sizeCell*6+3*sizeAlum, Y,X+sizeCell*7.0+3*sizeAlum,Y+sizeCell); //bin5
  mFT0CFrame->AddBin(X+sizeCell*7+3*sizeAlum, Y,X+sizeCell*8.0+3*sizeAlum,Y+sizeCell);  //bin6
  mFT0CFrame->AddBin(X+sizeCell*8+4*sizeAlum, Y,X+sizeCell*9.0+4*sizeAlum,Y+sizeCell);  //bin7
  mFT0CFrame->AddBin(X+sizeCell*9+4*sizeAlum, Y,X+sizeCell*10.0+4*sizeAlum,Y+sizeCell);  //bin8
  //row 2
  mFT0CFrame->AddBin(X+sizeCell*2+sizeAlum, Y+sizeCell,X+sizeCell*3.0+sizeAlum,Y+sizeCell*2);     //bin1
  mFT0CFrame->AddBin(X+sizeCell*3+sizeAlum, Y+sizeCell,X+sizeCell*4.0+sizeAlum,Y+sizeCell*2);     //bin2
  mFT0CFrame->AddBin(X+sizeCell*4+2*sizeAlum, Y+sizeCell,X+sizeCell*5.0+2*sizeAlum,Y+sizeCell*2); //bin3
  mFT0CFrame->AddBin(X+sizeCell*5+2*sizeAlum, Y+sizeCell,X+sizeCell*6.0+2*sizeAlum,Y+sizeCell*2); //bin4
  mFT0CFrame->AddBin(X+sizeCell*6+3*sizeAlum, Y+sizeCell,X+sizeCell*7.0+3*sizeAlum,Y+sizeCell*2); //bin5
  mFT0CFrame->AddBin(X+sizeCell*7+3*sizeAlum, Y+sizeCell,X+sizeCell*8.0+3*sizeAlum,Y+sizeCell*2);  //bin6
  mFT0CFrame->AddBin(X+sizeCell*8+4*sizeAlum, Y+sizeCell,X+sizeCell*9.0+4*sizeAlum,Y+sizeCell*2);  //bin7
  mFT0CFrame->AddBin(X+sizeCell*9+4*sizeAlum, Y+sizeCell,X+sizeCell*10.0+4*sizeAlum,Y+sizeCell*2);  //bin8
  //row 3
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*2+sizeAlum,X+sizeCell*1.0,Y+sizeCell*3+sizeAlum);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*2+sizeAlum,X+sizeCell*2.0,Y+sizeCell*3+sizeAlum);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*2+sizeAlum,X+sizeCell*3.0+sizeAlum,Y+sizeCell*3+sizeAlum);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*2+sizeAlum,X+sizeCell*4.0+sizeAlum,Y+sizeCell*3+sizeAlum);//bin4
  mFT0CFrame->AddBin(X+sizeCell*4.0+sizeAlum*2, Y+sizeCell*2+sizeAlum,X+sizeCell*5.0+sizeAlum*2,Y+sizeCell*3+sizeAlum);//bin5
  mFT0CFrame->AddBin(X+sizeCell*5.0+sizeAlum*2, Y+sizeCell*2+sizeAlum,X+sizeCell*6.0+sizeAlum*2,Y+sizeCell*3+sizeAlum);//bin6
  mFT0CFrame->AddBin(X+sizeCell*6.0+sizeAlum*3, Y+sizeCell*2+sizeAlum,X+sizeCell*7.0+sizeAlum*3,Y+sizeCell*3+sizeAlum);//bin7
  mFT0CFrame->AddBin(X+sizeCell*7.0+sizeAlum*3, Y+sizeCell*2+sizeAlum,X+sizeCell*8.0+sizeAlum*3,Y+sizeCell*3+sizeAlum);//bin8
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*2+sizeAlum,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*3+sizeAlum);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*2+sizeAlum,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*3+sizeAlum);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*2+sizeAlum,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*3+sizeAlum);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*2+sizeAlum,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*3+sizeAlum);//bin12
  //row 4
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*3+sizeAlum,X+sizeCell*1.0,Y+sizeCell*4+sizeAlum);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*3+sizeAlum,X+sizeCell*2.0,Y+sizeCell*4+sizeAlum);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*3+sizeAlum,X+sizeCell*3.0+sizeAlum,Y+sizeCell*4+sizeAlum);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*3+sizeAlum,X+sizeCell*4.0+sizeAlum,Y+sizeCell*4+sizeAlum);//bin4
  mFT0CFrame->AddBin(X+sizeCell*4.0+sizeAlum*2, Y+sizeCell*3+sizeAlum,X+sizeCell*5.0+sizeAlum*2,Y+sizeCell*4+sizeAlum);//bin5
  mFT0CFrame->AddBin(X+sizeCell*5.0+sizeAlum*2, Y+sizeCell*3+sizeAlum,X+sizeCell*6.0+sizeAlum*2,Y+sizeCell*4+sizeAlum);//bin6
  mFT0CFrame->AddBin(X+sizeCell*6.0+sizeAlum*3, Y+sizeCell*3+sizeAlum,X+sizeCell*7.0+sizeAlum*3,Y+sizeCell*4+sizeAlum);//bin7
  mFT0CFrame->AddBin(X+sizeCell*7.0+sizeAlum*3, Y+sizeCell*3+sizeAlum,X+sizeCell*8.0+sizeAlum*3,Y+sizeCell*4+sizeAlum);//bin8
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*3+sizeAlum,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*4+sizeAlum);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*3+sizeAlum,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*4+sizeAlum);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*3+sizeAlum,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*4+sizeAlum);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*3+sizeAlum,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*4+sizeAlum);//bin12
  //row 5
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*4+sizeAlum*2,X+sizeCell*1.0,Y+sizeCell*5+sizeAlum*2);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*4+sizeAlum*2,X+sizeCell*2.0,Y+sizeCell*5+sizeAlum*2);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*4+sizeAlum*2,X+sizeCell*3.0+sizeAlum,Y+sizeCell*5+sizeAlum*2);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*4+sizeAlum*2,X+sizeCell*4.0+sizeAlum,Y+sizeCell*5+sizeAlum*2);//bin4
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*4+sizeAlum*2,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*5+sizeAlum*2);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*4+sizeAlum*2,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*5+sizeAlum*2);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*4+sizeAlum*2,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*5+sizeAlum*2);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*4+sizeAlum*2,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*5+sizeAlum*2);//bin12
  //row 6
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*5+sizeAlum*2,X+sizeCell*1.0,Y+sizeCell*6+sizeAlum*2);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*5+sizeAlum*2,X+sizeCell*2.0,Y+sizeCell*6+sizeAlum*2);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*5+sizeAlum*2,X+sizeCell*3.0+sizeAlum,Y+sizeCell*6+sizeAlum*2);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*5+sizeAlum*2,X+sizeCell*4.0+sizeAlum,Y+sizeCell*6+sizeAlum*2);//bin4
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*5+sizeAlum*2,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*6+sizeAlum*2);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*5+sizeAlum*2,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*6+sizeAlum*2);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*5+sizeAlum*2,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*6+sizeAlum*2);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*5+sizeAlum*2,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*6+sizeAlum*2);//bin12
  //row 7
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*6+sizeAlum*3,X+sizeCell*1.0,Y+sizeCell*7+sizeAlum*3);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*6+sizeAlum*3,X+sizeCell*2.0,Y+sizeCell*7+sizeAlum*3);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*6+sizeAlum*3,X+sizeCell*3.0+sizeAlum,Y+sizeCell*7+sizeAlum*3);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*6+sizeAlum*3,X+sizeCell*4.0+sizeAlum,Y+sizeCell*7+sizeAlum*3);//bin4
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*6+sizeAlum*3,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*7+sizeAlum*3);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*6+sizeAlum*3,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*7+sizeAlum*3);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*6+sizeAlum*3,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*7+sizeAlum*3);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*6+sizeAlum*3,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*7+sizeAlum*3);//bin12
  //row 8
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*7+sizeAlum*3,X+sizeCell*1.0,Y+sizeCell*8+sizeAlum*3);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*7+sizeAlum*3,X+sizeCell*2.0,Y+sizeCell*8+sizeAlum*3);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*7+sizeAlum*3,X+sizeCell*3.0+sizeAlum,Y+sizeCell*8+sizeAlum*3);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*7+sizeAlum*3,X+sizeCell*4.0+sizeAlum,Y+sizeCell*8+sizeAlum*3);//bin4
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*7+sizeAlum*3,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*8+sizeAlum*3);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*7+sizeAlum*3,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*8+sizeAlum*3);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*7+sizeAlum*3,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*8+sizeAlum*3);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*7+sizeAlum*3,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*8+sizeAlum*3);//bin12
  //row 9
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*8+sizeAlum*4,X+sizeCell*1.0,Y+sizeCell*9+sizeAlum*4);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*8+sizeAlum*4,X+sizeCell*2.0,Y+sizeCell*9+sizeAlum*4);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*8+sizeAlum*4,X+sizeCell*3.0+sizeAlum,Y+sizeCell*9+sizeAlum*4);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*8+sizeAlum*4,X+sizeCell*4.0+sizeAlum,Y+sizeCell*9+sizeAlum*4);//bin4
  mFT0CFrame->AddBin(X+sizeCell*4.0+sizeAlum*2, Y+sizeCell*8+sizeAlum*4,X+sizeCell*5.0+sizeAlum*2,Y+sizeCell*9+sizeAlum*4);//bin5
  mFT0CFrame->AddBin(X+sizeCell*5.0+sizeAlum*2, Y+sizeCell*8+sizeAlum*4,X+sizeCell*6.0+sizeAlum*2,Y+sizeCell*9+sizeAlum*4);//bin6
  mFT0CFrame->AddBin(X+sizeCell*6.0+sizeAlum*3, Y+sizeCell*8+sizeAlum*4,X+sizeCell*7.0+sizeAlum*3,Y+sizeCell*9+sizeAlum*4);//bin7
  mFT0CFrame->AddBin(X+sizeCell*7.0+sizeAlum*3, Y+sizeCell*8+sizeAlum*4,X+sizeCell*8.0+sizeAlum*3,Y+sizeCell*9+sizeAlum*4);//bin8
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*8+sizeAlum*4,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*9+sizeAlum*4);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*8+sizeAlum*4,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*9+sizeAlum*4);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*8+sizeAlum*4,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*9+sizeAlum*4);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*8+sizeAlum*4,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*9+sizeAlum*4);//bin12
  //row 10
  mFT0CFrame->AddBin(X+sizeCell*0.0, Y+sizeCell*9+sizeAlum*4,X+sizeCell*1.0,Y+sizeCell*10+sizeAlum*4);//bin1
  mFT0CFrame->AddBin(X+sizeCell*1.0, Y+sizeCell*9+sizeAlum*4,X+sizeCell*2.0,Y+sizeCell*10+sizeAlum*4);//bin2
  mFT0CFrame->AddBin(X+sizeCell*2.0+sizeAlum, Y+sizeCell*9+sizeAlum*4,X+sizeCell*3.0+sizeAlum,Y+sizeCell*10+sizeAlum*4);//bin3
  mFT0CFrame->AddBin(X+sizeCell*3.0+sizeAlum, Y+sizeCell*9+sizeAlum*4,X+sizeCell*4.0+sizeAlum,Y+sizeCell*10+sizeAlum*4);//bin4
  mFT0CFrame->AddBin(X+sizeCell*4.0+sizeAlum*2, Y+sizeCell*9+sizeAlum*4,X+sizeCell*5.0+sizeAlum*2,Y+sizeCell*10+sizeAlum*4);//bin5
  mFT0CFrame->AddBin(X+sizeCell*5.0+sizeAlum*2, Y+sizeCell*9+sizeAlum*4,X+sizeCell*6.0+sizeAlum*2,Y+sizeCell*10+sizeAlum*4);//bin6
  mFT0CFrame->AddBin(X+sizeCell*6.0+sizeAlum*3, Y+sizeCell*9+sizeAlum*4,X+sizeCell*7.0+sizeAlum*3,Y+sizeCell*10+sizeAlum*4);//bin7
  mFT0CFrame->AddBin(X+sizeCell*7.0+sizeAlum*3, Y+sizeCell*9+sizeAlum*4,X+sizeCell*8.0+sizeAlum*3,Y+sizeCell*10+sizeAlum*4);//bin8
  mFT0CFrame->AddBin(X+sizeCell*8.0+sizeAlum*4, Y+sizeCell*9+sizeAlum*4,X+sizeCell*9.0+sizeAlum*4,Y+sizeCell*10+sizeAlum*4);//bin9
  mFT0CFrame->AddBin(X+sizeCell*9.0+sizeAlum*4, Y+sizeCell*9+sizeAlum*4,X+sizeCell*10.0+sizeAlum*4,Y+sizeCell*10+sizeAlum*4);//bin10
  mFT0CFrame->AddBin(X+sizeCell*10.0+sizeAlum*5, Y+sizeCell*9+sizeAlum*4,X+sizeCell*11.0+sizeAlum*5,Y+sizeCell*10+sizeAlum*4);//bin11
  mFT0CFrame->AddBin(X+sizeCell*11.0+sizeAlum*5, Y+sizeCell*9+sizeAlum*4,X+sizeCell*12.0+sizeAlum*5,Y+sizeCell*10+sizeAlum*4);//bin12
  //row11
  mFT0CFrame->AddBin(X+sizeCell*2+sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*3.0+sizeAlum,Y+sizeCell*11+sizeAlum*5);     //bin1
  mFT0CFrame->AddBin(X+sizeCell*3+sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*4.0+sizeAlum,Y+sizeCell*11+sizeAlum*5);     //bin2
  mFT0CFrame->AddBin(X+sizeCell*4+2*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*5.0+2*sizeAlum,Y+sizeCell*11+sizeAlum*5); //bin3
  mFT0CFrame->AddBin(X+sizeCell*5+2*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*6.0+2*sizeAlum,Y+sizeCell*11+sizeAlum*5); //bin4
  mFT0CFrame->AddBin(X+sizeCell*6+3*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*7.0+3*sizeAlum,Y+sizeCell*11+sizeAlum*5); //bin5
  mFT0CFrame->AddBin(X+sizeCell*7+3*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*8.0+3*sizeAlum,Y+sizeCell*11+sizeAlum*5);  //bin6
  mFT0CFrame->AddBin(X+sizeCell*8+4*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*9.0+4*sizeAlum,Y+sizeCell*11+sizeAlum*5);  //bin7
  mFT0CFrame->AddBin(X+sizeCell*9+4*sizeAlum, Y+sizeCell*10+sizeAlum*5,X+sizeCell*10.0+4*sizeAlum,Y+sizeCell*11+sizeAlum*5);  //bin8
  //row 12
  mFT0CFrame->AddBin(X+sizeCell*2+sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*3.0+sizeAlum,Y+sizeCell*12+sizeAlum*5);     //bin1
  mFT0CFrame->AddBin(X+sizeCell*3+sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*4.0+sizeAlum,Y+sizeCell*12+sizeAlum*5);     //bin2
  mFT0CFrame->AddBin(X+sizeCell*4+2*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*5.0+2*sizeAlum,Y+sizeCell*12+sizeAlum*5); //bin3
  mFT0CFrame->AddBin(X+sizeCell*5+2*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*6.0+2*sizeAlum,Y+sizeCell*12+sizeAlum*5); //bin4
  mFT0CFrame->AddBin(X+sizeCell*6+3*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*7.0+3*sizeAlum,Y+sizeCell*12+sizeAlum*5); //bin5
  mFT0CFrame->AddBin(X+sizeCell*7+3*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*8.0+3*sizeAlum,Y+sizeCell*12+sizeAlum*5);  //bin6
  mFT0CFrame->AddBin(X+sizeCell*8+4*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*9.0+4*sizeAlum,Y+sizeCell*12+sizeAlum*5);  //bin7
  mFT0CFrame->AddBin(X+sizeCell*9+4*sizeAlum, Y+sizeCell*11+sizeAlum*5,X+sizeCell*10.0+4*sizeAlum,Y+sizeCell*12+sizeAlum*5);  //bin8
  }

} // namespace o2::quality_control_modules::ft0
