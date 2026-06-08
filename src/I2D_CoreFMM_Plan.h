/*
 *  I2D_CoreFMM_Plan.h
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

extern double _THETA;
#define _FMMSILENT

#include "I2D_FMMTypes.h"

// #include "I2D_CoreFMM_AggressiveExpansion.h"
#include "VortexExpansions.h"
#include "hcfmm_box.h"
#include "hcfmm_boxBuilder_serial.h"
#ifndef _MRAG_TBB
#include "hcfmm_evaluator_serial.h"
#else
#include "hcfmm_evaluator_tbb.h"
#endif

class Plan {
public:
  typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,
                     _FMM_MAX_LEVEL_>
      tBox;

  virtual void addDirectInteraction(const tBox *const, int, bool) = 0;
  virtual void addIndirectInteraction(const tBox *const, int) = 0;
  virtual void merge_direct_intervals(int) = 0;
};
