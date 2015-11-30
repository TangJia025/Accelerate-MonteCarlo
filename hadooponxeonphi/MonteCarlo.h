//
// ./copyright
//
// INTEL CONFIDENTIAL 
//
// Copyright 2011 Intel Corporation All Rights Reserved.  
//
// The source code contained or described herein and all documents related to the 
// source code ("Material") are owned by Intel Corporation or its suppliers
// or licensors. Title to the Material remains with Intel Corporation or its suppliers 
// and licensors. The Material contains trade secrets and proprietary and confidential 
// information of Intel or its suppliers and licensors. The Material is protected by 
// worldwide copyright and trade secret laws and treaty provisions. No part of the 
// Material may be used, copied, reproduced, modified, published, uploaded, posted,
// transmitted, distributed, or disclosed in any way without Intel.s prior express 
// written permission.
//
// No license under any patent, copyright, trade secret or other intellectual property 
// right is granted to or conferred upon you by disclosure or delivery of the Materials, 
// either expressly, by implication, inducement, estoppel or otherwise. Any license under 
// such intellectual property rights must be express and approved by Intel in writing.
//
//
//
#include "omp.h"

#ifdef __INTEL_OFFLOAD
#pragma omp declare target
#endif
static const int        OPT_N = 512*64*60;
static const int       RAND_N = 262144;
static const int       RAND_BLOCK_LENGTH = 16*1024;
static const float   RISKFREE = 0.06;
static const float VOLATILITY = 0.10;
#ifdef __INTEL_OFFLOAD
#pragma omp end declare target
#endif

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif
EXTERNC double MonteCarlo(unsigned int);
