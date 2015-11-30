// ./copyright
//
// INTEL CONFIDENTIAL 
//
// Copyright 2013 Intel Corporation All Rights Reserved.  
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
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "MonteCarlo.h"
#include <immintrin.h>
#include "omp.h"
#include "math.h"
#define RANDSEED 123
#define SIMDALIGN 1024

#ifdef __INTEL_OFFLOAD
#define EXP expf
#else
#define EXP exp_sp_ep
#endif 

#ifdef __INTEL_OFFLOAD
#pragma omp declare target
#endif
static const float  RVV = RISKFREE-0.5f*VOLATILITY*VOLATILITY;
static const float  INV_RAND_N = 1.0f/RAND_N;
static const float  F_RAND_N = static_cast<float>(RAND_N);
static const float STDDEV_DENOM = 1 / (F_RAND_N * (F_RAND_N - 1.0f));
static const float CONFIDENCE_DENOM = 1 / sqrtf(F_RAND_N);
static const int NBLOCKS = RAND_N/RAND_BLOCK_LENGTH;
#ifdef __INTEL_OFFLOAD
#pragma omp end declare target
#endif

double second()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

inline float RandFloat(float low, float high, unsigned int * seed){
    float t = (float)rand_r(seed) / (float)RAND_MAX;
    return (1.0f - t) * low + t * high;
}

///////////////////////////////////////////////////////////////////////////////
// Polynomial approximation of cumulative normal distribution function
///////////////////////////////////////////////////////////////////////////////

double CND(double d){
    const double       A1 = 0.31938153;
    const double       A2 = -0.356563782;
    const double       A3 = 1.781477937;
    const double       A4 = -1.821255978;
    const double       A5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;

    double
        K = 1.0 / (1.0 + 0.2316419 * fabs(d));

    double
        cnd = RSQRT2PI * exp(- 0.5 * d * d) *
        (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));

    if(d > 0)
        cnd = 1.0 - cnd;

    return cnd;
}

void BlackScholesFormula(
    double& callResult,
    double Sf, //Stock price
    double Xf, //Option strike
    double Tf, //Option years
    double Rf, //Riskless rate
    double Vf  //Volatility rate
){
    double S = Sf, X = Xf, T = Tf, R = Rf, V = Vf;

    double sqrtT = sqrt(T);
    double    d1 = (log(S / X) + (R + 0.5 * V * V) * T) / (V * sqrtT);
    double    d2 = d1 - V * sqrtT;
    double CNDD1 = CND(d1);
    double CNDD2 = CND(d2);

    double expRT = exp(- R * T);
    callResult   = (S * CNDD1 - X * expRT * CNDD2);
}

#ifndef __INTEL_OFFLOAD
/*
//               INTEL CORPORATION PROPRIETARY INFORMATION
//  This software is supplied under the terms of a license agreement or
//  nondisclosure agreement with Intel Corporation and may not be copied
//  or disclosed except in accordance with the terms of that agreement.
//    Copyright (c) 1996-2013 Intel Corporation. All Rights Reserved.
//
*/





/*
 * Support for casting between various INT, FP types.
 * Note that these do no conversion of values, they just change
 * the type.
 */


    typedef          __int16  VINT16;
    typedef   signed __int16 VSINT16;
    typedef unsigned __int16 VUINT16;

    typedef          __int32  VINT32;
    typedef   signed __int32 VSINT32;
    typedef unsigned __int32 VUINT32;

    typedef          __int64  VINT64;
    typedef   signed __int64 VSINT64;
    typedef unsigned __int64 VUINT64;

#pragma warning( disable : 2571 )
typedef struct
{
        __declspec(align(64)) VUINT32 _sInvLn2[1];
        __declspec(align(64)) VUINT32 _sShifter[1];
        __declspec(align(64)) VUINT32 _sLn2[1];
        __declspec(align(64)) VUINT32 _iBias[1];
        __declspec(align(64)) VUINT32 _sPC0[1];
        __declspec(align(64)) VUINT32 _sPC1[1];
        __declspec(align(64)) VUINT32 _sPC2[1];
        __declspec(align(64)) VUINT32 _sPC3[1];
        __declspec(align(64)) VUINT32 _sPC0_fma[1];
        __declspec(align(64)) VUINT32 _sPC1_fma[1];
        __declspec(align(64)) VUINT32 _sPC2_fma[1];
        __declspec(align(64)) VUINT32 _sPC3_fma[1];
    __declspec(align(64)) VUINT32 _iAbsMask[1];
    __declspec(align(64)) VUINT32 _iDomainRange[1];
} sExp_Table_Type;




static __declspec(align(64)) const sExp_Table_Type vsexp_ep_data =
{
        { (VUINT32)(0x3FB8AA3Bu) }, /* _sInvLn2  */  //k=0
        { (VUINT32)(0x4b400000u) }, /* _sShifter */
        { (VUINT32)(0x3F317218u) }, /* _sLn2 */
        { (VUINT32)(0x0000007fu) }, /* _iBias */
        // here we approximate 2^x on [-0.5, 0.5], leading coef. is set to 1.
        // this is only beneficial in the absence of FMA instruction
        { (VUINT32)(0x3F800000u) },  /* _sPC0  */
        { (VUINT32)(0x3f317422u) },  /* _sPC1  */
        { (VUINT32)(0x3e77d66au) },  /* _sPC2  */
        { (VUINT32)(0x3d63582bu) },  /* _sPC3  */
        // these coeffs are used when FMA is available
        { (VUINT32)(0x3F800000u) },  /* _sPC0_fma  */
        { (VUINT32)(0x3F8003DEu) },  /* _sPC1_fma  */
        { (VUINT32)(0x3F00F2D6u) },  /* _sPC2_fma  */
        { (VUINT32)(0x3E2963ACu) },  /* _sPC3_fma  */

    { (VUINT32)(0x7fffffffu) },         /* _iAbsMask */
    { (VUINT32)(0x42aeac4fu) },         /* _iDomainRange */
}; /*sExp_Table*/


__forceinline
float exp_sp_ep (float x)
{
        /* temporary VLANG variable declarations */
        __declspec(align(64)) float va1;
        __declspec(align(64)) float vr1;
        __declspec(align(64)) VUINT32 vm;

        /* load arguments into temporary VLANG valiables */
        va1 = x;
        /* call VLANG core algorithm implementation */
        {

    __declspec(align(64)) float sN;         __declspec(align(64)) float sM;
    __declspec(align(64)) float sR;
    __declspec(align(64)) float sP;
    __declspec(align(64)) float s2N;
    __declspec(align(64)) VUINT32 iX;
    __declspec(align(64)) VUINT32 iAbsX;
    __declspec(align(64)) VUINT32 iRangeMask;
    __declspec(align(64)) VUINT32 iRes;
    __declspec(align(64)) VUINT32 iM;
    __declspec(align(64)) VUINT32 iP;

    /* Constants */
    __declspec(align(64)) float sInvLn2;       __declspec(align(64)) float sShifter;
    __declspec(align(64)) VUINT32 iBias;
    __declspec(align(64)) VUINT32 iAbsMask;
    __declspec(align(64)) VUINT32 iDomainRange;

        __declspec(align(64)) float sPC[4];
        __declspec(align(64)) float sPC_fma[4];
        __declspec(align(64)) float sLn2;

    sInvLn2 = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sInvLn2));
    sShifter = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sShifter));
        sLn2 = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sLn2));
    iBias = (*(const VUINT32 *)&vsexp_ep_data._iBias);
    iAbsMask = (*(const VUINT32 *)&vsexp_ep_data._iAbsMask);
    iDomainRange = (*(const VUINT32 *)&vsexp_ep_data._iDomainRange);
        sPC[0] = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sPC0));
        sPC[1] = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sPC1));
        sPC[2] = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sPC2));
        sPC[3] = _castu32_f32((*(const VUINT32 *)&vsexp_ep_data._sPC3));
    /* -------------------- Implementation --------------------- */
        // we save one operation in case FMA is not available
        sR = __fence( va1 * sInvLn2 );
        sM = __fence( sR + sShifter );
    sN = __fence( sM - sShifter );
    /* ...............Check for overflow\underflow ............. */
    iX = _castf32_u32(va1);
    iAbsX = __fence( iX & iAbsMask );
    iRangeMask = __fence( (VUINT32)( -(VSINT32)((VSINT32)iAbsX > (VSINT32)iDomainRange) ) );
    vm = 0; vm = iRangeMask;
    iM = _castf32_u32(sM);
    /* ............... 2^N ..................................... */
        iM = __fence( (VUINT32)(iM) << (23) );
    /* ................... R ................................... */
        // we save one operation in case FMA is not available
        // notice we initialized sR in above code under the same
        // ifdef guards.
        sR = __fence( sR - sN );

    /* ................... Polynomial .......................... */
    /* ---------------ULP=1.0-------------*/
        sP = __fence( __fence( sPC[3] * sR ) + sPC[2] );
        sP = __fence( __fence( sP * sR ) + sPC[1] );
        sP = __fence( __fence( sP * sR ) + sPC[0] );
    /* ................... Reconstruction ...................... */
            iP = _castf32_u32(sP);
            iRes = __fence( iM + iP );
            vr1 = _castu32_f32(iRes);
        }
	return vr1;
}
#endif

double MonteCarlo(unsigned int rid)
{
    float 
        *CallResultParallel,
        *CallConfidence,
        *StockPrice,
        *OptionStrike,
        *OptionYears;
    double
	sTime, eTime;
    int 
        mem_size, rand_size; 

    static int verbose = 1;

    mem_size = sizeof(float)*OPT_N;
    rand_size = sizeof(float)*RAND_N;

    printf("Monte Carlo European Option Pricing in Single Precision\n");	

    posix_memalign((void **)&CallResultParallel, SIMDALIGN, mem_size);
    posix_memalign((void **)&CallConfidence, SIMDALIGN, mem_size);
    posix_memalign((void **)&StockPrice, SIMDALIGN, mem_size);
    posix_memalign((void **)&OptionStrike, SIMDALIGN, mem_size);
    int res = posix_memalign((void **)&OptionYears, SIMDALIGN, mem_size);
    if (res != 0)
       printf("malloc failed\n");	

    if (verbose)
    {
        printf("...generating the input data.\n");
    }
    for(int i = 0; i < OPT_N; i++)
    {
        CallResultParallel[i] = 0.0;
        CallConfidence[i]= -1.0;
        StockPrice[i]    = RandFloat(5.0f, 50.0f, &rid);
        OptionStrike[i]  = RandFloat(10.0f, 25.0f, &rid);
        OptionYears[i]   = RandFloat(1.0f, 5.0f, &rid);
    }
    printf("Pricing %d options with path length of %d.\n", OPT_N, RAND_N);
    sTime = second();
#ifdef __INTEL_OFFLOAD
#pragma omp target device(0) map(to:StockPrice[0:OPT_N], OptionStrike[0:OPT_N], OptionYears[0:OPT_N])\
				map(from:CallResultParallel[0:OPT_N], CallConfidence[0:OPT_N]) 
#endif
{
    __attribute__((align(1024))) float random [RAND_N];

#pragma omp parallel for
    for (int k = 0; k < RAND_N; k++)
        random[k] = cdfnorminv((k+1.0)/(RAND_N+1.0));

#pragma omp parallel for
    for(int opt = 0; opt < OPT_N; opt++)
    {
        CallResultParallel[opt]     = 0;
        CallConfidence[opt] = 0;
    }

    for(int block = 0; block < NBLOCKS; ++block)
    {
        const float *randblock = random + block * RAND_BLOCK_LENGTH;
#pragma omp parallel for
    for(int opt = 0; opt < OPT_N; opt++) 
    {
        float VBySqrtT = VOLATILITY * sqrtf(OptionYears[opt]);
        float MuByT = RVV * OptionYears[opt];
        float Sval = StockPrice[opt];
        float Xval = OptionStrike[opt];
        float val = 0.0, val2 = 0.0;
#pragma vector aligned
#pragma simd reduction(+:val) reduction(+:val2)
#pragma unroll(4) 
        for(int pos = 0; pos < RAND_BLOCK_LENGTH; pos++)
        {
            float callValue  = Sval * EXP(MuByT + VBySqrtT * randblock[pos]) - Xval;
            callValue = (callValue > 0) ? callValue : 0;
            val  += callValue;
            val2 += callValue * callValue;
        }
        CallResultParallel[opt]  +=  val;
        CallConfidence[opt]      +=  val2;
     }
    }
#pragma omp parallel for
    for(int opt = 0; opt < OPT_N; opt++)
    {
        const float val      = CallResultParallel[opt];
        const float val2     = CallConfidence[opt];
        float exprt = expf(-RISKFREE *OptionYears[opt]);
        CallResultParallel[opt] = exprt * val * INV_RAND_N;
        float stdDev = sqrtf((F_RAND_N * val2 - val * val)* STDDEV_DENOM);
        CallConfidence[opt] = (exprt * 1.96f * stdDev * CONFIDENCE_DENOM);
    }
}
    eTime = second();
    printf("Completed in %8.4f seconds.\n",eTime-sTime );
    printf("Computation rate - %8.3f options per second.\n", OPT_N/(eTime-sTime));

    double
        delta, sum_delta, sum_ref, L1norm, sumReserve;
    double CallMaster;

    if (verbose)
    {

        sum_delta = 0;
        sum_ref   = 0;
        sumReserve = 0;

        for(int i = 0; i < OPT_N; i++)
        {
            BlackScholesFormula(CallMaster,
                (double) StockPrice[i], 
                (double) OptionStrike[i],
                (double) OptionYears[i],
                (double) RISKFREE,
                (double) VOLATILITY);
            delta = fabs(CallMaster - CallResultParallel[i]);
            sum_delta += delta;
            sum_ref   += fabs(CallMaster);
            if(delta > 1e-6) 
                sumReserve += CallConfidence[i] / delta;
        }
        sumReserve /= (double)OPT_N;
        L1norm = sum_delta / sum_ref;
        printf("L1 norm: %E\n", L1norm);
        printf("Average reserve: %f\n", sumReserve);

        printf("...freeing CPU memory.\n");
        printf((sumReserve > 1.0f) ? "PASSED\n" : "FAILED\n");
    }
    free(CallResultParallel);
    free(CallConfidence);
    free(StockPrice);
    free(OptionStrike);
    free(OptionYears);

    return sumReserve;
}
