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
#include <string.h>
#include <math.h>
#include "MonteCarlo.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    static int verbose = 1; 
    const int RAND_N = 1 << 18;

    if (argc != 1)
    {
       	printf("usage: MonteCarlo <verbose>  where verbose = 1 for validtating result, the default is not to validate result. \n");
        	exit(1);
    }
    if (argc == 1)
	verbose = 0;
    if (argc == 2)
	verbose = atoi(argv[1]);

    MonteCarlo();

    return 0;
}
