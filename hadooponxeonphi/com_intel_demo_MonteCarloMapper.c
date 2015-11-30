#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jni.h>
#include "com_intel_demo_MonteCarloMapper.h"

extern double MonteCarlo(int);

JNIEXPORT jdouble JNICALL Java_com_intel_demo_MonteCarloMapper_monteCarlo
  (JNIEnv *env, jobject obj, jint jid)
{
  int id = jid;

  unsigned int uid = (unsigned int)id;

  double ret = MonteCarlo(id);
  
  return ret;
 }
