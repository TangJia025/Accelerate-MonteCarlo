#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jni.h>
#include "MonteCarloDemo.h"
extern void MonteCarlo(void);

JNIEXPORT jstring JNICALL Java_MonteCarloDemo_nativeMonteCarlo(JNIEnv *env, jobject obj)
{
  int i;
  int ds_ret;
  char *newstring;

  MonteCarlo();
  jstring ret = 0;

  newstring = (char *)malloc(30);

  if(newstring == NULL) 
  {      
      return ret;
  }

  memset(newstring, 0, 30);  
  char * tmpstring = "Test program of JNI.\n";

  strcpy(newstring, tmpstring);

  ret = (*env)->NewStringUTF(env, newstring);

  free(newstring);

  return ret;
}
