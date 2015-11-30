package com.intel.demo;

import java.io.IOException;

import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MonteCarloMapper extends Mapper<LongWritable, Text, Text, DoubleWritable> {
  private static final Logger LOG = LoggerFactory.getLogger(MonteCarloMapper.class);
  static {
    LOG.info("env: " + System.getenv());
    LOG.info("java.library.path: " + System.getProperty("java.library.path"));
    System.loadLibrary("MonteCarloMapper");
  }
  protected void setup(Context context){
    LOG.info("java.library.path: " + System.getProperty("java.library.path"));
  }
  
  protected void map(LongWritable key, Text value, Context context) 
      throws IOException, InterruptedException{
    int rndId = Integer.parseInt(value.toString());
    double result = monteCarlo(rndId);
    System.out.println(result);
    context.write(new Text("MonteCarlo"), new DoubleWritable(result));
  }
  
  public native double monteCarlo(int rid);
  
  public static void main(String[] args){
    System.out.println(new MonteCarloMapper().monteCarlo(1));
  }
}
