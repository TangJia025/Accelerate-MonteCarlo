package com.intel.demo;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.UUID;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.thrift.generated.IllegalArgument;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.NLineInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Driver implements Tool{
  private static final Logger LOG = LoggerFactory.getLogger(Driver.class);
  
  private static final String JOB_NAME_PREFIX = "MONTECARLO-";
  private static final String MONTECARLO_FOLDER = "MONTECARLO";
  private static final String MONTECARLO_INPUT = "INPUT";
  
  private Configuration conf;
  
  public static void main(String[] args) throws Exception{
    ToolRunner.run(new Driver(), args);
  }
  
  public void start(int mapNum){
    try {
      FileSystem hdfs = FileSystem.get(conf);
      Path homeDir = hdfs.getHomeDirectory();
      LOG.info("HomeDir: " + homeDir.makeQualified(hdfs));
      hdfs.delete(new Path(homeDir, MONTECARLO_INPUT), true);
      Path inputPath = new Path(homeDir, MONTECARLO_FOLDER + "/" + MONTECARLO_INPUT);
      FSDataOutputStream outputStream = hdfs.create(inputPath);
      PrintWriter printWriter = new PrintWriter(outputStream);
      for(int i=1;i<=mapNum;i++){
        printWriter.println(i);
      }
      printWriter.close();
      
      conf.set("mapred.task.timeout", "0");
      conf.set("mapred.map.tasks.speculative.execution", "false");
      conf.set("mapred.reduce.tasks.speculative.execution", "false");
      Job job = new Job(conf, getJobName());
      job.setJarByClass(MonteCarloMapper.class);
      NLineInputFormat.setInputPaths(job, inputPath);
      NLineInputFormat.setNumLinesPerSplit(job, 1);
      job.setInputFormatClass(NLineInputFormat.class);
      job.setMapperClass(MonteCarloMapper.class);
      job.setReducerClass(MonteCarloReducer.class);
      job.setNumReduceTasks(1);
      hdfs.delete(new Path("/MonteCarlo"), true);
      FileOutputFormat.setOutputPath(job, new Path("/MonteCarlo"));
      job.setOutputKeyClass(Text.class);
      job.setOutputValueClass(DoubleWritable.class);
      
      job.waitForCompletion(true);
    } catch (IOException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
    }
  }
  
  public String getJobName(){
    return JOB_NAME_PREFIX + UUID.randomUUID().toString();
  }

  @Override
  public Configuration getConf() {
    return null;
  }

  @Override
  public void setConf(Configuration arg0) {
    conf = arg0;
  }

  @Override
  public int run(String[] remainingArgs) throws Exception {
    int mapNum = 0;
    try{
      if(remainingArgs.length<1){
        throw new IllegalArgument();
      }
      
      mapNum = Integer.parseInt(remainingArgs[0]);
    }catch(Exception e){
      System.err.println("Usage: Driver <mapNum>");
      System.exit(1);
    }
    
    start(mapNum);
    return 0;
  }
}
