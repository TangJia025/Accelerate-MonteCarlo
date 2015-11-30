package com.intel.demo;

import java.io.IOException;

import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

public class MonteCarloReducer extends Reducer<Text, DoubleWritable, NullWritable, DoubleWritable> {
  
  protected void setup(Context context) {
    
  }
  
  protected void reduce(Text key, Iterable<DoubleWritable> values, Context context) 
      throws IOException, InterruptedException{
    int counter = 0;
    double sum = 0;
    for(DoubleWritable writable:values){
      counter++;
      sum += writable.get();
    }
    
    context.write(NullWritable.get(), new DoubleWritable(sum/counter));
  }

}
