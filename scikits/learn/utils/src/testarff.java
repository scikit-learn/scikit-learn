import java.io.BufferedReader;
import java.io.FileReader;

import weka.core.converters.ConverterUtils.DataSource;
import weka.core.Instances;
import weka.core.Attribute;
import weka.core.AttributeStats;

class testarff {
        public static void main(String[] args)  {
                if (args.length < 1) {
                        System.out.println("usage is : testarff file.arff");
                        System.exit(-1);
                }
                String filename = args[0];

                try {
                        DataSource source = new DataSource(filename);
                        Instances data = source.getDataSet();
                        //System.out.println("data has "+ data.numAttributes() + " attributes");
                        System.out.println(data.numAttributes());

                        // Cache the stats of every attribute of the dataset
                        AttributeStats[] stats = new AttributeStats [data.numAttributes()];
                        for (int i = 0; i < data.numAttributes(); ++i) {
                                stats[i] = data.attributeStats(i);
                        }

                        //System.out.println(data.numInstances() + " instance");
                        System.out.println(data.numInstances());
                        // Get the name of every attribute of the dataset
                        for (int i = 0; i < data.numAttributes(); ++i) {
                                //System.out.println("Attribute " + i + " is " + data.attribute(i).name());
                                System.out.print(data.attribute(i).name());
                                if (stats[i].numericStats != null) {
                                        System.out.print("," + stats[i].numericStats.min);
                                        System.out.print("," + stats[i].numericStats.max);
                                        System.out.print("," + stats[i].numericStats.mean);
                                        System.out.println("," + stats[i].numericStats.stdDev);
                                } else {
                                        if (stats[i].nominalCounts.length > 0) {
                                                System.out.print(",{" + data.attribute(i).value(0));
                                        }
                                        for (int j = 1; j < stats[i].nominalCounts.length; ++j) {
                                                System.out.print("," + data.attribute(i).value(j));
                                        }
                                        System.out.print("}");
                                        System.out.println();
                                }
                        }
                } catch (Exception e) {
                        System.out.println("caught exception: was " + e);
                }
        }
};
