package weka.classifiers;

import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.RankUtilities;


/**
 * The abstract Ranker class used for ranking problems.
 * 
 *
 */
public abstract class AbstractRanker extends AbstractClassifier {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3887189196758367179L;

	@Override
	public void buildClassifier(Instances data) throws Exception{}
		
	/* 	
    @Override
	public double classifyInstance(Instance instance){
    	RankDecoder rd = new RankDecoder();
    	rd.setLabels(LabelList.getLabels());
    	double labels=0;
    	
    	labels = instance.value(instance.numValues());
		return labels;
	}*/
    /*
	public double[] distributionForInstance(Instance inst){
		 return null;
	 }
	*/
}
