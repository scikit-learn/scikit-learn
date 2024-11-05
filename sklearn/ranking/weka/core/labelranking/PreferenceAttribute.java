/**
 * The attribute class suited to ranking attributes.
 */

package weka.core.labelranking;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import weka.core.Attribute;
import weka.core.Copyable;
import weka.core.RevisionHandler;

public class PreferenceAttribute extends Attribute 
implements Copyable, Serializable, RevisionHandler {

	

	/** for serialization */
	static final long serialVersionUID = -742180568732916383L;

	/** Constant set for ranking attributes. */
	public static final int RANKING = 5;

	/** The keyword used to denote a ranking attribute */
	public final static String ARFF_ATTRIBUTE_RANKING = "ranking";


	public PreferenceAttribute(String attributeName, List<String> values) {
		super(attributeName, values);
	}


	public PreferenceAttribute(String attributeName, ArrayList<String> attributeValues, int size) {
		super(attributeName,attributeValues,size);
	}

	
	/**
	 * Altered value function from the Attribute class.
	 * Here inputs are decoded. Thus label rankings are returned.
	 * @param valIndex
	 * @return A string including the label ranking of one instance.
	 
	
	@Override
	public String value(int valIndex) {
		
        Object val;
        ArrayList<Object>m_Values = getM_values();
		RankDecoder rnkDec = new RankDecoder();
		rnkDec.setLabels(labels);				
		 
	    	  try{
	    		 val = rnkDec.toString(rnkDec.decode((Double)instances.get(valIndex)));
	    	  }
	    	  catch(Exception e){
	    		val = m_Values.get(valIndex);
	    	  }
	    	  
	    	  return (String) val;
	      
	}*/

}
