package weka.core.labelranking;

import java.util.HashMap;

import weka.core.DenseInstance;

/**
 * @author senge
 *
 */
public class PreferenceDenseInstance extends DenseInstance {

	private static final long serialVersionUID = 2170442456435172975L;

	/**
	 * Contains all the preference informations for the preference attributes.
	 * key: attribute index (attribute is supposed to be of type PreferenceAttribute 
	 */

	private HashMap<Integer, int[][]> m_preferences = null;
	private HashMap<Integer,int[]> m_predictions = null;
	
	//rankingString is used for the edit area inside of the PreprocessPanel.
	String rankingString = "";
	
	
	/**
	 * Creates a new PreferenceDenseInstance from an existing one.
	 * @param instance
	 */
	public PreferenceDenseInstance(DenseInstance instance) {
		super(instance);
		if(instance instanceof PreferenceDenseInstance) {
			this.m_preferences = (HashMap<Integer, int[][]>)
			((PreferenceDenseInstance)instance).m_preferences.clone();
			m_predictions = new HashMap<Integer,int[]>();
			setRankingString(((PreferenceDenseInstance)instance).getRankingString());
		}
	}

	/**
	 * Creates a preference instance with weight, values and preferences.
	 * @param weight
	 * @param values
	 * @param preferences
	 */
	public PreferenceDenseInstance(double weight, double[] values, HashMap<Integer, int[][]> preferences){
		super(weight, values);
		this.m_preferences = preferences;
		this.m_predictions = new HashMap<Integer,int[]>();
	}

	/**
	 * Returns the pairwise preferences for the given attribute.
	 * @param attributeIndex
	 * @return
	 */
	public int[][] getPairwisePreferences(int attributeIndex) {
		if(this.m_preferences.containsKey(attributeIndex)) {
			return m_preferences.get(attributeIndex);
		} else {
			throw new RuntimeException("There is no preference attribute for the index: " + attributeIndex);
		}
	}

	/**
	 * Returns 1, iff label1 is preferred to label2,
	 * returns -1, iff label2 is preferred to label1 and
	 * return 0, if the labels are unrelated.
	 * @param attributeIndex
	 * @param label1
	 * @param label2
	 * @return
	 */
	public int preference(int attributeIndex, int label1, int label2) {
		if(this.m_preferences.containsKey(attributeIndex)) {
			int[][] prefs = m_preferences.get(attributeIndex);
			for (int i = 0; i < prefs.length; i++) {
				if(prefs[i][0] == label1 && prefs[i][1] == label2) {
					return 1;
				} else if(prefs[i][0] == label2 && prefs[i][1] == label1) {
					return -1;
				}
			}
			return 0;
		} else {
			throw new RuntimeException("There is no preference attribute for the index: " + attributeIndex);
		}
	}

	/**
	 * Returns a HashMap with a ranking for a specific PreferenceDenseInstance.
	 * @return
	 */
	public HashMap<Integer,int[][]> getHashMap(){
		return m_preferences;
	}
	
	/**
	 * Setting the HashMap with stored preferences in
	 * adjacent matrix format.
	 * @param input
	 */
	public void setHashMap(HashMap<Integer,int[][]> input){
		m_preferences = input;
	}
	
	/**
	 * Returns a HashMap with a predicted ranking for a specific PreferenceDenseInstance.
	 * @return
	 */
	public HashMap<Integer,int[]> getPrediction(){
		return m_predictions;
	}
	
	/**
	 * Inserts a prediction in a given HashMap.
	 * @param i
	 * @param rnk
	 */
	public void insertPrediction (Integer i,int[] rnk){
		m_predictions.put(i, rnk);
	}
	
	/**
	 * Overwrites the HashMap with predicted ranking of a PreferenceDenseInstance.
	 * @param input
	 */
	public void setPredictions(HashMap<Integer,int[]> input){
		m_predictions = input;
	}

	
	@Override
	public Object copy() {
		PreferenceDenseInstance result = new PreferenceDenseInstance(this);
	    result.m_Dataset = m_Dataset;
	    result.m_predictions = m_predictions;
	    return result;
	}
	
	/**
	 * Setting the ranking String used in the edit view of the PreprocessPanel.
	 * @param input
	 */
	public void setRankingString(String input){
		rankingString = input;
	}
	
	/**
	 * Returns the rankingString to the edit field inside of the PreprocessPanel.
	 * @return
	 */
	public String getRankingString(){
		return rankingString;
	}


}
