/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    BVDecomposeSegCVSub.java
 *    Copyright (C) 2003 Paul Conilione
 *
 *    Based on the class: BVDecompose.java by Len Trigg (1999)
 */


/*
 *    DEDICATION
 *
 *    Paul Conilione would like to express his deep gratitude and appreciation
 *    to his Chinese Buddhist Taoist Master Sifu Chow Yuk Nen for the abilities
 *    and insight that he has been taught, which have allowed him to program in
 *    a clear and efficient manner.
 *
 *    Master Sifu Chow Yuk Nen's Teachings are unique and precious. They are
 *    applicable to any field of human endeavour. Through his unique and powerful
 *    ability to skilfully apply Chinese Buddhist Teachings, people have achieved
 *    success in; Computing, chemical engineering, business, accounting, philosophy
 *    and more.
 *
 */

package weka.classifiers;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.labelranking.PreferenceAttribute;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * This class performs Bias-Variance decomposion on any classifier using the sub-sampled cross-validation procedure as specified in (1).<br/>
 * The Kohavi and Wolpert definition of bias and variance is specified in (2).<br/>
 * The Webb definition of bias and variance is specified in (3).<br/>
 * <br/>
 * Geoffrey I. Webb, Paul Conilione (2002). Estimating bias and variance from data. School of Computer Science and Software Engineering, Victoria, Australia.<br/>
 * <br/>
 * Ron Kohavi, David H. Wolpert: Bias Plus Variance Decomposition for Zero-One Loss Functions. In: Machine Learning: Proceedings of the Thirteenth International Conference, 275-283, 1996.<br/>
 * <br/>
 * Geoffrey I. Webb (2000). MultiBoosting: A Technique for Combining Boosting and Wagging. Machine Learning. 40(2):159-196.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;misc{Webb2002,
 *    address = {School of Computer Science and Software Engineering, Victoria, Australia},
 *    author = {Geoffrey I. Webb and Paul Conilione},
 *    institution = {Monash University},
 *    title = {Estimating bias and variance from data},
 *    year = {2002},
 *    PDF = {http://www.csse.monash.edu.au/\~webb/Files/WebbConilione04.pdf}
 * }
 *
 * &#64;inproceedings{Kohavi1996,
 *    author = {Ron Kohavi and David H. Wolpert},
 *    booktitle = {Machine Learning: Proceedings of the Thirteenth International Conference},
 *    editor = {Lorenza Saitta},
 *    pages = {275-283},
 *    publisher = {Morgan Kaufmann},
 *    title = {Bias Plus Variance Decomposition for Zero-One Loss Functions},
 *    year = {1996},
 *    PS = {http://robotics.stanford.edu/\~ronnyk/biasVar.ps}
 * }
 *
 * &#64;article{Webb2000,
 *    author = {Geoffrey I. Webb},
 *    journal = {Machine Learning},
 *    number = {2},
 *    pages = {159-196},
 *    title = {MultiBoosting: A Technique for Combining Boosting and Wagging},
 *    volume = {40},
 *    year = {2000}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 *
 * <pre> -c &lt;class index&gt;
 *  The index of the class attribute.
 *  (default last)</pre>
 *
 * <pre> -D
 *  Turn on debugging output.</pre>
 *
 * <pre> -l &lt;num&gt;
 *  The number of times each instance is classified.
 *  (default 10)</pre>
 *
 * <pre> -p &lt;proportion of objects in common&gt;
 *  The average proportion of instances common between any two training sets</pre>
 *
 * <pre> -s &lt;seed&gt;
 *  The random number seed used.</pre>
 *
 * <pre> -t &lt;name of arff file&gt;
 *  The name of the arff file used for the decomposition.</pre>
 *
 * <pre> -T &lt;number of instances in training set&gt;
 *  The number of instances in the training set.</pre>
 *
 * <pre> -W &lt;classifier class name&gt;
 *  Full class name of the learner used in the decomposition.
 *  eg: weka.classifiers.bayes.NaiveBayes</pre>
 *
 * <pre>
 * Options specific to learner weka.classifiers.rules.ZeroR:
 * </pre>
 *
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 *
 <!-- options-end -->
 *
 * Options after -- are passed to the designated sub-learner. <p>
 *
 * @author Paul Conilione (paulc4321@yahoo.com.au)
 * @version $Revision: 6041 $
 */
public class BVDecomposeSegCVSub
    implements OptionHandler, TechnicalInformationHandler, RevisionHandler {

    /** Debugging mode, gives extra output if true. */
    protected boolean m_Debug;

    /** An instantiated base classifier used for getting and testing options. */
    protected Classifier m_Classifier = new weka.classifiers.rules.ZeroR();

    /** The options to be passed to the base classifier. */
    protected String [] m_ClassifierOptions;

    /** The number of times an instance is classified*/
    protected int m_ClassifyIterations;

    /** The name of the data file used for the decomposition */
    protected String m_DataFileName;

    /** The index of the class attribute */
    protected int m_ClassIndex = -1;

    /** The random number seed */
    protected int m_Seed = 1;

    /** The calculated Kohavi & Wolpert bias (squared) */
    protected double m_KWBias;

    /** The calculated Kohavi & Wolpert variance */
    protected double m_KWVariance;

    /** The calculated Kohavi & Wolpert sigma */
    protected double m_KWSigma;

    /** The calculated Webb bias */
    protected double m_WBias;

    /** The calculated Webb variance */
    protected double m_WVariance;

    /** The error rate */
    protected double m_Error;

    /** The training set size */
    protected int m_TrainSize;

    /** Proportion of instances common between any two training sets. */
    protected double m_P;

    /**
     * Returns a string describing this object
     * @return a description of the classifier suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
      return
          "This class performs Bias-Variance decomposion on any classifier using the "
        + "sub-sampled cross-validation procedure as specified in (1).\n"
        + "The Kohavi and Wolpert definition of bias and variance is specified in (2).\n"
        + "The Webb definition of bias and variance is specified in (3).\n\n"
        + getTechnicalInformation().toString();
    }

    /**
     * Returns an instance of a TechnicalInformation object, containing
     * detailed information about the technical background of this class,
     * e.g., paper reference or book this class is based on.
     *
     * @return the technical information about this class
     */
    public TechnicalInformation getTechnicalInformation() {
      TechnicalInformation 	result;
      TechnicalInformation 	additional;

      result = new TechnicalInformation(Type.MISC);
      result.setValue(Field.AUTHOR, "Geoffrey I. Webb and Paul Conilione");
      result.setValue(Field.YEAR, "2002");
      result.setValue(Field.TITLE, "Estimating bias and variance from data");
      result.setValue(Field.INSTITUTION, "Monash University");
      result.setValue(Field.ADDRESS, "School of Computer Science and Software Engineering, Victoria, Australia");
      result.setValue(Field.PDF, "http://www.csse.monash.edu.au/~webb/Files/WebbConilione04.pdf");

      additional = result.add(Type.INPROCEEDINGS);
      additional.setValue(Field.AUTHOR, "Ron Kohavi and David H. Wolpert");
      additional.setValue(Field.YEAR, "1996");
      additional.setValue(Field.TITLE, "Bias Plus Variance Decomposition for Zero-One Loss Functions");
      additional.setValue(Field.BOOKTITLE, "Machine Learning: Proceedings of the Thirteenth International Conference");
      additional.setValue(Field.PUBLISHER, "Morgan Kaufmann");
      additional.setValue(Field.EDITOR, "Lorenza Saitta");
      additional.setValue(Field.PAGES, "275-283");
      additional.setValue(Field.PS, "http://robotics.stanford.edu/~ronnyk/biasVar.ps");

      additional = result.add(Type.ARTICLE);
      additional.setValue(Field.AUTHOR, "Geoffrey I. Webb");
      additional.setValue(Field.YEAR, "2000");
      additional.setValue(Field.TITLE, "MultiBoosting: A Technique for Combining Boosting and Wagging");
      additional.setValue(Field.JOURNAL, "Machine Learning");
      additional.setValue(Field.VOLUME, "40");
      additional.setValue(Field.NUMBER, "2");
      additional.setValue(Field.PAGES, "159-196");

      return result;
    }

    /**
     * Returns an enumeration describing the available options.
     *
     * @return an enumeration of all the available options.
     */
    public Enumeration listOptions() {

        Vector newVector = new Vector(8);

        newVector.addElement(new Option(
        "\tThe index of the class attribute.\n"+
        "\t(default last)",
        "c", 1, "-c <class index>"));
        newVector.addElement(new Option(
        "\tTurn on debugging output.",
        "D", 0, "-D"));
        newVector.addElement(new Option(
        "\tThe number of times each instance is classified.\n"
        +"\t(default 10)",
        "l", 1, "-l <num>"));
        newVector.addElement(new Option(
        "\tThe average proportion of instances common between any two training sets",
        "p", 1, "-p <proportion of objects in common>"));
        newVector.addElement(new Option(
        "\tThe random number seed used.",
        "s", 1, "-s <seed>"));
        newVector.addElement(new Option(
        "\tThe name of the arff file used for the decomposition.",
        "t", 1, "-t <name of arff file>"));
        newVector.addElement(new Option(
        "\tThe number of instances in the training set.",
        "T", 1, "-T <number of instances in training set>"));
        newVector.addElement(new Option(
        "\tFull class name of the learner used in the decomposition.\n"
        +"\teg: weka.classifiers.bayes.NaiveBayes",
        "W", 1, "-W <classifier class name>"));

        if ((m_Classifier != null) &&
        (m_Classifier instanceof OptionHandler)) {
            newVector.addElement(new Option(
            "",
            "", 0, "\nOptions specific to learner "
            + m_Classifier.getClass().getName()
            + ":"));
            Enumeration enu = ((OptionHandler)m_Classifier).listOptions();
            while (enu.hasMoreElements()) {
                newVector.addElement(enu.nextElement());
            }
        }
        return newVector.elements();
    }


    /**
     * Sets the OptionHandler's options using the given list. All options
     * will be set (or reset) during this call (i.e. incremental setting
     * of options is not possible). <p/>
     *
     <!-- options-start -->
     * Valid options are: <p/>
     *
     * <pre> -c &lt;class index&gt;
     *  The index of the class attribute.
     *  (default last)</pre>
     *
     * <pre> -D
     *  Turn on debugging output.</pre>
     *
     * <pre> -l &lt;num&gt;
     *  The number of times each instance is classified.
     *  (default 10)</pre>
     *
     * <pre> -p &lt;proportion of objects in common&gt;
     *  The average proportion of instances common between any two training sets</pre>
     *
     * <pre> -s &lt;seed&gt;
     *  The random number seed used.</pre>
     *
     * <pre> -t &lt;name of arff file&gt;
     *  The name of the arff file used for the decomposition.</pre>
     *
     * <pre> -T &lt;number of instances in training set&gt;
     *  The number of instances in the training set.</pre>
     *
     * <pre> -W &lt;classifier class name&gt;
     *  Full class name of the learner used in the decomposition.
     *  eg: weka.classifiers.bayes.NaiveBayes</pre>
     *
     * <pre>
     * Options specific to learner weka.classifiers.rules.ZeroR:
     * </pre>
     *
     * <pre> -D
     *  If set, classifier is run in debug mode and
     *  may output additional info to the console</pre>
     *
     <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {
        setDebug(Utils.getFlag('D', options));

        String classIndex = Utils.getOption('c', options);
        if (classIndex.length() != 0) {
            if (classIndex.toLowerCase().equals("last")) {
                setClassIndex(0);
            } else if (classIndex.toLowerCase().equals("first")) {
                setClassIndex(1);
            } else {
                setClassIndex(Integer.parseInt(classIndex));
            }
        } else {
            setClassIndex(0);
        }

        String classifyIterations = Utils.getOption('l', options);
        if (classifyIterations.length() != 0) {
            setClassifyIterations(Integer.parseInt(classifyIterations));
        } else {
            setClassifyIterations(10);
        }

        String prob = Utils.getOption('p', options);
        if (prob.length() != 0) {
            setP( Double.parseDouble(prob));
        } else {
            setP(-1);
        }
        //throw new Exception("A proportion must be specified" + " with a -p option.");

        String seedString = Utils.getOption('s', options);
        if (seedString.length() != 0) {
            setSeed(Integer.parseInt(seedString));
        } else {
            setSeed(1);
        }

        String dataFile = Utils.getOption('t', options);
        if (dataFile.length() != 0) {
            setDataFileName(dataFile);
        } else {
            throw new Exception("An arff file must be specified"
            + " with the -t option.");
        }

        String trainSize = Utils.getOption('T', options);
        if (trainSize.length() != 0) {
            setTrainSize(Integer.parseInt(trainSize));
        } else {
            setTrainSize(-1);
        }
        //throw new Exception("A training set size must be specified" + " with a -T option.");

        String classifierName = Utils.getOption('W', options);
        if (classifierName.length() != 0) {
            setClassifier(AbstractClassifier.forName(classifierName, Utils.partitionOptions(options)));
        } else {
            throw new Exception("A learner must be specified with the -W option.");
        }
    }

    /**
     * Gets the current settings of the CheckClassifier.
     *
     * @return an array of strings suitable for passing to setOptions
     */
    public String [] getOptions() {

        String [] classifierOptions = new String [0];
        if ((m_Classifier != null) &&
        (m_Classifier instanceof OptionHandler)) {
            classifierOptions = ((OptionHandler)m_Classifier).getOptions();
        }
        String [] options = new String [classifierOptions.length + 14];
        int current = 0;
        if (getDebug()) {
            options[current++] = "-D";
        }
        options[current++] = "-c"; options[current++] = "" + getClassIndex();
        options[current++] = "-l"; options[current++] = "" + getClassifyIterations();
        options[current++] = "-p"; options[current++] = "" + getP();
        options[current++] = "-s"; options[current++] = "" + getSeed();
        if (getDataFileName() != null) {
            options[current++] = "-t"; options[current++] = "" + getDataFileName();
        }
        options[current++] = "-T"; options[current++] = "" + getTrainSize();
        if (getClassifier() != null) {
            options[current++] = "-W";
            options[current++] = getClassifier().getClass().getName();
        }

        options[current++] = "--";
        System.arraycopy(classifierOptions, 0, options, current,
        classifierOptions.length);
        current += classifierOptions.length;
        while (current < options.length) {
            options[current++] = "";
        }
        return options;
    }

    /**
     * Set the classifiers being analysed
     *
     * @param newClassifier the Classifier to use.
     */
    public void setClassifier(Classifier newClassifier) {

        m_Classifier = newClassifier;
    }

    /**
     * Gets the name of the classifier being analysed
     *
     * @return the classifier being analysed.
     */
    public Classifier getClassifier() {

        return m_Classifier;
    }

    /**
     * Sets debugging mode
     *
     * @param debug true if debug output should be printed
     */
    public void setDebug(boolean debug) {

        m_Debug = debug;
    }

    /**
     * Gets whether debugging is turned on
     *
     * @return true if debugging output is on
     */
    public boolean getDebug() {

        return m_Debug;
    }


    /**
     * Sets the random number seed
     *
     * @param seed the random number seed
     */
    public void setSeed(int seed) {

        m_Seed = seed;
    }

    /**
     * Gets the random number seed
     *
     * @return the random number seed
     */
    public int getSeed() {

        return m_Seed;
    }

    /**
     * Sets the number of times an instance is classified
     *
     * @param classifyIterations number of times an instance is classified
     */
    public void setClassifyIterations(int classifyIterations) {

        m_ClassifyIterations = classifyIterations;
    }

    /**
     * Gets the number of times an instance is classified
     *
     * @return the maximum number of times an instance is classified
     */
    public int getClassifyIterations() {

        return m_ClassifyIterations;
    }

    /**
     * Sets the name of the dataset file.
     *
     * @param dataFileName name of dataset file.
     */
    public void setDataFileName(String dataFileName) {

        m_DataFileName = dataFileName;
    }

    /**
     * Get the name of the data file used for the decomposition
     *
     * @return the name of the data file
     */
    public String getDataFileName() {

        return m_DataFileName;
    }

    /**
     * Get the index (starting from 1) of the attribute used as the class.
     *
     * @return the index of the class attribute
     */
    public int getClassIndex() {

        return m_ClassIndex + 1;
    }

    /**
     * Sets index of attribute to discretize on
     *
     * @param classIndex the index (starting from 1) of the class attribute
     */
    public void setClassIndex(int classIndex) {

        m_ClassIndex = classIndex - 1;
    }

    /**
     * Get the calculated bias squared according to the Kohavi and Wolpert definition
     *
     * @return the bias squared
     */
    public double getKWBias() {

        return m_KWBias;
    }

    /**
     * Get the calculated bias according to the Webb definition
     *
     * @return the bias
     *
     */
    public double getWBias() {

        return m_WBias;
    }


    /**
     * Get the calculated variance according to the Kohavi and Wolpert definition
     *
     * @return the variance
     */
    public double getKWVariance() {

        return m_KWVariance;
    }

    /**
     * Get the calculated variance according to the Webb definition
     *
     * @return the variance according to Webb
     *
     */
    public double getWVariance() {

        return m_WVariance;
    }

    /**
     * Get the calculated sigma according to the Kohavi and Wolpert definition
     *
     * @return the sigma
     *
     */
    public double getKWSigma() {

        return m_KWSigma;
    }

    /**
     * Set the training size.
     *
     * @param size the size of the training set
     *
     */
    public void setTrainSize(int size) {

        m_TrainSize = size;
    }

    /**
     * Get the training size
     *
     * @return the size of the training set
     *
     */
    public int getTrainSize() {

        return m_TrainSize;
    }

    /**
     * Set the proportion of instances that are common between two training sets
     * used to train a classifier.
     *
     * @param proportion the proportion of instances that are common between training
     * sets.
     *
     */
    public void setP(double proportion) {

        m_P = proportion;
    }

    /**
     * Get the proportion of instances that are common between two training sets.
     *
     * @return the proportion
     *
     */
    public double getP() {

        return m_P;
    }

    /**
     * Get the calculated error rate
     *
     * @return the error rate
     */
    public double getError() {

        return m_Error;
    }

    /**
     * Carry out the bias-variance decomposition using the sub-sampled cross-validation method.
     *
     * @throws Exception if the decomposition couldn't be carried out
     */
    public void decompose() throws Exception {

        Reader dataReader;
        Instances data;

        int tps; // training pool size, size of segment E.
        int k; // number of folds in segment E.
        int q; // number of segments of size tps.

        dataReader = new BufferedReader(new FileReader(m_DataFileName)); //open file
        data = new Instances(dataReader); // encapsulate in wrapper class called weka.Instances()

        if (m_ClassIndex < 0) {
            data.setClassIndex(data.numAttributes() - 1);
        } else {
            data.setClassIndex(m_ClassIndex);
        }

        if (data.classAttribute().type() != Attribute.NOMINAL && data.classAttribute().type() != PreferenceAttribute.RANKING) {
            throw new Exception("Class attribute must be nominal");
        }
        int numClasses = data.numClasses();

        data.deleteWithMissingClass();
        if ( data.checkForStringAttributes() ) {
            throw new Exception("Can't handle string attributes!");
        }

        // Dataset size must be greater than 2
        if ( data.numInstances() <= 2 ){
            throw new Exception("Dataset size must be greater than 2.");
        }

        if ( m_TrainSize == -1 ){ // default value
            m_TrainSize = (int) Math.floor( (double) data.numInstances() / 2.0 );
        }else  if ( m_TrainSize < 0 || m_TrainSize >= data.numInstances() - 1 ) {  // Check if 0 < training Size < D - 1
            throw new Exception("Training set size of "+m_TrainSize+" is invalid.");
        }

        if ( m_P == -1 ){ // default value
            m_P = (double) m_TrainSize / ( (double)data.numInstances() - 1 );
        }else if (  m_P < ( m_TrainSize / ( (double)data.numInstances() - 1 ) ) || m_P >= 1.0  ) { //Check if p is in range: m/(|D|-1) <= p < 1.0
            throw new Exception("Proportion is not in range: "+ (m_TrainSize / ((double) data.numInstances() - 1 )) +" <= p < 1.0 ");
        }

        //roundup tps from double to integer
        tps = (int) Math.ceil( ((double)m_TrainSize / (double)m_P) + 1 );
        k = (int) Math.ceil( tps / (tps - (double) m_TrainSize));

        // number of folds cannot be more than the number of instances in the training pool
        if ( k > tps ) {
            throw new Exception("The required number of folds is too many."
            + "Change p or the size of the training set.");
        }

        // calculate the number of segments, round down.
        q = (int) Math.floor( (double) data.numInstances() / (double)tps );

        //create confusion matrix, columns = number of instances in data set, as all will be used,  by rows = number of classes.
        double [][] instanceProbs = new double [data.numInstances()][numClasses];
        int [][] foldIndex = new int [ k ][ 2 ];
        Vector segmentList = new Vector(q + 1);

        //Set random seed
        Random random = new Random(m_Seed);

        data.randomize(random);

        //create index arrays for different segments

        int currentDataIndex = 0;

        for( int count = 1; count <= (q + 1); count++ ){
            if( count > q){
                int [] segmentIndex = new int [ (data.numInstances() - (q * tps)) ];
                for(int index = 0; index < segmentIndex.length; index++, currentDataIndex++){

                    segmentIndex[index] = currentDataIndex;
                }
                segmentList.add(segmentIndex);
            } else {
                int [] segmentIndex = new int [ tps ];

                for(int index = 0; index < segmentIndex.length; index++, currentDataIndex++){
                    segmentIndex[index] = currentDataIndex;
                }
                segmentList.add(segmentIndex);
            }
        }

        int remainder = tps % k; // remainder is used to determine when to shrink the fold size by 1.

        //foldSize = ROUNDUP( tps / k ) (round up, eg 3 -> 3,  3.3->4)
        int foldSize = (int) Math.ceil( (double)tps /(double) k); //roundup fold size double to integer
        int index = 0;
        int currentIndex;

        for( int count = 0; count < k; count ++){
            if( remainder != 0 && count == remainder ){
                foldSize -= 1;
            }
            foldIndex[count][0] = index;
            foldIndex[count][1] = foldSize;
            index += foldSize;
        }

        for( int l = 0; l < m_ClassifyIterations; l++) {

            for(int i = 1; i <= q; i++){

                int [] currentSegment = (int[]) segmentList.get(i - 1);

                randomize(currentSegment, random);

                //CROSS FOLD VALIDATION for current Segment
                for( int j = 1; j <= k; j++){

                    Instances TP = null;
                    for(int foldNum = 1; foldNum <= k; foldNum++){
                        if( foldNum != j){

                            int startFoldIndex = foldIndex[ foldNum - 1 ][ 0 ]; //start index
                            foldSize = foldIndex[ foldNum - 1 ][ 1 ];
                            int endFoldIndex = startFoldIndex + foldSize - 1;

                            for(int currentFoldIndex = startFoldIndex; currentFoldIndex <= endFoldIndex; currentFoldIndex++){

                                if( TP == null ){
                                    TP = new Instances(data, currentSegment[ currentFoldIndex ], 1);
                                }else{
                                    TP.add( data.instance( currentSegment[ currentFoldIndex ] ) );
                                }
                            }
                        }
                    }

                    TP.randomize(random);

                    if( getTrainSize() > TP.numInstances() ){
                        throw new Exception("The training set size of " + getTrainSize() + ", is greater than the training pool "
                        + TP.numInstances() );
                    }

                    Instances train = new Instances(TP, 0, m_TrainSize);

                    Classifier current = AbstractClassifier.makeCopy(m_Classifier);
                    current.buildClassifier(train); // create a clssifier using the instances in train.

                    int currentTestIndex = foldIndex[ j - 1 ][ 0 ]; //start index
                    int testFoldSize = foldIndex[ j - 1 ][ 1 ]; //size
                    int endTestIndex = currentTestIndex + testFoldSize - 1;

                    while( currentTestIndex <= endTestIndex ){

                        Instance testInst = data.instance( currentSegment[currentTestIndex] );
                        int pred = (int)current.classifyInstance( testInst );


                        if(pred != testInst.classValue()) {
                            m_Error++; // add 1 to mis-classifications.
                        }
                        instanceProbs[ currentSegment[ currentTestIndex ] ][ pred ]++;
                        currentTestIndex++;
                    }

                    if( i == 1 && j == 1){
                        int[] segmentElast = (int[])segmentList.lastElement();
                        for( currentIndex = 0; currentIndex < segmentElast.length; currentIndex++){
                            Instance testInst = data.instance( segmentElast[currentIndex] );
                            int pred = (int)current.classifyInstance( testInst );
                            if(pred != testInst.classValue()) {
                                m_Error++; // add 1 to mis-classifications.
                            }

                            instanceProbs[ segmentElast[ currentIndex ] ][ pred ]++;
                        }
                    }
                }
            }
        }

        m_Error /= (double)( m_ClassifyIterations * data.numInstances() );

        m_KWBias = 0.0;
        m_KWVariance = 0.0;
        m_KWSigma = 0.0;

        m_WBias = 0.0;
        m_WVariance = 0.0;

        for (int i = 0; i < data.numInstances(); i++) {

            Instance current = data.instance( i );

            double [] predProbs = instanceProbs[ i ];
            double pActual, pPred;
            double bsum = 0, vsum = 0, ssum = 0;
            double wBSum = 0, wVSum = 0;

            Vector centralTendencies = findCentralTendencies( predProbs );

            if( centralTendencies == null ){
                throw new Exception("Central tendency was null.");
            }

            for (int j = 0; j < numClasses; j++) {
                pActual = (current.classValue() == j) ? 1 : 0;
                pPred = predProbs[j] / m_ClassifyIterations;
                bsum += (pActual - pPred) * (pActual - pPred) - pPred * (1 - pPred) / (m_ClassifyIterations - 1);
                vsum += pPred * pPred;
                ssum += pActual * pActual;
            }

            m_KWBias += bsum;
            m_KWVariance += (1 - vsum);
            m_KWSigma += (1 - ssum);

            for( int count = 0; count < centralTendencies.size(); count++ ) {

                int wB = 0, wV = 0;
                int centralTendency = ((Integer)centralTendencies.get(count)).intValue();

                // For a single instance xi, find the bias and variance.
                for (int j = 0; j < numClasses; j++) {

                    //Webb definition
                    if( j != (int)current.classValue() && j == centralTendency ) {
                        wB += predProbs[j];
                    }
                    if( j != (int)current.classValue() && j != centralTendency ) {
                        wV += predProbs[j];
                    }

                }
                wBSum += (double) wB;
                wVSum += (double) wV;
            }

            // calculate bais by dividing bSum by the number of central tendencies and
            // total number of instances. (effectively finding the average and dividing
            // by the number of instances to get the nominalised probability).

            m_WBias += ( wBSum / ((double) ( centralTendencies.size() * m_ClassifyIterations )));
            // calculate variance by dividing vSum by the total number of interations
            m_WVariance += ( wVSum / ((double) ( centralTendencies.size() * m_ClassifyIterations )));

        }

        m_KWBias /= (2.0 * (double) data.numInstances());
        m_KWVariance /= (2.0 * (double) data.numInstances());
        m_KWSigma /= (2.0 * (double) data.numInstances());

        // bias = bias / number of data instances
        m_WBias /= (double) data.numInstances();
        // variance = variance / number of data instances.
        m_WVariance /= (double) data.numInstances();

        if (m_Debug) {
            System.err.println("Decomposition finished");
        }

    }

    /** Finds the central tendency, given the classifications for an instance.
     *
     * Where the central tendency is defined as the class that was most commonly
     * selected for a given instance.<p>
     *
     * For example, instance 'x' may be classified out of 3 classes y = {1, 2, 3},
     * so if x is classified 10 times, and is classified as follows, '1' = 2 times, '2' = 5 times
     * and '3' = 3 times. Then the central tendency is '2'. <p>
     *
     * However, it is important to note that this method returns a list of all classes
     * that have the highest number of classifications.
     *
     * In cases where there are several classes with the largest number of classifications, then
     * all of these classes are returned. For example if 'x' is classified '1' = 4 times,
     * '2' = 4 times and '3' = 2 times. Then '1' and '2' are returned.<p>
     *
     * @param predProbs the array of classifications for a single instance.
     *
     * @return a Vector containing Integer objects which store the class(s) which
     * are the central tendency.
     */
    public Vector findCentralTendencies(double[] predProbs) {

        int centralTValue = 0;
        int currentValue = 0;
        //array to store the list of classes the have the greatest number of classifictions.
        Vector centralTClasses;

        centralTClasses = new Vector(); //create an array with size of the number of classes.

        // Go through array, finding the central tendency.
        for( int i = 0; i < predProbs.length; i++) {
            currentValue = (int) predProbs[i];
            // if current value is greater than the central tendency value then
            // clear vector and add new class to vector array.
            if( currentValue > centralTValue) {
                centralTClasses.clear();
                centralTClasses.addElement( new Integer(i) );
                centralTValue = currentValue;
            } else if( currentValue != 0 && currentValue == centralTValue) {
                centralTClasses.addElement( new Integer(i) );
            }
        }
        //return all classes that have the greatest number of classifications.
        if( centralTValue != 0){
            return centralTClasses;
        } else {
            return null;
        }

    }

    /**
     * Returns description of the bias-variance decomposition results.
     *
     * @return the bias-variance decomposition results as a string
     */
    public String toString() {

        String result = "\nBias-Variance Decomposition Segmentation, Cross Validation\n" +
        "with subsampling.\n";

        if (getClassifier() == null) {
            return "Invalid setup";
        }

        result += "\nClassifier    : " + getClassifier().getClass().getName();
        if (getClassifier() instanceof OptionHandler) {
            result += Utils.joinOptions(((OptionHandler)m_Classifier).getOptions());
        }
        result += "\nData File     : " + getDataFileName();
        result += "\nClass Index   : ";
        if (getClassIndex() == 0) {
            result += "last";
        } else {
            result += getClassIndex();
        }
        result += "\nIterations    : " + getClassifyIterations();
        result += "\np             : " + getP();
        result += "\nTraining Size : " + getTrainSize();
        result += "\nSeed          : " + getSeed();

        result += "\n\nDefinition   : " +"Kohavi and Wolpert";
        result += "\nError         :" + Utils.doubleToString(getError(), 4);
        result += "\nBias^2        :" + Utils.doubleToString(getKWBias(), 4);
        result += "\nVariance      :" + Utils.doubleToString(getKWVariance(), 4);
        result += "\nSigma^2       :" + Utils.doubleToString(getKWSigma(), 4);

        result += "\n\nDefinition   : " +"Webb";
        result += "\nError         :" + Utils.doubleToString(getError(), 4);
        result += "\nBias          :" + Utils.doubleToString(getWBias(), 4);
        result += "\nVariance      :" + Utils.doubleToString(getWVariance(), 4);

        return result;
    }

    /**
     * Returns the revision string.
     *
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6041 $");
    }

    /**
     * Test method for this class
     *
     * @param args the command line arguments
     */
    public static void main(String [] args) {

        try {
            BVDecomposeSegCVSub bvd = new BVDecomposeSegCVSub();

            try {
                bvd.setOptions(args);
                Utils.checkForRemainingOptions(args);
            } catch (Exception ex) {
                String result = ex.getMessage() + "\nBVDecompose Options:\n\n";
                Enumeration enu = bvd.listOptions();
                while (enu.hasMoreElements()) {
                    Option option = (Option) enu.nextElement();
                    result += option.synopsis() + "\n" + option.description() + "\n";
                }
                throw new Exception(result);
            }

            bvd.decompose();

            System.out.println(bvd.toString());

        } catch (Exception ex) {
            System.err.println(ex.getMessage());
        }

    }

    /**
     * Accepts an array of ints and randomises the values in the array, using the
     * random seed.
     *
     *@param index is the array of integers
     *@param random is the Random seed.
     */
    public final void randomize(int[] index, Random random) {
        for( int j = index.length - 1; j > 0; j-- ){
            int k = random.nextInt( j + 1 );
            int temp = index[j];
            index[j] = index[k];
            index[k] = temp;
        }
    }
}
