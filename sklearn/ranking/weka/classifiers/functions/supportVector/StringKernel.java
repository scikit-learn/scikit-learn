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
 * StringKernel.java
 * 
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.functions.supportVector;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Implementation of the subsequence kernel (SSK) as described in [1] and of the subsequence kernel with lambda pruning (SSK-LP) as described in [2].<br/>
 * <br/>
 * For more information, see<br/>
 * <br/>
 * Huma Lodhi, Craig Saunders, John Shawe-Taylor, Nello Cristianini, Christopher J. C. H. Watkins (2002). Text Classification using String Kernels. Journal of Machine Learning Research. 2:419-444.<br/>
 * <br/>
 * F. Kleedorfer, A. Seewald (2005). Implementation of a String Kernel for WEKA. Wien, Austria.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Lodhi2002,
 *    author = {Huma Lodhi and Craig Saunders and John Shawe-Taylor and Nello Cristianini and Christopher J. C. H. Watkins},
 *    journal = {Journal of Machine Learning Research},
 *    pages = {419-444},
 *    title = {Text Classification using String Kernels},
 *    volume = {2},
 *    year = {2002},
 *    HTTP = {http://www.jmlr.org/papers/v2/lodhi02a.html}
 * }
 * 
 * &#64;techreport{Kleedorfer2005,
 *    address = {Wien, Austria},
 *    author = {F. Kleedorfer and A. Seewald},
 *    institution = {Oesterreichisches Forschungsinstitut fuer Artificial Intelligence},
 *    number = {TR-2005-13},
 *    title = {Implementation of a String Kernel for WEKA},
 *    year = {2005}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Enables debugging output (if available) to be printed.
 *  (default: off)</pre>
 * 
 * <pre> -no-checks
 *  Turns off all checks - use with caution!
 *  (default: checks on)</pre>
 * 
 * <pre> -P &lt;0|1&gt;
 *  The pruning method to use:
 *  0 = No pruning
 *  1 = Lambda pruning
 *  (default: 0)</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The size of the cache (a prime number).
 *  (default: 250007)</pre>
 * 
 * <pre> -IC &lt;num&gt;
 *  The size of the internal cache (a prime number).
 *  (default: 200003)</pre>
 * 
 * <pre> -L &lt;num&gt;
 *  The lambda constant. Penalizes non-continuous subsequence
 *  matches. Must be in (0,1).
 *  (default: 0.5)</pre>
 * 
 * <pre> -ssl &lt;num&gt;
 *  The length of the subsequence.
 *  (default: 3)</pre>
 * 
 * <pre> -ssl-max &lt;num&gt;
 *  The maximum length of the subsequence.
 *  (default: 9)</pre>
 * 
 * <pre> -N
 *  Use normalization.
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 * 
 * <h1>Theory</h1>
 * <h2>Overview</h2>
 * The algorithm computes a measure of similarity between two texts based on
 * the number and form of their common subsequences, which need not be
 * contiguous. This method can be parametrized by specifying the subsequence
 * length k, the penalty factor lambda, which penalizes non-contiguous matches,
 * and optional 'lambda pruning', which takes maxLambdaExponent,
 * <code>m</code>, as parameter. Lambda pruning causes very 'stretched'
 * substring matches not to be counted, thus speeding up the computation. The
 * functionality of SSK and SSK-LP is explained in the following using simple
 * examples.
 * 
 * <h2>Explanation &amp; Examples</h2>
 * for all of the following examples, we assume these parameter values: 
 *<pre> 
 *k=2
 *lambda=0.5
 *m=8 (for SSK-LP examples)
 *</pre>
 * 
 * <h3>SSK</h3>
 * 
 * <h4>Example 1</h4>
 * 
 * <pre>
 *SSK(2,"ab","axb")=0.5^5 = 0,03125
 *</pre>
 * There is one subsequence of the length of 2 that both strings have in
 * common, "ab".  The result of SSK is computed by raising lambda to the power
 * of L, where L is the length of the subsequence match in the one string plus
 * the length of the subsequence match in the other, in our case:
 * <pre>
 *&nbsp;  ab    axb
 *L= 2  +   3 = 5
 * </pre>
 * hence, the kernel yields 0.5^5 = 0,03125
 *
 * <h4>Example 2</h4>
 * <pre>
 *SSK(2,"ab","abb")=0.5^5 + 0.5^4 = 0,09375
 *</pre>
 * Here, we also have one subsequence of the length of 2 that both strings have
 * in common, "ab".  The result of SSK is actually computed by summing over all
 * values computed for each occurrence of a common subsequence match. In this
 * example, there are two possible cases:
 * <pre>
 *ab    abb
 *--    --  L=4
 *--    - - L=5
 * </pre>
 * we have two matches, one of the length of 2+2=4, one of the length of 2+3=5, 
 * so we get the result 0.5^5 + 0.5^4 = 0,09375.
 *
 * <h3>SSK-LP</h3>
 * Without lambda pruning, the string kernel finds *all* common subsequences of
 * the given length, whereas with lambda pruning, common subsequence matches
 * that are too much stretched in both strings are not taken into account. It
 * is argued that the value yielded for such a common subsequence is too low
 * (<code>lambda ^(length[match_in_s] + length[match_in_t]</code>) . Tests have
 * shown that a tremendous speedup can be achieved using this technique while
 * suffering from very little quality loss. <br>
 * Lambda pruning is parametrized by the maximum lambda exponent. It is a good
 * idea to choose that value to be about 3 or 4 times the subsequence length as
 * a rule of thumb. YMMV.
 *
 * <h4>Example 3</h4>
 * Without lambda pruning, one common subsequence, 
 * "AB" would be found in the following two strings. (With k=2)  
 * <pre>
 *SSK(2,"ab","axb")=0.5^14 = 0,00006103515625
 *</pre>
 * lambda pruning allows for the control of the match length. So, if m 
 * (the maximum lambda exponent) is e.g. 8, these two strings would 
 * yield a kernel value of 0:
 * <pre>
 *with lambda pruning:    SSK-LP(2,8,"AxxxxxxxxxB","AyB")= 0
 *without lambda pruning: SSK(2,"AxxxxxxxxxB","AyB")= 0.5^14 = 0,00006103515625  
 *</pre>
 * This is because the exponent for lambda (=the length of the subsequence
 * match) would be 14, which is &gt; 8. In Contrast, the next result is
 * &gt; 0
 *<pre>
 *m=8
 *SSK-LP(2,8,"AxxB","AyyB")=0.5^8 = 0,00390625
 *</pre>
 * because the lambda exponent would be 8, which is just accepted by lambda
 * pruning.
 *
 * <h3>Normalization</h3>
 * When the string kernel is used for its main purpose, as the kernel of a
 * support vector machine, it is not normalized.  The normalized kernel can be
 * switched on by -F (feature space normalization) but is much slower.  Like
 * most unnormalized kernels, K(x,x) is not a fixed value, see the next
 * example.
 *
 * <h4>Example 4</h4> 
 *<pre>
 *SSK(2,"ab","ab")=0.5^4 = 0.0625
 *SSK(2,"AxxxxxxxxxB","AxxxxxxxxxB") = 12.761724710464478
 *</pre>
 * SSK is evaluated twice, each time for two identical strings. A good measure
 * of similarity would produce the same value in both cases, which should  
 * indicate the same level of similarity. The value of the normalized SSK would
 * be 1.0 in both cases. So for the purpose of computing string similarity the
 * normalized kernel should be used. For SVM the unnormalized kernel is usually
 * sufficient.
 *
 * <h2>Complexity of SSK and SSK-LP</h2>
 * The time complexity of this method (without lambda pruning and with an
 * infinitely large cache) is<br> 
 * <pre>O(k*|s|*|t|)</pre>
 * Lambda Pruning has a complexity (without caching) of<br> 
 * <pre>O(m*binom(m,k)^2*(|s|+n)*|t|)</pre> <br>  
 * <pre>
 *k...          subsequence length (ssl)
 *s,t...        strings
 *|s|...        length of string s
 *binom(x,y)... binomial coefficient (x!/[(x-y)!y!])
 *m...          maxLambdaExponent (ssl-max)
 *</pre>
 * 
 * Keep in mind that execution time can increase fast for long strings 
 * and big values for k, especially if you don't use lambda pruning.
 * With lambda pruning, computation is usually so fast that switching
 * on the cache leads to slower computation because of setup costs. Therefore
 * caching is switched off for lambda pruning.
 * <br>
 * <br>
 * For details and qualitative experiments about SSK, see [1] <br>
 * For details about lambda pruning and performance comparison of SSK 
 * and SSK-LP (SSK with lambda pruning), see [2]   
 * Note that the complexity estimation in [2] assumes no caching of
 * intermediate results, which has been implemented in the meantime and
 * greatly improves the speed of the SSK without lambda pruning.
 *<br>
 *
 *<h1>Notes for usage within Weka</h1>
 * Only instances of the following form can be processed using string kernels:
 * <pre>
 *+----------+-------------+---------------+
 *|attribute#|     0       |       1       |
 *+----------+-------------+---------------+
 *| content  | [text data] | [class label] |
 *+----------------------------------------+
 * ... or ...
 *+----------+---------------+-------------+
 *|attribute#|     0         |     1       |
 *+----------+---------------+-------------+
 *| content  | [class label] | [text data] |
 *+----------------------------------------+
 *</pre>
 *
 * @author Florian Kleedorfer (kleedorfer@austria.fm)
 * @author Alexander K. Seewald (alex@seewald.at)
 * @version $Revision: 5450 $
 */
public class StringKernel 
  extends Kernel
  implements TechnicalInformationHandler {
  
  /** for serialization */
  private static final long serialVersionUID = -4902954211202690123L;

  /** The size of the cache (a prime number) */
  private int m_cacheSize = 250007;

  /** The size of the internal cache for intermediate results (a prime number) */
  private int m_internalCacheSize = 200003;

  /** The attribute number of the string attribute */
  private int m_strAttr;

  /** Kernel cache (i.e., cache for kernel evaluations) */
  private double[] m_storage;
  private long[] m_keys;

  /** Counts the number of kernel evaluations. */
  private int m_kernelEvals;

  /** The number of instance in the dataset */
  private int m_numInsts;

  /** Pruning method: No Pruning */
  public final static int PRUNING_NONE = 0;
  /** Pruning method: Lambda See [2] for details. */
  public final static int PRUNING_LAMBDA = 1;
  /** Pruning methods */
  public static final Tag [] TAGS_PRUNING = {
    new Tag(PRUNING_NONE, "No pruning"),
    new Tag(PRUNING_LAMBDA, "Lambda pruning"),
  };
  
  /** the pruning method */
  protected int m_PruningMethod = PRUNING_NONE;

  /** the decay factor that penalizes non-continuous substring matches. See [1]
   * for details. */
  protected double m_lambda = 0.5;

  /** The substring length */
  private int m_subsequenceLength = 3;

  /** The maximum substring length for lambda pruning */
  private int m_maxSubsequenceLength = 9;

  /** powers of lambda are prepared prior to kernel evaluations.
   * all powers between 0 and this value are precalculated */
  protected static final int MAX_POWER_OF_LAMBDA = 10000;

  /** the precalculated powers of lambda */
  protected double[] m_powersOflambda = null;

  /** flag for switching normalization on or off. This defaults to false and
   * can be turned on by the switch for feature space normalization in SMO
   */
  private boolean m_normalize = false;

  /** private cache for intermediate results */	
  private int maxCache; // is set in unnormalizedKernel(s1,s2)
  private double[] cachekh;
  private int[] cachekhK;
  private double[] cachekh2;
  private int[] cachekh2K;
  /** cached indexes for private cache */
  private int m_multX;
  private int m_multY;
  private int m_multZ;
  private int m_multZZ;

  private boolean m_useRecursionCache = true;

  /**
   * default constructor
   */
  public StringKernel() {
    super();
  }
  
  /**
   * creates a new StringKernel object. Initializes the kernel cache and the
   * 'lambda cache', i.e. the precalculated powers of lambda from lambda^2 to
   * lambda^MAX_POWER_OF_LAMBDA
   *
   * @param data		the dataset to use
   * @param cacheSize		the size of the cache
   * @param subsequenceLength	the subsequence length
   * @param lambda		the lambda value
   * @param debug		whether to output debug information
   * @throws Exception		if something goes wrong
   */
  public StringKernel(Instances data, int cacheSize, int subsequenceLength,
    double lambda, boolean debug) throws Exception {

    setDebug(debug);
    setCacheSize(cacheSize);
    setInternalCacheSize(200003);
    setSubsequenceLength(subsequenceLength);
    setMaxSubsequenceLength(-1);
    setLambda(lambda);
    
    buildKernel(data);
  }
  
  /**
   * Returns a string describing the kernel
   * 
   * @return a description suitable for displaying in the
   *         explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Implementation of the subsequence kernel (SSK) as described in [1] "
      + "and of the subsequence kernel with lambda pruning (SSK-LP) as "
      + "described in [2].\n\n"
      + "For more information, see\n\n"
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
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Huma Lodhi and Craig Saunders and John Shawe-Taylor and Nello Cristianini and Christopher J. C. H. Watkins");
    result.setValue(Field.YEAR, "2002");
    result.setValue(Field.TITLE, "Text Classification using String Kernels");
    result.setValue(Field.JOURNAL, "Journal of Machine Learning Research");
    result.setValue(Field.VOLUME, "2");
    result.setValue(Field.PAGES, "419-444");
    result.setValue(Field.HTTP, "http://www.jmlr.org/papers/v2/lodhi02a.html");

    additional = result.add(Type.TECHREPORT);
    additional.setValue(Field.AUTHOR, "F. Kleedorfer and A. Seewald");
    additional.setValue(Field.YEAR, "2005");
    additional.setValue(Field.TITLE, "Implementation of a String Kernel for WEKA");
    additional.setValue(Field.INSTITUTION, "Oesterreichisches Forschungsinstitut fuer Artificial Intelligence");
    additional.setValue(Field.ADDRESS, "Wien, Austria");
    additional.setValue(Field.NUMBER, "TR-2005-13");
    
    return result;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector		result;
    Enumeration		en;
    String		desc;
    String		param;
    int			i;
    SelectedTag		tag;
    
    result = new Vector();

    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    desc  = "";
    param = "";
    for (i = 0; i < TAGS_PRUNING.length; i++) {
      if (i > 0)
	param += "|";
      tag = new SelectedTag(TAGS_PRUNING[i].getID(), TAGS_PRUNING);
      param += "" + tag.getSelectedTag().getID();
      desc  +=   "\t" + tag.getSelectedTag().getID() 
      	       + " = " + tag.getSelectedTag().getReadable()
      	       + "\n";
    }

    result.addElement(new Option(
	"\tThe pruning method to use:\n"
	+ desc
	+ "\t(default: " + PRUNING_NONE + ")",
	"P", 1, "-P <" + param + ">"));

    result.addElement(new Option(
	"\tThe size of the cache (a prime number).\n"
	+ "\t(default: 250007)",
	"C", 1, "-C <num>"));

    result.addElement(new Option(
	"\tThe size of the internal cache (a prime number).\n"
	+ "\t(default: 200003)",
	"IC", 1, "-IC <num>"));

    result.addElement(new Option(
	"\tThe lambda constant. Penalizes non-continuous subsequence\n"
	+ "\tmatches. Must be in (0,1).\n"
	+ "\t(default: 0.5)",
	"L", 1, "-L <num>"));

    result.addElement(new Option(
	"\tThe length of the subsequence.\n"
	+ "\t(default: 3)",
	"ssl", 1, "-ssl <num>"));

    result.addElement(new Option(
	"\tThe maximum length of the subsequence.\n"
	+ "\t(default: 9)",
	"ssl-max", 1, "-ssl-max <num>"));

    result.addElement(new Option(
	"\tUse normalization.\n"
	+ "\t(default: no)",
	"N", 0, "-N"));

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Enables debugging output (if available) to be printed.
   *  (default: off)</pre>
   * 
   * <pre> -no-checks
   *  Turns off all checks - use with caution!
   *  (default: checks on)</pre>
   * 
   * <pre> -P &lt;0|1&gt;
   *  The pruning method to use:
   *  0 = No pruning
   *  1 = Lambda pruning
   *  (default: 0)</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The size of the cache (a prime number).
   *  (default: 250007)</pre>
   * 
   * <pre> -IC &lt;num&gt;
   *  The size of the internal cache (a prime number).
   *  (default: 200003)</pre>
   * 
   * <pre> -L &lt;num&gt;
   *  The lambda constant. Penalizes non-continuous subsequence
   *  matches. Must be in (0,1).
   *  (default: 0.5)</pre>
   * 
   * <pre> -ssl &lt;num&gt;
   *  The length of the subsequence.
   *  (default: 3)</pre>
   * 
   * <pre> -ssl-max &lt;num&gt;
   *  The maximum length of the subsequence.
   *  (default: 9)</pre>
   * 
   * <pre> -N
   *  Use normalization.
   *  (default: no)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('P', options);
    if (tmpStr.length() != 0)
      setPruningMethod(
	  new SelectedTag(Integer.parseInt(tmpStr), TAGS_PRUNING));
    else
      setPruningMethod(
	  new SelectedTag(PRUNING_NONE, TAGS_PRUNING));

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setCacheSize(Integer.parseInt(tmpStr));
    else
      setCacheSize(250007);
    
    tmpStr = Utils.getOption("IC", options);
    if (tmpStr.length() != 0)
      setInternalCacheSize(Integer.parseInt(tmpStr));
    else
      setInternalCacheSize(200003);
    
    tmpStr = Utils.getOption('L', options);
    if (tmpStr.length() != 0)
      setLambda(Double.parseDouble(tmpStr));
    else
      setLambda(0.5);
    
    tmpStr = Utils.getOption("ssl", options);
    if (tmpStr.length() != 0)
      setSubsequenceLength(Integer.parseInt(tmpStr));
    else
      setSubsequenceLength(3);
    
    tmpStr = Utils.getOption("ssl-max", options);
    if (tmpStr.length() != 0)
      setMaxSubsequenceLength(Integer.parseInt(tmpStr));
    else
      setMaxSubsequenceLength(9);

    setUseNormalization(Utils.getFlag('N', options));

    if (getMaxSubsequenceLength()<2*getSubsequenceLength()) {
      throw new IllegalArgumentException("Lambda Pruning forbids even contiguous substring matches! " +
      "Use a bigger value for ssl-max (at least 2*ssl).");
    }
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Kernel.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector    result;
    String[]  options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-P");
    result.add("" + m_PruningMethod);

    result.add("-C");
    result.add("" + getCacheSize());

    result.add("-IC");
    result.add("" + getInternalCacheSize());

    result.add("-L");
    result.add("" + getLambda());

    result.add("-ssl");
    result.add("" + getSubsequenceLength());

    result.add("-ssl-max");
    result.add("" + getMaxSubsequenceLength());

    if (getUseNormalization())
      result.add("-L");

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String pruningMethodTipText() {
    return "The pruning method.";
  }

  /**
   * Sets the method used to for pruning. 
   *
   * @param value 	the pruning method to use.
   */
  public void setPruningMethod(SelectedTag value) {
    if (value.getTags() == TAGS_PRUNING)
      m_PruningMethod = value.getSelectedTag().getID();
  }

  /**
   * Gets the method used for pruning. 
   *
   * @return 		the pruning method to use.
   */
  public SelectedTag getPruningMethod() {
    return new SelectedTag(m_PruningMethod, TAGS_PRUNING);
  }


  /**
   * Sets the size of the cache to use (a prime number)
   * 
   * @param value	the size of the cache
   */
  public void setCacheSize(int value) {
    if (value >= 0) {
      m_cacheSize = value;
      clean();
    }
    else {
      System.out.println(
	  "Cache size cannot be smaller than 0 (provided: " + value + ")!");
    }
  }
  
  /**
   * Gets the size of the cache
   * 
   * @return 		the cache size
   */
  public int getCacheSize() {
    return m_cacheSize;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String cacheSizeTipText() {
    return "The size of the cache (a prime number).";
  }

  /** 
   * sets the size of the internal cache for intermediate results. Memory
   * consumption is about 16x this amount in bytes. Only use when lambda
   * pruning is switched off.
   *
   * @param value	the size of the internal cache
   */
  public void setInternalCacheSize(int value) {
    if (value >= 0) {
      m_internalCacheSize = value;
      clean();
    } else {
      System.out.println(
	  "Cache size cannot be smaller than 0 (provided: " + value + ")!");
    }
  }
  
  /**
   * Gets the size of the internal cache
   * 
   * @return 		the cache size
   */
  public int getInternalCacheSize() {
    return m_internalCacheSize;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String internalCacheSizeTipText() {
    return "The size of the internal cache (a prime number).";
  }

  /**
   * Sets the length of the subsequence.
   * 
   * @param value	the length
   */
  public void setSubsequenceLength(int value) {
    m_subsequenceLength = value;
  }
  
  /**
   * Returns the length of the subsequence
   * 
   * @return		the length
   */
  public int getSubsequenceLength() {
    return m_subsequenceLength;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String subsequenceLengthTipText() {
    return "The subsequence length.";
  }

  /**
   * Sets the maximum length of the subsequence.
   * 
   * @param value	the maximum length
   */
  public void setMaxSubsequenceLength(int value) {
    m_maxSubsequenceLength = value;
  }
  
  /**
   * Returns the maximum length of the subsequence
   * 
   * @return		the maximum length
   */
  public int getMaxSubsequenceLength() {
    return m_maxSubsequenceLength;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String maxSubsequenceLengthTipText() {
    return "The maximum subsequence length (theta in the paper)";
  }
  
  /**
   * Sets the lambda constant used in the string kernel
   * 
   * @param value	the lambda value to use
   */
  public void setLambda(double value) {
    m_lambda = value;
  }
  
  /**
   * Gets the lambda constant used in the string kernel
   * 
   * @return		the current lambda constant
   */  
  public double getLambda() {
    return m_lambda;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String lambdaTipText(){
    return "Penalizes non-continuous subsequence matches, from (0,1)";
  }

  /**
   * Sets whether to use normalization.
   * Each time this value is changed, the kernel cache is cleared.
   *
   * @param value	whether to use normalization
   */
  public void setUseNormalization(boolean value) {
    if  (value != m_normalize) 
      clean();
    
    m_normalize = value;
  }
  
  /**
   * Returns whether normalization is used.
   * 
   * @return		true if normalization is used
   */
  public boolean getUseNormalization() {
    return m_normalize;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String useNormalizationTipText(){
    return "Whether to use normalization.";
  }

  /**
   * Computes the result of the kernel function for two instances.
   * If id1 == -1, eval use inst1 instead of an instance in the dataset.
   *
   * @param id1 the index of the first instance in the dataset
   * @param id2 the index of the second instance in the dataset
   * @param inst1 the instance corresponding to id1 (used if id1 == -1)
   * @return the result of the kernel function
   * @throws Exception if something goes wrong
   */
  public double eval(int id1, int id2, Instance inst1) throws Exception {
    if (m_Debug && id1>-1 && id2>-1) {
      System.err.println("\nEvaluation of string kernel for");
      System.err.println(m_data.instance(id1).stringValue(m_strAttr));
      System.err.println("and");
      System.err.println(m_data.instance(id2).stringValue(m_strAttr));
    }

    //the normalized kernel returns 1 for comparison of 
    //two identical strings
    if (id1 == id2 && m_normalize) 
      return 1.0;

    double result = 0;
    long key = -1;
    int location = -1;

    // we can only cache if we know the indexes
    if ((id1 >= 0) && (m_keys != null)) {
      if (id1 > id2) {
        key = (long)id1 * m_numInsts + id2;
      } else {
        key = (long)id2 * m_numInsts + id1;
      }
      if (key < 0) {
        throw new Exception("Cache overflow detected!");
      }
      location = (int)(key % m_keys.length);
      if (m_keys[location] == (key + 1)) {
        if (m_Debug) 
          System.err.println("result (cached): " + m_storage[location]);
        return m_storage[location];
      }
    }

    m_kernelEvals++;
    long start = System.currentTimeMillis();

    Instance inst2 = m_data.instance(id2);
    char[] s1 = inst1.stringValue(m_strAttr).toCharArray();
    char[] s2 = inst2.stringValue(m_strAttr).toCharArray();

    // prevent the kernel from returning NaN
    if (s1.length == 0 || s2.length == 0) return 0;

    if (m_normalize) {
      result = normalizedKernel(s1,s2);
    } else {
      result = unnormalizedKernel(s1, s2);
    }

    if (m_Debug) {
      long duration = System.currentTimeMillis() - start;
      System.err.println("result: " + result);
      System.err.println("evaluation time:" + duration +"\n");
    }

    // store result in cache
    if (key != -1){
      m_storage[location] = result;
      m_keys[location] = (key + 1);
    }
    return result;
  }

  /**
   * Frees the memory used by the kernel.
   * (Useful with kernels which use cache.)
   * This function is called when the training is done.
   * i.e. after that, eval will be called with id1 == -1.
   */
  public void clean() {
    m_storage = null;
    m_keys = null;
  }

  /**
   * Returns the number of kernel evaluation performed.
   *
   * @return the number of kernel evaluation performed.
   */
  public int numEvals() {
    return m_kernelEvals;
  }

  /**
   * Returns the number of dot product cache hits.
   *
   * @return the number of dot product cache hits, or -1 if not supported by
   * this kernel.
   */
  public int numCacheHits() {
    // TODO: implement!
    return -1;
  }
  
  /**
   * evaluates the normalized kernel between s and t. See [1] for details about
   * the normalized SSK. 
   *
   * @param s first input string
   * @param t second input string
   * @return a double indicating their distance, or similarity
   */
  public double normalizedKernel(char[] s, char[] t){
    double k1 = unnormalizedKernel(s, s);
    double k2 = unnormalizedKernel(t, t);
    double normTerm = Math.sqrt( k1*k2 );
    return unnormalizedKernel(s, t) / normTerm;
  }

  /**
   * evaluates the unnormalized kernel between s and t. See [1] for details
   * about the unnormalized SSK.
   *
   * @param s first input string
   * @param t second input string
   * @return a double indicating their distance, or similarity
   */
  public double unnormalizedKernel(char[] s, char[] t){
    if (t.length > s.length) {
      //swap because the algorithm is faster if s is
      //the longer string
      char[] buf = s;
      s = t;
      t = buf;
    }
    if (m_PruningMethod == PRUNING_NONE) {
      m_multX=(s.length+1)*(t.length+1); 
      m_multY=(t.length+1); 
      m_multZ=1;
      maxCache = m_internalCacheSize;
      if (maxCache==0) {
	maxCache=(m_subsequenceLength+1)*m_multX;
      }
      else if ((m_subsequenceLength+1)*m_multX<maxCache) { 
	maxCache=(m_subsequenceLength+1)*m_multX; 
      }
      m_useRecursionCache=true;
      cachekhK = new int[maxCache];
      cachekh2K = new int[maxCache];
      cachekh = new double[maxCache];
      cachekh2 = new double[maxCache];
    } else if (m_PruningMethod == PRUNING_LAMBDA) {
      maxCache=0; 
      m_useRecursionCache=false;
    }

    double res;
    if (m_PruningMethod == PRUNING_LAMBDA) {
      res = kernelLP(
              m_subsequenceLength,s,s.length-1,t,t.length-1,
              m_maxSubsequenceLength);
    } else {
      res = kernel(
              m_subsequenceLength,s,s.length-1, t, t.length-1);
    }
    cachekh = null;
    cachekhK = null;
    cachekh2 = null;
    cachekh2K = null;
    
    return res;
  }

  /**
   * Recursion-ending function that is called at the end of each 
   * recursion branch.
   * 
   * @param n
   * @return
   */
  protected double getReturnValue(int n){
    if (n == 0) return 1; else return 0;
  }

  /**
   * the kernel function (Kn). This function performs the outer loop
   * character-wise over the first input string s. For each character
   * encountered, a recursion branch is started that identifies all
   * subsequences in t starting with that character. <br> See [1] for details
   * but note that this code is optimized and may be hard to recognize.
   * 
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @return a double indicating the distance or similarity between s and t, 
   * according to and depending on the initial value for n.
   */
  protected double kernel(int n, char[] s,int endIndexS, char[] t, 
    int endIndexT) {

    //normal recursion ending case
    if (Math.min(endIndexS+1,endIndexT+1) < n) return getReturnValue(n);

    //accumulate all recursion results in one:
    double result = 0;

    //the tail-recursive function defined in [1] is turned into a
    //loop here, preventing stack overflows.
    //skim s from back to front
    for (int iS=endIndexS; iS > n-2; iS--) {
      double buf = 0;
      //let the current character in s be x 
      char x = s[iS];
      // iterate over all occurrences of x in t
      for (int j=0; j <= endIndexT; j++) {
        if (t[j] == x){
          //this is a match for the current character, hence
          //1. use previous chars in both strings (iS-1, j-1)
          //2. decrement the remainingMatchLength (n-1)
          //and start a recursion branch for these parameters
          buf += kernelHelper(n-1,s,iS-1, t, j-1);
        }
      }
      //ok, all occurrences of x in t have been found
      //multiply the result with lambda^2
      //  (one lambda for x, and the other for all matches of x in t)
      result += buf * m_powersOflambda[2];
    }
    return result;
  }

  /**
   * The kernel helper function, called K' in [1] and [2].
   *
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @return a partial result for K
   */
  protected double kernelHelper (int n, char[] s,int endIndexS, char[] t, 
    int endIndexT) {

    //recursion ends if the current subsequence has maximal length, 
    //which is the case here
    if (n <= 0 ) {
      return getReturnValue(n);
    }

    //recursion ends, too, if the current subsequence is shorter than
    //maximal length, but there is no chance that it will reach maximal length.
    //in this case, normally 0 is returned, but the EXPERIMENTAL 
    //minSubsequenceLength feature allows shorter subsequence matches
    //also to contribute
    if (Math.min(endIndexS+1,endIndexT+1) < n) {
      return getReturnValue(n);
    }
    int adr = 0;
    if (m_useRecursionCache) {
      adr=m_multX*n+m_multY*endIndexS+m_multZ*endIndexT;
      if ( cachekhK[adr % maxCache] == adr+1) return cachekh[adr % maxCache];
    }

    //the tail-recursive function defined in [1] is turned into a
    //loop here, preventing stack overflows.
    //loop over s, nearly from the start (skip the first n-1 characters)
    //and only up until endIndexS, and recursively apply K''. Thus, every
    //character between n-1 and endIndexS in s is counted once as 
    //being part of the subsequence match and once just as a gap. 
    //In both cases lambda is multiplied with the result.
    double result = 0;
    /*
       for (int iS = n-1; iS <= endIndexS;iS++) {
       result *= m_lambda;
       result += kernelHelper2(n,s,iS, t, endIndexT);
       }
       if (m_useRecursionCache) {
       cachekhK[adr % maxCache]=adr+1; cachekh[adr % maxCache]=result;
       }
       return result;
       */
    /* ^^^ again, above code segment does not store some intermediate results... */
    result = m_lambda*kernelHelper(n,s,endIndexS-1,t,endIndexT)
             + kernelHelper2(n,s,endIndexS,t,endIndexT);
    if (m_useRecursionCache) {
      cachekhK[adr % maxCache]=adr+1; cachekh[adr % maxCache]=result;
    }
    return result;
  }

  /** 
   * helper function for the evaluation of the kernel K'' see section
   * 'Efficient Computation of SSK' in [1]
   *
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @return a partial result for K'
   */
  protected double kernelHelper2(int n, char[] s, int endIndexS, char[] t, 
    int endIndexT) {

    //recursion ends if one of the indices in both strings is <0
    if (endIndexS <0 || endIndexT <0) {
      return getReturnValue(n);
    }

    int adr = 0;
    if (m_useRecursionCache) {
      adr=m_multX*n+m_multY*endIndexS+m_multZ*endIndexT;
      if ( cachekh2K[adr % maxCache] == adr+1) return cachekh2[adr % maxCache];
    }

    //spot the last character in s, we'll need it
    char x = s[endIndexS];

    //recurse if the last characters of s and t, x (and y) are identical.
    //which is an easy case: just add up two recursions, 
    // 1. one that counts x and y as a part of the subsequence match
    //	 -> n, endIndexS and endIndexT are decremented for next recursion level
    // 	 -> lambda^2 is multiplied with the result to account for the length
    //      of 2 that has been added to the length of the subsequence match
    //      by accepting x and y.
    // 2. one that counts y as a gap in the match 
    //   -> only endIndexT is decremented for next recursion level
    // 	 -> lambda is multiplied with the result to account for the length
    //      of 1 that has been added to the length of the subsequence match
    //		by omitting y.
    if (x == t[endIndexT]) {
      double ret =  m_lambda * (kernelHelper2(n,s,endIndexS, t, endIndexT-1)
          + m_lambda * kernelHelper(n-1,s,endIndexS-1, t, endIndexT-1));
      if (m_useRecursionCache) {
        cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret;
      }
      return ret;
    } else {
      double ret = m_lambda*kernelHelper2(n,s,endIndexS,t,endIndexT-1);
      if (m_useRecursionCache) {
        cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret;
      }
      return ret;
    }

    //look for x in t from back to front. 
    //this is actually an optimization from [1] that spares unneccessary
    //recursions iff
    //x is actually found in t, but not at the last position.
    /*
       int i;
       int threshold = n>0?n-1:0;
       for (i=endIndexT-1; i >= threshold;i--) {
       if (x == t[i]) {
       double ret=getPowerOfLambda(endIndexT-i) * kernelHelper2(n,s,endIndexS, t, i);
       if (m_useRecursionCache) {
       cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret;
       }
       return ret;
       }
       }
       */		
    //end the recursion if x is not found in t.
    /*        double ret = getReturnValue(n);
              if (m_useRecursionCache) {
              cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret;
              }
              return ret;*/
  }

  /**
   * the kernel function K explained in [1] using lambda pruning, explained in
   * [2].  An additional parameter is introduced, which denotes the maximum
   * length of a subsequence match. This allows for the control of how relaxed
   * the subsequence matches are. <br>
   *
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @param remainingMatchLength actually the initial value for
   * maxLambdaExponent
   * @return a double indicating the distance or similarity between s and t, 
   * according to and depending on the initial value for n. 
   */
  protected double kernelLP(int n, char[] s, int endIndexS,char[] t,
      int endIndexT,int remainingMatchLength) {
    //see code docs in kernel()
    if (Math.min(endIndexS+1,endIndexT +1) < n) {
      return getReturnValue(n);
    }
    //lambda pruning check 
    //stops recursion if the match is so long that the resulting
    //power of lambda is smaller than minLambda
    //if lambda pruning is not used, the remainingMatchLength is < 0
    //and this check never stops the recursion
    if (remainingMatchLength == 0) return getReturnValue(n);
    double result = 0;
    //see code docs in kernel()
    for (int iS =endIndexS; iS > n-2; iS--) {
      double buf = 0;
      char x = s[iS];
      for (int j=0; j <= endIndexT; j++) {
        if (t[j] == x){
          //both t[j] and x are considered part of the subsequence match, hence
          //subtract 2 from the remainingMatchLength
          buf += kernelHelperLP(n-1,s,iS-1,t,j-1,remainingMatchLength-2);
        }
      }
      result += buf * m_powersOflambda[2];
    }
    return result;
  }

  /**
   * helper function for the evaluation of the kernel (K'n) using lambda pruning
   *
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @param remainingMatchLength the number of characters that may still be
   * used 
   * for matching (i.e. gaps + matches in both strings)
   * @return a partial result for K 
   */
  protected double kernelHelperLP (int n, char[] s, int endIndexS,char[] t,
    int endIndexT,int remainingMatchLength) {

    //see code docs in kernelHelper()
    if (n == 0) {
      return getReturnValue(n);      

    }
    //see code docs in kernelHelper()
    if (Math.min(endIndexS+1,endIndexT +1) < n) {;
      return getReturnValue(n);
    }

    //lambda pruning check
    //stops recursion if the match is so long that the resulting
    //power of lambda is smaller than minLambda
    //if lambda pruning is not used, the remainingMatchLength is < 0
    //and this check never stops the recursion
    if (remainingMatchLength < 2*n) { 
      return getReturnValue(n);
    }
    int adr=0;
    if (m_useRecursionCache) {
      adr = m_multX*n+m_multY*endIndexS+m_multZ*endIndexT 
            + m_multZZ * remainingMatchLength;
      if (cachekh2K[adr % maxCache]==adr+1) { 
        return cachekh2[adr % maxCache]; 
      }
    }

    int rml = 0; //counts the remaining match length
    double result = 0;
    //see code docs in kernelHelper()
    //difference to implementation in kernelHelper:
    //*)choose different starting point, which is found counting
    //the maximal remaining match length from endIndexS.
    //*)keep track of the remaining match length, rml, which is
    //  incremented each loop
    for (int iS = (endIndexS-remainingMatchLength); iS <= endIndexS;iS++) {
      result *= m_lambda;
      result += kernelHelper2LP(n,s,iS, t, endIndexT,rml++);
    }

    if (m_useRecursionCache && endIndexS >= 0 && endIndexT >= 0 && n >= 0) {
      cachekhK[adr % maxCache]=adr+1; cachekh[adr % maxCache]=result; 
    }
    return result;
  }

  /**
   * helper function for the evaluation of the kernel (K''n) using lambda
   * pruning
   *
   * @param n the current length of the matching subsequence
   * @param s first string, as a char array
   * @param t second string, as a char array
   * @param endIndexS the portion of s currently regarded is s[1:endIndexS]
   * @param endIndexT the portion of t currently regarded is t[1:endIndexT]
   * @param remainingMatchLength the number of characters that may still be
   * used 
   * for matching (i.e. gaps + matches in both strings)
   * @return a partial result for K' 
   */
  protected double kernelHelper2LP(int n, char[] s, int endIndexS,char[] t,
    int endIndexT,int remainingMatchLength) {

    //lambda pruning check
    //stops recursion if the match is so long that the resulting
    //power of lambda is smaller than minLambda
    //if lambda pruning is not used, the remainingMatchLength is < 0
    //and this check never stops the recursion
    //if (remainingMatchLength <= 0) return 0;
    if (remainingMatchLength < 2*n) return getReturnValue(n);

    //see code docs in kernelHelper2()
    if (endIndexS <0 || endIndexT <0) return getReturnValue(n);
    int adr=0;
    if (m_useRecursionCache){
      adr = m_multX*n+m_multY*endIndexS+m_multZ*endIndexT
            + m_multZZ * remainingMatchLength;
      if (cachekh2K[adr % maxCache]==adr+1) { 
        return cachekh2[adr % maxCache]; 
      }
    }

    char x = s[endIndexS];
    if (x == t[endIndexT]) {
      double ret = 
        m_lambda 
        * (kernelHelper2LP(n,s,endIndexS,t,endIndexT-1,remainingMatchLength-1)
        + m_lambda 
        * kernelHelperLP(n-1,s,endIndexS-1,t,endIndexT-1,remainingMatchLength-2));
      if (m_useRecursionCache && endIndexS >= 0 && endIndexT >= 0 && n >= 0) {
        cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret; }
      return ret;
    }

    //see code docs in kernelHelper()
    //differences to implementation in kernelHelper():
    //*) choose a different ending point for the loop
    //   based on the remaining match length
    int i;
    int minIndex = endIndexT - remainingMatchLength;
    if (minIndex < 0) minIndex = 0;
    for (i=endIndexT; i >= minIndex;i--) {
      if (x == t[i]) {
        int skipLength = endIndexT -i;
        double ret = getPowerOfLambda(skipLength) *
          kernelHelper2LP(n,s,endIndexS,t,i,remainingMatchLength-skipLength);
        if (m_useRecursionCache && endIndexS >= 0 && endIndexT >= 0 && n >= 0) { 
          cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret; 
        }
        return ret;
      }
    }
    double ret = getReturnValue(n);
    if (m_useRecursionCache && endIndexS >= 0 && endIndexT >= 0 && n >= 0) { 
      cachekh2K[adr % maxCache]=adr+1; cachekh2[adr % maxCache]=ret; 
    }
    return ret;
  }

  /**
   * precalculates small powers of lambda to speed up the kernel evaluation
   *
   * @return		the powers
   */
  private double[] calculatePowersOfLambda(){
    double[] powers = new double[MAX_POWER_OF_LAMBDA+1];
    powers[0] = 1.0;
    double val = 1.0;
    for (int i = 1; i<=MAX_POWER_OF_LAMBDA;i++) {
      val *= m_lambda;
      powers[i] = val;
    }
    return powers;
  }

  /**
   * retrieves a power of lambda from the lambda cache or calculates it
   * directly
   *
   * @param exponent	the exponent to calculate
   * @return 		the exponent-th power of lambda
   */
  private double getPowerOfLambda(int exponent){
    if (exponent > MAX_POWER_OF_LAMBDA) 
      return Math.pow(m_lambda,exponent);
    
    if (exponent < 0) 
      throw new IllegalArgumentException(
          "only positive powers of lambda may be computed");

    return m_powersOflambda[exponent];
  }

  /**
   * initializes variables etc.
   * 
   * @param data	the data to use
   */
  protected void initVars(Instances data) {
    super.initVars(data);

    m_kernelEvals    = 0;
    // take the first string attribute
    m_strAttr        = -1;
    for (int i = 0; i < data.numAttributes(); i++) {
      if (i == data.classIndex())
	continue;
      if (data.attribute(i).type() == Attribute.STRING) {
	m_strAttr = i;
	break;
      }
    }
    m_numInsts       = m_data.numInstances();
    m_storage        = new double[m_cacheSize];
    m_keys           = new long[m_cacheSize];
    m_powersOflambda = calculatePowersOfLambda();
  }

  /** 
   * Returns the Capabilities of this kernel.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    
    result.enable(Capability.STRING_ATTRIBUTES);
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }
  
  /**
   * builds the kernel with the given data.
   * 
   * @param data	the data to base the kernel on
   * @throws Exception	if something goes wrong, e.g., the data does not
   * 			consist of one string attribute and the class
   */
  public void buildKernel(Instances data) throws Exception {
    super.buildKernel(data);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5450 $");
  }
}

