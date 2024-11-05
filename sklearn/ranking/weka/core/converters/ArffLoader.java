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
 *    ArffLoader.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.net.URL;
import java.text.ParseException;
import java.util.ArrayList;

/**
 <!-- globalinfo-start -->
 * Reads a source that is in arff (attribute relation file format) format.
 * <p/>
 <!-- globalinfo-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6339 $
 * @see Loader
 */
public class ArffLoader 
  extends AbstractFileLoader 
  implements BatchConverter, IncrementalConverter, URLSourcedLoader {

  /** for serialization */
  static final long serialVersionUID = 2726929550544048587L;
  
  /** the file extension */
  public static String FILE_EXTENSION = Instances.FILE_EXTENSION;
  public static String FILE_EXTENSION_COMPRESSED = FILE_EXTENSION + ".gz";

  /** the url */
  protected String m_URL = "http://";

  /** The reader for the source file. */
  protected transient Reader m_sourceReader = null;

  /** The parser for the ARFF file */
  protected transient ArffReader m_ArffReader = null;
  
  /**
   * Reads data from an ARFF file, either in incremental or batch mode. <p/>
   * 
   * Typical code for batch usage:
   * <pre>
   * BufferedReader reader = new BufferedReader(new FileReader("/some/where/file.arff"));
   * ArffReader arff = new ArffReader(reader);
   * Instances data = arff.getData();
   * data.setClassIndex(data.numAttributes() - 1);
   * </pre>
   * 
   * Typical code for incremental usage:
   * <pre>
   * BufferedReader reader = new BufferedReader(new FileReader("/some/where/file.arff"));
   * ArffReader arff = new ArffReader(reader, 1000);
   * Instances data = arff.getStructure();
   * data.setClassIndex(data.numAttributes() - 1);
   * Instance inst;
   * while ((inst = arff.readInstance(data)) != null) {
   *   data.add(inst);
   * }
   * </pre>
   * 
   * @author  Eibe Frank (eibe@cs.waikato.ac.nz)
   * @author  Len Trigg (trigg@cs.waikato.ac.nz)
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 6339 $
   */
  public static class ArffReader
    implements RevisionHandler {

    /** the tokenizer for reading the stream */
    protected StreamTokenizer m_Tokenizer;
    
    /** Buffer of values for sparse instance */
    protected double[] m_ValueBuffer;

    /** Buffer of indices for sparse instance */
    protected int[] m_IndicesBuffer;

    /** the actual data */
    protected Instances m_Data;

    /** the number of lines read so far */
    protected int m_Lines;
    
    /**
     * Reads the data completely from the reader. The data can be accessed
     * via the <code>getData()</code> method.
     * 
     * @param reader		the reader to use
     * @throws IOException	if something goes wrong
     * @see			#getData()
     */
    public ArffReader(Reader reader) throws IOException {
      m_Tokenizer = new StreamTokenizer(reader);
      initTokenizer();

      readHeader(1000);
      initBuffers();
      
      Instance inst;
      while ((inst = readInstance(m_Data)) != null) {
        m_Data.add(inst);
      };
      
      compactify();
    }
    
    /**
     * Reads only the header and reserves the specified space for instances.
     * Further instances can be read via <code>readInstance()</code>.
     * 
     * @param reader			the reader to use
     * @param capacity 			the capacity of the new dataset 
     * @throws IOException		if something goes wrong
     * @throws IllegalArgumentException	if capacity is negative
     * @see				#getStructure()
     * @see				#readInstance(Instances)
     */
    public ArffReader(Reader reader, int capacity) throws IOException {
      if (capacity < 0)
	throw new IllegalArgumentException("Capacity has to be positive!");

      m_Tokenizer = new StreamTokenizer(reader);
      initTokenizer();

      readHeader(capacity);
      initBuffers();
    }
    
    /**
     * Reads the data without header according to the specified template.
     * The data can be accessed via the <code>getData()</code> method.
     * 
     * @param reader		the reader to use
     * @param template		the template header
     * @param lines		the lines read so far
     * @throws IOException	if something goes wrong
     * @see			#getData()
     */
    public ArffReader(Reader reader, Instances template, int lines) throws IOException {
      this(reader, template, lines, 100);

      Instance inst;
      while ((inst = readInstance(m_Data)) != null) {
        m_Data.add(inst);
      };

      compactify();
    }
    
    /**
     * Initializes the reader without reading the header according to the 
     * specified template. The data must be read via the 
     * <code>readInstance()</code> method.
     * 
     * @param reader		the reader to use
     * @param template		the template header
     * @param lines		the lines read so far
     * @param capacity 		the capacity of the new dataset 
     * @throws IOException	if something goes wrong
     * @see			#getData()
     */
    public ArffReader(Reader reader, Instances template, int lines, int capacity) throws IOException {
      m_Lines     = lines;
      m_Tokenizer = new StreamTokenizer(reader);
      initTokenizer();

      m_Data = new Instances(template, capacity);
      initBuffers();
    }

    /**
     * initializes the buffers for sparse instances to be read
     * 
     * @see			#m_ValueBuffer
     * @see			#m_IndicesBuffer
     */
    protected void initBuffers() {
      m_ValueBuffer = new double[m_Data.numAttributes()];
      m_IndicesBuffer = new int[m_Data.numAttributes()];
    }
    
    /**
     * compactifies the data
     */
    protected void compactify() {
      if (m_Data != null)
        m_Data.compactify();
    }
    
    /**
     * Throws error message with line number and last token read.
     *
     * @param msg 		the error message to be thrown
     * @throws IOException 	containing the error message
     */
    protected void errorMessage(String msg) throws IOException {
      String str = msg + ", read " + m_Tokenizer.toString();
      if (m_Lines > 0) {
	int line = Integer.parseInt(str.replaceAll(".* line ", ""));
	str = str.replaceAll(" line .*", " line " + (m_Lines + line - 1));
      }
      throw new IOException(str);
    }

    /**
     * returns the current line number
     * 
     * @return			the current line number
     */
    public int getLineNo() {
      return m_Lines + m_Tokenizer.lineno();
    }
    
    /**
     * Gets next token, skipping empty lines.
     *
     * @throws IOException 	if reading the next token fails
     */
    protected void getFirstToken() throws IOException {
      while (m_Tokenizer.nextToken() == StreamTokenizer.TT_EOL) {};
      
      if ((m_Tokenizer.ttype == '\'') ||
  	(m_Tokenizer.ttype == '"')) {
        m_Tokenizer.ttype = StreamTokenizer.TT_WORD;
      } else if ((m_Tokenizer.ttype == StreamTokenizer.TT_WORD) &&
  	       (m_Tokenizer.sval.equals("?"))){
        m_Tokenizer.ttype = '?';
      }
    }

    /**
     * Gets index, checking for a premature and of line.
     *
     * @throws IOException 	if it finds a premature end of line
     */
    protected void getIndex() throws IOException {
      if (m_Tokenizer.nextToken() == StreamTokenizer.TT_EOL) {
        errorMessage("premature end of line");
      }
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        errorMessage("premature end of file");
      }
    }
    
    /**
     * Gets token and checks if its end of line.
     *
     * @param endOfFileOk 	whether EOF is OK
     * @throws IOException 	if it doesn't find an end of line
     */
    protected void getLastToken(boolean endOfFileOk) throws IOException {
      if ((m_Tokenizer.nextToken() != StreamTokenizer.TT_EOL) &&
  	((m_Tokenizer.ttype != StreamTokenizer.TT_EOF) || !endOfFileOk)) {
        errorMessage("end of line expected");
      }
    }

    /**
     * Gets the value of an instance's weight (if one exists)
     *
     * @return the value of the instance's weight, or NaN
     * if no weight has been supplied in the file
     */
    protected double getInstanceWeight() throws IOException {
      double weight = Double.NaN;
      m_Tokenizer.nextToken();
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOL || 
          m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        return weight;
      }
      // see if we can read an instance weight
      //      m_Tokenizer.pushBack();
      if (m_Tokenizer.ttype == '{') {
        m_Tokenizer.nextToken();
        String weightS = m_Tokenizer.sval;
        // try to parse weight as a double
        try {
          weight = Double.parseDouble(weightS);
        } catch (NumberFormatException e) {
          // quietly ignore
          return weight;
        }
        // see if we have the closing brace
        m_Tokenizer.nextToken();
        if (m_Tokenizer.ttype != '}') {
          errorMessage("Problem reading instance weight");
        }
      }
      return weight;
    }

    /**
     * Gets next token, checking for a premature and of line.
     *
     * @throws IOException 	if it finds a premature end of line
     */
    protected void getNextToken() throws IOException {
      if (m_Tokenizer.nextToken() == StreamTokenizer.TT_EOL) {
        errorMessage("premature end of line");
      }
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        errorMessage("premature end of file");
      } else if ((m_Tokenizer.ttype == '\'') ||
  	       (m_Tokenizer.ttype == '"')) {
        m_Tokenizer.ttype = StreamTokenizer.TT_WORD;
      } else if ((m_Tokenizer.ttype == StreamTokenizer.TT_WORD) &&
  	       (m_Tokenizer.sval.equals("?"))){
        m_Tokenizer.ttype = '?';
      }
    }
  	
    /**
     * Initializes the StreamTokenizer used for reading the ARFF file.
     */
    protected void initTokenizer(){
      m_Tokenizer.resetSyntax();         
      m_Tokenizer.whitespaceChars(0, ' ');    
      m_Tokenizer.wordChars(' '+1,'\u00FF');
      m_Tokenizer.whitespaceChars(',',',');
      m_Tokenizer.commentChar('%');
      m_Tokenizer.quoteChar('"');
      m_Tokenizer.quoteChar('\'');
      m_Tokenizer.ordinaryChar('{');
      m_Tokenizer.ordinaryChar('}');
      m_Tokenizer.eolIsSignificant(true);
    }
    
    /**
     * Reads a single instance using the tokenizer and returns it. 
     *
     * @param structure 	the dataset header information, will get updated 
     * 				in case of string or relational attributes
     * @return 			null if end of file has been reached
     * @throws IOException 	if the information is not read 
     * successfully
     */ 
    public Instance readInstance(Instances structure) throws IOException {
      return readInstance(structure, true);
    }
    
    /**
     * Reads a single instance using the tokenizer and returns it. 
     *
     * @param structure 	the dataset header information, will get updated 
     * 				in case of string or relational attributes
     * @param flag 		if method should test for carriage return after 
     * 				each instance
     * @return 			null if end of file has been reached
     * @throws IOException 	if the information is not read 
     * successfully
     */ 
    public Instance readInstance(Instances structure, boolean flag) throws IOException {
      return getInstance(structure, flag);
    }
    
    /**
     * Reads a single instance using the tokenizer and returns it. 
     *
     * @param structure 	the dataset header information, will get updated 
     * 				in case of string or relational attributes
     * @param flag 		if method should test for carriage return after 
     * 				each instance
     * @return 			null if end of file has been reached
     * @throws IOException 	if the information is not read 
     * 				successfully
     */ 
    protected Instance getInstance(Instances structure, boolean flag) throws IOException {
      m_Data = structure;
      
      // Check if any attributes have been declared.
      if (m_Data.numAttributes() == 0) {
        errorMessage("no header information available");
      }

      // Check if end of file reached.
      getFirstToken();
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        return null;
      }
      
      // Parse instance
      if (m_Tokenizer.ttype == '{') {
        return getInstanceSparse(flag);
      } else {
        return getInstanceFull(flag);
      }
    }

    /**
     * Reads a single instance using the tokenizer and returns it.
     *
     * @param flag 		if method should test for carriage return after 
     * 				each instance
     * @return 			null if end of file has been reached
     * @throws IOException 	if the information is not read 
     * 				successfully
     */ 
    protected Instance getInstanceSparse(boolean flag) throws IOException {
      int valIndex, numValues = 0, maxIndex = -1;
      
      // Get values
      do {
        // Get index
        getIndex();
        if (m_Tokenizer.ttype == '}') {
  	break;
        }
   
        // Is index valid?
        try{
  	m_IndicesBuffer[numValues] = Integer.valueOf(m_Tokenizer.sval).intValue();
        } catch (NumberFormatException e) {
  	errorMessage("index number expected");
        }
        if (m_IndicesBuffer[numValues] <= maxIndex) {
  	errorMessage("indices have to be ordered");
        }
        if ((m_IndicesBuffer[numValues] < 0) || 
  	  (m_IndicesBuffer[numValues] >= m_Data.numAttributes())) {
  	errorMessage("index out of bounds");
        }
        maxIndex = m_IndicesBuffer[numValues];

        // Get value;
        getNextToken();

        // Check if value is missing.
        if  (m_Tokenizer.ttype == '?') {
          m_ValueBuffer[numValues] = Utils.missingValue();
        } else {

  	// Check if token is valid.
  	if (m_Tokenizer.ttype != StreamTokenizer.TT_WORD) {
  	  errorMessage("not a valid value");
  	}
          switch (m_Data.attribute(m_IndicesBuffer[numValues]).type()) {
            case Attribute.NOMINAL:
              // Check if value appears in header.
              valIndex = 
                m_Data.attribute(m_IndicesBuffer[numValues]).indexOfValue(m_Tokenizer.sval);
              if (valIndex == -1) {
                errorMessage("nominal value not declared in header");
              }
              m_ValueBuffer[numValues] = (double)valIndex;
              break;
  	case Attribute.NUMERIC:
  	  // Check if value is really a number.
  	  try{
  	    m_ValueBuffer[numValues] = Double.valueOf(m_Tokenizer.sval).
  	      doubleValue();
  	  } catch (NumberFormatException e) {
  	    errorMessage("number expected");
  	  }
            break;
  	case Attribute.STRING:
  	  m_ValueBuffer[numValues] = 
  	    m_Data.attribute(m_IndicesBuffer[numValues]).addStringValue(m_Tokenizer.sval);
            break;
          case Attribute.DATE:
            try {
              m_ValueBuffer[numValues] = 
                m_Data.attribute(m_IndicesBuffer[numValues]).parseDate(m_Tokenizer.sval);
            } catch (ParseException e) {
              errorMessage("unparseable date: " + m_Tokenizer.sval);
            }
            break;
          case Attribute.RELATIONAL:
            try {
              ArffReader arff = new ArffReader(new StringReader(m_Tokenizer.sval), m_Data.attribute(m_IndicesBuffer[numValues]).relation(), 0);
              Instances data = arff.getData();
              m_ValueBuffer[numValues] = m_Data.attribute(m_IndicesBuffer[numValues]).addRelation(data);
            }
            catch (Exception e) {
              throw new IOException(e.toString() + " of line " + getLineNo());
            }
            break;
          default:
            errorMessage("unknown attribute type in column " + m_IndicesBuffer[numValues]);
  	}
        }
        numValues++;
      } while (true);

      double weight = 1.0;
      if (flag) {
        // check for an instance weight
        weight = getInstanceWeight();
        if (!Double.isNaN(weight)) {
          getLastToken(true);
        } else {
          weight = 1.0;
        }        
      }
        
      // Add instance to dataset
      double[] tempValues = new double[numValues];
      int[] tempIndices = new int[numValues];
      System.arraycopy(m_ValueBuffer, 0, tempValues, 0, numValues);
      System.arraycopy(m_IndicesBuffer, 0, tempIndices, 0, numValues);
      Instance inst = new SparseInstance(weight, tempValues, tempIndices, m_Data.numAttributes());
      inst.setDataset(m_Data);
      
      return inst;
    }

    /**
     * Reads a single instance using the tokenizer and returns it.
     *
     * @param flag 		if method should test for carriage return after 
     * 				each instance
     * @return 			null if end of file has been reached
     * @throws IOException 	if the information is not read 
     * 				successfully
     */ 
    protected Instance getInstanceFull(boolean flag) throws IOException {
      double[] instance = new double[m_Data.numAttributes()];
      int index;
      
      // Get values for all attributes.
      for (int i = 0; i < m_Data.numAttributes(); i++){
        // Get next token
        if (i > 0) {
  	getNextToken();
        }
              
        // Check if value is missing.
        if  (m_Tokenizer.ttype == '?') {
  	instance[i] = Utils.missingValue();
        } else {

  	// Check if token is valid.
  	if (m_Tokenizer.ttype != StreamTokenizer.TT_WORD) {
  	  errorMessage("not a valid value");
  	}
          switch (m_Data.attribute(i).type()) {
          case Attribute.NOMINAL:
  	  // Check if value appears in header.
  	  index = m_Data.attribute(i).indexOfValue(m_Tokenizer.sval);
  	  if (index == -1) {
  	    errorMessage("nominal value not declared in header");
  	  }
  	  instance[i] = (double)index;
            break;
  	case Attribute.NUMERIC:
  	  // Check if value is really a number.
  	  try{
  	    instance[i] = Double.valueOf(m_Tokenizer.sval).
  	      doubleValue();
  	  } catch (NumberFormatException e) {
  	    errorMessage("number expected");
  	  }
            break;
  	case Attribute.STRING:
  	  instance[i] = m_Data.attribute(i).addStringValue(m_Tokenizer.sval);
            break;
          case Attribute.DATE:
            try {
              instance[i] = m_Data.attribute(i).parseDate(m_Tokenizer.sval);
            } catch (ParseException e) {
              errorMessage("unparseable date: " + m_Tokenizer.sval);
            }
            break;
          case Attribute.RELATIONAL:
            try {
              ArffReader arff = new ArffReader(new StringReader(m_Tokenizer.sval), m_Data.attribute(i).relation(), 0);
              Instances data = arff.getData();
              instance[i] = m_Data.attribute(i).addRelation(data);
            }
            catch (Exception e) {
              throw new IOException(e.toString() + " of line " + getLineNo());
            }
            break;
          default:
            errorMessage("unknown attribute type in column " + i);
  	}
        }
      }
      
      double weight = 1.0;
      if (flag) {
        // check for an instance weight
        weight = getInstanceWeight();
        if (!Double.isNaN(weight)) {
          getLastToken(true);
        } else {
          weight = 1.0;
        }
      }
        
      // Add instance to dataset
      Instance inst = new DenseInstance(weight, instance);
      inst.setDataset(m_Data);
      
      return inst;
    }

    /**
     * Reads and stores header of an ARFF file.
     *
     * @param capacity 		the number of instances to reserve in the data 
     * 				structure
     * @throws IOException 	if the information is not read 
     * 				successfully
     */ 
    protected void readHeader(int capacity) throws IOException {
      m_Lines = 0;
      String relationName = "";
      
      // Get name of relation.
      getFirstToken();
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        errorMessage("premature end of file");
      }
      if (Instances.ARFF_RELATION.equalsIgnoreCase(m_Tokenizer.sval)) {
        getNextToken();
        relationName = m_Tokenizer.sval;
        getLastToken(false);
      } else {
        errorMessage("keyword " + Instances.ARFF_RELATION + " expected");
      }

      // Create vectors to hold information temporarily.
      ArrayList<Attribute> attributes = new ArrayList<Attribute>();
   
      // Get attribute declarations.
      getFirstToken();
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
        errorMessage("premature end of file");
      }

      while (Attribute.ARFF_ATTRIBUTE.equalsIgnoreCase(m_Tokenizer.sval)) {
        attributes = parseAttribute(attributes);
      }

      // Check if data part follows. We can't easily check for EOL.
      if (!Instances.ARFF_DATA.equalsIgnoreCase(m_Tokenizer.sval)) {
        errorMessage("keyword " + Instances.ARFF_DATA + " expected");
      }
      
      // Check if any attributes have been declared.
      if (attributes.size() == 0) {
        errorMessage("no attributes declared");
      }
      
      m_Data = new Instances(relationName, attributes, capacity);
    }

    /**
     * Parses the attribute declaration.
     *
     * @param attributes 		the current attributes vector
     * @return 			the new attributes vector
     * @throws IOException 	if the information is not read 
     * 				successfully
     */
    protected ArrayList<Attribute> parseAttribute(ArrayList<Attribute> attributes) throws IOException {
      String attributeName;
      ArrayList<String> attributeValues;

      // Get attribute name.
      getNextToken();
      attributeName = m_Tokenizer.sval;
      getNextToken();
      
      // Check if attribute is nominal.
      if (m_Tokenizer.ttype == StreamTokenizer.TT_WORD) {
        
        // Attribute is real, integer, or string.
        if (m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_REAL) ||
            m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_INTEGER) ||
            m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_NUMERIC)) {
          attributes.add(new Attribute(attributeName, attributes.size()));
          readTillEOL();
        } else if (m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_STRING)) {
          attributes.add(new Attribute(attributeName, (ArrayList<String>)null,
                attributes.size()));
          readTillEOL();
        } else if (m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_DATE)) {
          String format = null;
          if (m_Tokenizer.nextToken() != StreamTokenizer.TT_EOL) {
            if ((m_Tokenizer.ttype != StreamTokenizer.TT_WORD) &&
                (m_Tokenizer.ttype != '\'') &&
                (m_Tokenizer.ttype != '\"')) {
              errorMessage("not a valid date format");
            }
            format = m_Tokenizer.sval;
            readTillEOL();
          } else {
            m_Tokenizer.pushBack();
          }
          attributes.add(new Attribute(attributeName, format, attributes.size()));
          
        } else if (m_Tokenizer.sval.equalsIgnoreCase(Attribute.ARFF_ATTRIBUTE_RELATIONAL)) {
          readTillEOL();
          
          // Read attributes for subrelation
          // First, save current set of attributes
          ArrayList<Attribute> atts = attributes;
          attributes = new ArrayList<Attribute>();
          
          // Now, read attributes until we hit end of declaration of relational value
          getFirstToken();
          if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF) {
            errorMessage("premature end of file");
          }
          do {
            if (Attribute.ARFF_ATTRIBUTE.equalsIgnoreCase(m_Tokenizer.sval)) {
              attributes = parseAttribute(attributes);
            } else if (Attribute.ARFF_END_SUBRELATION.equalsIgnoreCase(m_Tokenizer.sval)) {
              getNextToken();
              if (!attributeName.equalsIgnoreCase(m_Tokenizer.sval)) {
                errorMessage("declaration of subrelation " + attributeName + 
                      " must be terminated by " + "@end " + attributeName);
              }
              break;
            } else {
              errorMessage("declaration of subrelation " + attributeName + 
                    " must be terminated by " + "@end " + attributeName);
            }
          } while (true);
          
          // Make relation and restore original set of attributes
          Instances relation = new Instances(attributeName, attributes, 0);
          attributes = atts;
          attributes.add(new Attribute(attributeName, relation, attributes.size()));
        } else {
          errorMessage("no valid attribute type or invalid "+
                "enumeration");
        }
      } else {
        
        // Attribute is nominal.
        attributeValues = new ArrayList<String>();
        m_Tokenizer.pushBack();
        
        // Get values for nominal attribute.
        if (m_Tokenizer.nextToken() != '{') {
          errorMessage("{ expected at beginning of enumeration");
        }
        while (m_Tokenizer.nextToken() != '}') {
          if (m_Tokenizer.ttype == StreamTokenizer.TT_EOL) {
            errorMessage("} expected at end of enumeration");
          } else {
            attributeValues.add(m_Tokenizer.sval);
          }
        }
        attributes.add(new Attribute(attributeName, attributeValues,
              attributes.size()));
      }
      getLastToken(false);
      getFirstToken();
      if (m_Tokenizer.ttype == StreamTokenizer.TT_EOF)
        errorMessage("premature end of file");
      
      return attributes;
    }

    /**
     * Reads and skips all tokens before next end of line token.
     *
     * @throws IOException 	in case something goes wrong
     */
    protected void readTillEOL() throws IOException {
      while (m_Tokenizer.nextToken() != StreamTokenizer.TT_EOL) {};
      
      m_Tokenizer.pushBack();
    }

    /**
     * Returns the header format
     * 
     * @return			the header format
     */
    public Instances getStructure() {
      return new Instances(m_Data, 0);
    }
    
    /**
     * Returns the data that was read
     * 
     * @return			the data
     */
    public Instances getData() {
      return m_Data;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6339 $");
    }
  }

  /**
   * Returns a string describing this Loader
   * @return a description of the Loader suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Reads a source that is in arff (attribute relation file format) "
      +"format. ";
  }

  /**
   * Get the file extension used for arff files
   *
   * @return the file extension
   */
  public String getFileExtension() {
    return FILE_EXTENSION;
  }

  /**
   * Gets all the file extensions used for this type of file
   *
   * @return the file extensions
   */
  public String[] getFileExtensions() {
    return new String[]{FILE_EXTENSION, FILE_EXTENSION_COMPRESSED};
  }

  /**
   * Returns a description of the file type.
   *
   * @return a short file description
   */
  public String getFileDescription() {
    return "Arff data files";
  }

  /**
   * Resets the Loader ready to read a new data set or the
   * same data set again.
   * 
   * @throws IOException if something goes wrong
   */
  public void reset() throws IOException {
    m_structure = null;
    setRetrieval(NONE);
    
    if (m_File != null && !(new File(m_File).isDirectory())) {
      setFile(new File(m_File));
    } else if (m_URL != null && !m_URL.equals("http://")) {
      setURL(m_URL);
    }
  }

  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied url.
   *
   * @param url the source url.
   * @throws IOException if an error occurs
   */
  public void setSource(URL url) throws IOException {
    m_structure = null;
    setRetrieval(NONE);
    
    setSource(url.openStream());

    m_URL = url.toString();
    // make sure that the file is null so that any calls to
    // reset() work properly
    m_File = null;
  }
  

  /**
   * get the File specified as the source
   *
   * @return the source file
   */
  public File retrieveFile() {
    return new File(m_File);
  }

  /**
   * sets the source File
   *
   * @param file the source file
   * @throws IOException if an error occurs
   */
  public void setFile(File file) throws IOException {
    m_File = file.getPath();
    setSource(file);
  }

  /**
   * Set the url to load from
   *
   * @param url the url to load from
   * @throws IOException if the url can't be set.
   */
  public void setURL(String url) throws IOException {
    m_URL = url;
    setSource(new URL(url));
  }

  /**
   * Return the current url
   *
   * @return the current url
   */
  public String retrieveURL() {
    return m_URL;
  }

  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied InputStream.
   *
   * @param in the source InputStream.
   * @throws IOException always thrown.
   */
  public void setSource(InputStream in) throws IOException {
    m_File = (new File(System.getProperty("user.dir"))).getAbsolutePath();
    m_URL = "http://";

    m_sourceReader = new BufferedReader(new InputStreamReader(in));
  }

  /**
   * Determines and returns (if possible) the structure (internally the 
   * header) of the data set as an empty set of instances.
   *
   * @return the structure of the data set as an empty set of Instances
   * @throws IOException if an error occurs
   */
  public Instances getStructure() throws IOException {

    if (m_sourceReader == null) {
      throw new IOException("No source has been specified");
    }

    if (m_structure == null) {
      try {
	m_ArffReader = new ArffReader(m_sourceReader, 1);
	m_structure  = m_ArffReader.getStructure();
      } catch (Exception ex) {
	throw new IOException("Unable to determine structure as arff (Reason: " + ex.toString() + ").");
      }
    }

    return new Instances(m_structure, 0);
  }

  /**
   * Return the full data set. If the structure hasn't yet been determined
   * by a call to getStructure then method should do so before processing
   * the rest of the data set.
   *
   * @return the structure of the data set as an empty set of Instances
   * @throws IOException if there is no source or parsing fails
   */
  public Instances getDataSet() throws IOException {

    if (m_sourceReader == null) {
      throw new IOException("No source has been specified");
    }
    if (getRetrieval() == INCREMENTAL) {
      throw new IOException("Cannot mix getting Instances in both incremental and batch modes");
    }
    setRetrieval(BATCH);
    if (m_structure == null) {
      getStructure();
    }

    // Read all instances
    Instance inst;
    while ((inst = m_ArffReader.readInstance(m_structure)) != null)
      m_structure.add(inst);
    
    Instances readIn = new Instances(m_structure);

    // close the stream
    m_sourceReader.close();
    
    return readIn;
  }

  /**
   * Read the data set incrementally---get the next instance in the data 
   * set or returns null if there are no
   * more instances to get. If the structure hasn't yet been 
   * determined by a call to getStructure then method should do so before
   * returning the next instance in the data set.
   *
   * @param structure the dataset header information, will get updated in 
   * case of string or relational attributes
   * @return the next instance in the data set as an Instance object or null
   * if there are no more instances to be read
   * @throws IOException if there is an error during parsing
   */
  public Instance getNextInstance(Instances structure) throws IOException {

    m_structure = structure;

    if (getRetrieval() == BATCH) {
      throw new IOException("Cannot mix getting Instances in both incremental and batch modes");
    }
    setRetrieval(INCREMENTAL);

    Instance current = m_ArffReader.readInstance(m_structure);
    if (current == null) {
      try {
        // close the stream
        m_sourceReader.close();
        //        reset();
      } catch (Exception ex) {
        ex.printStackTrace();
      }
    }
    return current;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6339 $");
  }

  /**
   * Main method.
   *
   * @param args should contain the name of an input file.
   */
  public static void main(String [] args) {
    runFileLoader(new ArffLoader(), args);
  }
}
