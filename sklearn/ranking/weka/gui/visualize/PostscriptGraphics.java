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
 *    PostscriptGraphics.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import java.awt.image.renderable.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.font.*;
import java.io.*;
import java.util.*;
import java.awt.image.*;
import java.text.*;


/** 
 * The PostscriptGraphics class extends the Graphics2D class to 
 * produce an encapsulated postscript file rather than on-screen display.
 * <p>
 * Currently only a small (but useful) subset of Graphics methods have been 
 * implemented. 
 * To handle the ability to Clone a Graphics object, the graphics state of the 
 * eps is set from the graphics state of the local PostscriptGraphics before output.
 * To use, create a PostscriptGraphics object, and pass it to the PaintComponent
 * method of a JComponent.
 * <p>
 * If necessary additional font replacements can be inserted, since some fonts 
 * might be displayed incorrectly.
 *
 * @see #addPSFontReplacement(String, String)
 * @see #m_PSFontReplacement
 * @author Dale Fletcher (dale@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.5 $
 */

public class PostscriptGraphics extends Graphics2D {
  
  /**
   * This inner class is used to maintain the graphics states of the PostScript
   * file and graphics context.
   */
  private class GraphicsState {
    /** The current pen color */
    protected Color m_currentColor;
    
    /** The current Font */
    protected Font m_currentFont;
    
    /** The current Stroke (not yet used) */ 
    protected Stroke m_currentStroke;
    
    /** x,y Translation */
    protected int m_xOffset;
    protected int m_yOffset;
    
    /** the scale factors */
    protected double m_xScale;
    protected double m_yScale;
    
    /**
     * Create a new GraphicsState with default values.
     */
    GraphicsState(){
      m_currentColor  = Color.white;
      m_currentFont   = new Font ("Courier", Font.PLAIN, 11);
      m_currentStroke = new BasicStroke();
      m_xOffset       = 0;
      m_yOffset       = 0;
      m_xScale        = 1.0;
      m_yScale        = 1.0;
    }
    
    /**
     * Create a new cloned GraphicsState
     *
     * @param copy The GraphicsState to clone
     */
    GraphicsState(GraphicsState copy){
      m_currentColor  = copy.m_currentColor;
      m_currentFont   = copy.m_currentFont;
      m_currentStroke = copy.m_currentStroke;
      m_xOffset       = copy.m_xOffset;
      m_yOffset       = copy.m_yOffset;
      m_xScale        = copy.m_xScale;
      m_yScale        = copy.m_yScale;
    }
    
    /* Stroke Methods */
    protected Stroke getStroke(){
      return m_currentStroke;
    }
    
    protected void setStroke(Stroke s){
      m_currentStroke = s;
    }
    
    /* Font Methods */
    protected Font getFont(){
      return m_currentFont;
    }
    
    protected void setFont(Font f){
      m_currentFont = f;
    }
    
    /* Color Methods */
    protected Color getColor(){
      return m_currentColor;
    }
    
    protected void setColor(Color c){
      m_currentColor = c;
    }
    
    /* Translation methods */
    protected void setXOffset(int xo){
      m_xOffset = xo;
    }
    
    protected void setYOffset(int yo){
      m_yOffset = yo;
    }
    
    protected int getXOffset(){
      return m_xOffset;
    }
    
    protected int getYOffset(){
      return m_yOffset;
    }
    
    protected void setXScale(double x){
      m_xScale = x;
    }
    
    protected void setYScale(double y){
      m_yScale = y;
    }
    
    protected double getXScale(){
      return m_xScale;
    }
    
    protected double getYScale(){
      return m_yScale;
    }
  }
  
  /** The bounding box of the output */
  protected Rectangle m_extent;
  
  /** The output file */
  protected PrintStream m_printstream;
  
  /** The current global PostScript graphics state for all cloned objects */
  protected GraphicsState m_psGraphicsState;
  
  /** The current local graphics state for this PostscriptGraphics object */
  protected GraphicsState m_localGraphicsState;
  
  /** whether to print some debug information */
  protected final static boolean DEBUG = false;
  
  /** the font replacement */
  protected static Hashtable m_PSFontReplacement;
  
  /** output if we're in debug mode */
  static {
    if (DEBUG)
      System.err.println(PostscriptGraphics.class.getName() + ": DEBUG ON");
    
    // get font replacements
    m_PSFontReplacement = new Hashtable();
    m_PSFontReplacement.put("SansSerif.plain", "Helvetica.plain");   // SansSerif.plain is displayed as Courier in GV???
    m_PSFontReplacement.put("Dialog.plain", "Helvetica.plain");  // dialog is a Sans Serif font, but GV displays it as Courier???
    m_PSFontReplacement.put("Microsoft Sans Serif", "Helvetica.plain");  // MS Sans Serif is a Sans Serif font (hence the name!), but GV displays it as Courier???
    m_PSFontReplacement.put("MicrosoftSansSerif", "Helvetica.plain");  // MS Sans Serif is a Sans Serif font (hence the name!), but GV displays it as Courier???
  }
  
  /** 
   * Constructor
   * Creates a new PostscriptGraphics object, given dimensions and 
   * output file.
   *
   * @param width The width of eps in points.
   * @param height The height of eps in points.
   * @param os File to send postscript to.
   */
  public PostscriptGraphics(int width, int height, OutputStream os ){
    
    m_extent             = new Rectangle(0, 0, height, width);
    m_printstream        = new PrintStream(os);
    m_localGraphicsState = new GraphicsState();
    m_psGraphicsState    = new GraphicsState();
    
    Header();
  }
  
  /** 
   * Creates a new cloned PostscriptGraphics object.
   *
   * @param copy The PostscriptGraphics object to clone.
   */
  PostscriptGraphics(PostscriptGraphics copy){
    
    m_extent             = new Rectangle(copy.m_extent);
    m_printstream        = copy.m_printstream;
    m_localGraphicsState = new GraphicsState(copy.m_localGraphicsState); // create a local copy of the current state
    m_psGraphicsState    = copy.m_psGraphicsState; // link to global state of eps file
  }
  
  /**
   * Finalizes output file.
   */
  public void finished(){
    m_printstream.flush();
  }
  
  /**
   * Output postscript header to PrintStream, including helper macros.
   */
  private void Header(){
    m_printstream.println("%!PS-Adobe-3.0 EPSF-3.0");
    m_printstream.println("%%BoundingBox: 0 0 " + xScale(m_extent.width) + " " + yScale(m_extent.height));
    m_printstream.println("%%CreationDate: " + Calendar.getInstance().getTime());
    
    m_printstream.println("/Oval { % x y w h filled");
    m_printstream.println("gsave");
    m_printstream.println("/filled exch def /h exch def /w exch def /y exch def /x exch def");
    m_printstream.println("x w 2 div add y h 2 div sub translate");
    m_printstream.println("1 h w div scale");
    m_printstream.println("filled {0 0 moveto} if");
    m_printstream.println("0 0 w 2 div 0 360 arc");
    m_printstream.println("filled {closepath fill} {stroke} ifelse grestore} bind def");
    
    m_printstream.println("/Rect { % x y w h filled");
    m_printstream.println("/filled exch def /h exch def /w exch def /y exch def /x exch def");
    m_printstream.println("newpath ");
    m_printstream.println("x y moveto");    
    m_printstream.println("w 0 rlineto");
    m_printstream.println("0 h neg rlineto");
    m_printstream.println("w neg 0 rlineto");
    m_printstream.println("closepath");
    m_printstream.println("filled {fill} {stroke} ifelse} bind def");
    
    m_printstream.println("%%BeginProlog\n%%EndProlog");
    m_printstream.println("%%Page 1 1");
    setFont(null); // set to default
    setColor(null); // set to default
    setStroke(null); // set to default
  }
  
  /**
   * adds the PS font name to replace and its replacement in the replacement
   * hashtable
   * 
   * @param replace       the PS font name to replace
   * @param with          the PS font name to replace the font with 
   */
  public static void addPSFontReplacement(String replace, String with) {
    m_PSFontReplacement.put(replace, with);
  }
  
  /**
   * Convert Java Y coordinate (0 = top) to PostScript (0 = bottom)
   * Also apply current Translation
   * @param y Java Y coordinate
   * @return translated Y to postscript
   */
  private int yTransform(int y){
    return (m_extent.height - (m_localGraphicsState.getYOffset() + y));
  }
  
  /**
   * Apply current X Translation
   * @param x Java X coordinate
   * @return translated X to postscript
   */
  private int xTransform(int x){
    return (m_localGraphicsState.getXOffset() + x);
  }
  
  /**
   * scales the given number with the provided scale factor
   */
  private int doScale(int number, double factor) {
    return (int) StrictMath.round(number * factor);
  }
  
  /**
   * scales the given x value with current x scale factor
   */
  private int xScale(int x) {
    return doScale(x, m_localGraphicsState.getXScale());
  }
  
  /**
   * scales the given y value with current y scale factor
   */
  private int yScale(int y) {
    return doScale(y, m_localGraphicsState.getYScale());
  }
  
  /** Set the current eps graphics state to that of the local one
   */
  private void setStateToLocal(){
    setColor(this.getColor());
    setFont(this.getFont());
    setStroke(this.getStroke());
  }
  
  /**
   * returns a two hexadecimal representation of i, if shorter than 2 chars
   * then an additional "0" is put in front   
   */
  private String toHex(int i) {
    String      result;
    
    result = Integer.toHexString(i);
    if (result.length() < 2)
      result = "0" + result;
    
    return result;
  }

  /***** overridden Graphics methods *****/  
  
  /**
   * Draw a filled rectangle with the background color.
   *
   * @param x starting x coord
   * @param y starting y coord
   * @param width rectangle width
   * @param height rectangle height
   */
  public void clearRect(int x, int y, int width, int height) {
    setStateToLocal();
    Color saveColor = getColor();
    setColor(Color.white); // background color for page
    m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " " + xScale(width) + " " + yScale(height) + " true Rect");
    setColor(saveColor);
  }
  
  /**
   * Not implemented
   */
  public void clipRect(int x, int y, int width, int height) {}
  
  /**
   * Not implemented
   */
  public void copyArea(int x, int y, int width, int height, int dx, int dy) {}
  
  /**
   * Clone a PostscriptGraphics object
   */  
  public Graphics create() {
    if (DEBUG)
      m_printstream.println("%create");
    PostscriptGraphics psg = new PostscriptGraphics(this);
    return(psg);
  }
  
  /**
   * Not implemented
   */
  public void dispose(){}
  
  /**
   * Draw an outlined rectangle with 3D effect in current pen color.
   * (Current implementation: draw simple outlined rectangle)
   *
   * @param x starting x coord
   * @param y starting y coord
   * @param width rectangle width
   * @param height rectangle height
   * @param raised True: appear raised, False: appear etched
   */
  public void draw3DRect(int x, int y, int width, int height, boolean raised){
    drawRect(x,y,width,height);
  }
  
  /**
   * Not implemented
   */
  public void drawArc(int x, int y, int width, int height, int startAngle, int arcAngle){}
  
  /**
   * simply calls drawString(String,int,int)
   * 
   * @see #drawString(String,int,int)
   */
  public void drawBytes(byte[] data, int offset, int length, int x, int y) {
    drawString(new String(data, offset, length), x, y);
  }
  
  /**
   * simply calls drawString(String,int,int)
   * 
   * @see #drawString(String,int,int)
   */
  public void drawChars(char[] data, int offset, int length, int x, int y) {
    drawString(new String(data, offset, length), x, y);
  }
  
  /**
   * calls drawImage(Image,int,int,int,int,Color,ImageObserver)
   * 
   * @see #drawImage(Image,int,int,int,int,Color,ImageObserver)
   */
  public boolean drawImage(Image img, int x, int y, Color bgcolor, ImageObserver observer){
    return drawImage(img, x, y, img.getWidth(observer), img.getHeight(observer), bgcolor, observer);
  }
  
  /**
   * calls drawImage(Image,int,int,Color,ImageObserver) with Color.WHITE as 
   * background color
   * 
   * @see #drawImage(Image,int,int,Color,ImageObserver)
   * @see Color#WHITE
   */
  public boolean drawImage(Image img, int x, int y, ImageObserver observer){
    return drawImage(img, x, y, Color.WHITE, observer);
  }
  
  /**
   * PS see http://astronomy.swin.edu.au/~pbourke/geomformats/postscript/
   * Java http://show.docjava.com:8086/book/cgij/doc/ip/graphics/SimpleImageFrame.java.html
   */
  public boolean drawImage(Image img, int x, int y, int width, int height, Color bgcolor, ImageObserver observer){
    try {
      // get data from image
      int[] pixels = new int[width * height];
      PixelGrabber grabber = new PixelGrabber(img, 0, 0, width, height, pixels, 0, width);
      grabber.grabPixels();
      ColorModel model = ColorModel.getRGBdefault();
      
      // print data to ps
      m_printstream.println("gsave");
      m_printstream.println(xTransform(xScale(x)) + " " + (yTransform(yScale(y)) - yScale(height)) + " translate");
      m_printstream.println(xScale(width) + " " + yScale(height) + " scale");
      m_printstream.println(width + " " + height + " " + "8" + " [" + width + " 0 0 " + (-height) + " 0 " + height + "]");
      m_printstream.println("{<");

      int index;
      for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
          index = i * width + j;
          m_printstream.print(toHex(model.getRed(pixels[index])));
          m_printstream.print(toHex(model.getGreen(pixels[index])));
          m_printstream.print(toHex(model.getBlue(pixels[index])));
        }
        m_printstream.println();
      }
      
      m_printstream.println(">}");
      m_printstream.println("false 3 colorimage");
      m_printstream.println("grestore");
      return true;
    }
    catch (Exception e) {
      e.printStackTrace();
      return false;
    }
  }
  
  /**
   * calls drawImage(Image,int,int,int,int,Color,ImageObserver) with the color 
   * WHITE as background
   * 
   * @see #drawImage(Image,int,int,int,int,Color,ImageObserver)
   * @see Color#WHITE
   */
  public boolean drawImage(Image img, int x, int y, int width, int height, ImageObserver observer){
    return drawImage(img, x, y, width, height, Color.WHITE, observer);
  }
  
  /**
   * Not implemented
   */
  public boolean drawImage(Image img, int dx1, int dy1, int dx2, int dy2, int sx1, int sy1, int sx2, int sy2, Color bgcolor, ImageObserver  observer){
    return false;
  }
  
  /**
   * calls drawImage(Image,int,int,int,int,int,int,int,int,Color,ImageObserver)
   * with Color.WHITE as background color
   * 
   * @see #drawImage(Image,int,int,int,int,int,int,int,int,Color,ImageObserver)
   */
  public boolean drawImage(Image img, int dx1, int dy1, int dx2, int dy2, int sx1, int sy1, int sx2, int sy2, ImageObserver observer){
    return drawImage(img, dx1, dy1, dx2, dy2, sx1, sy1, sx2, sy2, Color.WHITE, observer);
  }
  
  
  /**
   * Draw a line in current pen color.
   *
   * @param x1 starting x coord
   * @param y1 starting y coord
   * @param x2 ending x coord
   * @param y2 ending y coord
   */
  public void drawLine(int x1, int y1, int x2, int y2){
    setStateToLocal();
    m_printstream.println(xTransform(xScale(x1)) + " " + yTransform(yScale(y1)) + " moveto " + xTransform(xScale(x2)) + " " + yTransform(yScale(y2)) + " lineto stroke");
  }
  
  /**
   * Draw an Oval outline in current pen color.
   *
   * @param x x-axis center of oval
   * @param y y-axis center of oval
   * @param width oval width
   * @param height oval height
   */
  public void drawOval(int x, int y, int width, int height){
    setStateToLocal();
    m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " " + xScale(width) + " " + yScale(height) + " false Oval");
  }
  
  /**
   * Not implemented
   */
  public void drawPolygon(int[] xPoints, int[] yPoints, int nPoints){}
  
  /**
   * Not implemented
   */
  public void drawPolyline(int[] xPoints, int[] yPoints, int nPoints){}
  
  /**
   * Draw an outlined rectangle in current pen color.
   *
   * @param x starting x coord
   * @param y starting y coord
   * @param width rectangle width
   * @param height rectangle height
   */
  public void drawRect(int x, int y, int width, int height){   
    setStateToLocal();
    m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " " + xScale(width) + " " + yScale(height) + " false Rect");   
  }
  
  /**
   * Not implemented
   */
  public void drawRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight){}
  
  /**
   * Not implemented
   */
  public void drawString(AttributedCharacterIterator iterator, int x, int y){}
  
  /**
   * Escapes brackets in the string with backslashes.
   * 
   * @param s the string to escape
   * @return the escaped string
   */
  protected String escape(String s) {
    StringBuffer	result;
    int			i;
    
    result = new StringBuffer();
    
    for (i = 0; i < s.length(); i++) {
      if ( (s.charAt(i) == '(') || (s.charAt(i) == ')') )
	result.append('\\');
      result.append(s.charAt(i));
    }
    
    return result.toString();
  }
  
  /**
   * Draw text in current pen color.
   *
   * @param str Text to output
   * @param x starting x coord
   * @param y starting y coord
   */
  public void drawString(String str, int x, int y){
    setStateToLocal();
    m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " moveto" + " (" + escape(str) + ") show stroke");
  }
  
  /**
   * Draw a filled rectangle with 3D effect in current pen color.
   * (Current implementation: draw simple filled rectangle)
   *
   * @param x starting x coord
   * @param y starting y coord
   * @param width rectangle width
   * @param height rectangle height
   * @param raised True: appear raised, False: appear etched
   */
  public void fill3DRect(int x, int y, int width, int height, boolean raised){
    fillRect(x, y, width, height);
  }
  
  /**
   * Not implemented
   */
  public void fillArc(int x, int y, int width, int height, int startAngle, int arcAngle){}
  
  /**
   * Draw a filled Oval in current pen color.
   *
   * @param x x-axis center of oval
   * @param y y-axis center of oval
   * @param width oval width
   * @param height oval height
   */
  public void fillOval(int x, int y, int width, int height){
    setStateToLocal();
    m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " " + xScale(width) + " " + yScale(height) + " true Oval");
  }
  
  /**
   * Not implemented
   */
  public void fillPolygon(int[] xPoints, int[] yPoints, int nPoints){}
  
  /**
   * Not implemented
   */
  public void fillPolygon(Polygon p){}
  
  /**
   * Draw a filled rectangle in current pen color.
   *
   * @param x starting x coord
   * @param y starting y coord
   * @param width rectangle width
   * @param height rectangle height
   */
  
  public void fillRect(int x, int y, int width, int height){
    if (width == m_extent.width && height == m_extent.height) {
      clearRect(x, y, width, height); // if we're painting the entire background, just make it white
    } else {
      if (DEBUG)
        m_printstream.println("% fillRect");
      setStateToLocal();
      m_printstream.println(xTransform(xScale(x)) + " " + yTransform(yScale(y)) + " " + xScale(width) + " " + yScale(height) + " true Rect");
    }
  }
  
  /**
   * Not implemented
   */
  public void fillRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight){}
  
  /**
   * Not implemented
   */
  public void finalize(){}
  
  /**
   * Not implemented
   */
  public Shape getClip(){
    return(null);
  }
  
  /**
   * This returns the full current drawing area
   * @return full drawing area
   */
  public Rectangle getClipBounds(){
    return(new Rectangle(0, 0, m_extent.width, m_extent.height));
  }
  
  /**
   * This returns the full current drawing area
   * @return full drawing area
   */
  public Rectangle getClipBounds(Rectangle r) {
    r.setBounds(0, 0, m_extent.width, m_extent.height);
    return r;
  }
  
  /**
   * Not implemented
   */
  public Rectangle getClipRect() {return null;}
  
  /**
   * Get current pen color.
   *
   * @return current pen color.
   */
  public Color getColor(){
    return (m_localGraphicsState.getColor());
  }
  
  /**
   * Get current font.
   *
   * @return current font.
   */
  public Font getFont(){
    return (m_localGraphicsState.getFont());
  }
  
  /**
   * Get Font metrics
   *
   * @param f Font 
   * @return Font metrics.
   */
  public FontMetrics getFontMetrics(Font f){
    return(Toolkit.getDefaultToolkit().getFontMetrics(f));  
    
  }
  
  /**
   * Not implemented
   */
  public void setClip(int x, int y, int width, int height) {}
  
  /**
   * Not implemented
   */
  public void setClip(Shape clip){}
  
  /**
   * Set current pen color. Default to black if null.
   *
   * @param c new pen color.
   */
  public void setColor(Color c){
    if (c != null){
      m_localGraphicsState.setColor(c);
      if (m_psGraphicsState.getColor().equals(c)) {
        return;
      }
      m_psGraphicsState.setColor(c);
    } else {
      m_localGraphicsState.setColor(Color.black);
      m_psGraphicsState.setColor(getColor());
    }
    m_printstream.print(getColor().getRed()/255.0);
    m_printstream.print(" ");
    m_printstream.print(getColor().getGreen()/255.0);
    m_printstream.print(" ");
    m_printstream.print(getColor().getBlue()/255.0);
    m_printstream.println(" setrgbcolor");
  }
  
  /**
   * replaces the font (PS name) if necessary and returns the new name
   */
  private static String replacePSFont(String font) {
    String      result;
    
    result = font;
    
    // do we have to replace it? -> same style, size
    if (m_PSFontReplacement.containsKey(font)) {
      result = m_PSFontReplacement.get(font).toString();
      if (DEBUG)
        System.out.println("switched font from '" + font + "' to '" + result +  "'");
    }
    
    return result;
  }
  
  /**
   * Set current font. Default to Plain Courier 11 if null.
   *
   * @param font new font.
   */
  public void setFont(Font font){
    
    if (font != null){
      m_localGraphicsState.setFont(font);
      if (   font.getName().equals(m_psGraphicsState.getFont().getName())
          && (m_psGraphicsState.getFont().getStyle() == font.getStyle())
          && (m_psGraphicsState.getFont().getSize() == yScale(font.getSize())))
        return;
      m_psGraphicsState.setFont(new Font(font.getName(), font.getStyle(), yScale(getFont().getSize())));
    } 
    else {
      m_localGraphicsState.setFont(new Font ("Courier", Font.PLAIN, 11));
      m_psGraphicsState.setFont(getFont());
    }
    
    m_printstream.println("/(" + replacePSFont(getFont().getPSName()) + ")" + " findfont");
    m_printstream.println(yScale(getFont().getSize()) + " scalefont setfont");        
  }
  
  /**
   * Not implemented
   */
  public void setPaintMode(){}
  
  /**
   * Not implemented
   */
  public void setXORMode(Color c1){}
  
  /**
   * Translates the origin of the graphics context to the point (x, y) in the 
   * current coordinate system. Modifies this graphics context so that its new 
   * origin corresponds to the point (x, y) in this graphics context's original 
   * coordinate system. All coordinates used in subsequent rendering operations 
   * on this graphics context will be relative to this new origin.
   * 
   * @param x the x coordinate.
   * @param y the y coordinate.
   */
  public void translate(int x, int y){
    if (DEBUG)
      System.out.println("translate with x = " + x + " and y = " + y);
    m_localGraphicsState.setXOffset(m_localGraphicsState.getXOffset() + xScale(x));
    m_localGraphicsState.setYOffset(m_localGraphicsState.getYOffset() + yScale(y));
    m_psGraphicsState.setXOffset(m_psGraphicsState.getXOffset() + xScale(x));
    m_psGraphicsState.setYOffset(m_psGraphicsState.getYOffset() + yScale(y));
  }
  /***** END overridden Graphics methods *****/
  
  /***** START overridden Graphics2D methods *****/
  
  public FontRenderContext getFontRenderContext(){
    return (new FontRenderContext(null,true,true));
  }
  public void clip(Shape s){}
  public Stroke getStroke(){
    return(m_localGraphicsState.getStroke());
  }
  
  public Color getBackground(){
    return(Color.white);
  }
  public void setBackground(Color c){}
  public Composite getComposite(){
    return(AlphaComposite.getInstance(AlphaComposite.SRC));
  }
  public Paint getPaint(){
    return((Paint) (new Color(getColor().getRed(),getColor().getGreen(),getColor().getBlue())));
  }
  public AffineTransform getTransform(){
    return(new AffineTransform());
  }
  public void setTransform(AffineTransform at) {}
  public void transform(AffineTransform at) {}
  public void shear(double d1, double d2){}
  public void scale(double d1, double d2) {
    m_localGraphicsState.setXScale(d1);
    m_localGraphicsState.setYScale(d2);
    if (DEBUG)
      System.err.println("d1 = " + d1 + ", d2 = " + d2);
  }
  public void rotate(double d1, double d2, double d3){}
  public void rotate(double d1){}
  public void translate(double d1, double d2) {}
  public RenderingHints getRenderingHints(){
    return(new RenderingHints(null));
  }
  public void addRenderingHints(Map m){}
  public void setRenderingHints(Map m){}
  public Object getRenderingHint(RenderingHints.Key key){
    return(null);
  }
  public void setRenderingHint(RenderingHints.Key key, Object o){}
  public void setStroke(Stroke s){
    if (s != null){
      m_localGraphicsState.setStroke(s);
      if (s.equals(m_psGraphicsState.getStroke())) {
        return;
      }
      m_psGraphicsState.setStroke(s); 
    } else {
      m_localGraphicsState.setStroke(new BasicStroke());
      m_psGraphicsState.setStroke(getStroke());
    }
    // ouput postscript here to set stroke.
  }
  public void setPaint(Paint p){
  }
  public void setComposite(Composite c){}
  public GraphicsConfiguration getDeviceConfiguration(){
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    GraphicsDevice gd = ge.getDefaultScreenDevice();
    return(gd.getDefaultConfiguration()); 
  }
  public boolean hit(Rectangle r, Shape s, boolean onstroke){
    return(false);
  }
  public void fill(Shape s){}
  public void drawGlyphVector(GlyphVector gv, float f1, float f2){} 
  public void drawString(AttributedCharacterIterator aci, float f1, float f2){}
  public void drawString(String str, float x, float y){
    drawString(str,(int)x, (int)y);
  }
  public void drawRenderableImage(RenderableImage ri, AffineTransform at){}
  public void drawRenderedImage(RenderedImage ri, AffineTransform af){}
  public void drawImage(BufferedImage bi, BufferedImageOp bio, int i1, int i2){}
  public boolean drawImage(Image im, AffineTransform at, ImageObserver io){
    return(false);
  }
  public void draw(Shape s){}
  /***** END *****/
}
