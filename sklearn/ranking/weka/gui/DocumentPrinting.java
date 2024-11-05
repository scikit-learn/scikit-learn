/*
 * DocumentPrinting.java (original classname: PrintMe)
 * (C) Jan Michael Soan (http://it.toolbox.com/wiki/index.php/How_to_print_in_Java)
 * (C) 2009 University of Waikato, Hamilton NewZealand
 */

package weka.gui;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import javax.swing.JTextPane;
import javax.swing.text.Document;
import javax.swing.text.View;

/**
 * DocumentPrinting is a class that lets you print documents on
 * the fly for free ;) Printing in JDK 1.2 - 1.5 is hard.
 * With this, you just simply call it in your application
 * and add your text component like JTextPane:
 * <pre>
 * new DocumentPrinting().print(YourJTextComponent);
 * </pre>
 * Reminder: Just make sure there is a text on your component ;P
 * <pre>
 * Author : Jan Michael Soan
 * WebSite: http://www.jmsoan.com
 * Date   : 04/17/2004 
 * Time   : 2:20 PM 
 * </pre>
 * 
 * Found on <a href="http://it.toolbox.com/wiki/index.php/How_to_print_in_Java" target="_blank">Toolbox</a>
 * (<a href="http://www.toolbox.com/TermsofUse.aspx" target="_blank">Terms of Use</a>).
 * 
 * @author Jan Michael Soan (<a href="http://it.toolbox.com/wiki/index.php/How_to_print_in_Java" target="_blank">original code</a>)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 */
public class DocumentPrinting
  implements Printable {

  /** the current page. */
  protected int m_CurrentPage = -1;
  
  /** the JTextPane which content is to be printed. */
  protected JTextPane m_PrintPane; 
  
  /** the page end. */
  protected double m_PageEndY = 0;
  
  /** the page start. */
  protected double m_PageStartY = 0;
  
  /** whether to scale the width to fit. */
  protected boolean m_ScaleWidthToFit = true; 
  
  /** the pageformat. */
  protected PageFormat m_PageFormat;
  
  /** the printer job. */
  protected PrinterJob m_PrinterJob;

  /**
   * Initializes the printing.
   */
  public DocumentPrinting() {
    m_PrintPane  = new JTextPane();
    m_PageFormat = new PageFormat();
    m_PrinterJob = PrinterJob.getPrinterJob();
  }

  /**
   * Brings up the page dialog.
   */
  protected void pageDialog() {
    m_PageFormat = m_PrinterJob.pageDialog(m_PageFormat);
  }

  /**
   * Prints the page.
   * 
   * @param graphics	the graphics context
   * @param pageFormat	the format of the page
   * @param pageIndex	the page index
   * @return		either NO_SUCH_PAGE or PAGE_EXISTS
   * @see		Printable#NO_SUCH_PAGE
   * @see		Printable#PAGE_EXISTS
   */
  public int print(Graphics graphics, PageFormat pageFormat, int pageIndex) {
    double scale = 1.0;
    Graphics2D graphics2D;
    View rootView;

    graphics2D = (Graphics2D) graphics;
    m_PrintPane.setSize((int) pageFormat.getImageableWidth(),Integer.MAX_VALUE);
    m_PrintPane.validate();

    rootView = m_PrintPane.getUI().getRootView(m_PrintPane);

    if ((m_ScaleWidthToFit) && (m_PrintPane.getMinimumSize().getWidth() > pageFormat.getImageableWidth())) {
      scale = pageFormat.getImageableWidth()/
      m_PrintPane.getMinimumSize().getWidth();
      graphics2D.scale(scale,scale);
    }

    graphics2D.setClip(
	(int) (pageFormat.getImageableX()/scale),
	(int) (pageFormat.getImageableY()/scale),
	(int) (pageFormat.getImageableWidth()/scale),
	(int) (pageFormat.getImageableHeight()/scale));

    if (pageIndex > m_CurrentPage) {
      m_CurrentPage = pageIndex;
      m_PageStartY += m_PageEndY;
      m_PageEndY = graphics2D.getClipBounds().getHeight();
    }

    graphics2D.translate(
	graphics2D.getClipBounds().getX(),
	graphics2D.getClipBounds().getY());
    Rectangle allocation = new Rectangle(
	0,
	(int) -m_PageStartY,
	(int) (m_PrintPane.getMinimumSize().getWidth()),
	(int) (m_PrintPane.getPreferredSize().getHeight()));

    if (printView(graphics2D,allocation,rootView)) {
      return Printable.PAGE_EXISTS;
    }
    else {
      m_PageStartY = 0;
      m_PageEndY = 0;
      m_CurrentPage = -1;
      return Printable.NO_SUCH_PAGE;
    }
  }

  /**
   * Prints the document in the JTextPane.
   * 
   * @param pane	the document to print
   */
  public void print(JTextPane pane) {
    setDocument(pane);
    printDialog();
  }

  /**
   * Shows the print dialog.
   */
  public void printDialog() {
    if (m_PrinterJob.printDialog()) {
      m_PrinterJob.setPrintable(this,m_PageFormat);
      try {
	m_PrinterJob.print();
      }
      catch (PrinterException printerException) {
	m_PageStartY = 0;
	m_PageEndY = 0;
	m_CurrentPage = -1;
	System.out.println("Error Printing Document");
      }
    }
  }

  /**
   * Shows a print view.
   * 
   * @param graphics2D	the graphics context
   * @param allocation	the allocation
   * @param view	the view
   * @return		true if the page exists
   */
  protected boolean printView(Graphics2D graphics2D, Shape allocation, View view) {
    boolean pageExists = false;
    Rectangle clipRectangle = graphics2D.getClipBounds();
    Shape childAllocation;
    View childView;

    if (view.getViewCount() > 0) {
      for (int i = 0; i < view.getViewCount(); i++) {
	childAllocation = view.getChildAllocation(i,allocation);
	if (childAllocation != null) {
	  childView = view.getView(i);
	  if (printView(graphics2D,childAllocation,childView)) {
	    pageExists = true;
	  }
	}
      }
    } else {
      if (allocation.getBounds().getMaxY() >= clipRectangle.getY()) {
	pageExists = true;
	if ((allocation.getBounds().getHeight() > clipRectangle.getHeight()) &&
	    (allocation.intersects(clipRectangle))) {
	  view.paint(graphics2D,allocation);
	} else {
	  if (allocation.getBounds().getY() >= clipRectangle.getY()) {
	    if (allocation.getBounds().getMaxY() <= clipRectangle.getMaxY()) {
	      view.paint(graphics2D,allocation);
	    } else {
	      if (allocation.getBounds().getY() < m_PageEndY) {
		m_PageEndY = allocation.getBounds().getY();
	      }
	    }
	  }
	}
      }
    }
    return pageExists;
  }

  /**
   * Sets the content type.
   * 
   * @param type	the content type
   */
  public void setContentType(String type) {
    m_PrintPane.setContentType(type);
  }

  /**
   * Returns the document to print.
   * 
   * @return		the document or null
   */
  public Document getDocument() {
    if (m_PrintPane != null) 
      return m_PrintPane.getDocument();
    else 
      return null;
  }

  /**
   * Sets the document from the given JTextPane.
   * 
   * @param pane	the JTextPane to get the document from
   */
  public void setDocument(JTextPane pane) { 
    m_PrintPane = new JTextPane();
    setDocument(pane.getContentType(), pane.getDocument());
  }

  /**
   * Sets the document and the according content type.
   * 
   * @param type	the content type
   * @param document	the document to print
   */
  public void setDocument(String type, Document document) {
    setContentType(type);
    m_PrintPane.setDocument(document);
  }

  /**
   * Sets whether to scale the width to fit.
   * 
   * @param scaleWidth	if true then the width will be scaled
   */
  public void setScaleWidthToFit(boolean scaleWidth) {
    m_ScaleWidthToFit = scaleWidth;
  }

  /**
   * Returns whether the width is to be scaled.
   * 
   * @return		true if scaled
   */
  public boolean getScaleWidthToFit() {
    return m_ScaleWidthToFit;
  }
} 