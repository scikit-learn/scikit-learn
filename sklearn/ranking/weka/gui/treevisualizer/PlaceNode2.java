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
 *    PlaceNode2.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

import java.util.*;

/**
 * This class will place the Nodes of a tree. <p>
 * 
 * It will place these nodes so that they fall at evenly below their parent.
 * It will then go through and look for places where nodes fall on the wrong 
 * side of other nodes
 * when it finds one it will trace back up the tree to find the first common 
 * sibling group these two nodes have
 * And it will adjust the spacing between these two siblings so that the two 
 * nodes no longer overlap.
 * This is nasty to calculate with , and takes a while with the current 
 * algorithm I am using to do this.<p>
 *
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class PlaceNode2 implements NodePlace {
  /** The space each row will take up. */
  private double m_yRatio;

  /** An array that lists the groups and information about them. */
  private Group[] m_groups;

  /** An array that lists the levels and information about them. */
  private Level[] m_levels;

  /** The Number of groups the tree has */
  private int m_groupNum;

  /** The number of levels the group tree has */
  private int m_levelNum;
  
  /** 
   * The Funtion to call to have the nodes arranged.
   *
   * @param r The top node of the tree to arrange.
   */
  public void place(Node r) {
    //note i might write count groups into the node class as well as
    //it may be useful too;
    
    m_groupNum = Node.getGCount(r,0); //i could swap over to the node class 
    //group count,but this works os i'm not gonna
    m_groups = new Group[m_groupNum];
        
    for (int noa = 0;noa < m_groupNum;noa++) {
      m_groups[noa] = new Group();
      m_groups[noa].m_gap = 3;
      m_groups[noa].m_start = -1;
    }
    
    groupBuild(r);
    m_levelNum = Node.getHeight(r,0);
    m_yRatio = 1 / (double)(m_levelNum + 1);
    
    m_levels = new Level[m_levelNum];
    
    for (int noa = 0;noa < m_levelNum;noa++) {
      m_levels[noa] = new Level();
    }
    r.setTop(m_yRatio);
    yPlacer();
    r.setCenter(0);
    xPlacer(0);


    //ok now i just have to untangle then scale down
    //note instead of starting with coords between 1 and 0 i will
    //use ints then scale them down 
    //i will scale them down either by all relative to the largest
    //line or by each line individually

    untangle2();
    
    scaleByMax();
    //scaleByInd();
  }


  /*
  private void thinner()
  {
    //what this function does is it retains the symmetry of the
     // parent node about the children but the children are no longer evenly
      //spaced this stops children from being pushed too far to the sides
      //,note this algorithm may need the method altered as it may 
     // require heavy optimisation to go at any decent speed   
  
    Node r,s;
    Edge e;
    double parent_x;
    for (int noa = group_num - 1;noa >= 0;noa--)
      {
	Vector shifts = new Vector(20,10);
	shifts.addElement(0);
	int g_num = 0;//this is the offset from groups.m_start to get the right 1
	r = groups[noa].m_p;
	parent_x = r.getCenter();
	for (int nob = 1;(e = r.getChild(nob)) != null;nob++)
	  {
	    double margin;
	    s = e.getTarget();
	    margin = s_getCenter - r.getChild(nob - 1).getTarget().getCenter-1
	             - shift.elementAt(nob-1);
	    if (margin > 0)
	      {
		margin = check_down(s,g_num,margin);
		if (margin > 0)
		  {
		    shift.addElement(-margin);
		  }
		else
		  {
		    shift.addElement(0);
		  }
	      }
	    else
	      {
		shift.addElement(0);
	      }
	    if (s.getChild(0) != null)
	      {
		g_num++;
	      }
	  }
      }
  }


  private double check_down(Node r,int gn,double m)
  {
    //note i need to know where the children of the 
    //other changers are to properly overlap check
    //to do this i think the best way is to go up the other group
    //parents line and see if it goes through the current group
    //this means to save time i need to know the level that is being 
    //worked with along with the group
    
    Edge e;
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++)
      {
	
      }
  }
*/


  /**
   * This will set initial places for the x coord of the nodes.
   * @param start The `number for the first group to start on (I think).
   */
  private void xPlacer(int start) {
    //this can be one of a few x_placers (the first)
    //it will work by placing 1 space inbetween each node
    //ie the first at 0 the second at 1 and so on
    //then it will add to this value the place of the parent 
    //node - half of the size
    //i will break this up into several functions
    //first the gap setter;
    //then the shifter
    //it will require a vector shift function added to the node class
    //i will write an additional shifter for the untangler 
    //for its particular situation

    Node r;
    Edge e;
    if (m_groupNum > 0) {
      m_groups[0].m_p.setCenter(0);
      for (int noa = start;noa < m_groupNum;noa++) {
	int nob,alter =0;
	double c = m_groups[noa].m_gap;
	r = m_groups[noa].m_p;
	for (nob = 0;(e = r.getChild(nob)) != null;nob++) {
	  if (e.getTarget().getParent(0) == e) {
	    e.getTarget().setCenter(nob * c);
	  }
	  else {
	    alter++;
	  }
	}
	m_groups[noa].m_size = (nob - 1 - alter) * c;
	xShift(noa);
      }
    }
  }

  /**
   * This will shift a group of nodes to be aligned under their parent.
   * @param n The group number to shift
   */
  private void xShift(int n) {
    Edge e;
    Node r = m_groups[n].m_p;
    double h = m_groups[n].m_size / 2;
    double c = m_groups[n].m_p.getCenter();
    double m = c - h;
    m_groups[n].m_left = m;
    m_groups[n].m_right = c + h;
    
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
      if (e.getTarget().getParent(0) == e) {
	e.getTarget().adjustCenter(m);
      }
    }
  }

  /**
   * This scales all the x values to be between 0 and 1.
   */
  private void scaleByMax() {
    //ammendment to what i may have commented before
    //this takes the lowest x and highest x  and uses that as the scaling
    //factor
    double l_x = 5000,h_x = -5000;
    for (int noa = 0;noa < m_groupNum;noa++) {
      if (l_x > m_groups[noa].m_left) {
	l_x = m_groups[noa].m_left;
      }

      if (h_x < m_groups[noa].m_right) {
	h_x = m_groups[noa].m_right;
      }
    }
    
    Edge e;
    Node r,s;
    double m_scale = h_x - l_x + 1;
    if (m_groupNum > 0) {
      r = m_groups[0].m_p;
      r.setCenter((r.getCenter() - l_x) / m_scale);
      //System.out.println("from scaler " + l_x + " " + m_scale);
      for (int noa = 0; noa < m_groupNum;noa++) {
	r = m_groups[noa].m_p;
	for (int nob = 0;(e = r.getChild(nob)) != null;nob++) {
	  s = e.getTarget();
	  if (s.getParent(0) == e) {
	    s.setCenter((s.getCenter() - l_x) / m_scale);
	  }
	}
      }
    }
  }
  
  /**
   * This scales the x values to between 0 and 1 for each individual line
   * rather than doing them all at once.
   */
  private void scaleByInd() {
    //ammendment to what i may have commented before
    //this takes the lowest x and highest x  on each line and uses that for 
    //the line in question
    double l_x,h_x;

    Edge e;
    Node r,s;
    r = m_groups[0].m_p;
    r.setCenter(.5);
    double m_scale;
    for (int noa = 0;noa < m_levelNum;noa++) {
      l_x = m_groups[m_levels[noa].m_start].m_left;
      h_x = m_groups[m_levels[noa].m_end].m_right;
      m_scale = h_x - l_x + 1;
      for (int nob = m_levels[noa].m_start; nob <= m_levels[noa].m_end;nob++) {
	r = m_groups[nob].m_p;
	for (int noc = 0;(e = r.getChild(noc)) != null;noc++) {
	  s = e.getTarget();
	  if (s.getParent(0) == e) {
	    s.setCenter((s.getCenter() - l_x) / m_scale);
	  }
	}
      }
    }
  }
  
  /**
   * This untangles the nodes so that they will will fall on the correct
   * side of the other nodes along their row.
   */
  private void untangle2() {
    Ease a;
    Edge e;
    Node r,nf = null,ns = null,mark;
    int l = 0,times = 0;
    int f,s,tf = 0,ts = 0,pf,ps;
    while ((a = overlap(l)) != null) {
      times++;
      //System.out.println("from untang 2 " + group_num);
      f = a.m_place;
      s = a.m_place + 1;
      while (f != s) {
	a.m_lev--;
	tf = f;
	ts = s;
	f = m_groups[f].m_pg;
	s = m_groups[s].m_pg;
      }
      l = a.m_lev;
      pf = 0;
      ps = 0;
      r = m_groups[f].m_p;
      mark = m_groups[tf].m_p;
      nf = null;
      ns = null;
      for (int noa = 0; nf != mark;noa++) {
	pf++;
	nf = r.getChild(noa).getTarget();
      }
      mark = m_groups[ts].m_p;
      for (int noa = pf; ns != mark;noa++) {
	ps++; //the number of gaps between the two nodes
	ns = r.getChild(noa).getTarget();
      }
      //m_groups[f].gap =
      //              Math.ceil((a.amount / (double)ps) + m_groups[f].gap);
      //note for this method i do not need the group gap ,but i will leave
      //it for the other methods;
      Vector o_pos = new Vector(20,10);
      for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
	if (e.getTarget().getParent(0) == e) {
	  Double tem = new Double(e.getTarget().getCenter());
	  o_pos.addElement(tem);
	}
      }

      pf--;
      double inc = a.m_amount / (double)ps;
      for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
	ns = e.getTarget();
	if (ns.getParent(0) == e) {
	  if (noa > pf + ps) {
	    ns.adjustCenter(a.m_amount);
	  }
	  else if (noa > pf) {
	    ns.adjustCenter(inc * (double)(noa - pf));
	  }
	}
      }

      nf = r.getChild(0).getTarget();
      inc = ns.getCenter() - nf.getCenter();
      m_groups[f].m_size = inc;
      m_groups[f].m_left = r.getCenter() - inc / 2; 
      m_groups[f].m_right = m_groups[f].m_left + inc;
      inc = m_groups[f].m_left - nf.getCenter();

      double shift;
      int g_num = 0;
      for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
	ns = e.getTarget();
	if (ns.getParent(0) == e) {
	  ns.adjustCenter(inc);
	  shift = ns.getCenter() - 
	    ((Double)o_pos.elementAt(noa)).doubleValue();
	  if (ns.getChild(0) != null) {
	    moveSubtree(m_groups[f].m_start + g_num,shift);
	    g_num++;
	  }
	}
	//ns.adjustCenter(-shift);
      }
      //zero_offset(r);
      
      //x_placer(f);
    }
  }


  /**
   * This will recursively shift a sub there to be centered about
   * a particular value.
   * @param n The first group in the sub tree.
   * @param o The point to start shifting the subtree.
   */
  private void moveSubtree(int n,double o) {
    Edge e;
    Node r = m_groups[n].m_p;
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
      if (e.getTarget().getParent(0) == e) {
	e.getTarget().adjustCenter(o);
      }
    }
    m_groups[n].m_left += o;
    m_groups[n].m_right += o;
    if (m_groups[n].m_start != -1) {
      for (int noa = m_groups[n].m_start;noa <= m_groups[n].m_end;noa++) {
	moveSubtree(noa,o);
      }
    }
  }


  /**
   * This will untangle the nodes in the tree so that they fall on the
   * correct side of each other.
   */
  private void untangle() {
    Ease a;
    Edge e;
    Node r,nf = null,ns = null,mark;
    int l = 0,times = 0;
    int f,s,tf = 0,ts = 0,pf,ps;
    while ((a = overlap(l)) != null) {
      times++;
      //System.out.println(group_num);
      f = a.m_place;
      s = a.m_place + 1;
      while (f != s) {
	a.m_lev--;
	tf = f;
	ts = s;
	f = m_groups[f].m_pg;
	s = m_groups[s].m_pg;
      }
      l = a.m_lev;
      pf = 0;
      ps = 0;
      r = m_groups[f].m_p;
      mark = m_groups[tf].m_p;
      nf = null;
      ns = null;
      for (int noa = 0; nf != mark;noa++) {
	pf++;
	nf = r.getChild(noa).getTarget();
      }
      mark = m_groups[ts].m_p;
      for (int noa = pf; ns != mark;noa++) {
	ps++; //the number of gaps between the two nodes
	ns = r.getChild(noa).getTarget();
      }
      m_groups[f].m_gap =
	Math.ceil((a.m_amount / (double)ps) + m_groups[f].m_gap);
      
      xPlacer(f);
    }
  }
  
  /**
   * This will find an overlap and then return information about that overlap
   * @param l The level to start on.
   * @return null if there was no overlap , otherwise an object containing
   * the group number that overlaps (only need one) how much they overlap by,
   * and the level they overlap on.
   */
  private Ease overlap(int l) {
    Ease a = new Ease();
    for (int noa = l;noa < m_levelNum;noa++) {
      for (int nob = m_levels[noa].m_start;nob < m_levels[noa].m_end;nob++) {
	a.m_amount = m_groups[nob].m_right - m_groups[nob+1].m_left + 2;
	//System.out.println(m_groups[nob].m_right + " + " + 
	//	       m_groups[nob+1].m_left + " = " + a.amount);
	if (a.m_amount >= 0) {
	  a.m_amount++;
	  a.m_lev = noa;
	  a.m_place = nob;
	  return a;
	}
      }
    }
    return null;
  }
  
  /* private int count_m_groups(Node r,int l)
  {
    Edge e;
    if (r.getChild(0) != null)
      {
	l++;
      }
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++)
      {
	l = count_groups(e.getTarget(),l);
      }

    return l;
  }
  */

  /**
   * This function sets up the height of each node, and also fills the
   * levels array with information about what the start and end groups on that
   * level are.
   */
  private void yPlacer() {
    //note this places the y height and sets up the levels array
    double changer = m_yRatio;
    int lev_place = 0;
    if (m_groupNum > 0) {
      m_groups[0].m_p.setTop(m_yRatio);
      m_levels[0].m_start = 0;
      
      for (int noa = 0;noa < m_groupNum;noa++) {
	if (m_groups[noa].m_p.getTop() != changer) {
	  m_levels[lev_place].m_end = noa - 1;
	  lev_place++;
	  m_levels[lev_place].m_start = noa;
	  changer = m_groups[noa].m_p.getTop();
	}
	nodeY(m_groups[noa].m_p);
      }
      m_levels[lev_place].m_end = m_groupNum - 1;
    }
  }

  /**
   * This will set all of the children node of a particular node to their
   * height.
   * @param r The parent node of the children to set their height. 
   */
  private void nodeY(Node r) {
    Edge e;
    double h = r.getTop() + m_yRatio;
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
      if (e.getTarget().getParent(0) == e) {
	e.getTarget().setTop(h);
	if (!e.getTarget().getVisible()) {
	  //System.out.println("oh bugger");
	}
      }
    }
  }
  
  /**
   * This starts to create the information about the sibling groups.
   * As more groups are created the for loop in this will check those groups
   * for lower groups.
   * @param r The top node.
   */
  private void groupBuild(Node r) {
    if (m_groupNum > 0) {
      m_groupNum = 0;
      m_groups[0].m_p = r;
      m_groupNum++;
      //note i need to count up the num of groups first
      //woe is me
      for (int noa = 0;noa < m_groupNum ;noa++) {
	groupFind(m_groups[noa].m_p,noa);
      }
    }
  }
  
  /**
   * This is called to build the rest of the grouping information.
   * @param r The parent of the group.
   * @param pg The number for the parents group.
   */
  private void groupFind(Node r,int pg) {
    Edge e;
    boolean first = true;
    for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
      if (e.getTarget().getParent(0) == e) {
	if (e.getTarget().getChild(0) != null && e.getTarget().getCVisible()) {
	  if (first) {
	    m_groups[pg].m_start = m_groupNum;
	    first = false;
	  }
	  m_groups[pg].m_end = m_groupNum;
	  m_groups[m_groupNum].m_p = e.getTarget();
	  m_groups[m_groupNum].m_pg = pg;
	  m_groups[m_groupNum].m_id = m_groupNum; //just in case I ever need
	  //this info
	  m_groupNum++;
	}
      }
    }
  }
  

  //note these three classes are only to help organise the data and are
  //inter related between each other and this placer class
  //so don't mess with them or try to use them somewhere else
  //(because that would be a mistake and I would pity you)
  

  /**
   * Inner class for containing the level data.
   */
  private class Level {
    /** The number for the group on the left of this level. */
    public int m_start;
    /** The number for the group on the right of this level. */
    public int m_end;

    /** These two params would appear to not be used. */
    public int m_left;
    public int m_right;
  }

  /**
   * Inner class for containing the grouping data.
   */
  private class Group {
    /** The parent node of this group. */
    public Node m_p;

    /** The group number for the parent of this group. */
    public int m_pg;

    /** The gap size for the distance between the nodes in this group. */
    public double m_gap;

    /** The leftmost position of this group. */
    public double m_left;

    /** The rightmost position of this group. */
    public double m_right;

    /** The size of this group. */
    public double m_size;

    /** The start node of this group. */
    public int m_start;

    /** The end node of this group. */
    public int m_end;

    /** The group number for this group. (may not be used!?). */
    public int m_id;
  }
  
  /**
   * An inner class used to report information about any tangles found. 
   */
  private class Ease {
    /** The number of the group on the left of the tangle. */
    public int m_place;
    /** The distance they were tangled. */
    public double m_amount;
    /** The level on which they were tangled. */
    public int m_lev;
  }
}














