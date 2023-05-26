package ch.epfl.bii.ij2command.TirfTracking;
import java.awt.Color;

import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import java.io.Serializable;

public class Spot implements Serializable {
	private static final long serialVersionUID = 1L;
	public int x;
	public int y;
	public int t;
	// to build a graph if needed
	// "next" become a "spot object" in this class
	Spot next = null;
	private Color color;
	private boolean realTraj;

	// assign a random color to the spot
	public Spot(int x, int y, int t) {
		// "this" refer to the current object. 
		// For user to have access to object attributes 
		this.x = x;
		this.y = y;
		this.t = t;
		this.realTraj = true;
		// generate a random color
		color = Color.getHSBColor((float)Math.random(), 1f, 1f);
		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 255); // no transparency
	}

	// compute the distance between two spots
	public double distance(Spot spot) {
		double dx = x - spot.x;
		double dy = y - spot.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	// draw the spot and the link to the next spot
	public void draw(Overlay overlay, int radius) {
		double xp = x + 0.5;
		double yp = y + 0.5;
		// OvalRoi is a class in ImageJ, which is a circle
		OvalRoi roi = new OvalRoi(xp - radius, yp - radius, 2 * radius, 2 * radius);
		// display roi in the t+1 frame (t is one-based)
		roi.setPosition(t+1); // display roi in one frame
		// display the circle outline
		roi.setStrokeColor(new Color(255, 255, 0, 150));
		roi.setStrokeWidth(0.3);
		overlay.add(roi);
		// draw trajectory
		if (next != null) {
			Line line = new Line(x, y, next.x, next.y);
			line.setStrokeColor(color);
			line.setStrokeWidth(0.5);
			overlay.add(line);
		}
		
		//TextRoi text = new TextRoi(x, y-10, "" + value);
		//text.setPosition(t+1);
		//overlay.add(text);
	}
	
	// usage: current.link(next);
	public void link(Spot a) {
		// if a is null, do nothing and return
		if (a == null)
			return;
		a.next = this; // take the spot "current" as the next spot of "a" a->current
		a.color = this.color;
	}

	public String toString() {
		return "(" + x + ", " + y + ", " + t + ")";
	}
}