package ch.epfl.bii.ij2command;
import java.awt.Color;

import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;

public class Spot {
	public int x;
	public int y;
	public int t;
	// to build a graph if needed
	private Spot next = null;
	private Color color;

	// assign a random color to the spot
	public Spot(int x, int y, int t) {
		this.x = x;
		this.y = y;
		this.t = t;
		color = Color.getHSBColor((float)Math.random(), 1f, 1f);
		color = new Color(color.getRed(), color.getGreen(), color.getBlue(), 120);
	}

	// compute the distance between two spots
	public double distance(Spot spot) {
		double dx = x - spot.x;
		double dy = y - spot.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	// draw the spot and the link to the next spot
	public void draw(Overlay overlay) {
		double xp = x + 0.5;
		double yp = y + 0.5;
		int radius = 5;
		OvalRoi roi = new OvalRoi(xp - radius, yp - radius, 2 * radius, 2 * radius);
		// display roi in the t+1 frame (t is one-based)
		roi.setPosition(t+1); // display roi in one frqme
		roi.setStrokeColor(new Color(255, 0, 0, 120));
		roi.setStrokeWidth(1);
		overlay.add(roi);
		if (next != null) {
			Line line = new Line(x, y, next.x, next.y);
			line.setStrokeColor(color);
			line.setStrokeWidth(2);
			overlay.add(line);
		}
		
		//TextRoi text = new TextRoi(x, y-10, "" + value);
		//text.setPosition(t+1);
		//overlay.add(text);
	}
	
	public void link(Spot a) {
		if (a == null)
			return;
		a.next = this;
		a.color = this.color;
	}
	
	public String toString() {
		return "(" + x + ", " + y + ", " + t + ")";
	}
}