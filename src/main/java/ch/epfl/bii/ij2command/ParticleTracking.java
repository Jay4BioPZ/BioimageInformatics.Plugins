package ch.epfl.bii.ij2command;

import java.util.ArrayList;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.ImageCalculator;
import ij.process.ImageProcessor;

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>Homework Tracking Template")
public class ParticleTracking implements Command {
	// the template for the homework
	public void run() {
		double sigma = 5;
		double threshold = 10;
		double distance_max = 5;
		ImagePlus imp = IJ.getImage();
		int nt = imp.getNFrames();
		ImagePlus dog = dog(imp, sigma);
		// create dynamic array
		// use localMax to find the spot and store them in an array of Spots
		ArrayList<Spot> localmax[] = localMax(dog);
		// select the spots that are above the threshold after localMax
		ArrayList<Spot> spots[] = filter(dog, localmax, threshold);
		// connect the spots in time
		// PART THAT YOU HAVE TO IMPLEMENT AND MODIFY
		for (int t = 0; t < nt - 1; t++) {
			for (Spot current : spots[t]) {
				for (Spot next : spots[t+1]) {
					// if the distance between the two spots is less than 3 pixels, link them
					// a naive method
					// i.e., if dist. > 3 pixels, there is a breakage in the track
					if (current.distance(next) < distance_max)
						current.link(next);
				}
			}
		}
		Overlay overlay = new Overlay();
		draw(overlay, spots);
		imp.setOverlay(overlay);
	}


	private void draw(Overlay overlay, ArrayList<Spot> spots[]) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++)
			for (Spot spot : spots[t])
				spot.draw(overlay);
	}

	// the dog filter with sigma and sigma*sqrt(2)
	private ImagePlus dog(ImagePlus imp, double sigma) {
		ImagePlus g1 = imp.duplicate();
		ImagePlus g2 = imp.duplicate();
		IJ.run(g1, "Gaussian Blur...", "sigma=" + sigma + " stack");
		IJ.run(g2, "Gaussian Blur...", "sigma=" + (Math.sqrt(2) * sigma) + " stack");
		ImagePlus dog = ImageCalculator.run(g1, g2, "Subtract create stack");
		dog.show();
		return dog;
	}

	private ArrayList<Spot>[] filter(ImagePlus dog, ArrayList<Spot> spots[], double threshold) {
		int nt = spots.length;
		ArrayList<Spot> out[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			out[t] = new Spots();
			for (Spot spot : spots[t]) {
				dog.setPosition(1, 1, t + 1);
				double value = dog.getProcessor().getPixelValue(spot.x, spot.y);
				if (value > threshold)
					out[t].add(spot);
			}
		}
		return out;
	}

	public Spots[] localMax(ImagePlus imp) {
		int nt = imp.getNFrames();
		int nx = imp.getWidth();
		int ny = imp.getHeight();
		Spots spots[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			imp.setPosition(1, 1, t + 1);
			ImageProcessor ip = imp.getProcessor();
			spots[t] = new Spots();
			// loop over the pixels
			for (int x = 1; x < nx - 1; x++) {
				for (int y = 1; y < ny - 1; y++) {
					double v = ip.getPixelValue(x, y);
					double max = -1;
					// loop over the 3x3 neighborhood
					for (int k = -1; k <= 1; k++)
						for (int l = -1; l <= 1; l++)
							max = Math.max(max, ip.getPixelValue(x + k, y + l));
					// if the pixel is the maximum of the neighborhood, add it to the list
							if (v == max)
						spots[t].add(new Spot(x, y, t));
				}
			}
		}
		// return a list of spots
		return spots;
	}

}
