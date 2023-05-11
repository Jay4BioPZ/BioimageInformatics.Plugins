package ch.epfl.bii.ij2command.TirfTracking;

/*
 *  BIO-410 Bioimage Informatics Miniproject - TIRF protein tracking
 *  Author: Yung-Cheng Chiang
 *  Latest update: 2023/04/29
 *  
 *  Abstract:
 */

import java.util.ArrayList;

import org.apache.commons.io.filefilter.TrueFileFilter;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.ImageCalculator;
import ij.process.ImageProcessor;

// Workflow
// 1. Load the image
// 2. Apply the DoG filter
// 3. Find the local maxima
// 4. Filter the local maxima
// 5. Link the local maxima
//  5.1. Loop over the frames
//  5.2. Loop over the spots in current frame and next frame
//  5.3. Use a GOOD COST FUNCTION to determine if the two spots are linked
// 6. Display the result

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>TIRFTracking")
public class TirfTracking implements Command {
	
	// pass parameter from the ParticleTrackingTest()
	@Parameter
	private double lambda; // lambda ~ 0 is the best !!!
	
	public void run() {
		// IMAGE PROCESSING VARIABLES
		double sigma = 1; // dog filter
		int bd = 5; // localMax kernel size would be 2*bd+1
		// base on the line profile, a threshold 1000 seems ok
		double threshold = 2000; // localMax thresholding (input is a 16 bit image
		int radius = 2; // spot circle contour
		
		// general parameters
		ImagePlus imp = IJ.getImage();
		int nt = imp.getNFrames();
		System.out.println("Number of frames: " + nt);
		int xmax = imp.getWidth();
		int ymax = imp.getHeight();
		
		// finding potential offspring
		// set zero (assume tirf protein won't divide, so there is no offspring)
		double distance_max = 0;
		
		// cost function
		double fmax = findMax(imp, nt);
		double dmax = Math.sqrt(xmax*xmax + ymax*ymax);
		double intdiff = 0;
		double temp_cost = 1;
		
		// start the workflow
		ImagePlus dog = dog(imp, sigma);
		// create dynamic array
		// use localMax to find the spot and store them in an array of Spots
		ArrayList<Spot> localmax[] = localMax(dog, bd);
		// select the spots that are above the threshold after localMax
		ArrayList<Spot> spots[] = filter(dog, localmax, threshold);
		// connect the spots in positive time order
		// PART THAT YOU HAVE TO IMPLEMENT AND MODIFY
		for (int t = 0; t < nt - 1; t++) {
			System.out.println("Frame " + (t+1));
			for (Spot current : spots[t]) { // for each spot in the current frame
				// for method 2
				double min_cost = 1;
				int min_ind = 0;
				int ind = 0;
				for (Spot next : spots[t+1]) { // for each spot in the next frame
					// COST FUNCTION
					// DISTANCE and INTENSITY
					// calculate the cost for every spot next in frame t+1
					// link the spot with the lowest cost
					intdiff = findIntdiff(imp, t, current, next);
					// compute the cost
					temp_cost = getDISINT(current, next, dmax, fmax, lambda, intdiff);
					if (temp_cost <= min_cost) {
						min_cost = temp_cost; 
						min_ind = ind;
					}
					ind++;
					// identify possible offspring
					if (current.distance(next) < distance_max) {
						current.offspring(next);
					}
				}
				// link the spot with the minimum cost
				Spot good_next = spots[t+1].get(min_ind);
				current.link(good_next);
			}
		}
		// bring back to the first frame
		Overlay overlay = new Overlay();
		draw(overlay, spots, radius);
		imp.setOverlay(overlay);
		//System.out.println("Workflow finished");
	}


	private void draw(Overlay overlay, ArrayList<Spot> spots[], int radius) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++) {
			// print the time point
			System.out.println("overlaying time point: " + (t+1));
			for (Spot spot : spots[t]) {
				// call the draw function in the Spot class
				spot.draw(overlay, radius);
			}
		}
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
					// set a threshold for spotting
					out[t].add(spot);
			}
		}
		return out;
	}

	public Spots[] localMax(ImagePlus imp, int bd) {
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
					// loop over the bd neighborhood
					for (int k = -bd; k <= bd; k++)
						for (int l = -bd; l <= bd; l++)
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

	// compute the max of the image
	public double findMax(ImagePlus imp, int nt) {
		double fmax = 0;
		for (int t = 0; t < nt; t++) {
			imp.setSlice(t+1);
			ImageProcessor ip = imp.getProcessor();
			double max = ip.getMax();
			if (max > fmax)
				fmax = max;
		}
		return fmax;
	}
	
	// get the coordinates of the image
	public double findIntdiff(ImagePlus imp, int t, Spot current, Spot next) {
		imp.setSlice(t+1);
		ImageProcessor pimp1 = imp.getProcessor();
		double v1 = pimp1.getPixelValue(current.x, current.y);
		imp.setSlice(t+2);
		ImageProcessor pimp2 = imp.getProcessor();
		double v2 = pimp2.getPixelValue(next.x, next.y);
		double intdiff = Math.abs(v1-v2);
		return intdiff;
	}

	// compute the DISINT cost function
	public double getDISINT(Spot current, Spot next, double dmax, double fmax, double lambda, double intdiff) {
		double cost = 0;
		double dist = current.distance(next);
		// do something
		cost = (1-lambda)*(dist/dmax)+lambda*(intdiff)/fmax;
		return cost;
	}

}
