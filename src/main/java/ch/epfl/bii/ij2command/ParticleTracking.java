package ch.epfl.bii.ij2command;

import java.util.ArrayList;

import org.apache.commons.io.filefilter.TrueFileFilter;
import org.scijava.command.Command;
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
//  5.4. Using a GOOD COST FUNCTION to determine if the two spots are linked
// 6. Display the result

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>Homework Tracking Template")
public class ParticleTracking implements Command {
	// the template for the homework
	public void run() {
		// for model selection
		boolean method1 = false;
		boolean method2 = true;
		double lambda = 0.3;
		
		// for dog and localmax
		double sigma = 5;
		double threshold = 10;
		
		// for general parameters
		ImagePlus imp = IJ.getImage();
		int nt = imp.getNFrames();
		int xmax = imp.getWidth();
		int ymax = imp.getHeight();
		
		// for method 1
		double distance_max = 5;
		
		// for method 2
		double fmax = findMax(imp, nt);
		double dmax = Math.sqrt(xmax*xmax + ymax*ymax);
		double intdiff = 0;
		double temp_cost = 1;
		
		// start the workflow
		ImagePlus dog = dog(imp, sigma);
		// create dynamic array
		// use localMax to find the spot and store them in an array of Spots
		ArrayList<Spot> localmax[] = localMax(dog);
		// select the spots that are above the threshold after localMax
		ArrayList<Spot> spots[] = filter(dog, localmax, threshold);
		// connect the spots in positive time order
		// PART THAT YOU HAVE TO IMPLEMENT AND MODIFY
		for (int t = 0; t < nt - 1; t++) {
			System.out.println("Frame " + (t+1));
			int cur_ind = 0;
			for (Spot current : spots[t]) { // for each spot in the current frame
				// for method 2
				double min_cost = 1;
				int min_ind = 0;
				int ind = 0;
				for (Spot next : spots[t+1]) { // for each spot in the next frame
					// COST FUNCTION
					// 1. THRESHOLDING
					// if the distance between the two spots is less than distance_max, link them
					// i.e., if dist. > 3 pixels, there is a breakage in the track
					if (method1 == true) {
						if (current.distance(next) < distance_max) {
							current.link(next);
							min_ind = ind;
							// System.out.println("Frame " + (t+1) + " add link between cur " + cur_ind + "(" + current.x + ", " + current.y + ") and next " + min_ind + "(" + next.x + ", " + next.y + ")");
							break;
						}
						ind++;
					}
					// 2. DISTANCE and INTENSITY
					// calculate the cost for every spot next in frame t+1
					// link the spot with the lowest cost
					if (method2 == true) {
						intdiff = findIntdiff(imp, t, current, next);
						// compute the cost
						temp_cost = getDISINT(current, next, dmax, fmax, lambda, intdiff);
						if (temp_cost <= min_cost) {
							min_cost = temp_cost; 
							min_ind = ind;
						}
						ind++;
					}
				}
				// link the spot with the minimum cost
				//System.out.println("time point: " + (t+1) + " find min cost " + min_cost + " for cur " + cur_ind + " and min_ind " + min_ind);
				// for method 2
				if (method2 == true) {
					Spot good_next = spots[t+1].get(min_ind);
					//System.out.println("Frame " + (t+1) + " add link between cur " + cur_ind + "(" + current.x + ", " + current.y + ") and next " + min_ind + "(" + good_next.x + ", " + good_next.y + ")");
					current.link(good_next);
				}
				cur_ind++;
			}
		}
		// bring back to the first frame
		Overlay overlay = new Overlay();
		draw(overlay, spots);
		imp.setOverlay(overlay);
		//System.out.println("Workflow finished");
	}


	private void draw(Overlay overlay, ArrayList<Spot> spots[]) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++) {
			// print the time point
			//System.out.println("overlaying time point: " + (t+1));
			for (Spot spot : spots[t]) {
				// call the draw function in the Spot class
				spot.draw(overlay);
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
