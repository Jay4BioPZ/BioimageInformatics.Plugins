
package ch.epfl.bii.ij2command.TirfTracking;

/*
 *  BIO-410 Bioimage Informatics Miniproject - TIRF protein tracking
 *  Author: Yung-Cheng Chiang
 *  Latest update: 2023/04/29
 *  
 *  Abstract:
 */

import java.util.ArrayList;
// for object serialization
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.FileInputStream;
import java.io.ObjectInputStream;

import org.apache.commons.io.filefilter.TrueFileFilter;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.CurveFitter;
import ij.plugin.ImageCalculator;
import ij.process.ImageProcessor;
import ij.plugin.frame.RoiManager;
import ij.measure.Measurements;

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>TIRFTracking")
public class TirfTracking implements Command {
	// WORKFLOW
	/* 1. Detection
	 * 	1.1 Hyperstack reordering: convert from z-stack to time stack
	 * 	1.2 Exponential correction (skipped)
	 *  1.3 DoG filter with sigma1 = 1, sigma2 = 6
	 *  1.4 Median filter (skipped)
	 *  1.5 Gaussian blur time axis with sigma = 1 + TopHat filter with radius 2
	 *  1.6 LocalMax with a dynamic threshold: threstime*std of the frame + avoiding multi-maximal spots within the filter
	 *  
	 * 2. Linking (pseudo-code)
	 * 	for current particle in frame t
	 * 		for next particle in frame from t+1 to t+5 (allow frame skipping)
	 * 			compute the cost function 
	 * 			if the cost is the minimal, and distance between current and next doesn't exceed long_link
	 * 				link between current and next
	 * 				as long as the link forms, go ahead to the next current particle (break the inner loop
	 *  draw the linkage between frames
	 */
	
	public void run() {
		// generate a generic dialog
		GenericDialog gd = new GenericDialog("TIRF particle tracking");
		gd.addMessage("The plugin performs TIRF particle tracking and trajectories analysis.\nPlease enter the following parameters to proceed.");
		gd.addStringField("Image stack path (.tif): ", "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Project/data/");
		//gd.addNumericField("Threstime (>=5): ", 5);
		gd.addNumericField("Lambda (0-1): ", 0.3);
		gd.addCheckbox("Diffusing particle", false);
		gd.addMessage("Authors: BII2023 Group 8");
		gd.showDialog();
		if(gd.wasCanceled()) {
			return;
		}
		
		String imagePath = (String)gd.getNextString();
		//double threstime = (double)gd.getNextNumber();
		double lambda = (double)gd.getNextNumber();
		boolean diffuse = (boolean)gd.getNextBoolean();
		
		// IMAGE PROCESSING VARIABLES
		double sigma1, sigma2;
		if (diffuse == true) {
			sigma1 = 1;
			sigma2 = 6;
		}
		else {
			sigma1 = 1;
			sigma2 = Math.sqrt(2)*sigma1;
		}
		int bd = 7; // localMax kernel size would be 2*bd+1
		int radius = 1; // spot circle contour
		
		// DETECTION
		// general parameters
		ImagePlus imp = IJ.openImage(imagePath);
		imp.show();
		int nz = imp.getNSlices();
		int xmax = imp.getWidth();
		int ymax = imp.getHeight();
		
		// start the workflow
		// STEP D1. Hyperstack reordering
		imp.setDimensions(1, 1, nz);
		int nt = imp.getNFrames();
		IJ.log("The dataset have " + nt + " time frames");
		
		// STEP D2. Exponential correction
		//IJ.log("Launch exponential correction...");
		if (diffuse == false) {
			imp = imagecorrection(imp, nt, xmax, ymax);
		}
		
		// STEP D3. DoG
		IJ.log("Start particle detection...");
		ImagePlus imp_dog = dog(imp, sigma1, sigma2);
		
		// STEP D4. Median filter (skipped)
		
		// STEP D5. Gaussian blur 3D
		IJ.run(imp_dog, "Gaussian Blur 3D...", "x=0 y=0 z=1");
		IJ.run(imp_dog, "Top Hat...", "radius=2 stack");
		
		// STEP D6. LocalMax with a dynamic threshold
		//ArrayList<Spot> localmax[] = localMax(imp_dog, bd);
		//ArrayList<Spot> spots[] = filter(imp_dog, localmax, threstime);
		ArrayList<Spot> spots[] = thresDetect(imp_dog, nt, xmax, ymax);
		
		// LINKING
		IJ.log("Start particle linking...");
		// cost function
		double fmax = findMax(imp, nt);
		double dmax = Math.sqrt(xmax*xmax + ymax*ymax);
		double intdiff = 0;
		double temp_cost = 1;
		double long_link;
		
		for (int t = 0; t < nt - 1; t++) {
			int indc = 0;
			for (Spot current : spots[t]) { // for each spot in the current frame
				// for method 2
				double min_cost = 1;
				int min_ind = 0;
				int min_ts = 1;
				boolean findspot = false;
				// find candidates in next "few frames"
				for (int ts = 1; ts < Math.min(nt-t, 6); ts++) {
					// estimate the local density of spots and compute a good long_link
					long_link = estimateLocalDens(current, spots, t, ts, xmax, ymax, diffuse);
					int indn = 0;
					for (Spot next : spots[t+ts]) { // for each spot in the next frame
						intdiff = findIntdiff(imp, t, current, next);
						// compute the cost
						temp_cost = getDISINT(current, next, dmax, fmax, lambda, intdiff);
						if (temp_cost <= min_cost && current.distance(next) <= long_link) {
							min_cost = temp_cost;
							min_ind = indn;
							min_ts = ts;
							findspot = true;
						}
						indn++;
					}
					if (findspot == true) {
						break;
					}
				}
				// only when proper candidate is found then form link
				// otherwise just cut the trajectories
				if (findspot == true) {
					int nextt = t+min_ts;
					System.out.println("Link (" + t + ", " + indc + ") to (" + nextt + ", " + min_ind + ")");
					Spot good_next = spots[t+min_ts].get(min_ind);
					current.link(good_next);
				}
				indc++;
			}
		}
		//ArrayList<Spot> spots_filt = removeBadTraj(spots, nt); // remove single spot trajectory by changing the "realTraj" attribute in spot to false
		// overlaying the trajectories to the image stack
		IJ.log("Visualizing trajectories (it takes several minutes)...");
		Overlay overlay = new Overlay();
		draw(overlay, spots, radius);
		imp_dog.setOverlay(overlay);
		IJ.log("Finalizing the visualization");
	}

	
	// STEP D2. exponential correction
	public ImagePlus imagecorrection(ImagePlus imp, int nt, int nx, int ny) {
		//System.out.println("Step D1. Exponential Detection...");
		// plot the z-intensity profile before processing
		IJ.run(imp, "Select All", "");
		IJ.run(imp, "Plot Z-axis Profile", "");
		IJ.run(imp, "32-bit", "");
		ImagePlus out = imp.duplicate();
		out.show();
		
		// Get parameters
		double[] params = expgetparameter(imp, nt);
		// double A = params[0];
		double tau = params[1];
		double C = params[2];
		for(int t=0; t<nt; t++) {
			imp.setSlice(t+1);
			out.setSlice(t+1);
			ImageProcessor iin = imp.getProcessor();
			ImageProcessor iout = out.getProcessor();
			for (int x = 0; x < nx; x++) {
				for (int y = 0; y < ny; y++) {
					double new_v = iin.getPixelValue(x, y);
					new_v = (new_v-C)*Math.exp(t/tau)+C;
					iout.putPixelValue(x, y, new_v);
				}
			}	
		}
		out.setSlice(1);
		out.setTitle("out");
		// plot the z-intensity profile after processing
		IJ.run(out, "Select All", "");
		IJ.run(out, "Plot Z-axis Profile", "");
		return out;
	}
	
	// STEP D2. fitting with imageJ CurveFitter
	public double[] exponentialfit2(int[] x, double[] y) {
		double[] dx = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			dx[i] = (double) x[i]; 
		}
		CurveFitter fit = new CurveFitter(dx, y);
		fit.doFit(CurveFitter.EXP_WITH_OFFSET, true);
		double[] coefficients = fit.getParams();
		return coefficients;
	}
	
	// STEP D2. use three time point to evaluate the parameter of the exponential model
	public double[] expgetparameter(ImagePlus imp, int nt) {
		int[] ind = {0,nt/3,2*nt/3, nt-1};
		double[] mean_stack = new double[ind.length];
		for (int i = 0; i < ind.length; i++) {
			imp.setSlice(ind[i]+1); // t is one-based
			double mean = imp.getStatistics().mean;
			//System.out.println(mean);
			mean_stack[i] = mean;
		}
		double[] params = exponentialfit2(ind, mean_stack);
		double A = params[0];
		double tau = Math.abs(1/params[1]);
		double C = params[2];
		// System.out.println("A: "+A+"; B: "+tau+"; C: "+C);
		return new double[]{A, tau, C};
	}
	
	private void draw(Overlay overlay, ArrayList<Spot> spots[], int radius) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++) {
			// print the time point
			//System.out.println("overlaying time point: " + (t+1));
			for (Spot spot : spots[t]) {
				// call the draw function in the Spot class
				spot.draw(overlay, radius);
			}
		}
	}

	// the dog filter with sigma and sigma*sqrt(2)
	
	// STEP D3. DoG filter
	private ImagePlus dog(ImagePlus imp, double sigma1, double sigma2) {
		//System.out.println("Step D2. DoG filtering...");
		ImagePlus g1 = imp.duplicate();
		ImagePlus g2 = imp.duplicate();
		if (sigma1 != 0) {
			IJ.run(g1, "Gaussian Blur...", "sigma=" + sigma1 + " stack");
		}
		IJ.run(g2, "Gaussian Blur...", "sigma=" + (Math.sqrt(2) * sigma2) + " stack");
		ImagePlus dog = ImageCalculator.run(g1, g2, "Subtract create stack");
		dog.show();
		return dog;
	}
	

	private ArrayList<Spot>[] filter(ImagePlus imp, ArrayList<Spot> spots[], double threstime) {
		int nt = spots.length;
		ArrayList<Spot> out[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			out[t] = new Spots();
			// set threshold to be 2 times std of the image
			imp.setPosition(1, 1, t + 1);
			double threshold = threstime*imp.getStatistics().stdDev;
			for (Spot spot : spots[t]) {
				double value = imp.getProcessor().getPixelValue(spot.x, spot.y);
				if (value > threshold)
					// set a threshold for spotting
					out[t].add(spot);
			}
		}
		return out;
	}

	public Spots[] thresDetect(ImagePlus imp, int nt, int xmax, int ymax) {
		RoiManager rm = new RoiManager();
		int options = Measurements.CENTER_OF_MASS;
		ImagePlus imp_bn = imp.duplicate();
		IJ.run(imp_bn, "Make Binary", "method=Default background=Default calculate black");
		//IJ.run(imp_bn, "Watershed", "stack");
		//IJ.run(imp_bn, "Analyze Particles...", "size=7-Infinity show=Masks clear stack");
		//IJ.run(imp_bn, "Options...", "iterations=1 count=1 black do=Skeletonize stack");
		IJ.run(imp_bn, "Analyze Particles...", "size=10-Infinity display clear overlay add stack");
		rm.runCommand(imp_bn, "Show None");
		Spots spots[] = new Spots[nt];
		Roi[] roiArray = rm.getRoisAsArray();
		int count = 1;
		int dcount = 0;
		int lframe = -1;
		for (Roi roi : roiArray) {
			int cd = count/1000;
			if (cd>dcount) {
				dcount++;
				//IJ.log("Detected rois: " + count);
			}
			int frame = roi.getPosition() - 1;
			if (frame >= lframe) {
				lframe++;
				spots[frame] = new Spots();
			}
			imp_bn.setRoi(roi);
			int x = (int) Math.round(imp_bn.getStatistics(options).xCenterOfMass);
			int y = (int) Math.round(imp_bn.getStatistics(options).yCenterOfMass);
			//IJ.log("Adding: " + x +", " + y + ", " + frame);
			spots[frame].add(new Spot(x, y, frame));
			count++;
		}
		return spots;
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
					for (int k = -bd; k <= bd; k++) {
						for (int l = -bd; l <= bd; l++) {
							max = Math.max(max, ip.getPixelValue(x + k, y + l));
						}
					}
					// count the number of max pixels, if larger then 1, do not generate spot
					int count = 0;
					for (int k = -bd; k <= bd; k++) {
						for (int l = -bd; l <= bd; l++) {
							if (ip.getPixelValue(x+k, y+l) == max)
								count ++;
						}
					}
					// if the pixel is the "only" maximum in the neighborhood, add it to the list
					if (v == max && count < 2) {
						// skip it if not the only largest!
						spots[t].add(new Spot(x, y, t));
					}
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
	
	// estimate local density and the long_link
	public double estimateLocalDens(Spot current, ArrayList<Spot>[] spots, int t, int ts, int xmax, int ymax, boolean diffuse) {
		double long_link;
		int ccount = 0;
		int ncount = 0;
		int bound;
		if (diffuse == true) {
			bound = 8;
		}
		else {
			bound = 4;
		}
		
		int cx = current.x;
		int cy = current.y;
		for (Spot spot: spots[t]) {
			if (spot.x >= cx - bound && spot.x <= cx + bound && spot.y >= cy - bound && spot.y <= cy + bound
                    && spot.x >= 0 && spot.x <= xmax && spot.y >= 0 && spot.y <= ymax
                    && spot.x != cx && spot.y != cy) {
                ccount++;
            }
		}
		for (Spot spot: spots[t+ts]) {
			if (spot.x >= cx - bound && spot.x <= cx + bound && spot.y >= cy - bound && spot.y <= cy + bound
                    && spot.x >= 0 && spot.x <= xmax && spot.y >= 0 && spot.y <= ymax) {
                ncount++;
            }
		}
		if (ccount == 0) {
			if (ncount == 0) {
				long_link = 0;
			}
			else {
				long_link = (double)bound/(ncount*ts);
			}
		}
		else {
			long_link = (double)bound/(ccount+ncount*ts);
		}
		return long_link;
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
	
	// remove single point trajectories
//	public ArrayList<Spot> removeBadTraj(ArrayList<Spot> spots[], int nt) {
//		ArrayList<Spot> out[] = new Spots[nt];
//		for (int t = 0; t < nt; t++) {
//			out[t] = new Spots();
//			for (Spot current : spots[t]) {
//				if(current)
//			}
//		}
//		return out;
//	}
	
	// write serialized file
	public void writeSpots(ArrayList<Spot> spots[]) {
		// save the spots in the current directory
		try{
			FileOutputStream fos = new FileOutputStream("C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Project/data/spots.ser");
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(spots);
			oos.flush();
			oos.close();
		} catch (Exception ioe) {
			ioe.printStackTrace();
		}
	}

	
}
