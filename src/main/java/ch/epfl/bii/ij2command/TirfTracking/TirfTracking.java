
package ch.epfl.bii.ij2command.TirfTracking;

/*
 *  BIO-410 Bioimage Informatics Miniproject - TIRF protein tracking
 */

import java.util.ArrayList;

import org.apache.commons.io.filefilter.TrueFileFilter;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import java.nio.file.Path;
import java.nio.file.Paths;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.plugin.ImageCalculator;
import ij.process.ImageProcessor;
import ij.plugin.frame.RoiManager;
import ij.measure.Measurements;

// for object serialization
import java.io.*;

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>TIRFTracking")
public class TirfTracking implements Command {
	
	public void run() {
		// generate a generic dialog
		GenericDialog gd = new GenericDialog("TIRF particle tracking");
		gd.addMessage("TIRF particle tracking and trajectories analysis. Please enter the following parameters to proceed.\n ");
		gd.setInsets(0, 50, 0);
		gd.addFileField("Image stack path (.tif): ", "C:\\EPFL\\2022-2023 Microtechnique MA4\\Bioimage informatics\\l_4360-crop_outer_long_traj.tif");
		gd.addNumericField("Lambda (0-1): ", 0.2);
		gd.addNumericField("Gamma (0-1): ", 0.3);
		gd.setInsets(0, 180, 0);
		gd.addCheckbox("Exponential correction", false);
		gd.setInsets(0, 180, 0);
		gd.addCheckbox("Diffusing particle", false);
		gd.addMessage("\nAuthors: Bioimage Informatics Group 8. Y.Chiang, C.Liu, K.Aydin (2023.05)");
		gd.showDialog();
		if(gd.wasCanceled()) {
			return;
		}
		
		String imagePath = (String)gd.getNextString();
		Path path = Paths.get(imagePath);
        Path folderPath = path.getParent();
        String folder = folderPath.toString();
        
		double lambda = (double)gd.getNextNumber();
		double gamma = (double)gd.getNextNumber();
		boolean expcor = (boolean)gd.getNextBoolean();
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
		IJ.run("Coordinates...", "left=0 right="+xmax+" top=0 bottom="+ymax+" front=0 back="+nz);
		imp.setDimensions(1, 1, nz);
		int nt = imp.getNFrames();
		IJ.log("The dataset have " + nt + " time frames");
		// STEP D2. Exponential correction
		//IJ.log("Launch exponential correction...");
		if (expcor == true) {
			imp = imagecorrection(imp, nt, xmax, ymax);
		}
		// STEP D3. DoG
		IJ.log("Start particle detection...");
		ImagePlus imp_dog = dog(imp, sigma1, sigma2);
		// STEP D4. Gaussian blur 3D
		IJ.run(imp_dog, "Gaussian Blur 3D...", "x=0 y=0 z=1");
		// STEP D5. TopHat filtering
		IJ.run(imp_dog, "Top Hat...", "radius=2 stack");
		
		ImagePlus imp_dog_traj = imp_dog.duplicate();
		ImagePlus imp_dog_arrw = imp_dog.duplicate();
		
		// STEP D6. Detect spots either with local max or center of mass method
		ArrayList<Spot> spots[];
		if (diffuse == true) {
			spots = thresDetect(imp_dog, nt, xmax, ymax);
		}
		else {
			int bd = 7;
			int threstime = 8;
			ArrayList<Spot> localmax[] = localMax(imp_dog, bd);
			spots = localMaxfilter(imp_dog, localmax, threstime);
		}
		
		// STEP D7. Prepare the orientation term images
		ImagePlus imp_future = foreseeDirection(imp);
		imp_future.show();
		
		// LINKING
		IJ.log("Start particle linking...");
		// cost function
		double f1max = findMax(imp, nt);
		//IJ.log("f1max "+f1max);
		double f2max = findMax(imp_future, nt);
		//IJ.log("f2max "+f2max);
		double dmax = Math.sqrt(xmax*xmax + ymax*ymax);
		//IJ.log("dmax "+dmax);
		double intdiff = 0;
		double dirvalue = 0;
		double temp_cost = 1;
		double long_link;
		
		for (int t = 0; t < nt - 1; t++) {
			int indc = 0;
			for (Spot current : spots[t]) {
				double min_cost = 1;
				int min_ind = 0;
				int min_ts = 1;
				boolean findspot = false;
				// find candidates in next "few frames"
				for (int ts = 1; ts < Math.min(nt-t, 6); ts++) {
					// STEP L1. Estimate the local density of spots and compute a good long_link
					long_link = estimateLocalDens(current, spots, t, ts, xmax, ymax, diffuse);
					//System.out.println(long_link);
					int indn = 0;
					for (Spot next : spots[t+ts]) {
						// STEP L2. Compute each cost terms
						intdiff = findIntdiff(imp, t, current, next);
						dirvalue = findDirValue(imp_future, t, next);
						// STEP L3. Compute the cost and locate the best candidate
						temp_cost = getCOST(current, next, dmax, f1max, lambda, intdiff, f2max, gamma, dirvalue);
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
				// STEP L4. Form links
				if (findspot == true) {
					int nextt = t+min_ts;
					System.out.println("Link (" + t + ", " + indc + ") to (" + nextt + ", " + min_ind + ")");
					Spot good_next = spots[t+min_ts].get(min_ind);
					current.link(good_next);
				}
				else {
					//System.out.println("No link found, long_link");
				}
				indc++;
			}
		}
		
		IJ.log("Tracing individual trajectories...");
		ArrayList<ArrayList<int[]>> Arrayofall = new ArrayList<>();
		
		// TRAJECTORY ANALYSIS
		// STEP T1. Convert Spots to a list of trajectories for later parameter computation
		// Get through all the frames.
		// In an order of finding one spot and tracking all of its trajectory points, then switch to the next beginning spot.
		// Add Arrayofone into Arrayofall.
		
		for (int t = 0; t < nt-1; t++) {
		    for (Spot current : spots[t]) {
		        if (current.track == false) {
		        	Spot first = current;
		        	while(first.prev != null) {
		        		first = first.prev;
		        	}
		        	assignDirections(first);
		            ArrayList<int[]> Arrayofone = new ArrayList<>();
		            Arrayofone = findNext(first, Arrayofone);
		            Arrayofall.add(Arrayofone);
		        }
		    }
		}
		
		
		// STEP T2. Characterize the trajectories with diffusion coefficient and speed
		double[] diffCoeffs = diffusionCalculator(Arrayofall);
		double[] speeds = speedCalculator(Arrayofall);
		
		IJ.log("Saving the data...");
		
		// STEP T3. Save results as csv
		try {
			saveCoordinates(Arrayofall, speeds, diffCoeffs, folder);
		} catch (Exception e) {
			System.out.println("Error writing to the file");
		}
	
		IJ.log("Saved the data");
		
		IJ.log("Plotting...");
		
		// STEP T5. Plot histograms
		drawHistogram(diffCoeffs, "Diffusion Coefficients", "Pixel^2/s");
		
		drawHistogram(speeds, "Speeds", "Pixel/s");
		IJ.log("Plots done");
		
		IJ.log("Finalizing the visualization... (takes tens of minutes for raw images)");
		// Create overlay of trajectories
		Overlay overlay_traj = new Overlay();
		draw(overlay_traj, spots, radius);
		// Overlay direction arrows of the long trajectories (displacement > 5 pixel)
		Overlay overlay_arrw = new Overlay();
		drawDirections(overlay_arrw, spots, radius);
		imp_dog_traj.setOverlay(overlay_traj);
		imp_dog_arrw.setOverlay(overlay_arrw);
		imp_dog_traj.show();
		imp_dog_arrw.show();
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
	
	// STPE T2. Draw direction arrows
	private void drawDirections(Overlay overlay, ArrayList<Spot> spots[], int radius) {
		int nt = spots.length;
		for (int t = 0; t < nt; t++) {
			// print the time point
			//System.out.println("overlaying time point: " + (t+1));
			for (Spot spot : spots[t]) {
				// call the draw function in the Spot class
				spot.drawDirections(overlay, radius);
			}
		}
	}

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

	// STEP D6a. Center of mass detection
	public Spots[] thresDetect(ImagePlus imp, int nt, int xmax, int ymax) {
		RoiManager rm = new RoiManager();
		int options = Measurements.CENTER_OF_MASS;
		ImagePlus imp_bn = imp.duplicate();
		IJ.run(imp_bn, "Make Binary", "method=Default background=Default calculate black");
		//IJ.run(imp_bn, "Watershed", "stack");
		//IJ.run(imp_bn, "Analyze Particles...", "size=7-Infinity show=Masks clear stack");
		//IJ.run(imp_bn, "Options...", "iterations=1 count=1 black do=Skeletonize stack");
		IJ.run(imp_bn, "Analyze Particles...", "size=3-Infinity display clear overlay add stack");
		rm.runCommand(imp_bn, "Show None");
		Spots spots[] = new Spots[nt];
		Roi[] roiArray = rm.getRoisAsArray();
//		int num_roi = roiArray.length;
//		if (num_roi >= 200*nt) {
//			System.out.println("Too many particles ("+num_roi/nt+"/frame) detected. Try the 'Noisy data' option.");
//			System.exit(1);
//		}
//		else {
//			IJ.log("Particles ("+num_roi/nt+"/frame) are detected.");
//		}
		int count = 1;
		int dcount = 0;
		int lframe = -1;
		for (Roi roi : roiArray) {
			int cd = count/2000;
			if (cd>dcount) {
				dcount++;
				IJ.log("Detected rois: " + count);
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

	// Step D6b. LocalMax
	private ArrayList<Spot>[] localMaxfilter(ImagePlus imp, ArrayList<Spot> spots[], double threstime) {
		int nt = spots.length;
		int count = 1;
		int dcount = 0;
		ArrayList<Spot> out[] = new Spots[nt];
		for (int t = 0; t < nt; t++) {
			out[t] = new Spots();
			// set threshold to be 2 times std of the image
			imp.setPosition(1, 1, t + 1);
			double threshold = threstime*imp.getStatistics().stdDev;
			for (Spot spot : spots[t]) {
				double value = imp.getProcessor().getPixelValue(spot.x, spot.y);
				if (value > threshold) {
					// set a threshold for spotting
					int cd = count/2000;
					if (cd>dcount) {
						dcount++;
						IJ.log("Detected rois: " + count);	
					}
					count++;
					if (count >= 200*nt) {
						System.out.println("Too many particles ("+count/nt+"/frame) detected. Check your input.");
						System.exit(1);
					}
					out[t].add(spot);
				}
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
	
	// STEP D7. Preparation of the orentation map
	// convolution using a [0,0,0,0,0,0,3,3,4,4,5] filter on time axis
	public ImagePlus foreseeDirection(ImagePlus imp) {
		ImagePlus future = imp.duplicate();
		IJ.run(future, "Reslice [/]...", "output=1.000 start=Top avoid");
		IJ.run(future, "Convolve...", "text1=0\n0\n0\n0\n0\n0\n3\n3\n4\n4\n5\n normalize stack");
		IJ.run(future, "Reslice [/]...", "output=1.000 start=Top avoid");
		IJ.run(future, "Invert", "stack");
		IJ.run(future, "8-bit", "");
		return future;
	}
	
	// STEP L1. Estimate local density and give a good long_link
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
		//long_link = 5;
		return long_link;
	}
	

	
	// STEP L2. Find the intensity term for specific pixel
	public double findIntdiff(ImagePlus imp, int t, Spot current, Spot next) {
		imp.setSlice(t+1);
		ImageProcessor pimp1 = imp.getProcessor();
		double v1 = pimp1.getPixelValue(current.x, current.y);
		//IJ.log("intdiff v1 "+v1);
		imp.setSlice(t+2);
		ImageProcessor pimp2 = imp.getProcessor();
		double v2 = pimp2.getPixelValue(next.x, next.y);
		//IJ.log("intdiff v2 "+v2);
		double intdiff = Math.abs(v1-v2);
		return intdiff;
	}
	
	// STEP L2. Find the orientation term for specific pixel
	public double findDirValue(ImagePlus future, int t, Spot next) {
		future.setSlice(t+1);
		//IJ.log("DirValue x y: "+ next.x + next.y);
		ImageProcessor pimp = future.getProcessor();
		double value = pimp.getPixelValue(next.x, next.y);
		return value;
	}

	// STEP L3. Compute the cost function
	public double getCOST(Spot current, Spot next, double dmax, double f1max, double lambda, double intdiff, double f2max, double gamma, double dirvalue) {
		//IJ.log("intdiff "+intdiff+" dirvalue "+dirvalue);
		double cost;
		double dist = current.distance(next);
		// do something
		cost = (double)(1-lambda-gamma)*(dist/dmax) + lambda*(intdiff)/(f1max) + gamma*(dirvalue)/(f2max);
		return cost;
	}
	
	// STEP T1. Recursively track the trajectory starting from specific points
	// Use recursion to go through all current points with next point.
	// Add coordinates into ArrayOfOne.
	public ArrayList<int[]> findNext(Spot P, ArrayList<int[]> ArrayOfOne) {
	    int[] coordinates = { P.x, P.y, P.t };
	    ArrayOfOne.add(coordinates);
	    P.track = true;

	    if (P.next == null) {
	        return ArrayOfOne;
	    } 
	    else {
	        ArrayOfOne = findNext(P.next, ArrayOfOne);
	    }

	    return ArrayOfOne;
	}
	
	// STEP T2. Compute diffusion parameter
	public double[] diffusionCalculator(ArrayList<ArrayList<int[]>> spotlist) {
		int spotNumber = spotlist.size();
		
		//diffusion coefficient for each detected particle
		double[] diffusionCoefficients = new double[spotNumber];
		
		for(int i = 0; i<spotNumber; i++) {
			
			int min_t = Integer.MAX_VALUE;
			int max_t = -1;
			
			int x_1=0;
			int x_2=0;
			int y_1=0;
			int y_2=0;
			
			double mse_distance = 0;
			
			for(int[] spotsnapshot : spotlist.get(i)) {
				
				//find the first and last occurence of the same spot
				if(spotsnapshot[2]<min_t) {
					min_t=spotsnapshot[2];
					x_1=spotsnapshot[0];
					y_1=spotsnapshot[1];
				} if(spotsnapshot[2]>max_t) {
					max_t=spotsnapshot[2];
					x_2=spotsnapshot[0];
					y_2=spotsnapshot[1];
				}
			}
			
			//calculate distance between start and the end
			mse_distance = Math.pow(x_1-x_2, 2) + Math.pow(y_1-y_2, 2);
			
			//calculate diffusion coefficient: d=4Dt
			if(max_t-min_t==0) {
				diffusionCoefficients[i]=0;
				continue;
			}
			
			diffusionCoefficients[i] = mse_distance/(4*(max_t-min_t));
		}
		
		
		double avg=0;
		
		for (int i = 0; i < diffusionCoefficients.length; i++) {
            avg += diffusionCoefficients[i];
        }
 
        double average = avg / diffusionCoefficients.length;
        
        IJ.log("Average diffusion coefficient is: " + average);
		
		return diffusionCoefficients;
	}
	
	// STEP T2. Compute speed parameter
	public double[] speedCalculator(ArrayList<ArrayList<int[]>> spotlist) {
		int spotNumber = spotlist.size();
		
		double[] speeds = new double[spotNumber];
		
		for(int i = 0; i<spotNumber; i++) {
			
			int min_t = Integer.MAX_VALUE;
			int max_t = -1;
			
			double distance = 0;
			
			for(int j=1; j<spotlist.get(i).size(); j++) {
				
				 int[] spotsnapshot = spotlist.get(i).get(j);
				 int[] previoussnapshot = spotlist.get(i).get(j-1);
				
				if(spotsnapshot[2]<min_t) {
					min_t=spotsnapshot[2];
				} if(spotsnapshot[2]>max_t) {
					max_t=spotsnapshot[2];
				}
				
				//calculates the distance advanced between the previous spot and this one
				//keeps adding
				distance += Math.sqrt(Math.pow(spotsnapshot[0]-previoussnapshot[0], 2) 
									+ Math.pow(spotsnapshot[1]-previoussnapshot[1], 2));
			}
			
			if(max_t-min_t==0) {
				//In case there is a single data for the spot
				speeds[i]=0;
				continue;
			}
			
			//linear speed is total distance divided by total time
			speeds[i] = distance/(max_t-min_t);
		}
		
		double avg=0;
		
		for (int i = 0; i < speeds.length; i++) {
            avg += speeds[i];
        }
 
        double average = avg / speeds.length;
        
        IJ.log("Average speed is: " + average);
		
		return speeds;
	}
	
	// STEP T5. Compute direction vector
	public void assignDirections(Spot spot) {
		int x_1 = spot.x;
		int y_1 = spot.y;
		
		Spot last=spot;
		
		while(last.next != null) {
			last = last.next;
		}
		
		int x_2 = last.x;
		int y_2 = last.y;
		
		double distance = spot.distance(last);
		
		if(distance<5) return;
		
		while(spot.next != null) {
			spot.direction = new int[] {x_1-x_2, y_1-y_2};
			spot = spot.next;
		}
	}
	
	// STEP T4. Plot histogram
	public void drawHistogram(double[] values, String title, String xTitle) {
		Plot p = new Plot(title, xTitle, "Counts");
		
		p.addHistogram(values);
		p.show();
	}
	
	// STEP T3. Save trajectory coordinates and parameters in csv file
	//col1: spot number, col2: x, col3: y, col4: t
	public void saveCoordinates(ArrayList<ArrayList<int[]>> spotlist, double[] speeds, double[] diffusionCoeffs, String folder) throws IOException {
		String filename = folder + "\\Spots.csv";
		File csvFile = new File(filename);
		System.out.println(csvFile.getAbsolutePath());
		FileWriter writer = new FileWriter(csvFile);
		
		int spotnumber = 0;
		
		for (ArrayList<int[]> spot : spotlist) {
			
			++spotnumber;
			
			for(int[] spotsnapshot : spot) {
				
				StringBuilder line = new StringBuilder();
				
				if(spotsnapshot.length < 3) continue;
				
				line.append(String.valueOf(spotnumber));
				line.append(',');
			    
			    for (int i = 0; i < spotsnapshot.length; i++) {
			    	
			        line.append(String.valueOf(spotsnapshot[i]));
			        line.append(',');
			    }
			    
		        line.append(String.valueOf(speeds[spotnumber]));
		        line.append(',');
		        line.append(String.valueOf(diffusionCoeffs[spotnumber]));
			    line.append("\n");
			    
			    writer.write(line.toString());
			}
			
		}
		
		writer.close();
	}
}
