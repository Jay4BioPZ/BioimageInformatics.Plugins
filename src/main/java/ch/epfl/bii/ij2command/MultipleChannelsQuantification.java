package ch.epfl.bii.ij2command;

import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.ImageCalculator;
import ij.process.ImageProcessor;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.gui.Plot;
import ij.gui.PlotWindow;

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>MultipleChannelsQuantification")
public class MultipleChannelsQuantification implements Command {
	@Parameter
	private ImagePlus image1; // nucleus
	
	@Parameter
	private ImagePlus image2; // cytoplasm
	
	@Override
	public void run() {
		// duplicate the input
		ImagePlus nimp = image1.duplicate();
		ImagePlus cimp = image2.duplicate();
		// return the outline of the nucleus
		ImagePlus nmask = findnuclearmask(nimp);
		nmask.show();
		cytomeasure(cimp, nmask);
		// pseudo code
		/* blue channel
		 * preprocessing -> thresholding -> segmentation -> create the periphery band
		 * 
		 * Preprocessing, to remove background and shot noise
		 * - DoG: G7-G10
		 * - Make Binary - Huang method
		 * - Binary operation: closing -> fill holes
		 * x Dilate 4 - original photo -> obtain a mask of 5-pixel band around the nucleus
		 * - Binary find outline
		 * 
		 * Overlay band mask as ROI to cyto channel
		 * x Convert image from 32-bit to 8-bit
		 * - Analyze > Analyze Particles > click add to manager
		 * 		IJ.run(imp, "Analyze Particles...", "add stack");
		 * 
		 * green channel
		 * - Go through every ROI, convert it into a band
		 * 		// Go to image cyto (maybe duplicate it)
		 * 		nR = rm.count;
		 * 		t = 1;
		 * 		For ROI r: (start from zero)
		 * 		{
		 * 			// Select ROI, get time information
		 * 			rm.select(r);
		 * 			rName = rm.getName()
		 * 			t = //get the time variable from the name string
		 * 			
		 * 			// Make a band of 5 pixels and measure the mean
		 * 			IJ.run("Make Band...", "band=5"); 
		 * 			rm.runCommand(imp, "Update");
		 * 			rm.runCommand(imp,"Measure");
		 * 			
		 * 			// (macro) getResult("Mean", i)
		 * 			// Plot the dot onto the time, intensity plot
		 * 			
		 * 		}
		 */
	}
	
	public ImagePlus findnuclearmask(ImagePlus imp) {
		ImagePlus Gaussian_blur_1 = imp.duplicate();
		ImagePlus Gaussian_blur_2 = imp.duplicate();
		// Perform DoG filtering
		double sigma1 = 7;
		//double sigma2 = Math.sqrt(2)*sigma1;
		double sigma2 = 12;
		IJ.run(Gaussian_blur_1, "Gaussian Blur...", "sigma=" + sigma1 + " stack");
		IJ.run(Gaussian_blur_2, "Gaussian Blur...", "sigma=" + sigma2 + " stack");
		ImagePlus DoG = ImageCalculator.run(Gaussian_blur_1, Gaussian_blur_2, "Subtract create 32-bit stack");
		//DoG.show(); // 32-bit
	
		// Make Binary (8-bit) + binary operations
		IJ.run(DoG, "Make Binary", "method=Huang background=Dark calculate black");
		IJ.run(DoG, "Options...", "iterations=1 count=1 black do=Close stack");
		IJ.run(DoG, "Options...", "iterations=1 count=1 black do=[Fill Holes] stack");
		
		// To make nucleus outline a bit bigger, not overlapping with the nucleus
		IJ.run(DoG, "Options...", "iterations=1 count=1 black do=Dilate stack");
		
		// Return outline of the nucleus
		IJ.run(DoG, "Options...", "iterations=2 count=1 black do=Outline stack");
		return DoG;
	}
	
	public void cytomeasure(ImagePlus cimp, ImagePlus nmask) {
		// Create roi manager and measurement table option
		RoiManager rm = new RoiManager();
		int options = Measurements.MEAN;
		
		// Create a new scatter plot
		Plot scatterPlot = new Plot("Results", "Time (frame)", "Intensity (A.U.)");
		
		int sizethres = 50;
		// Analyze particle, put outline mask as ROI to ROI manager
		// Set a threshold for size to get rid of small noisy outlines
		IJ.run(nmask, "Analyze Particles...", "size=" + sizethres + "-Infinity add stack");
		rm.runCommand(nmask,"Show None");
		// Go through and update ROIs from outline to band
		// During the update, calculate the mean intensity of fluorescence and plot it
		int nR = rm.getCount();
		//System.out.println("nR: " + nR);
		// time frame for x-axis matching
		double t = 1;
		double[] tlist = new double[nR];
		double[] meanlist = new double[nR];
		
		// Get ROI list
		Roi[] roiArray = rm.getRoisAsArray();
		// Loop through all ROI
		for (Roi roi : roiArray) {
			// add individual ROI onto the cyto image
			cimp.setRoi(roi);
			// get time information, update t if necessary
			String rName = roi.getName();
			String rframe = rName.substring(0, 4).replaceFirst("^0+(?!$)", "");
			int rt = Integer.parseInt(rframe);
			if (rt != t) {
				t = t+1;
			}
			// update ROI from outline to circular band
			IJ.run(cimp, "Make Band...", "band=5");
			Roi modifiedRoi = cimp.getRoi();
			int roiIndex = rm.getRoiIndex(roi);
			rm.setRoi(modifiedRoi, roiIndex);
			tlist[roiIndex] = t;
			meanlist[roiIndex] = cimp.getStatistics(options).mean;
		}
		
		//cimp.show();
		
		// add data points in tlist and meanlist to scatterplot
		scatterPlot.addPoints(tlist, meanlist, PlotWindow.CIRCLE);
		scatterPlot.show();
	}
	
	public static void main(final String... args) throws Exception {
		
	}
}
