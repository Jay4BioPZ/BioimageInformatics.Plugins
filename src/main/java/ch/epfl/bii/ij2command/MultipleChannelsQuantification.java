package ch.epfl.bii.ij2command;

/*
 *  BIO-410 Bioimage Informatics Homework C
 *  Author: Yung-Cheng Chiang
 *  Latest update: 2023/03/27
 *  
 *  Abstract:
 *  The plugin computes the gene expression (in form of fluorescence) of the cytoplasm 
 *  at the very close periphery of every nucleus in a sequence of images. The input should includes two stacks, 
 *  one for nucleus and the other for cytoplasm. The plugin first generates a mask of nucleus 
 *  and then overlays ROIs onto the cytoplasm image. The measured values would be plotted into 
 *  a time vs. intensity scatter plot.
 *  
 */

import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.ImageCalculator;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.measure.Measurements;
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
		// time frame for x-axis matching
		// double t = 1;
		double[] tlist = new double[nR];
		double[] meanlist = new double[nR];
		
		// Get ROI list
		Roi[] roiArray = rm.getRoisAsArray(); //将由make particles生成的、在roi manager内的rois存成一个array
		// Loop through all ROI
		// 对于每个在array内的单个roi, 进行以下操作
		for (Roi roi : roiArray) {
			// select ROI, get the frame number and rm index of this ROI
			// cimp是ImagePlus图档
			cimp.setRoi(roi); // 把它投放到cimp图上，会看到某个frame的某个位置圈出了黄色的outline
			int frame = roi.getPosition() - 1; // 记录这个roi所在的time frame，以便等一下画图
			int roiIndex = rm.getRoiIndex(roi); // 记录这个roi在roi manager里面所处的顺序，以便等一下画图
			// update ROI from outline to circular band
			IJ.run(cimp, "Make Band...", "band=5"); // 将这个选好的outline转换成5pixel的band
			Roi modifiedRoi = cimp.getRoi(); // 选取band
			rm.setRoi(modifiedRoi, roiIndex); // 用band取代outline
			// add the measurement to tlist and meanlist for later plotting
			tlist[roiIndex] = frame; // 把frame(时间, 也就是scatterplot的x-axis)存到tlist里面
			meanlist[roiIndex] = cimp.getStatistics(options).mean; // 测量band位置的mean, 存到meanlist里面
		}
		// 全部roi都改成band，并且测量、储存完数据之后再一并画scatterplot
		// composite image was created manually after the analysis
		// add data points in tlist and meanlist to scatterplot
		scatterPlot.addPoints(tlist, meanlist, PlotWindow.CIRCLE);
		scatterPlot.show();
	}
	
	public static void main(final String... args) throws Exception {
		
	}
}
