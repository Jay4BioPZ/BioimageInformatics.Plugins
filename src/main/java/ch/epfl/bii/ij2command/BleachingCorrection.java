package ch.epfl.bii.ij2command;
/*
 *  BIO-410 Bioimage Informatics Homework B
 *  Author: Yung-Cheng Chiang
 *  Latest update: 2023/03/18
 *  
 *  Abstract:
 *  The purpose of this code is to perform bleaching correction of a time-seried image using ij.measure.CurveFitter
 *  
 *  model: y = A*exp(-t/tau) + C
 *  
 */


import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.measure.CurveFitter;

@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>BleachCorrection")
public class BleachingCorrection implements Command {
	
	// global variables
	ImagePlus imp = IJ.getImage();
	int nt = imp.getNFrames();
	int nx = imp.getWidth();
	int ny = imp.getHeight();
	
	// the function that is implemented when clicking the plugin
	// the interface Command have a run() method
	@Override
	public void run() {
		imagecorrection();
	}
	
	public void imagecorrection() {
		// plot the z-intensity profile before processing
		IJ.run(imp, "Select All", "");
		IJ.run(imp, "Plot Z-axis Profile", "");
		IJ.run(imp, "32-bit", "");
		IJ.run(imp, "Make Montage...", "columns=8 rows=7 scale=0.25");
		IJ.run("Enhance Contrast", "saturated=0.35");
		ImagePlus out = imp.duplicate();
		out.show();
		
		// Get parameters
		double[] params = getparameter(imp, nt);
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
					new_v = (new_v-C)*Math.exp(t/tau);
					iout.putPixelValue(x, y, new_v);
				}
			}	
		}
		out.setSlice(1);
		out.setTitle("out");
		// plot the z-intensity profile after processing
		IJ.run(out, "Select All", "");
		IJ.run(out, "Plot Z-axis Profile", "");
		IJ.run(out, "Make Montage...", "columns=8 rows=7 scale=0.25");
		IJ.run("Enhance Contrast", "saturated=0.35");
	}
	
	// fitting with imageJ CurveFitter
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
	
	public double[] getparameter(ImagePlus imp, int nt) {
		// find the mean of the 1, nt/2, nt frame
		double[] mean_stack = new double[3];
		int[] ind = {0,nt/2,nt-1};
		for (int i = 0; i < 3; i++) {
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
	
	// man() is the entry point of a java code. 
	// intentionally empty
	public static void main(final String... args) throws Exception {
		
	}
}