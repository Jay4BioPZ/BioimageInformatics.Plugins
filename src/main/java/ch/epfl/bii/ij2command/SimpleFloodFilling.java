package ch.epfl.bii.ij2command;

import net.imagej.ImageJ;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.platform.PlatformService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.process.FloodFiller;
import ij.process.ImageProcessor;

import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 *  This example illustrates how to create an ImageJ 2 {@link Command} plugin.
 * The pom file of this project is customized for the course BioImageInformatics BIO-410 @EPFL
 * https://edu.epfl.ch/coursebook/fr/bioimage-informatics-BIO-410
 * 
 * <p>
 * The code here is performing a simple recursive flood filling algorithm (4-connections). 
 * The command can be tested in the java DummyCommandTest class.
 * 
 * </p>
 */

// decorator 
@Plugin(type = Command.class, menuPath = "Plugins>BII 2023>FloodFiller")
public class SimpleFloodFilling implements Command {

    ImageProcessor outIP;
    
    public void run() {
    	ImagePlus imp=WindowManager.getCurrentImage();
    	if (imp==null) {IJ.showMessage("Input image needed");return;}
    	ImageProcessor inputIP=imp.getProcessor();
    	outIP=inputIP.duplicate();
    	int width=outIP.getWidth();
    	int height=outIP.getHeight();
    	int count=1;
    	for (int nx=0;nx<width;nx++) {
    		for (int ny=0;ny<height;ny++) {
    			if (outIP.getPixel(nx, ny)==255) {
    				checkConnection(nx,ny,count,inputIP);
    				count++;
    			}
    		}
    	}
    	
    	ImagePlus out=new ImagePlus("Test",outIP);
    	out.show();
    	
    	

    }
    
    ImageProcessor checkConnection(int x,int y,int value, ImageProcessor ip) {
    	
    	int valueSeed=(int) outIP.getValue(x, y);
    	outIP.putPixel(x,y,value);
    	if (outIP.getValue(x+1,y)==valueSeed) {
    		checkConnection(x+1,y,value,ip);
    	}
    	if (outIP.getValue(x-1,y)==valueSeed) {
    		checkConnection(x-1,y,value,ip);
    	}
    	if (outIP.getValue(x,y+1)==valueSeed) {
    		checkConnection(x,y+1,value,ip);
    	}
    	if (outIP.getValue(x,y-1)==valueSeed) {
    		checkConnection(x,y-1,value,ip);
    	}
    	return ip;
    }

    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    public static void main(final String... args) throws Exception {
        
    }
    
}
