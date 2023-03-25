import ch.epfl.bii.ij2command.SimpleFloodFilling;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import net.imagej.ImageJ;


public class SimpleFloodFillingTest {

    public static void main(String... args ) throws Exception {
    	// create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        IJ.run("Blobs (25K)");
        ImagePlus imp=WindowManager.getCurrentImage();
        IJ.setAutoThreshold(imp, "Default");
        IJ.run(imp, "Convert to Mask", "");
        ij.command().run(SimpleFloodFilling.class, true);
    }
}