import ch.epfl.bii.ij2command.TirfTracking.TirfTracking;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;

// explanation of the code was written in ParticleTracking.java

public class TirfTrackingTest {
    public static void main(String... args ) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        // process on a background corrected image, remember to rearrange the hyperstack to frames instead of slices
        String imagePath = "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Project/data/l_4360_slice_test_subtract_time_projection_short.tif";
		ImagePlus imp = IJ.openImage(imagePath);
		imp.show();
		double lambda = 0; // lambda ~ 0 is the best !!!
        ij.command().run(TirfTracking.class, true, "lambda", lambda);
    }
}