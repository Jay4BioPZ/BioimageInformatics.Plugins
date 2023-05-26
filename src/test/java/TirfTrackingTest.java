import ch.epfl.bii.ij2command.TirfTracking.TirfTracking;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;

// explanation of the code was written in ParticleTracking.java

public class TirfTrackingTest {
    public static void main(String... args ) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(TirfTracking.class, true);
    }
}