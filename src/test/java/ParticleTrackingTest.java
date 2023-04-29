import ch.epfl.bii.ij2command.ParticleTracking;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;

// explanation of the code was written in ParticleTracking.java

public class ParticleTrackingTest {
    public static void main(String... args ) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        String imagePath = "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Homework/Homework D-Tracking/homework.tif";
		ImagePlus imp = IJ.openImage(imagePath);
		imp.show();
		double lambda = 0; // lambda ~ 0 is the best !!!
        ij.command().run(ParticleTracking.class, true, "lambda", lambda);
    }
}