import ch.epfl.bii.ij2command.ParticleTracking;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;

public class ParticleTrackingTest {
    public static void main(String... args ) throws Exception {
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        String imagePath = "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Homework/Homework D-Tracking/homework.tif";
		ImagePlus imp = IJ.openImage(imagePath);
		imp.show();
        ij.command().run(ParticleTracking.class, true);
    }
}