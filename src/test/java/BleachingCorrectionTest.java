import ch.epfl.bii.ij2command.BleachingCorrection;
import ij.ImagePlus;
import net.imagej.ImageJ;
import ij.io.Opener;


public class BleachingCorrectionTest {

    public static void main(String... args ) throws Exception {
    	// create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        
        // show imageJ window
        ij.ui().showUI();
        
        // Open image
        String imagePath = "C:\\Users\\asus\\Desktop\\EPFL\\Course\\2023 Spring\\Bioimage informatics\\Project\\data\\l_4360_reorder.tif";
        Opener opener = new Opener();
        ImagePlus image = opener.openImage(imagePath);
        image.show();
        ij.command().run(BleachingCorrection.class, true);
    }
}