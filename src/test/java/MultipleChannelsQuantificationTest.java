import ch.epfl.bii.ij2command.MultipleChannelsQuantification;
import ij.ImagePlus;
import ij.io.Opener;
import net.imagej.ImageJ;

public class MultipleChannelsQuantificationTest {
	public static void main(String... args ) throws Exception {
    	// create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        
        // show imageJ window
        ij.ui().showUI();
        
        // Open 2 images
        String imagePath1 = "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Homework/Homework C-Quantification 2 channels/exp2/nucleus.tif";
        String imagePath2 = "C:/Users/asus/Desktop/EPFL/Course/2023 Spring/Bioimage informatics/Homework/Homework C-Quantification 2 channels/exp2/cyto.tif";
        
        Opener opener = new Opener();
        ImagePlus image1 = opener.openImage(imagePath1);
        ImagePlus image2 = opener.openImage(imagePath2);
        image1.show();
        image2.show();
        
        // Pass two images into the Class
        ij.command().run(MultipleChannelsQuantification.class, true, "image1", image1, "image2", image2);
    }
}
