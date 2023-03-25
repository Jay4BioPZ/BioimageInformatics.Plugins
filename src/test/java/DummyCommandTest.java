import net.imagej.ImageJ;
import ch.epfl.bii.ij2command.DummyCommand;

public class DummyCommandTest {

    public static void main(String... args ) throws Exception {
    	 // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();

        ij.command().run(DummyCommand.class, true);
    }
}