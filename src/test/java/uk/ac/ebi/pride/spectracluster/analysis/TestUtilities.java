package uk.ac.ebi.pride.spectracluster.analysis;

import java.io.File;
import java.net.URI;

/**
 * Created by jg on 14.07.14.
 */
public class TestUtilities {
    private TestUtilities() {

    }

    public static File getTestfile() throws Exception {
        URI testFileUri = TestUtilities.class.getClassLoader().getResource("testfile.clustering").toURI();
        File testFile = new File(testFileUri);

        return testFile;
    }
}
