package uk.ac.ebi.pride.spectracluster.analysis.util;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 05.01.15.
 */
public class TestSpectrumExtractor {
    private static File testFile;

    @Before
    public void setUp() {
        testFile = new File(TestSpectrumExtractor.class.getClassLoader().getResource("test_spectra.clustering").getPath());
    }

    @Test
    public void testExtractSpectra() throws Exception {
        File tmpDirectory = new File("/tmp");

        // don't run the test if the environment isn't right
        if (!tmpDirectory.exists() || !tmpDirectory.isDirectory() || !tmpDirectory.canWrite())
            return;

        SpectrumExtractor spectrumExtractor = new SpectrumExtractor(tmpDirectory.getAbsolutePath());

        ClusteringFileReader clusteringFileReader = new ClusteringFileReader(testFile);
        List<ICluster> clusters = clusteringFileReader.readAllClusters();

        clusteringFileReader = new ClusteringFileReader(testFile);
        List<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>();
        listeners.add(spectrumExtractor);

        clusteringFileReader.readClustersIteratively(listeners);
        
        for (ICluster cluster : clusters) {
            File mgfFilePath = new File(tmpDirectory, cluster.getId() + ".mgf");
            
            Assert.assertTrue(mgfFilePath.exists());
            mgfFilePath.delete();
        }
    }
}
