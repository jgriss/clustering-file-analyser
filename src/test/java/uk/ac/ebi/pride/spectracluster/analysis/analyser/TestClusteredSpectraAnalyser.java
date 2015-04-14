package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.analysis.TestUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceReader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 13.04.15.
 */
public class TestClusteredSpectraAnalyser {
    private File testFile;
    private IClusterSourceReader reader;

    @Before
    public void setUp() {
        try {
            testFile = TestUtilities.getTestfile();
            reader = new ClusteringFileReader(testFile);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testAnalysis() throws Exception {
        ClusteredSpectraAnalyser clusteredSpectraAnalyser = new ClusteredSpectraAnalyser();

        List<IClusterSourceListener> analysers = new ArrayList<IClusterSourceListener>(1);
        analysers.add(clusteredSpectraAnalyser);

        reader.readClustersIteratively(analysers);

        String resultString = clusteredSpectraAnalyser.getAnalysisResultString();

        Assert.assertEquals("min_cluster_size\tclustered_spectra\ttotal_cluster_size\tn_clusters\n" +
                "1\t7234\t8927\t960\n" +
                "50\t3320\t4292\t36\n" +
                "3\t6410\t8013\t375\n" +
                "5\t6096\t7609\t267\n" +
                "10\t5617\t6970\t170\n", resultString);
    }
}
