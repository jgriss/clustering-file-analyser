package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.analysis.TestUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceReader;

import java.io.File;
import java.io.StringWriter;
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

        StringWriter writer = new StringWriter();
        clusteredSpectraAnalyser.setWriter(writer);

        List<IClusterSourceListener> analysers = new ArrayList<IClusterSourceListener>(1);
        analysers.add(clusteredSpectraAnalyser);

        reader.readClustersIteratively(analysers);

        clusteredSpectraAnalyser.completeResultFile();
        writer.close();
        String resultString = writer.toString();

        Assert.assertEquals("min_cluster_size\tclustered_spectra\ttotal_cluster_size\tn_clusters\tmixed_clusters\tincorrectly_assigned_spectra\ttotal_spectra\n" +
                "1\t7234\t8927\t960\t322\t4203\t8927\n" +
                "50\t3320\t4292\t36\t35\t2557\t4292\n" +
                "3\t6410\t8013\t375\t299\t4180\t8013\n" +
                "5\t6096\t7609\t267\t230\t4055\t7609\n" +
                "10\t5617\t6970\t170\t154\t3770\t6970\n", resultString);
    }
}
