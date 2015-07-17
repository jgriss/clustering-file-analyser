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
 * Created by jg on 15.07.14.
 */
public class ClusterDuplicationAnalyserTest {
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
        ClusterDuplicationAnalyser analyzer = new ClusterDuplicationAnalyser();

        StringWriter writer = new StringWriter();
        analyzer.setWriter(writer);

        List<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(analyzer);

        reader.readClustersIteratively(listeners);

        analyzer.completeResultFile();
        writer.close();
        String result = writer.toString();

        Assert.assertEquals("count\tfrequency\n" +
                "1\t5990\n" +
                "2\t872\n" +
                "3\t302\n" +
                "4\t63\n" +
                "5\t7\n", result);

        Assert.assertEquals(1, (int) analyzer.getSpecIdCounts().get("PXD000637;PRIDE_Exp_Complete_Ac_33209.xml;spectrum=1410"));
    }
}
