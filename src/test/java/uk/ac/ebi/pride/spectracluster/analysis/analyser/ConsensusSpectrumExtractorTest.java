package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.analysis.TestUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.StringWriter;
import java.util.ArrayList;

/**
 * Created by jg on 18.07.15.
 */
public class ConsensusSpectrumExtractorTest {
    private File testFile;
    private IClusterSourceReader reader;
    private String referenceMgfString;

    @Before
    public void setUp() {
        try {
            testFile = new File(ConsensusSpectrumExtractorTest.class.getClassLoader().getResource("clusteringBin1709.clustering").toURI());
            reader = new ClusteringFileReader(testFile);

            BufferedReader reader = new BufferedReader(new InputStreamReader(ConsensusSpectrumExtractorTest.class.getClassLoader().getResourceAsStream("clusteringBin1709.clustering.mgf")));
            String line;
            StringBuilder mgfString = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                mgfString.append(line + "\n");
            }
            reader.close();
            referenceMgfString = mgfString.toString();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testConsensusSpectrumExtraction() throws Exception {
        AbstractClusteringSourceAnalyser analyser = new ConsensusSpectrumExtractor();

        StringWriter writer = new StringWriter();
        analyser.setWriter(writer);

        ArrayList<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(analyser);

        reader.readClustersIteratively(listeners);

        analyser.completeResultFile();
        writer.close();

        Assert.assertEquals(referenceMgfString, writer.toString());
    }
}
