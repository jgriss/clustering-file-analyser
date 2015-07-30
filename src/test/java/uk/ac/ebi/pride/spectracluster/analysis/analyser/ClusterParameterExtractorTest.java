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

/**
 * Created by jg on 14.07.14.
 */
public class ClusterParameterExtractorTest {
    private File testFile;
    private File unidentifiedTestfile;
    private IClusterSourceReader reader;

    @Before
    public void setUp() {
        try {
            testFile = TestUtilities.getTestfile();
            unidentifiedTestfile = new File(ClusterParameterExtractorTest.class.getClassLoader().getResource("unidentified.clustering").toURI());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testParameterExtractor() throws Exception {
        reader = new ClusteringFileReader(testFile);
        ClusterParameterExtractor parameterExtractor = new ClusterParameterExtractor();

        StringWriter writer = new StringWriter();
        parameterExtractor.setWriter(writer);

        ArrayList<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(parameterExtractor);

        reader.readClustersIteratively(listeners);

        parameterExtractor.completeResultFile();
        writer.close();

        String resultString = writer.toString();

        String[] lines = resultString.split("\n");

        Assert.assertEquals(961, lines.length);
        Assert.assertEquals("null\t305.000\t1.000\t2\t2\t0\t1.000\t1.000\t0.000\tYIAHLPAK:2\tYIAHLPAK\t2\t\tnull\t0\t1\t1\t", lines[1]);
        Assert.assertEquals("null\t305.010\t1.000\t1\t1\t0\t1.000\t1.000\t0.000\tKNYGK:1\tKNYGK\t1\t\tnull\t0\t1\t1\t", lines[11]);
    }

    @Test
    public void testUnidentifiedExtraction() throws Exception {
        reader = new ClusteringFileReader(unidentifiedTestfile);
        ClusterParameterExtractor parameterExtractor = new ClusterParameterExtractor();

        StringWriter writer = new StringWriter();
        parameterExtractor.setWriter(writer);

        ArrayList<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>(1);
        listeners.add(parameterExtractor);

        reader.readClustersIteratively(listeners);

        parameterExtractor.completeResultFile();
        writer.close();

        String resultString = writer.toString();

        Assert.assertEquals("id\tprecursor_mz\tprecursor_intensity\tsize\tidentified_spec_count\tunidentified_spec_count\tmax_ratio\tmax_il_ratio\tprecursor_mz_range\tsequences\tmax_sequence\tmax_sequence_count\tmax_sequence_mods\tsecond_max_sequence\tsecond_max_sequence_count\tproject_count\tassay_count\tspecies\n" +
                "b1d45cf2-a1ae-4084-bd8c-66ca54f255d1\t700.000\t1.000\t2\t2\t0\t1.000\t1.000\t0.000\tRQSSLTFQSSDPEHVR:2\tRQSSLTFQSSDPEHVR\t2\t0-MOD:01499,4-MOD:00696\tnull\t0\t1\t1\t10090\n" +
                "ccd02911-81b7-4583-9101-58a88cab1f04\t700.000\t1.000\t81\t51\t30\t0.961\t0.902\t0.005\tIMDPNIVGNEHYDVAR:1,EVYMGNVIQGGEGQAPTR:1,NVMLLPVGSADDGAHSQNEK:46,VIQFVCGNGLVCETMEEAR:3\tNVMLLPVGSADDGAHSQNEK\t46\t3-MOD:00719\tVIQFVCGNGLVCETMEEAR\t3\t12\t56\t10090,9606\n" +
                "61326dc1-6a3a-49b7-b668-5edd0dae5136\t700.000\t1.000\t5\t0\t5\t0.000\t0.000\t0.001\t\tnull\t0\t\tnull\t0\t1\t2\t338654\n" +
                "66c1ebcb-dffc-43b6-a862-3807a7d9847d\t700.000\t1.000\t4\t2\t2\t0.500\t0.500\t0.001\tLGVQDLFNSSK:1,VFVQKEILDKFTEEVVK:1\tLGVQDLFNSSK\t1\t\tVFVQKEILDKFTEEVVK\t1\t2\t4\t9606\n" +
                "268360b6-92d5-4976-8a15-a1b1456188aa\t700.000\t1.000\t3\t3\t0\t1.000\t1.000\t0.002\tIYEFPETDDEEENKLVK:3\tIYEFPETDDEEENKLVK\t3\t\tnull\t0\t1\t2\t9606\n" +
                "45fa0c9d-a5c0-48a9-b6c3-f596b7c65da6\t700.000\t1.000\t88\t0\t88\t0.000\t0.000\t0.040\t\tnull\t0\t\tnull\t0\t1\t5\t9606\n", resultString);
    }
}
