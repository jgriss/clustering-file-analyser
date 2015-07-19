package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
import java.io.StringWriter;
import java.util.List;

/**
 * Created by jg on 10.07.15.
 */
public class ProjectAnalyserTest {
    private File testFile;

    @Before
    public void setUp() {
        testFile = new File(ProjectAnalyserTest.class.getClassLoader().getResource("clusteringBin1709.clustering").getFile());
    }

    @Test
    public void testProjectAnalysis() throws Exception {
        ProjectAnalyser analyser = new ProjectAnalyser();

        StringWriter writer = new StringWriter();
        analyser.setWriter(writer);

        ClusteringFileReader clusteringFileReader = new ClusteringFileReader(testFile);
        List<ICluster> clusters = clusteringFileReader.readAllClusters();

        for (ICluster c : clusters) {
            analyser.onNewClusterRead(c);
        }

        analyser.completeResultFile();
        writer.close();

        Assert.assertEquals("project\tcorrect_ids\tincorrect_ids\tcorrect_contam\tincorrect_contam\tcorrect_psms\tincorrect_psms\tcorrect_contam_psms\tincorrect_contam_psms\n" +
                "PRD000073\t4\t0\t0\t0\t30\t0\t0\t0\n" +
                "PRD000280\t5\t0\t0\t0\t43\t0\t0\t0\n" +
                "PRD000748\t1\t0\t0\t0\t10\t0\t0\t0\n" +
                "PRD000349\t5\t0\t0\t0\t72\t0\t0\t0\n" +
                "PRD000019\t2\t0\t0\t0\t2\t0\t0\t0\n" +
                "PXD001910\t1\t0\t0\t0\t3\t0\t0\t0\n" +
                "PRD000545\t1\t0\t0\t0\t11\t0\t0\t0\n" +
                "PRD000478\t4\t0\t0\t0\t48\t0\t0\t0\n" +
                "PRD000044\t2\t0\t0\t0\t16\t0\t0\t0\n" +
                "PRD000397\t1\t0\t0\t0\t11\t0\t0\t0\n" +
                "PRD000562\t1\t0\t0\t0\t14\t0\t0\t0\n" +
                "PRD000097\t5\t0\t0\t0\t101\t0\t0\t0\n", writer.toString());
    }
}
