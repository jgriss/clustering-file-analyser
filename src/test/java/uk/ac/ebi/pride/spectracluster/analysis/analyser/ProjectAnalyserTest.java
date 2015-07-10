package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
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

        ClusteringFileReader clusteringFileReader = new ClusteringFileReader(testFile);
        List<ICluster> clusters = clusteringFileReader.readAllClusters();

        for (ICluster c : clusters) {
            analyser.onNewClusterRead(c);
        }

        Assert.assertEquals("project\tcorrect_ids\tincorrect_ids\tcorrect_contam\tincorrect_contam\n" +
                "PRD000073\t4\t0\t0\t0\n" +
                "PRD000280\t5\t0\t0\t0\n" +
                "PRD000748\t1\t0\t0\t0\n" +
                "PRD000349\t5\t0\t0\t0\n" +
                "PRD000019\t2\t0\t0\t0\n" +
                "PXD001910\t1\t0\t0\t0\n" +
                "PRD000545\t1\t0\t0\t0\n" +
                "PRD000478\t4\t0\t0\t0\n" +
                "PRD000044\t2\t0\t0\t0\n" +
                "PRD000397\t1\t0\t0\t0\n" +
                "PRD000562\t1\t0\t0\t0\n" +
                "PRD000097\t5\t0\t0\t0\n", analyser.getAnalysisResultString());
    }
}
