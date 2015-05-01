package uk.ac.ebi.pride.spectracluster.analysis.io;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.File;
import java.net.URI;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 01.05.15.
 */
public class TestMsClusterReader implements IClusterSourceListener {
    private File testFile;
    private int iterativeReads = 0;

    @Before
    public void setUp() throws Exception {
        URI testFileUri = TestMsClusterReader.class.getClassLoader().getResource("test.clust").toURI();
        testFile = new File(testFileUri);
    }

    @Test
    public void testReadAllClustes() throws Exception {
        MsClusterFileReader reader = new MsClusterFileReader(testFile);
        List<ICluster> clusters = reader.readAllClusters();

        Assert.assertEquals(20000, clusters.size());

        ICluster firstCluster = clusters.get(0);
        Assert.assertEquals(469, firstCluster.getPsmCount());
        Assert.assertEquals(469, firstCluster.getSpecCount());
        Assert.assertEquals(null, firstCluster.getMaxSequence());
        Assert.assertEquals(2.1880493F, firstCluster.getSpectrumPrecursorMzRange());
    }

    @Test
    public void readClustersIteratively() throws Exception {
        MsClusterFileReader reader = new MsClusterFileReader(testFile);
        List<IClusterSourceListener> listeners = new ArrayList<IClusterSourceListener>();
        listeners.add(this);
        iterativeReads = 0;
        reader.readClustersIteratively(listeners);

        Assert.assertEquals(20000, iterativeReads);
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        iterativeReads++;
    }
}
