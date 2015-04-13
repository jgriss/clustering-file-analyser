package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

/**
 * Created by jg on 06.11.14.
 */
abstract public class AbstractClusteringSourceAnalyser implements IClusteringSourceAnalyser {
    private int minClusterSize = Integer.MIN_VALUE;
    private int maxClusterSize = Integer.MAX_VALUE;
    private float minClusterRatio = Float.MIN_VALUE;
    private float maxClusterRatio = Float.MAX_VALUE;
    private float minPrecursorMz = Float.MIN_VALUE;
    private float maxPrecursorMz = Float.MAX_VALUE;

    public int getMinClusterSize() {
        return minClusterSize;
    }

    public void setMinClusterSize(int minClusterSize) {
        this.minClusterSize = minClusterSize;
    }

    public int getMaxClusterSize() {
        return maxClusterSize;
    }

    public void setMaxClusterSize(int maxClusterSize) {
        this.maxClusterSize = maxClusterSize;
    }

    public float getMinClusterRatio() {
        return minClusterRatio;
    }

    public void setMinClusterRatio(float minClusterRatio) {
        this.minClusterRatio = minClusterRatio;
    }

    public float getMaxClusterRatio() {
        return maxClusterRatio;
    }

    public void setMaxClusterRatio(float maxClusterRatio) {
        this.maxClusterRatio = maxClusterRatio;
    }

    public float getMinPrecursorMz() {
        return minPrecursorMz;
    }

    public void setMinPrecursorMz(float minPrecursorMz) {
        this.minPrecursorMz = minPrecursorMz;
    }

    public float getMaxPrecursorMz() {
        return maxPrecursorMz;
    }

    public void setMaxPrecursorMz(float maxPrecursorMz) {
        this.maxPrecursorMz = maxPrecursorMz;
    }

    protected boolean ignoreCluster(ICluster cluster) {
        if (cluster.getSpecCount() > maxClusterSize)
            return true;

        if (cluster.getSpecCount() < minClusterSize)
            return true;

        if (cluster.getAvPrecursorMz() > maxPrecursorMz)
            return true;

        if (cluster.getAvPrecursorMz() < minPrecursorMz)
            return true;

        if (cluster.getMaxRatio() > maxClusterRatio)
            return true;

        if (cluster.getMaxRatio() < minClusterRatio)
            return true;

        return false;
    }
}
