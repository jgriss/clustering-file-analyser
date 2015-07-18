package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.analysis.util.ModificationMapper;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IModification;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.io.Writer;
import java.util.Set;

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
    private Set<String> modifications;
    private Set<String> species;
    protected Writer writer;
    private boolean hasWritternHeader = false;

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

        // test for modifications
        if (modifications != null && modifications.size() > 0) {
            boolean clusterHasMod = false;
            for (String modification : modifications) {
                clusterHasMod = clusterHasMod || hasPrimarySequenceMod(cluster, modification);
                if (clusterHasMod)
                    break;
            }

            // ignore the cluster if it doesn't have any of the mods
            if (!clusterHasMod)
                return true;
        }

        // test for species
        if (species != null && species.size() > 0) {
            boolean clusterHasSpecies = false;
            for (String s : species) {
                clusterHasSpecies = clusterHasSpecies || hasPrimarySpecies(cluster, s);
                if (clusterHasSpecies)
                    break;
            }

            // ignore the cluster if it doesn't have any of the mods
            if (!clusterHasSpecies)
                return true;
        }

        return false;
    }

    protected boolean hasPrimarySpecies(ICluster cluster, String species) {
        ClusterUtilities clusterUtils = new ClusterUtilities(cluster);

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                if (!clusterUtils.getMaxSequence().equals(psm.getSequence()))
                    continue;

                String[] specRefSpecies = specRef.getSpecies().split(",");
                for (String clusterSpecies : specRefSpecies) {
                    if (species.equals(clusterSpecies))
                        return true;
                }

                break; // only process every specRef once
            }
        }

        return false;
    }

    public static boolean hasPrimarySequenceMod(ICluster cluster, String modification) {
        Set<String> modificationAccessions = ModificationMapper.getInstance().getAccessionsForMapping(modification);

        // get the most common sequence
        ClusterUtilities clusterUtils = new ClusterUtilities(cluster);
        String maxSequence = clusterUtils.getMaxSequence();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                // only apply to the most common sequence
                if (!maxSequence.equals(psm.getSequence()))
                    continue;

                for (IModification mod : psm.getModifications()) {
                    if (modificationAccessions.contains(mod.getAccession()))
                        return true;
                }
            }
        }

        return false;
    }

    public Writer getWriter() {
        return writer;
    }

    @Override
    public final void onNewClusterRead(ICluster newCluster) {
        if (ignoreCluster(newCluster))
            return;

        try {
            if (!hasWritternHeader && writer != null) {
                writer.write(getResultFileHeader());
                hasWritternHeader = true;
            }

            processClusterInternally(newCluster);
        }
        catch(Exception e) {
            throw new IllegalStateException(e);        }
    }

    abstract protected void processClusterInternally(ICluster newCluster) throws Exception;

    abstract protected String getResultFileHeader();

    public void setWriter(Writer writer) {
        this.writer = writer;
    }

    public Set<String> getModifications() {
        return modifications;
    }

    public void setModifications(Set<String> modifications) {
        this.modifications = modifications;
    }

    public Set<String> getSpecies() {
        return species;
    }

    public void setSpecies(Set<String> species) {
        this.species = species;
    }
}
