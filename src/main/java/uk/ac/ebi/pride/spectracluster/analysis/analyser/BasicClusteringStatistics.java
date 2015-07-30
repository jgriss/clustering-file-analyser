package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by jg on 12.07.14.
 */
public class BasicClusteringStatistics extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".basic_statistics.txt";
    public static String DESCRIPTION = "Generates basic statistics about a .clustering file like # clusters.";
    public static final float LARGE_PRECURSOR_MZ_RANGE = 1.5F;

    private int allIdentifiedClusters = 0;
    private int unidentifiedClusters = 0;
    private int mixedClusters = 0;
    private float nClusters = 0;
    private int nSingleClusters = 0;
    private float averageRatio = 0;
    private float averageClusterSize = 0;
    private int minSize = Integer.MAX_VALUE;
    private int maxSize = 0;
    private float minRatio = Float.MAX_VALUE;
    private float maxRatio = 0;
    private int stableClusters = 0;
    private float maxPrecursorMzRange = 0F;
    private long totalNumberOfSpectra = 0;
    private long totalNumberOfMismatchedSpectra = 0;
    private int cleanClusters = 0;
    private Map<String, Integer> speciesCounts = new HashMap<String, Integer>();

    /**
     * Number of clusters with a large precursor m/z range
     */
    private int nLargePrecursorMzRange = 0;
    /**
     * used to calculate the average precursor m/z range
     * after the run.
     */
    private float totalPrecursorMzRange = 0F;


    @Override
    public String getFileEnding() {
        return FILE_ENDING;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

    @Override
    public void processClusterInternally(ICluster newCluster) {
        nClusters++;

        if (minSize > newCluster.getSpecCount())
            minSize = newCluster.getSpecCount();
        if (maxSize < newCluster.getSpecCount())
            maxSize = newCluster.getSpecCount();

        // ignore clusters with only 1 spectrum for any other parameter
        if (newCluster.getSpecCount() < 2) {
            nSingleClusters++;
            return;
        }

        if (newCluster.getIdentifiedSpecCount() < 1)
            unidentifiedClusters++;
        else if (newCluster.getUnidentifiedSpecCount() < 1)
            allIdentifiedClusters++;
        else
            mixedClusters++;

        totalNumberOfSpectra += newCluster.getSpecCount();
        averageRatio = averageRatio / nClusters * (nClusters - 1) + newCluster.getMaxRatio() / nClusters;
        averageClusterSize = averageClusterSize / nClusters * (nClusters - 1) + newCluster.getSpecCount() / nClusters;
        totalPrecursorMzRange += newCluster.getSpectrumPrecursorMzRange();
        if (maxPrecursorMzRange < newCluster.getSpectrumPrecursorMzRange())
            maxPrecursorMzRange = newCluster.getSpectrumPrecursorMzRange();
        if (newCluster.getSpectrumPrecursorMzRange() >= LARGE_PRECURSOR_MZ_RANGE)
            nLargePrecursorMzRange++;

        if (minRatio > newCluster.getMaxRatio())
            minRatio = newCluster.getMaxRatio();
        if (maxRatio < newCluster.getMaxRatio())
            maxRatio = newCluster.getMaxRatio();

        processSpecies(newCluster);

        ClusterUtilities clusterUtilities = new ClusterUtilities(newCluster);

        if (clusterUtilities.getMaxILAngosticRatio() == 1) {
            cleanClusters++;
        }
        else {
            int mismatchedSpectra = Math.round  ( (1 - clusterUtilities.getMaxILAngosticRatio()) * newCluster.getSpecCount());
            totalNumberOfMismatchedSpectra += mismatchedSpectra;
        }

        if (clusterUtilities.isStable())
            stableClusters++;
    }

    private void processSpecies(ICluster newCluster) {
        for (ISpectrumReference spectrumReference : newCluster.getSpectrumReferences()) {
            String speciesString = spectrumReference.getSpecies();

            if (speciesString == null)
                continue;

            String[] speciesFields = speciesString.split(",");

            Set<String> uniqueSpecies = new HashSet<String>();
            for (String speciesField : speciesFields)
                uniqueSpecies.add(speciesField);

            for (String species : uniqueSpecies) {
                if (!speciesCounts.containsKey(species))
                    speciesCounts.put(species, 1);
                else
                    speciesCounts.put(species, speciesCounts.get(species) + 1);
            }
        }
    }

    @Override
    protected String getResultFileHeader() {
        return "";
    }

    @Override
    public void completeResultFile() throws Exception {
        String resultString = String.format("Number of clusters: %.0f (%d with 1 spec)\n" +
                        "All identified clusters: %d\n" +
                        "All unidentified clusters: %d\n" +
                        "Mixed clusters: %d\n" +
                        "Average maximum ratio: %.3f\n" +
                        "Average cluster size: %.3f\n" +
                        "Minimum size: %d\nMaximum size: %d\n" +
                        "Minimum ratio: %.3f\nMaximum ratio: %.3f\n" +
                        "Stable clusters: %d\n" +
                        "Average precursor m/z range: %.2f\n" +
                        "Max. precursor m/z range: %.2f\n" +
                        "Clusters with precursor m/z range > %.1f: %d\n" +
                        "Mismatched spectra: %.2f%%\n" +
                        "Clean clusters: %.2f%%\n",
                nClusters, nSingleClusters, allIdentifiedClusters, unidentifiedClusters, mixedClusters, averageRatio, averageClusterSize, minSize, maxSize,
                minRatio, maxRatio, stableClusters, totalPrecursorMzRange / (nClusters - nSingleClusters),
                maxPrecursorMzRange, LARGE_PRECURSOR_MZ_RANGE, nLargePrecursorMzRange,
                (float) (totalNumberOfMismatchedSpectra * 100 / totalNumberOfSpectra), (float) (cleanClusters * 100 / nClusters));

        for (String species : speciesCounts.keySet()) {
            resultString += "Species " + species + ": " + speciesCounts.get(species) + "\n";
        }

        writer.write(resultString);
    }

    @Override
    public void reset() {
        nClusters = 0;
        nSingleClusters = 0;
        averageRatio = 0;
        averageClusterSize = 0;
        minSize = Integer.MAX_VALUE;
        maxSize = 0;
        minRatio = Float.MAX_VALUE;
        maxRatio = 0;
        stableClusters = 0;
        totalPrecursorMzRange = 0;
        maxPrecursorMzRange = 0;
        nLargePrecursorMzRange = 0;
        totalNumberOfSpectra = 0;
        totalNumberOfMismatchedSpectra = 0;
        cleanClusters = 0;
    }
}
