package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.*;

/**
 * Created by jg on 13.04.15.
 */
public class ClusteredSpectraAnalyser extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".clustered_spectra.tsv";
    public static String DESCRIPTION = "Extracts the unique # of clustered spectra in clusters >= 1, 3, 5, and 10, 50 spectra.";

    Map<Integer, Set<String>> clusteredSpectraCounts;
    Map<Integer, Long> totalClusterSizes;
    Map<Integer, Integer> clusterCounts;

    Map<Integer, Integer> mixedClusters;
    Map<Integer, Integer> incorrectlyAssignedSpectra;
    Map<Integer, Integer> totalSpectra;

    public final int[] MIN_CLUSTER_SIZES = {1, 3, 5, 10, 50};

    public ClusteredSpectraAnalyser() {
        reset();
    }

    @Override
    public String getAnalysisResultString() {
        StringBuilder resultString = new StringBuilder(
                "min_cluster_size\tclustered_spectra\ttotal_cluster_size\tn_clusters\t" +
                        "mixed_clusters\tincorrectly_assigned_spectra\ttotal_spectra\n");

        for (Integer minSize : clusteredSpectraCounts.keySet()) {
            resultString.append(minSize)
                    .append("\t")
                    .append(clusteredSpectraCounts.get(minSize).size())
                    .append("\t")
                    .append(totalClusterSizes.get(minSize))
                    .append("\t")
                    .append(clusterCounts.get(minSize))
                    .append("\t")
                    .append(mixedClusters.get(minSize))
                    .append("\t")
                    .append(incorrectlyAssignedSpectra.get(minSize))
                    .append("\t")
                    .append(totalSpectra.get(minSize))
                    .append("\n");
        }

        return resultString.toString();
    }

    @Override
    public void reset() {
        this.clusteredSpectraCounts = new HashMap<Integer, Set<String>>();
        this.totalClusterSizes = new HashMap<Integer, Long>();
        this.clusterCounts = new HashMap<Integer, Integer>();

        mixedClusters = new HashMap<Integer, Integer>();
        incorrectlyAssignedSpectra = new HashMap<Integer, Integer>();
        totalSpectra = new HashMap<Integer, Integer>();
    }

    @Override
    public String getFileEnding() {
        return FILE_ENDING;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        for (int minSize : MIN_CLUSTER_SIZES) {
            if (newCluster.getSpecCount() < minSize)
                continue;

            // add the size to the total cluster size
            if (!totalClusterSizes.containsKey(minSize))
                totalClusterSizes.put(minSize, 0L);
            totalClusterSizes.put(minSize, totalClusterSizes.get(minSize) + newCluster.getSpecCount());

            // add to total number of clusters
            if (!clusterCounts.containsKey(minSize))
                clusterCounts.put(minSize, 1);
            else
                clusterCounts.put(minSize, clusterCounts.get(minSize) + 1);

            // store the spec ids in the clusters
            if (!clusteredSpectraCounts.containsKey(minSize))
                clusteredSpectraCounts.put(minSize, new HashSet<String>());

            for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
                clusteredSpectraCounts.get(minSize).add(specRef.getSpectrumId());
            }

            String mostCommonSequence = extractMostCommonSequence(newCluster);

            boolean isMixedCluster = false;

            // simply count non-fitting spectra, only taking the first PSM into consideration
            for (ISpectrumReference spectrumReference : newCluster.getSpectrumReferences()) {
                String sequence = extractSpectrumReferenceSequence(spectrumReference);

                if (!mostCommonSequence.equals(sequence)) {
                    isMixedCluster = true;

                    if (!incorrectlyAssignedSpectra.containsKey(minSize))
                        incorrectlyAssignedSpectra.put(minSize, 1);
                    else
                        incorrectlyAssignedSpectra.put(minSize, incorrectlyAssignedSpectra.get(minSize) + 1);
                }

                if (!totalSpectra.containsKey(minSize))
                    totalSpectra.put(minSize, newCluster.getSpecCount());
                else
                    totalSpectra.put(minSize, totalSpectra.get(minSize) + newCluster.getSpecCount());
            }

            if (isMixedCluster) {
                if (!mixedClusters.containsKey(minSize))
                    mixedClusters.put(minSize, 1);
                else
                    mixedClusters.put(minSize, mixedClusters.get(minSize) + 1);
            }
        }
    }

    private String extractSpectrumReferenceSequence(ISpectrumReference spectrumReference) {
        String id = spectrumReference.getSpectrumId();
        int startPos = id.indexOf("splib_sequence=") + 15;

        if (startPos < 15)
            return spectrumReference.getPSMs().get(0).getSequence().toUpperCase().replaceAll("[^A-Z]", "");

        int endPos = id.indexOf(",", startPos);

        // get the sequence
        String sequence = id.substring(startPos, endPos);
        sequence = sequence.toUpperCase().replaceAll("[^A-Z]", "");

        return sequence;
    }

    private String extractMostCommonSequence(ICluster cluster) {
        // test whether the cluster contains specially identified spectra
        if (cluster.getSpectrumReferences().get(0).getSpectrumId().contains("splib_sequence=")) {
            Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();

            for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
                String sequence = extractSpectrumReferenceSequence(spectrumReference);

                if (!sequenceCounts.containsKey(sequence))
                    sequenceCounts.put(sequence, 1);
                else
                    sequenceCounts.put(sequence, sequenceCounts.get(sequence) + 1);
            }

            int maxCount = Collections.max(sequenceCounts.values());
            String maxSequence = null;

            for (String sequence : sequenceCounts.keySet()) {
                if (sequenceCounts.get(sequence) == maxCount) {
                    maxSequence = sequence;
                    break;
                }
            }

            return maxSequence;
        }
        else {
            return cluster.getMaxSequence().toUpperCase().replaceAll("[^A-Z]", "");
        }
    }
}
