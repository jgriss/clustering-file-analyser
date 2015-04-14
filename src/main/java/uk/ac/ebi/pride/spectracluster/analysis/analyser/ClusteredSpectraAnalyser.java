package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by jg on 13.04.15.
 */
public class ClusteredSpectraAnalyser extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".clustered_spectra.tsv";
    public static String DESCRIPTION = "Extracts the unique # of clustered spectra in clusters >= 1, 3, 5, and 10, 50 spectra.";

    Map<Integer, Set<String>> clusteredSpectraCounts;
    Map<Integer, Long> totalClusterSizes;
    Map<Integer, Integer> clusterCounts;

    public final int[] MIN_CLUSTER_SIZES = {1, 3, 5, 10, 50};

    public ClusteredSpectraAnalyser() {
        reset();
    }

    @Override
    public String getAnalysisResultString() {
        StringBuilder resultString = new StringBuilder("min_cluster_size\tclustered_spectra\ttotal_cluster_size\tn_clusters\n");

        for (Integer minSize : clusteredSpectraCounts.keySet()) {
            resultString.append(minSize)
                    .append("\t")
                    .append(clusteredSpectraCounts.get(minSize).size())
                    .append("\t")
                    .append(totalClusterSizes.get(minSize))
                    .append("\t")
                    .append(clusterCounts.get(minSize))
                    .append("\n");
        }

        return resultString.toString();
    }

    @Override
    public void reset() {
        this.clusteredSpectraCounts = new HashMap<Integer, Set<String>>();
        this.totalClusterSizes = new HashMap<Integer, Long>();
        this.clusterCounts = new HashMap<Integer, Integer>();
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
                clusterCounts.put(minSize, 0);
            clusterCounts.put(minSize, clusterCounts.get(minSize) + 1);

            // store the spec ids in the clusters
            if (!clusteredSpectraCounts.containsKey(minSize))
                clusteredSpectraCounts.put(minSize, new HashSet<String>());

            for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
                clusteredSpectraCounts.get(minSize).add(specRef.getSpectrumId());
            }
        }
    }
}
