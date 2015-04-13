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

    public final int[] MIN_CLUSTER_SIZES = {1, 3, 5, 10, 50};

    public ClusteredSpectraAnalyser() {
        this.clusteredSpectraCounts = new HashMap<Integer, Set<String>>();
    }

    @Override
    public String getAnalysisResultString() {
        StringBuilder resultString = new StringBuilder("min_cluster_size\tclustered_spectra\n");

        for (Integer minSize : clusteredSpectraCounts.keySet()) {
            resultString.append(minSize)
                    .append("\t")
                    .append(clusteredSpectraCounts.get(minSize).size())
                    .append("\n");
        }

        return resultString.toString();
    }

    @Override
    public void reset() {
        clusteredSpectraCounts = new HashMap<Integer, Set<String>>();
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

            if (!clusteredSpectraCounts.containsKey(minSize))
                clusteredSpectraCounts.put(minSize, new HashSet<String>());

            for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
                clusteredSpectraCounts.get(minSize).add(specRef.getSpectrumId());
            }
        }
    }
}
