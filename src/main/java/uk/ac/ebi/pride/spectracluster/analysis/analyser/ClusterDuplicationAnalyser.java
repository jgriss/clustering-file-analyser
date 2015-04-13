package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jg on 15.07.14.
 *
 * Analysis in how many cluster a spectrum was observed in.
 */
public class ClusterDuplicationAnalyser extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".duplication_analysis.tsv";
    public static String DESCRIPTION = "Counts in how many clusters a certain spectrum was found.";

    Map<String, Integer> specIdCounts = new HashMap<String, Integer>();
    public static final char DELIMINATOR = '\t';
    public static final String TABLE_HEADER_COMPLETE = "spec_id" + DELIMINATOR + "cluster_count" + "\n";
    public static final String TABLE_HEADER_SUMMARY = "count" + DELIMINATOR + "frequency" + "\n";

    @Override
    public String getFileEnding() {
        return FILE_ENDING;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

    @Override
    public String getAnalysisResultString() {
        Map<Integer, Integer> countFrequency = new HashMap<Integer, Integer>();

        for (Integer count : specIdCounts.values()) {
            if (!countFrequency.containsKey(count)) {
                countFrequency.put(count, 0);
            }

            countFrequency.put(count, countFrequency.get(count) + 1);
        }

        // build the table
        StringBuilder stringBuilder = new StringBuilder(TABLE_HEADER_SUMMARY);

        for (Integer count : countFrequency.keySet()) {
            stringBuilder.append(String.format("%d%c%d\n",
                    count, DELIMINATOR,
                    countFrequency.get(count)));
        }

        return stringBuilder.toString();
    }

    /**
     * Returns the analysis result as a Map with the spectrum
     * id as key and the respective number of clusters the spectrum
     * was observed in as value.
     * @return
     */
    public Map<String, Integer> getSpecIdCounts() {
        return Collections.unmodifiableMap(specIdCounts);
    }

    /**
     * Returns a tab delimited table (as String) containing the
     * spectrum id as the first column and the respective number
     * of clusters it was found in as the second.
     * @return
     */
    public String getSpecIdCountsAsTable() {
        // build the string
        StringBuilder stringBuilder = new StringBuilder(TABLE_HEADER_COMPLETE);

        for (String sequence : specIdCounts.keySet()) {
            stringBuilder.append(String.format("%s%c%d\n",
                    sequence, DELIMINATOR,
                    specIdCounts.get(sequence)));
        }

        return stringBuilder.toString();
    }

    @Override
    public void reset() {
        specIdCounts.clear();
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        if (ignoreCluster(newCluster))
            return;

        for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
            String id = specRef.getSpectrumId();

            if (!specIdCounts.containsKey(id)) {
                specIdCounts.put(id, 0);
            }

            specIdCounts.put(id, specIdCounts.get(id) + 1);
        }
    }
}
