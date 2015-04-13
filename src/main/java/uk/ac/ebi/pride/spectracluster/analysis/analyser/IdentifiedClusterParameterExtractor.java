package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.*;

/**
 * Created by jg on 14.07.14.
 *
 * Extracts basic properties of every cluster (size,
 * max ratio, precursor_mz, precursor_intens) and returns
 * them as a tab-delimited table.
 * This class can only be run on clusters that were based
 * on the validation dataset (ie. search using SpectraST or
 * Pepitome)
 */
public class IdentifiedClusterParameterExtractor extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".id.cluster_parameters.tsv";
    public static String DESCRIPTION = "Extracts the basic properties of every cluster in the file and returns them in a TAB delimited table.";

    StringBuffer resultStringBuffer = new StringBuffer();
    public final char DELIMINATOR = '\t';
    public final String TABLE_HEADER=
                    "id" + DELIMINATOR +
                    "precursor_mz" + DELIMINATOR +
                    "precursor_intensity" + DELIMINATOR +
                    "size" + DELIMINATOR +
                    "max_ratio" + DELIMINATOR +
                    "max_il_ratio" + DELIMINATOR +
                    "precursor_mz_range" + DELIMINATOR +
                    "sequences" + DELIMINATOR +
                    "max_sequence" + DELIMINATOR +
                    "max_sequence_count" + DELIMINATOR +
                    "max_sequence_av_score" + DELIMINATOR +
                    "max_sequence_max_score" + DELIMINATOR +
                    "second_max_sequence" + DELIMINATOR +
                    "second_max_sequence_count" + DELIMINATOR +
                    "second_max_sequence_av_score" + DELIMINATOR +
                    "second_max_sequence_max_score" + DELIMINATOR +
                    "species" + "\n";


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
        return TABLE_HEADER + resultStringBuffer.toString();
    }

    @Override
    public void reset() {
        resultStringBuffer = new StringBuffer();
    }

    @Override
    public void onNewClusterRead(ICluster newCluster) {
        if (ignoreCluster(newCluster))
            return;

        List<SequenceCount> sequenceCounts = new ArrayList<SequenceCount>(getSequenceCounts(newCluster).values());
        Collections.sort(sequenceCounts);
        Collections.reverse(sequenceCounts);

        SequenceCount maxSequenceCount = sequenceCounts.get(0);
        SequenceCount secondSequenceCount = null;
        if (sequenceCounts.size() > 1)
            secondSequenceCount = sequenceCounts.get(1);

        // calculate the total size (including decoys)
        int size = 0;
        for (SequenceCount sequenceCount : sequenceCounts) {
            size += sequenceCount.getCount();
        }

        // build the sequence string
        StringBuilder sequenceString = new StringBuilder();
        for (SequenceCount sequenceCount : sequenceCounts) {
            if (sequenceString.length() > 0)
                sequenceString.append(",");

            sequenceString.append(sequenceCount.getSequence() + ":" + sequenceCount.getCount());
        }

        // get the species
        StringBuilder speciesString = new StringBuilder();
        for (String species : extractSpecies(newCluster)) {
            if (species == null)
                continue;

            if (speciesString.length() > 0) {
                speciesString.append(",");
            }
            speciesString.append(species);
        }

        // calculate the max ratio
        double maxRatio = (double) maxSequenceCount.getCount() / size;
        int maxIlAgnosticCount = getMaxIlAgnosticCount(sequenceCounts);

        // add the string representing the cluster to the result buffer
        resultStringBuffer.append(
                newCluster.getId() + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorMz()) + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorIntens()) + DELIMINATOR +
                size + DELIMINATOR +
                String.format("%.3f", maxRatio) + DELIMINATOR +
                String.format("%.3f", (float) maxIlAgnosticCount / size) + DELIMINATOR +
                String.format("%.3f", calculateMzRange(newCluster)) + DELIMINATOR +
                sequenceString.toString() + DELIMINATOR +
                maxSequenceCount.getSequence() + DELIMINATOR +
                maxSequenceCount.getCount() + DELIMINATOR +
                maxSequenceCount.getAvScore() + DELIMINATOR +
                maxSequenceCount.getMaxScore() + DELIMINATOR +
                (secondSequenceCount != null ? secondSequenceCount.getSequence() : "NA") + DELIMINATOR +
                (secondSequenceCount != null ? secondSequenceCount.getCount() : "NA") + DELIMINATOR +
                (secondSequenceCount != null ? secondSequenceCount.getAvScore() : "NA") + DELIMINATOR +
                (secondSequenceCount != null ? secondSequenceCount.getMaxScore() : "NA") + DELIMINATOR +
                speciesString.toString() + "\n"
        );
    }

    private double calculateMzRange(ICluster cluster) {
        double minMz = Double.MAX_VALUE, maxMz = 0;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            if (specRef.getPrecursorMz() > maxMz)
                maxMz = specRef.getPrecursorMz();
            if (specRef.getPrecursorMz() < minMz)
                minMz = specRef.getPrecursorMz();
        }

        return maxMz - minMz;
    }

    private int getMaxIlAgnosticCount(List<SequenceCount> sequenceCounts) {
        Map<String, Integer> sequenceCount = new HashMap<String, Integer>();

        for (SequenceCount sc : sequenceCounts) {
            String sequence = sc.getSequence().toUpperCase().replaceAll("I", "L");

            if (!sequenceCount.containsKey(sequence))
                sequenceCount.put(sequence, sc.getCount());
            else
                sequenceCount.put(sequence, sequenceCount.get(sequence) + sc.getCount());
        }

        return Collections.max(sequenceCount.values());
    }

    private Map<String, SequenceCount> getSequenceCounts(ICluster newCluster) {
        Map<String, SequenceCount> sequenceCounts = new HashMap<String, SequenceCount>();

        for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
            String id = specRef.getSpectrumId();

            // find the id
            int index = id.indexOf("splib_sequence=");

            if (index < 0)
                throw new IllegalStateException("Failed to extract identification from " + id);

            int endIndex = id.indexOf(",", index);

            if (endIndex < 0)
                throw new IllegalStateException("Failed to extract identification from " + id);

            String sequence = id.substring(index + 15, endIndex);

            // get the score
            int scoreIndex = id.indexOf("score=") + 6;
            if (scoreIndex < 6)
                scoreIndex = id.indexOf("mvh=") + 4;

            if (scoreIndex < 0)
                throw new IllegalStateException("Failed to extract score from " + id);

            int scoreEndIndex = id.indexOf(",", scoreIndex);
            if (scoreEndIndex < 0)
                throw new IllegalStateException("Failed to extract score from " + id);

            String scoreString = id.substring(scoreIndex, scoreEndIndex);

            double score;

            if (scoreString.equalsIgnoreCase("inf")) {
                score = 1.5;
            }
            else {
                score = Double.parseDouble(scoreString);
            }

            if (!sequenceCounts.containsKey(sequence))
                sequenceCounts.put(sequence, new SequenceCount(sequence, sequence.endsWith("_D")));

            sequenceCounts.get(sequence).addIdentification(score);
        }

        return sequenceCounts;
    }

    private Set<String> extractSpecies(ICluster cluster) {
        Set<String> species = new HashSet<String>();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            String s = specRef.getSpecies();

            if (s == null)
                continue;

            if (s.contains(",")) {
                String[] ss = s.split(",");
                for (String tmp : ss)
                    species.add(tmp);
            }
            else  {
                species.add(s);
            }
        }

        return species;
    }

    private class SequenceCount implements Comparable<SequenceCount> {
        private final String sequence;
        private final boolean isDecoy;
        private int count = 0;
        private double totalScore = 0;
        private double maxScore = 0;

        private SequenceCount(String sequence, boolean isDecoy) {
            this.sequence = sequence;
            this.isDecoy = isDecoy;
        }

        public void addIdentification(double score) {
            count++;
            totalScore += score;

            if (score > maxScore)
                maxScore = score;
        }

        public int getCount() {
            return count;
        }

        public double getAvScore() {
            return totalScore / count;
        }

        public double getMaxScore() {
            return maxScore;
        }

        @Override
        public int compareTo(SequenceCount o) {
            return Integer.compare(count, o.getCount());
        }

        public String getSequence() {
            return sequence;
        }

        public boolean isDecoy() {
            return isDecoy;
        }
    }
}
