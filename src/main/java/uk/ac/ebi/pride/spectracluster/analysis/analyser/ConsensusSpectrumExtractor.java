package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

/**
 * Created by jg on 18.07.15.
 */
public class ConsensusSpectrumExtractor extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".mgf";
    public static String DESCRIPTION = "Writes the clusters' consensus spectra to an MGF file";
    private int clusterCounter = 0;

    @Override
    protected void processClusterInternally(ICluster newCluster) throws Exception {
        ClusterUtilities clusterUtilities = new ClusterUtilities(newCluster);
        clusterCounter++;

        StringBuilder stringBuilder = new StringBuilder("BEGIN IONS\n");

        stringBuilder.append(String.format("TITLE=%s\n", (newCluster.getId() != null) ? newCluster.getId() : clusterUtilities.getMaxSequence(), clusterCounter));
        int charge = clusterUtilities.getCharge();
        stringBuilder.append(String.format("PEPMASS=%.3f\n", newCluster.getAvPrecursorMz()));
        stringBuilder.append(String.format("CHARGE=%d%c\n", Math.abs(charge), (charge > 0 ? '+' : '-')));
        stringBuilder.append(String.format("SEQUENCE=%s\n", clusterUtilities.getMaxSequence()));

        // add the peak list
        for (int i = 0; i < newCluster.getConsensusMzValues().size(); i++) {
            // ignore peaks with 0 m/z and 0 intensity
            if (newCluster.getConsensusMzValues().get(i) == 0) {
                System.out.println("Warning: Cluster " + newCluster.getId() + " contains empty peak (m/z).");
                continue;
            }
            if (newCluster.getConsensusIntensValues().get(i) == 0) {
                System.out.println("Warning: Cluster " + newCluster.getId() + " contains empty peak (intensity).");
                continue;
            }

            stringBuilder.append(newCluster.getConsensusMzValues().get(i)).append(" ").append(newCluster.getConsensusIntensValues().get(i)).append("\n");
        }

        stringBuilder.append("END IONS\n\n");

        writer.write(stringBuilder.toString());
    }

    @Override
    protected String getResultFileHeader() {
        return "";
    }

    @Override
    public void reset() {
        clusterCounter = 0;
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
    public void completeResultFile() throws Exception {

    }
}
