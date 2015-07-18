package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 14.07.14.
 *
 * Extracts basic properties of every cluster (size,
 * max ratio, precursor_mz, precursor_intens) and returns
 * them as a tab-delimited table.
 */
public class ClusterParameterExtractor extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".cluster_parameters.tsv";
    public static String DESCRIPTION = "Extracts the basic properties of every cluster in the file and returns them in a TAB delimited table.";

    private ClusterUtilities clusterUtilities = new ClusterUtilities();

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
                    "max_sequence_mods" + DELIMINATOR +
                    "second_max_sequence" + DELIMINATOR +
                    "second_max_sequence_count" + DELIMINATOR +
                    "project_count" + DELIMINATOR +
                    "assay_count" + DELIMINATOR +
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
    protected String getResultFileHeader() {
        return TABLE_HEADER;
    }

    @Override
    public void completeResultFile() throws Exception {

    }

    @Override
    public void reset() {

    }

    @Override
    protected void processClusterInternally(ICluster newCluster) throws Exception {
        clusterUtilities.processCluster(newCluster);

        // build the sequence string
        StringBuilder sequenceString = new StringBuilder();
        for (String sequence : clusterUtilities.getSequenceCounts().keySet()) {
            if (sequenceString.length() > 0)
                sequenceString.append(",");

            sequenceString.append(sequence + ":" + clusterUtilities.getSequenceCounts().get(sequence));
        }

        // get the species
        StringBuilder speciesString = new StringBuilder();
        for (String species : clusterUtilities.getSpecies()) {
            if (species == null)
                continue;

            if (speciesString.length() > 0) {
                speciesString.append(",");
            }
            speciesString.append(species);
        }

        // get the max sequence ptms
        String maxSequencePtms = extractMaxSequencePtms(clusterUtilities, newCluster);

        // add the string representing the cluster to the result buffer
       writer.write(
                newCluster.getId() + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorMz()) + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorIntens()) + DELIMINATOR +
                newCluster.getSpecCount() + DELIMINATOR +
                String.format("%.3f", newCluster.getMaxRatio()) + DELIMINATOR +
                String.format("%.3f", clusterUtilities.getMaxILAngosticRatio()) + DELIMINATOR +
                String.format("%.3f", clusterUtilities.getMzRange()) + DELIMINATOR +
                sequenceString.toString() + DELIMINATOR +
                clusterUtilities.getMaxSequence() + DELIMINATOR +
                clusterUtilities.getMaxSequenceCount() + DELIMINATOR +
                maxSequencePtms + DELIMINATOR +
                clusterUtilities.getSecondMaxSequence() + DELIMINATOR +
                clusterUtilities.getSecondMaxSequenceCount() + DELIMINATOR +
                clusterUtilities.getnProjects() + DELIMINATOR +
                clusterUtilities.getnAssays() + DELIMINATOR +
                speciesString.toString() + "\n"
        );
    }

    private String extractMaxSequencePtms(ClusterUtilities clusterUtilities, ICluster newCluster) {
        String maxSequence = clusterUtilities.getMaxSequence();
        List<IModification> modifications = new ArrayList<IModification>();


        for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                if (!maxSequence.equals(psm.getSequence()))
                    continue;

                if (psm.getModifications().size() > modifications.size()) {
                    modifications.clear();
                    modifications.addAll(psm.getModifications());
                }
            }
        }

        // create the modification string
        Collections.sort(modifications);  // sorts according to position

        StringBuilder modStringBuilder = new StringBuilder();

        for (IModification mod : modifications) {
            if (modStringBuilder.length() > 0)
                modStringBuilder.append(",");

            modStringBuilder.append(String.valueOf(mod.getPosition())).append("+").append(mod.getAccession());
        }

        return modStringBuilder.toString();
    }
}
