package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.analysis.util.CrapFastaFile;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by jg on 14.07.14.
 *
 * Extracts basic properties of every cluster (size,
 * max ratio, precursor_mz, precursor_intens) and returns
 * them as a tab-delimited table.
 */
public class ExtendedClusterParameterExtractor extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".ext_parameters.tsv";
    public static String DESCRIPTION = "Extracts the extended properties of every cluster in the file and returns them in a TAB delimited table.";

    private ClusterUtilities clusterUtilities = new ClusterUtilities();

    public final char DELIMINATOR = '\t';
    public final String TABLE_HEADER=
                    "id" + DELIMINATOR +
                    "precursor_mz" + DELIMINATOR +
                    "size" + DELIMINATOR +
                    "identified_spec_count" + DELIMINATOR +
                    "unidentified_spec_count" + DELIMINATOR +
                    "max_ratio" + DELIMINATOR +
                    "max_il_ratio" + DELIMINATOR +
                    "precursor_mz_range" + DELIMINATOR +
                    "sequences" + DELIMINATOR +
                    "max_sequence" + DELIMINATOR +
                    "max_sequence_count" + DELIMINATOR +
                    "max_sequence_projects" + DELIMINATOR +
                    "max_sequence_contaminant" + DELIMINATOR +
                    "second_max_sequence" + DELIMINATOR +
                    "second_max_sequence_count" + DELIMINATOR +
                    "second_max_sequence_projects" + DELIMINATOR +
                    "second_max_sequence_contaminant" + DELIMINATOR +
                    "third_max_sequence" + DELIMINATOR +
                    "third_max_sequence_count" + DELIMINATOR +
                    "third_max_sequence_projects" + DELIMINATOR +
                    "third_max_sequence_contaminant" + DELIMINATOR +
                    "is_contaminant" + DELIMINATOR +
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
        clusterUtilities = new ClusterUtilities(newCluster);

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

        // get the project strings
        String projectsMaxS = getProjectStringForSequence(newCluster, clusterUtilities.getMaxSequence());
        String projectSecondS = getProjectStringForSequence(newCluster, clusterUtilities.getSecondMaxSequence());
        String projectThridS = getProjectStringForSequence(newCluster, clusterUtilities.getThirdMaxSequence());

        String maxContam = CrapFastaFile.getInstance().getProteinAnnotation(clusterUtilities.getMaxSequence());
        String secContam = clusterUtilities.getSecondMaxSequence() != null ? CrapFastaFile.getInstance().getProteinAnnotation(clusterUtilities.getSecondMaxSequence()) : "";
        String thirdContam = clusterUtilities.getThirdMaxSequence() != null ? CrapFastaFile.getInstance().getProteinAnnotation(clusterUtilities.getThirdMaxSequence()): "";

        // add the string representing the cluster to the result buffer
        writer.write(
                newCluster.getId() + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorMz()) + DELIMINATOR +
                String.format("%.3f", newCluster.getAvPrecursorIntens()) + DELIMINATOR +
                newCluster.getSpecCount() + DELIMINATOR +
                newCluster.getIdentifiedSpecCount() + DELIMINATOR +
                newCluster.getUnidentifiedSpecCount() + DELIMINATOR +
                String.format("%.3f", newCluster.getMaxRatio()) + DELIMINATOR +
                String.format("%.3f", clusterUtilities.getMaxILAngosticRatio()) + DELIMINATOR +
                String.format("%.3f", clusterUtilities.getMzRange()) + DELIMINATOR +

                clusterUtilities.getMaxSequence() + DELIMINATOR +
                clusterUtilities.getMaxSequenceCount() + DELIMINATOR +
                projectsMaxS + DELIMINATOR +
                (maxContam != null ? maxContam : "") + DELIMINATOR +

                clusterUtilities.getSecondMaxSequence() + DELIMINATOR +
                clusterUtilities.getSecondMaxSequenceCount() + DELIMINATOR +
                projectSecondS + DELIMINATOR +
                (secContam != null ? secContam : "") + DELIMINATOR +

                clusterUtilities.getThirdMaxSequence() + DELIMINATOR +
                clusterUtilities.getThirdMaxSequenceCount() + DELIMINATOR +
                projectThridS + DELIMINATOR +
                (thirdContam != null ? thirdContam : "") + DELIMINATOR +

                clusterUtilities.getnProjects() + DELIMINATOR +
                clusterUtilities.getnAssays() + DELIMINATOR +
                speciesString.toString() + "\n"
        );
    }

    private String getProjectStringForSequence(ICluster cluster, String sequence) {
        Set<String> projects = new HashSet<String>();

        if (sequence == null)
            return "";

        String targetSequence = sequence.replace("I", "L");
        targetSequence = ClusterUtilities.cleanSequence(targetSequence);

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                String psmSequence = ClusterUtilities.cleanSequence( psm.getSequence().replace("I", "L") );

                if (targetSequence.equals(psmSequence)) {
                    // extract the project
                    int index = specRef.getSpectrumId().indexOf(";");
                    if (index < 0)
                        continue;

                    String projectId = specRef.getSpectrumId().substring(0, index);
                    projects.add(projectId);
                }
            }
        }

        // create the result String
        StringBuilder projectString = new StringBuilder();
        for (String projectId : projects) {
            if (projectString.length() > 0)
                projectString.append(",");

            projectString.append(projectId);
        }

        return projectString.toString();
    }
}
