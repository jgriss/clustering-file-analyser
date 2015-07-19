package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.analysis.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.analysis.util.CrapFastaFile;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;

import java.util.*;

/**
 * This analyser counts the number of correctly and incorrectly
 * identified spectra per project. Separate counts for known
 * contaminants are created.
 *
 * Created by jg on 10.07.15.
 */
public class ProjectAnalyser extends AbstractClusteringSourceAnalyser {
    public static String FILE_ENDING = ".project_analysis.tsv";
    public static String DESCRIPTION = "Estimates the number of correctly and incorrectly identified spectra per project.";

    private ClusterUtilities clusterUtilities = new ClusterUtilities();

    private Map<String, ProjectProperties> projectPropertiesMap = new HashMap<String, ProjectProperties>();

    @Override
    protected String getResultFileHeader() {
        return "project\tcorrect_ids\tincorrect_ids\tcorrect_contam\tincorrect_contam\t" +
                "correct_psms\tincorrect_psms\tcorrect_contam_psms\tincorrect_contam_psms\n";
    }

    @Override
    public void completeResultFile() throws Exception {
        StringBuilder stringBuilder = new StringBuilder();

        for (String project : projectPropertiesMap.keySet()) {
            ProjectProperties properties = projectPropertiesMap.get(project);

            stringBuilder.append(project)
                    .append("\t")
                    .append(String.valueOf(properties.getCorrectIds()))
                    .append("\t")
                    .append(String.valueOf(properties.getIncorrectIds()))
                    .append("\t")
                    .append(String.valueOf(properties.getCorrectContam()))
                    .append("\t")
                    .append(String.valueOf(properties.getIncorrectContam()))
                    .append("\t")
                    .append(String.valueOf(properties.getCorrectPSMs()))
                    .append("\t")
                    .append(String.valueOf(properties.getIncorrectPSMs()))
                    .append("\t")
                    .append(String.valueOf(properties.getCorrectContamPsms()))
                    .append("\t")
                    .append(String.valueOf(properties.getIncorrectContamPsms()))
                    .append("\n");
        }

        writer.write(stringBuilder.toString());
    }

    @Override
    public void reset() {
        projectPropertiesMap = new HashMap<String, ProjectProperties>();
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
    protected void processClusterInternally(ICluster newCluster) throws Exception {
        clusterUtilities.processCluster(newCluster);

        // only process reliable clusters
        if (!clusterUtilities.isStable())
            return;

        Map<String, Integer> sequenceCounts = clusterUtilities.getSequenceCounts();
        sequenceCounts = makeSequenceCountsIlAgnostic(sequenceCounts);

        // check whether it's a decoy cluster
        String protein = CrapFastaFile.getInstance().getProteinAnnotation(clusterUtilities.getMaxSequence());
        boolean isContaminant = protein != null;

        // get the correct and incorrect projects
        String maxSequence = clusterUtilities.getMaxSequence().replaceAll("I", "L");
        Set<String> correctProjects = getProjectsForSequence(newCluster, maxSequence);
        Set<String> incorrectProjects = new HashSet<String>();

        for (String sequence : sequenceCounts.keySet()) {
            sequence = ClusterUtilities.cleanSequence(sequence).replaceAll("I", "L");

            if (sequence.equals(maxSequence))
                continue;

            incorrectProjects.addAll(getProjectsForSequence(newCluster, sequence));
        }

        // remove potential ambiguous results
        incorrectProjects.removeAll(correctProjects);

        // save the results
        for (String project : correctProjects) {
            ProjectProperties projectProperties;

            if (!projectPropertiesMap.containsKey(project))
                projectPropertiesMap.put(project, new ProjectProperties());

            projectProperties = projectPropertiesMap.get(project);

            projectProperties.incrementCorrectIds();
            if (isContaminant)
                projectProperties.incrementCorrectContam();
        }

        for (String project : incorrectProjects) {
            ProjectProperties projectProperties;

            if (projectPropertiesMap.containsKey(project))
                projectProperties = projectPropertiesMap.get(project);
            else
                projectProperties = new ProjectProperties();

            projectProperties.incrementIncorrectIds();
            if (isContaminant)
                projectProperties.incrementIncorrectContam();
        }

        // save everything again, only this time per PSM
        String correctSequence = clusterUtilities.getMaxSequence().replaceAll("I", "L");
        for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
            boolean isCorrect = false;

            for (IPeptideSpectrumMatch psm : specRef.getPSMs()) {
                if (psm.getSequence().replaceAll("I", "L").equals(correctSequence)) {
                    isCorrect = true;
                    break;
                }
            }

            // get the project
            int index = specRef.getSpectrumId().indexOf(";");
            if (index < 0)
                continue;

            String projectId = specRef.getSpectrumId().substring(0, index);

            // the project should already exist
            if (!projectPropertiesMap.containsKey(projectId))
                projectPropertiesMap.put(projectId, new ProjectProperties());

            ProjectProperties projectProperties = projectPropertiesMap.get(projectId);

            if (isCorrect) {
                projectProperties.incrementCorrectPSMs();
                if (isContaminant)
                    projectProperties.incrementCorrectContamPSMs();
            }
            else {
                projectProperties.incrementIncorrectPSMs();
                if (isContaminant)
                    projectProperties.incrementIncorrectContamPSMs();
            }
        }
    }

    private Map<String, Integer> makeSequenceCountsIlAgnostic(Map<String, Integer> sequenceCounts) {
        Map<String, Integer> ilAngosticSequenceCounts = new HashMap<String, Integer>();

        for (String sequence : sequenceCounts.keySet()) {
            String ilAgnostic = sequence.replaceAll("I", "L");

            if (ilAngosticSequenceCounts.containsKey(ilAgnostic))
                ilAngosticSequenceCounts.put(ilAgnostic, ilAngosticSequenceCounts.get(ilAgnostic) + sequenceCounts.get(sequence));
            else
                ilAngosticSequenceCounts.put(ilAgnostic, sequenceCounts.get(sequence));
        }

        return ilAngosticSequenceCounts;
    }

    private Set<String> getProjectsForSequence(ICluster cluster, String sequence) {
        Set<String> projects = new HashSet<String>();

        if (sequence == null)
            return Collections.emptySet();

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

        return projects;
    }

    public class ProjectProperties {
        /**
         * These count the number of
         * clusters
         */
        private int correctIds = 0;
        private int incorrectIds = 0;
        private int correctContam = 0;
        private int incorrectContam = 0;

        /**
         * These count the number of actual spectra
         */
        private int correctPSMs = 0;
        private int incorrectPSMs = 0;
        private int correctContamPsms = 0;
        private int incorrectContamPsms = 0;

        public void incrementCorrectIds() {
            correctIds++;
        }

        public void incrementIncorrectIds() {
            incorrectIds++;
        }

        public void incrementCorrectContam() {
            correctContam++;
        }

        public void incrementIncorrectContam() {
            incorrectContam++;
        }

        public void incrementCorrectPSMs() {
            correctPSMs++;
        }

        public void incrementIncorrectPSMs() {
            incorrectPSMs++;
        }

        public void incrementCorrectContamPSMs() {
            correctContamPsms++;
        }

        public void incrementIncorrectContamPSMs() {
            incorrectContamPsms++;
        }

        public int getCorrectIds() {
            return correctIds;
        }

        public int getIncorrectIds() {
            return incorrectIds;
        }

        public int getCorrectContam() {
            return correctContam;
        }

        public int getIncorrectContam() {
            return incorrectContam;
        }

        public int getCorrectPSMs() {
            return correctPSMs;
        }

        public int getIncorrectPSMs() {
            return incorrectPSMs;
        }

        public int getCorrectContamPsms() {
            return correctContamPsms;
        }

        public int getIncorrectContamPsms() {
            return incorrectContamPsms;
        }
    }
}
