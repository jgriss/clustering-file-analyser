package uk.ac.ebi.pride.spectracluster.analysis.util;

import com.sun.org.apache.bcel.internal.generic.ILOAD;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.IPeptideSpectrumMatch;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.SequenceCount;

import java.util.*;

/**
 * Created by jg on 14.07.14.
 */
public class ClusterUtilities {
    private ICluster currentCluster;

    private String maxSequence;
    private float maxILAngosticRatio;
    private int maxSequenceCount;
    private int nProjects;
    private int nAssays;
    private double mzRange;
    private Set<String> species;
    private Map<String, Integer> sequenceCounts;
    private int charge = 0;

    private String secondMaxSequence;
    private int secondMaxSequenceCount;

    private String thirdMaxSequence;
    private int thirdMaxSequenceCount;

    public ClusterUtilities() {

    }

    public ClusterUtilities(ICluster cluster) {
        processCluster(cluster);
    }

    /**
     * Processed the passed cluster. This call overwrites
     * any results that may have been created before.
     * @param cluster
     */
    public void processCluster(ICluster cluster) {
        // this is only saved for potential future reference
        currentCluster = cluster;

        updateNumberOfProjects(cluster);

        // update the sequences
        updateSequences(cluster);
        updatePrecursorMzRange(cluster);
        updateSpecies(cluster);
        charge = calculateCharge(cluster);
    }

    private void updateSequences(ICluster cluster) {
        // update the most common sequence
        List<Object> maxSequenceProperties = getMaxSequence(cluster, Collections.EMPTY_SET);

        this.maxSequence = (String) maxSequenceProperties.get(0);
        this.maxSequenceCount = (Integer) maxSequenceProperties.get(1);
        this.maxILAngosticRatio = (float) this.maxSequenceCount / cluster.getSpectrumReferences().size();
        this.sequenceCounts = createSequenceCounts(cluster);

        // get the second most common sequence
        String maxIlAgnosticSequence = this.maxSequence.replaceAll("I", "L");
        Set<String> knownSequence = new HashSet<String>();
        knownSequence.add(maxIlAgnosticSequence);

        maxSequenceProperties = getMaxSequence(cluster, knownSequence);

        this.secondMaxSequence = (String) maxSequenceProperties.get(0);
        this.secondMaxSequenceCount = (Integer) maxSequenceProperties.get(1);

        // get the third max sequence
        if (secondMaxSequence != null) {
            knownSequence.add(secondMaxSequence.replace("I", "L"));
            maxSequenceProperties = getMaxSequence(cluster, knownSequence);

            this.thirdMaxSequence = (String) maxSequenceProperties.get(0);
            this.thirdMaxSequenceCount = (Integer) maxSequenceProperties.get(1);
        }
        else {
            this.thirdMaxSequence = null;
            this.thirdMaxSequenceCount = 0;
        }
    }

    private int calculateCharge(ICluster cluster) {
        // calculate the average charge
        int sumCharge = 0;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            sumCharge += specRef.getCharge();
        }

        float avCharge = (float) sumCharge / (float) cluster.getSpectrumReferences().size();
        int avChargeRounded = (int) (avCharge + 0.5);

        return avChargeRounded;
    }

    private void updateSpecies(ICluster cluster) {
        species = new HashSet<String>();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            species.add(specRef.getSpecies());
        }
    }

    public static String cleanSequence(String sequence) {
        if (sequence == null) {
            return null;
        }

        return sequence.toUpperCase().replaceAll("[^A-Z]", "");
    }

    private Map<String, Integer> createSequenceCounts(ICluster cluster) {
        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();

        for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
            IPeptideSpectrumMatch peptideSpectrumMatch = spectrumReference.getMostCommonPSM();
            String cleanSequence = cleanSequence(peptideSpectrumMatch.getSequence());

            if (!sequenceCounts.containsKey(cleanSequence))
                sequenceCounts.put(cleanSequence, 1);
            else
                sequenceCounts.put(cleanSequence, sequenceCounts.get(cleanSequence) + 1);
        }

        return sequenceCounts;
    }

    /**
     * Retruns the sequence and the max count when ignoring the sequences set in knownMaxSequence
     * @param cluster
     * @param knownMaxSequence
     * @return
     */
    private List<Object> getMaxSequence(ICluster cluster, Set<String> knownMaxSequence) {
        Map<String, Integer> sequenceCounts = new HashMap<String, Integer>();
        Map<String, String> ilCorrectedToOriginalSequence = new HashMap<String, String>();

        for (ISpectrumReference spectrumReference : cluster.getSpectrumReferences()) {
            IPeptideSpectrumMatch peptideSpectrumMatch = spectrumReference.getMostCommonPSM();
            String ilAgnosticSequence = peptideSpectrumMatch.getSequence().replaceAll("I", "L");
            ilAgnosticSequence = cleanSequence(ilAgnosticSequence);

            // ignore previously processed sequences
            if (knownMaxSequence.contains(ilAgnosticSequence))
                continue;

            if (!ilCorrectedToOriginalSequence.containsKey(ilAgnosticSequence)) {
                ilCorrectedToOriginalSequence.put(ilAgnosticSequence, peptideSpectrumMatch.getSequence());
            }

            if (!sequenceCounts.containsKey(ilAgnosticSequence)) {
                sequenceCounts.put(ilAgnosticSequence, 0);
            }

            sequenceCounts.put(ilAgnosticSequence, sequenceCounts.get(ilAgnosticSequence) + 1);
        }

        // get the max count
        int maxCount = 0;
        String maxSequence = null;

        for (String ilAgnosticSequence : sequenceCounts.keySet()) {
            int count = sequenceCounts.get(ilAgnosticSequence);

            if (count > maxCount) {
                maxSequence = ilCorrectedToOriginalSequence.get(ilAgnosticSequence);
                maxCount = count;
            }
        }

        List<Object> result = new ArrayList<Object>(2);
        result.add(maxSequence);
        result.add(new Integer(maxCount));

        return result;
    }

    private void updateNumberOfProjects(ICluster cluster) {
        Set<String> projects = new HashSet<String>();
        Set<String> assays = new HashSet<String>();

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            String id = specRef.getSpectrumId();

            String[] fields = id.split(";");

            // PXD000807;index=2376,splib_sequence=KYLYEIAR,score=0.639,peptideR2=,scoreR2=
            projects.add(fields[0]);
            assays.add((fields.length >= 3) ? fields[1] : "Unknown");
        }

        this.nProjects = projects.size();
        this.nAssays = assays.size();
    }

    private void updatePrecursorMzRange(ICluster cluster) {
        double minMZ = Float.MAX_VALUE, maxMz = 0;

        for (ISpectrumReference specRef : cluster.getSpectrumReferences()) {
            double mz = specRef.getPrecursorMz();

            if (mz < minMZ) {
                minMZ = mz;
            }
            if (mz > maxMz) {
                maxMz = mz;
            }
        }

        this.mzRange = maxMz - minMZ;
    }

    public ICluster getCurrentCluster() {
        return currentCluster;
    }

    public String getMaxSequence() {
        return maxSequence;
    }

    public float getMaxILAngosticRatio() {
        return maxILAngosticRatio;
    }

    public int getMaxSequenceCount() {
        return maxSequenceCount;
    }

    public int getnProjects() {
        return nProjects;
    }

    public int getnAssays() {
        return nAssays;
    }

    public double getMzRange() {
        return mzRange;
    }

    public Set<String> getSpecies() {
        return Collections.unmodifiableSet(species);
    }

    public String getSecondMaxSequence() {
        return secondMaxSequence;
    }

    public int getSecondMaxSequenceCount() {
        return secondMaxSequenceCount;
    }

    public Map<String, Integer> getSequenceCounts() {
        return Collections.unmodifiableMap(sequenceCounts);
    }

    public String getThirdMaxSequence() {
        return thirdMaxSequence;
    }

    public int getThirdMaxSequenceCount() {
        return thirdMaxSequenceCount;
    }

    public boolean isStable() {
        if (currentCluster.getSpectrumReferences().size() >= 10 & maxILAngosticRatio > 0.7)
            return true;

        return false;
    }

    public int getCharge() {
        return charge;
    }
}
