package uk.ac.ebi.pride.spectracluster.analysis.io;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ClusteringFileSpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.SequenceCount;

import java.util.*;

/**
 * Created by jg on 01.05.15.
 */
public class MsClusterCluster implements ICluster {
    private final String clusterId;
    private final int clusterSize;
    private final float precursorMz;
    private final int charge;
    private final List<MsClusterSpectrum> spectra;

    private Set<String> sequences;
    private List<ISpectrumReference> spectrumReferences;
    private List<SequenceCount> sequenceCounts;

    private float maxRatio = -1;
    private String maxSequence;
    private float mzRange = -1;

    public MsClusterCluster(String clusterId, int clusterSize, float precursorMz, int charge, List<MsClusterSpectrum> spectra) {
        this.clusterId = clusterId;
        this.clusterSize = clusterSize;
        this.precursorMz = precursorMz;
        this.charge = charge;
        this.spectra = spectra;
    }

    @Override
    public float getAvPrecursorMz() {
        return 0;
    }

    @Override
    public float getAvPrecursorIntens() {
        return 0;
    }

    @Override
    public Set<String> getSequences() {
        if (sequences == null) {
            sequences = new HashSet<String>();

            for (MsClusterSpectrum s : spectra) {
                if (s.getSequence() == null)
                    continue;

                sequences.add(s.getSequence());
            }
        }

        return sequences;
    }

    @Override
    public List<ISpectrumReference> getSpectrumReferences() {
        if (spectrumReferences == null) {
            spectrumReferences = new ArrayList<ISpectrumReference>(spectra.size());

            for (MsClusterSpectrum s : spectra) {
                spectrumReferences.add(new ClusteringFileSpectrumReference(
                        s.getSequence(), s.getCharge(), s.getPrecursorMz(), s.getDatasetIdx() + "_" +
                 s.getFileIdx() + "_" + s.getScanNumber(), s.getSimilarity(), "", ""));
            }
        }

        return spectrumReferences;
    }

    @Override
    public int getSpecCount() {
        return spectra.size();
    }

    @Override
    public int getPsmCount() {
        return spectra.size();
    }

    @Override
    public List<SequenceCount> getSequenceCounts() {
        if (sequenceCounts == null) {
            sequenceCounts = new ArrayList<SequenceCount>();
            Map<String, Integer> sequenceMap = new HashMap<String, Integer>();

            for (MsClusterSpectrum s : spectra) {
                if (s.getSequence().length() < 1)
                    continue;

                if (!sequenceMap.containsKey(s.getSequence()))
                    sequenceMap.put(s.getSequence(), 1);
                else
                    sequenceMap.put(s.getSequence(), sequenceMap.get(s.getSequence()) + 1);
            }

            for (String sequence : sequenceMap.keySet()) {
                sequenceCounts.add(new SequenceCount(sequence, sequenceMap.get(sequence)));
            }
        }

        return sequenceCounts;
    }

    @Override
    public Map<String, Integer> getPsmSequenceCounts() {
        Map<String, Integer> sequenceMap = new HashMap<String, Integer>();

        for (MsClusterSpectrum s : spectra) {
            if (s.getSequence().length() < 1)
                continue;

            if (!sequenceMap.containsKey(s.getSequence()))
                sequenceMap.put(s.getSequence(), 1);
            else
                sequenceMap.put(s.getSequence(), sequenceMap.get(s.getSequence()) + 1);
        }

        return sequenceMap;
    }

    @Override
    public float getMaxRatio() {
        if (maxRatio < 0) {
            List<SequenceCount> sequenceCounts1 = getSequenceCounts();

            int maxCount = 0;

            for (SequenceCount sq : sequenceCounts1) {
                if (sq.getCount() > maxCount)
                    maxCount = sq.getCount();
            }

            maxRatio = (float) maxCount / spectra.size();
        }

        return maxRatio;
    }

    @Override
    public String getMaxSequence() {
        if (maxSequence == null) {
            List<SequenceCount> sequenceCounts1 = getSequenceCounts();

            int maxCount = 0;

            for (SequenceCount sq : sequenceCounts1) {
                if (sq.getCount() > maxCount) {
                    maxCount = sq.getCount();
                    maxSequence = sq.getSequence();
                }
            }


        }

        return maxSequence;
    }

    @Override
    public float getSpectrumPrecursorMzRange() {
        if (mzRange < 0) {
            float maxMz = Float.MIN_VALUE;
            float minMz = Float.MAX_VALUE;

            for (MsClusterSpectrum s : spectra) {
                if (s.getPrecursorMz() > maxMz)
                    maxMz = s.getPrecursorMz();

                if (s.getPrecursorMz() < minMz)
                    minMz = s.getPrecursorMz();
            }

            mzRange = maxMz - minMz;
        }

        return mzRange;
    }

    @Override
    public List<Float> getConsensusMzValues() {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<Float> getConsensusIntensValues() {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getId() {
        return clusterId;
    }

    @Override
    public Set<String> getSpecies() {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getIdentifiedSpecCount() {
        return spectra.size();
    }

    @Override
    public int getUnidentifiedSpecCount() {
        return 0;
    }
}