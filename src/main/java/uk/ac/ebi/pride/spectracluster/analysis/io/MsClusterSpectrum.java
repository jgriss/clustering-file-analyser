package uk.ac.ebi.pride.spectracluster.analysis.io;

/**
 * Created by jg on 01.05.15.
 */
public class MsClusterSpectrum {
    private final int datasetIdx;
    private final int fileIdx;
    private final int scanNumber;
    private final float precursorMz;
    private final float similarity;
    private final float pValue;
    private final int charge;
    private final String sequence;

    public MsClusterSpectrum(int datasetIdx, int fileIdx, int scanNumber, float precursorMz, float similarity, float pValue, int charge, String sequence) {
        this.datasetIdx = datasetIdx;
        this.fileIdx = fileIdx;
        this.scanNumber = scanNumber;
        this.precursorMz = precursorMz;
        this.similarity = similarity;
        this.pValue = pValue;
        this.charge = charge;
        this.sequence = sequence;
    }

    public int getDatasetIdx() {
        return datasetIdx;
    }

    public int getFileIdx() {
        return fileIdx;
    }

    public int getScanNumber() {
        return scanNumber;
    }

    public float getPrecursorMz() {
        return precursorMz;
    }

    public float getSimilarity() {
        return similarity;
    }

    public float getpValue() {
        return pValue;
    }

    public int getCharge() {
        return charge;
    }

    public String getSequence() {
        return sequence;
    }
}

