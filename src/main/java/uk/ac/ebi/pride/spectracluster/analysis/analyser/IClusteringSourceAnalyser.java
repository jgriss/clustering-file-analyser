package uk.ac.ebi.pride.spectracluster.analysis.analyser;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;

/**
 * Created by jg on 12.07.14.
 */
public interface IClusteringSourceAnalyser extends IClusterSourceListener {
    public void reset();

    public String getFileEnding();

    public String getDescription();

    /**
     * Writes the file ending to the result
     * file.
     */
    public void completeResultFile() throws Exception ;
}
