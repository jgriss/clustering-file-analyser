package uk.ac.ebi.pride.spectracluster.analysis.util;

import org.apache.commons.io.FileUtils;
import uk.ac.ebi.pride.spectracluster.analysis.analyser.AbstractClusteringSourceAnalyser;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ClusteringFileSpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.PeakMzComparator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class extracts the spectra of a cluster and writes them to
 * a mgf file
 * Created by jg on 05.01.15.
 */
public class SpectrumExtractor extends AbstractClusteringSourceAnalyser {
    private final File spectraOutputDirectory;
    private static int nEmptySpectra = 0;
    private final static PeakMzComparator peakMzComparator = new PeakMzComparator();

    @Override
    protected String getResultFileHeader() {
        return "";
    }

    @Override
    public void completeResultFile() throws Exception {

    }

    @Override
    public void reset() {
        nEmptySpectra = 0;
    }

    @Override
    public String getFileEnding() {
        return ".mgf";
    }

    @Override
    public String getDescription() {
        return "Extracts all spectra for each cluster and writes them to a MGF file.";
    }

    public SpectrumExtractor(String spectraOutputDirectory) throws Exception {
        this.spectraOutputDirectory = new File(spectraOutputDirectory);

        if (!this.spectraOutputDirectory.exists())
            throw new Exception(spectraOutputDirectory + " does not exist");
        if (!this.spectraOutputDirectory.canWrite())
            throw new Exception("Cannot write to " + spectraOutputDirectory);
    }

    @Override
    protected void processClusterInternally(ICluster newCluster) throws Exception {
        try {
            // enable cluster filtering
            if (ignoreCluster(newCluster))
                return;

            String id = newCluster.getId();

            // create the file
            File tmpFile = File.createTempFile(id, ".mgf.tmp", spectraOutputDirectory);
            tmpFile.deleteOnExit();
            BufferedWriter writer = new BufferedWriter(new FileWriter(tmpFile));

            for (ISpectrumReference specRef : newCluster.getSpectrumReferences()) {
                // make sure the spectrum contains peaks
                if (!specRef.hasPeaks()) {
                    nEmptySpectra++;
                    continue;
                }

                // write the spectrum
                writer.write("BEGIN IONS\n");
                writer.write("TITLE=id=" + specRef.getSpectrumId() + ",sequence=" + specRef.getMostCommonPSM().getSequence() + ",taxid=" + specRef.getSpecies() + "\n");
                writer.write("PEPMASS=" + String.format("%.4f", specRef.getPrecursorMz()) + "\n");
                writer.write("CHARGE=" + specRef.getCharge() + "\n");

                // write the peaks
                List<ClusteringFileSpectrumReference.Peak> peaks = specRef.getPeaks();
                if (peaks == null)
                    throw new IllegalStateException("Spectrum without peaks encountered.");

                List<ClusteringFileSpectrumReference.Peak> sortedPeaks = new ArrayList<ClusteringFileSpectrumReference.Peak>(peaks);
                Collections.sort(sortedPeaks, peakMzComparator);

                for (ClusteringFileSpectrumReference.Peak peak : sortedPeaks) {
                    writer.write(String.format("%.3f %.3f\n", peak.getMz(), peak.getIntensity()));
                }

                writer.write("END IONS\n\n");
            }

            // done
            writer.close();

            // move the temporary file to the final one
            File finalFile = new File(spectraOutputDirectory, newCluster.getId() + ".mgf");
            FileUtils.copyFile(tmpFile, finalFile);
            tmpFile.delete();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
