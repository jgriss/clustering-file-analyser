package uk.ac.ebi.pride.spectracluster.analysis.io;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.IClusterSourceReader;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Reads clustering result files created by the MsCluster algorithm /
 * spectral achives algorithm.
 * Created by jg on 01.05.15.
 */
public class MsClusterFileReader implements IClusterSourceReader {
    private final File inputFile;
    private BufferedReader reader;
    private String currentLine;

    public MsClusterFileReader(File inputFile) throws FileNotFoundException {
        this.inputFile = inputFile;

        this.reader = new BufferedReader(new FileReader(this.inputFile));
    }

    @Override
    public List<ICluster> readAllClusters() throws Exception {
        List<ICluster> clusters = new ArrayList<ICluster>();

        ICluster cluster;

        while ((cluster = readNextCluster()) != null) {
            clusters.add(cluster);
        }

        return clusters;
    }

    @Override
    public boolean supportsReadAllClusters() {
        return true;
    }

    @Override
    public void readClustersIteratively(Collection<IClusterSourceListener> listeners) throws Exception {
        ICluster cluster;

        while ((cluster = readNextCluster()) != null) {
            for (IClusterSourceListener listener : listeners) {
                listener.onNewClusterRead(cluster);
            }
        }
    }

    private MsClusterCluster readNextCluster() throws Exception {
        if (currentLine == null) {
            currentLine = reader.readLine();
        }

        if (currentLine == null) {
            return null;
        }

        // get the cluster's details
        String[] clusterHeaderFields = currentLine.split("\t");

        if (clusterHeaderFields.length != 4 && clusterHeaderFields.length != 5) {
            throw new Exception("Cluster header line contains less than 4 and 5 fields: " + currentLine);
        }

        List<MsClusterSpectrum> currentSpectra = new ArrayList<MsClusterSpectrum>();

        while ((currentLine = reader.readLine()) != null) {
            if (currentLine.length() < 1)
                continue;

            String[] fields = currentLine.split("\t");

            // check if it's a cluster
            if (fields.length < 7) {
                break;
            }

            MsClusterSpectrum spectrum = new MsClusterSpectrum(
                    Integer.parseInt(fields[0]), Integer.parseInt(fields[1]), Integer.parseInt(fields[2]),
                    Float.parseFloat(fields[3]), Float.parseFloat(fields[4]), Float.parseFloat(fields[5]),
                    Integer.parseInt(fields[6]), fields.length > 7 ? fields[7] : ""
            );

            currentSpectra.add(spectrum);
        }

        // create the cluster
        return new MsClusterCluster(clusterHeaderFields[0], Integer.parseInt(clusterHeaderFields[1]),
                Float.parseFloat(clusterHeaderFields[2]), Integer.parseInt(clusterHeaderFields[3]), currentSpectra);
    }
}
