package uk.ac.ebi.pride.spectracluster.analysis.analyser;

/**
 * Created by jg on 15.07.14.
 */
public class AnalyserFactory {
    private AnalyserFactory() {

    }

    public static enum ANALYSERS {
        BASIC_CLUSTERING_STATISTICS("BasicClusteringAnalyzer"),
        CLUSTER_DUPLICATION_ANALYSER("ClusterDuplicationAnalyser"),
        CLUSTER_PARAMETER_EXTRACTOR("ClusterParameterExtractor"),
        EXTENDED_PARANETER_EXTRACTOR("ExtendedParameters"),
        ID_CLUSTER_PARAMETER_EXTRACTOR("IdClusterParameterExtractor"),
        CLUSTERED_SPECTRA_ANALYSER("ClusteredSpectra"),
        PROJECT_ANALYSER("ProjectAnalyser"),
        CONSENSUS_SPECTRUM("ConsensusSpectrum");

        private String name;
        private ANALYSERS(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public static ANALYSERS getAnalyserForString(String name) {
            for (ANALYSERS a : values()) {
                if (a.getName().trim().toLowerCase().equals(name.trim().toLowerCase()))
                    return a;
            }

            return null;
        }
    }

    public static AbstractClusteringSourceAnalyser getAnalyserForString(String name) {
        ANALYSERS analyser = ANALYSERS.getAnalyserForString(name);

        return getAnalyser(analyser);
    }

    public static AbstractClusteringSourceAnalyser getAnalyser(ANALYSERS analyser) {
        switch (analyser) {
            case BASIC_CLUSTERING_STATISTICS:
                return new BasicClusteringStatistics();
            case CLUSTER_DUPLICATION_ANALYSER:
                return new ClusterDuplicationAnalyser();
            case CLUSTER_PARAMETER_EXTRACTOR:
                return new ClusterParameterExtractor();
            case EXTENDED_PARANETER_EXTRACTOR:
                return new ExtendedClusterParameterExtractor();
            case ID_CLUSTER_PARAMETER_EXTRACTOR:
                return new IdentifiedClusterParameterExtractor();
            case CLUSTERED_SPECTRA_ANALYSER:
                return new ClusteredSpectraAnalyser();
            case PROJECT_ANALYSER:
                return new ProjectAnalyser();
            case CONSENSUS_SPECTRUM:
                return new ConsensusSpectrumExtractor();
            default:
                return null;
        }
    }
}
