package uk.ac.ebi.pride.spectracluster.analysis.cli;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

@SuppressWarnings("static-access")
public class CliOptions {

	public enum OPTIONS {
		LIST_ANALYSERS("list_analysers"),
        ANALYSER("analyser"),
        ALL_ANALYSERS("all_analysers"),
        OUTPUT_PATH("output_path"),
        CUMULATIVE_ANALYSIS("cumulative_analysis"),
        MIN_SIZE("min_size"),
        MAX_SIZE("max_size"),
        MIN_RATIO("min_ratio"),
        MAX_RATIO("max_ratio"),
        MIN_PRECURSOR("min_precursor"),
        MAX_PRECURSOR("max_precurosr"),
        EXTRACT_SPECTA("extract_spectra"),
        HELP("help");

		private String value;

		OPTIONS(String value) {
			this.value = value;
		}
		
		public String getValue() {
			return value;
		}

		@Override
		public String toString() {
			return value;
		}
	}

	private static final Options options = new Options();

	static {
        Option cumulative = OptionBuilder
                .withDescription("if set the analysis is not performed for each passed file but for all files together. In this case output_path MUST be set.")
                .create(OPTIONS.CUMULATIVE_ANALYSIS.getValue());
        options.addOption(cumulative);

        Option analyser = OptionBuilder
                .hasArg()
                .withArgName("ANALYSER")
                .withDescription("enables the defined analyser to run on the file.")
                .create(OPTIONS.ANALYSER.getValue());
        options.addOption(analyser);

        Option allAnalysers = OptionBuilder
                .withDescription("enables all available analysers.")
                .create(OPTIONS.ALL_ANALYSERS.getValue());
        options.addOption(allAnalysers);

        Option minSize = OptionBuilder
                .withDescription("minimum cluster size.")
                .hasArg()
                .create(OPTIONS.MIN_SIZE.getValue());
        options.addOption(minSize);

        Option maxSize = OptionBuilder
                .withDescription("maximum cluster size.")
                .hasArg()
                .create(OPTIONS.MAX_SIZE.getValue());
        options.addOption(maxSize);

        Option minRatio = OptionBuilder
                .withDescription("minimum ratio a cluster may have to still be processed.")
                .hasArg()
                .create(OPTIONS.MIN_RATIO.getValue());
        options.addOption(minRatio);

        Option maxRatio = OptionBuilder
                .withDescription("maximum ratio a cluster may have to still be processed.")
                .hasArg()
                .create(OPTIONS.MAX_RATIO.getValue());
        options.addOption(maxRatio);

        Option minPrecursor = OptionBuilder
                .withDescription("minimum (average) precursor m/z a cluster may have to be processed.")
                .hasArg()
                .create(OPTIONS.MIN_PRECURSOR.getValue());
        options.addOption(minPrecursor);

        Option maxPrecursor = OptionBuilder
                .withDescription("maximum (average) precursor m/z a cluster may have to be processed.")
                .hasArg()
                .create(OPTIONS.MIN_PRECURSOR.getValue());
        options.addOption(maxPrecursor);

        Option outputPath = OptionBuilder
                .hasArg()
                .withArgName("PATH")
                .withDescription("if specified the result files will be written to this path instead of the original files directory.")
                .create(OPTIONS.OUTPUT_PATH.getValue());
        options.addOption(outputPath);

        Option listAnalysers = OptionBuilder
                .withDescription("List all available analysers.")
                .create(OPTIONS.LIST_ANALYSERS.getValue());
        options.addOption(listAnalysers);

        Option extractSpectra = OptionBuilder
                .withDescription("extract the spectra from the clusters and write them to the specified directory")
                .withArgName("directory")
                .hasArg()
                .create(OPTIONS.EXTRACT_SPECTA.getValue());
        options.addOption(extractSpectra);

		Option help = OptionBuilder
                .withDescription("print this help.")
                .create(OPTIONS.HELP.getValue());
		options.addOption(help);
	}

	public static Options getOptions() {
		return options;
	}
}
