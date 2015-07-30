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
        MIN_IDENTIFIED_SPEC("min_identified_spectra"),
        MAX_IDENTIFIED_SPEC("max_identified_spectra"),
        MIN_UNIDENTIFIED_SPEC("min_unidentified_spectra"),
        MAX_UNIFENTIFIED_SPEC("max_unidentified_spectra"),
        MIN_PRECURSOR("min_precursor"),
        MAX_PRECURSOR("max_precurosr"),
        MODIFICATION("modification"),
        EXTRACT_SPECTA("extract_spectra"),
        LIST_MODIFICATIONS("list_modifications"),
        SPECIES("species"),
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

        Option minId = OptionBuilder
                .withDescription("minimum number of identified spectra")
                .hasArg()
                .create(OPTIONS.MIN_IDENTIFIED_SPEC.getValue());
        options.addOption(minId);

        Option maxId = OptionBuilder
                .withDescription("maximum number of identified spectra")
                .hasArg()
                .create(OPTIONS.MAX_IDENTIFIED_SPEC.getValue());
        options.addOption(maxId);

        Option minUnid = OptionBuilder
                .withDescription("minimum number of unidentified spectra")
                .hasArg()
                .create(OPTIONS.MIN_UNIDENTIFIED_SPEC.getValue());
        options.addOption(minUnid);

        Option maxUnid = OptionBuilder
                .withDescription("maximum number of unidentified spectra")
                .hasArg()
                .create(OPTIONS.MAX_UNIFENTIFIED_SPEC.getValue());
        options.addOption(maxUnid);

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

        Option modification = OptionBuilder
                .withArgName("modification name")
                .hasArg()
                .withDescription("only reports clusters where the most common sequence was identified with any of the specified modifications. This parameter may be set multiple times.")
                .create(OPTIONS.MODIFICATION.getValue());
        options.addOption(modification);

        Option listModifications = OptionBuilder
                .withDescription("lists the available modifications.")
                .create(OPTIONS.LIST_MODIFICATIONS.getValue());
        options.addOption(listModifications);

        Option species = OptionBuilder
               .withArgName("taxid")
               .hasArg()
                .withDescription("if set only cluster where the most common sequence was identified for the specied species are reported. The parameter may be specified multiple times.")
                .create(OPTIONS.SPECIES.getValue());
        options.addOption(species);

		Option help = OptionBuilder
                .withDescription("print this help.")
                .create(OPTIONS.HELP.getValue());
		options.addOption(help);
	}

	public static Options getOptions() {
		return options;
	}
}
