package uk.ac.ebi.pride.spectracluster.analysis.util;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by jg on 21.05.15.
 */
public final class CrapFastaFile {
    private Map<String, String> proteinSequences;

    private static CrapFastaFile instance = new CrapFastaFile();

    public static CrapFastaFile getInstance() {
        return instance;
    }

    private CrapFastaFile() {
        InputStream inputStream = CrapFastaFile.class.getClassLoader().getResourceAsStream("crap.fasta");
        try {
            this.proteinSequences = loadFastaFile(new BufferedReader(new InputStreamReader(inputStream)));
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    protected Map<String, String> loadFastaFile(BufferedReader reader) throws Exception {
        String line;
        String currentAnnotation = "UNKNOWN";
        StringBuilder currentSequence = new StringBuilder();
        Map<String, String> sequences = new HashMap<String, String>();

        while ((line = reader.readLine()) != null) {
            line = line.trim();

            // ignore empty lines
            if (line.length() < 1)
                continue;

            // ignore any comment line and store the previous protein sequence
            if (line.startsWith(">")) {
                if (currentSequence.length() > 0) {
                    sequences.put(currentAnnotation, currentSequence.toString());
                    currentSequence = new StringBuilder();
                }

                currentAnnotation = line.substring(1); // just ignore the ">"

                continue;
            }

            currentSequence.append(line);
        }

        if (currentSequence.length() > 0) {
            sequences.put(currentAnnotation, currentSequence.toString());
        }

        return sequences;
    }

    /**
     * Retruns the protein annotation for the given peptide or null in case
     * the peptide cannot be mapped to any protein.
     * @param peptideSequence
     * @return
     */
    public String getProteinAnnotation(String peptideSequence) {
        for (String annotation : proteinSequences.keySet()) {
            String proteinSequence = proteinSequences.get(annotation);

            if (proteinSequence.contains(peptideSequence))
                return annotation;
        }

        return null;
    }
}
