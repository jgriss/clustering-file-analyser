package uk.ac.ebi.pride.spectracluster.analysis.util;

import java.io.*;
import java.util.*;

/**
 * Created by jg on 17.07.15.
 */
public class ModificationMapper {
    private final static ModificationMapper instance = new ModificationMapper();

    public static ModificationMapper getInstance() {
        return instance;
    }

    private final Map<String, Set<String>> mappings;

    private ModificationMapper() {
        try {
            InputStream inputStream = ModificationMapper.class.getClassLoader().getResourceAsStream("ModificationMappings.txt");

            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

            String line;
            Map<String, Set<String>> fileMappings = new HashMap<String, Set<String>>();

            while ((line = reader.readLine()) != null) {
                // remove comments
                int commentIndex = line.indexOf('#');
                if (commentIndex >= 0)
                    line = line.substring(commentIndex + 1);

                // ignore empty lines
                line = line.trim();
                if (line.length() < 1)
                    continue;

                // process the mapping
                int index = line.indexOf(" => ");
                if (index < 0)
                    continue;

                String accessions = line.substring(0, index);
                String name = line.substring(index + " => ".length());

                Set<String> accessionSet = new HashSet<String>();
                String[] splitAccessions = accessions.split(",");
                for (String a : splitAccessions) {
                    accessionSet.add(a.trim());
                }

                fileMappings.put(name, accessionSet);
            }

            this.mappings = fileMappings;
        }
        catch (Exception e) {
            throw new IllegalStateException(e);
        }
    }

    public String getMappingForAccession(String accession) {
        for (String mapping : mappings.keySet()) {
            if (mappings.get(mapping).contains(accession))
                return mapping;
        }

        return null;
    }

    public Set<String> getAccessionsForMapping(String mapping) {
        return Collections.unmodifiableSet(mappings.get(mapping));
    }

    public Set<String> getAvailableModifications() {
        return mappings.keySet();
    }
}
