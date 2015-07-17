package uk.ac.ebi.pride.spectracluster.analysis.util;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by jg on 17.07.15.
 */
public class TestModificationMapper {
    @Test
    public void testModificationMapper() throws Exception {
        ModificationMapper mapper = ModificationMapper.getInstance();

        Assert.assertEquals("Carboxylation", mapper.getMappingForAccession("MOD:00123"));
        Assert.assertEquals(5, mapper.getAccessionsForMapping("Sulfo").size());
        Assert.assertTrue(mapper.getAccessionsForMapping("Sulfo").contains("MOD:00695"));
    }
}
