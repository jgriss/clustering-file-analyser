package uk.ac.ebi.pride.spectracluster.analysis.util;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by jg on 23.06.15.
 */
public class TestCrapFastaFile {
    @Test
    public void testIdentifyPeptide() throws Exception  {
        CrapFastaFile crapFastaFile = CrapFastaFile.getInstance();

        Assert.assertNull(crapFastaFile.getProteinAnnotation("ABC"));
        Assert.assertEquals("sp|AMYS_HUMAN|", crapFastaFile.getProteinAnnotation("GRVTEFKYGAKLGTVIRKWNGEKMSYLK"));
    }
}
