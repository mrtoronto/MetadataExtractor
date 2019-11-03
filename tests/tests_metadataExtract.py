""" Unit-testing for metadata extraction. Keep < 30sec if possible. 
    Age-specific functions moved to separate test .py """

import unittest, sys, os, re, warnings, requests

sys.path.append('./src/')
import setup_metadataExtract as util

gsmURL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='

class TestDataProcess(unittest.TestCase):
    
    def test_geoSampleCellCheck(self):
        """ Check: T/F on hand-picked cell examples. Detect cell lines vs types.
            Missed cell catch with 'incorrect' protocolEntries. """
        
        trueCell, sortCell, noCell = 'GSM937915', 'GSM1402452', 'GSM1418737'
        trueCellText = requests.get(gsmURL + trueCell).text
        sortCellText = requests.get(gsmURL + sortCell).text
        noCellText = requests.get(gsmURL + noCell).text
        
        self.assertTrue(util.geoSampleCellCheck(trueCellText))
        self.assertFalse(util.geoSampleCellCheck(sortCellText))
        self.assertFalse(util.geoSampleCellCheck(noCellText))

        cellDetectAlternative = 'cell lines?\:|cell types?\:'
        self.assertTrue(util.geoSampleCellCheck(trueCellText, 
            cellDetectChar = cellDetectAlternative))
        self.assertTrue(util.geoSampleCellCheck(sortCellText, 
            cellDetectChar = cellDetectAlternative))
        self.assertFalse(util.geoSampleCellCheck(noCellText, 
            cellDetectChar = cellDetectAlternative))
        
        protocolEntryAlternatives = ['Treatment']
        self.assertFalse(util.geoSampleCellCheck(trueCellText, 
            protocolEntries = protocolEntryAlternatives))
        self.assertFalse(util.geoSampleCellCheck(sortCellText, 
            protocolEntries = protocolEntryAlternatives))
        self.assertFalse(util.geoSampleCellCheck(noCellText, 
            protocolEntries = protocolEntryAlternatives))
       

    def test_extractGEOSampleInfo(self):
        """ Check: Starting from the sample ID, proper extractions (including gender) """

        samp1, samp2, samp3 = 'GSM400641', 'GSM937915', 'GSM1402452'
        samp4, samp5, samp6 = 'GSM1418737', 'GSM1282831', 'GSM1338039'
        fakeSamp = 'GSM9999999'

        with self.assertRaises(ValueError):
            util.extractGEOSampleInfo(fakeSamp)
        
        with self.assertRaises(ValueError):
            util.extractGEOSampleInfo(samp1, organism = 'musculo musculus')

        with self.assertRaises(AttributeError): 
            util.extractGEOSampleInfo(samp2, keepCells = False)

        samp1Dict = util.extractGEOSampleInfo(samp1)
        samp2Dict = util.extractGEOSampleInfo(samp2)
        samp3Dict = util.extractGEOSampleInfo(samp3, tryAgePMID = False)
        samp4Dict = util.extractGEOSampleInfo(samp4)
        samp5Dict = util.extractGEOSampleInfo(samp5)
        samp6Dict = util.extractGEOSampleInfo(samp6)
        
        self.assertAlmostEqual(4.28571428, samp1Dict['Age'])
        self.assertEqual('n/a', samp1Dict['Gender'])
        self.assertEqual('GSE16012', samp1Dict['Study'])
        
        self.assertEqual('n/a', samp2Dict['Age'])
        self.assertEqual('n/a', samp2Dict['Gender'])
        self.assertEqual(True, samp2Dict['Cells'])
        
        self.assertEqual('n/a', samp3Dict['Age'])
        self.assertEqual('n/a', samp3Dict['Gender'])
        self.assertEqual(False, samp3Dict['Cells'])
        self.assertEqual(samp3Dict['Flags']['Sort'], True)
        samp3Dict = util.extractGEOSampleInfo(samp3, tryAgePMID = False, flagSort = False)
        self.assertEqual(samp3Dict['Flags']['Sort'], False)
        
        self.assertEqual(13.5, samp4Dict['Age'])
        self.assertEqual('Female', samp4Dict['Gender'])
        self.assertEqual(False, samp4Dict['Cells'])

        self.assertEqual(16, samp5Dict['Age'])
        self.assertEqual('Male', samp5Dict['Gender'])
        self.assertEqual(False, samp5Dict['Cells'])

        self.assertAlmostEqual(9.42857142, samp6Dict['Age'])
        self.assertEqual('Female', samp6Dict['Gender'])
        self.assertEqual(False, samp6Dict['Cells'])


    def test_maxSectionMatch(self):
        """ Check: Proper url and text which contains section of interest 
            (e.g. methods). Altered return with different soupAttempts arg.
            nullReturn with either mismatched soupAttemps or reduced 
            possibleSections. 
            
            ** Add findSectionText ** 
            Check: Ensure example text found for multiple section tags. Return
                null with misspecification of soupDiv or invalid section tags.
                Mixed section tags will return first section type match"""
    
        links15452052 = ['https://iovs.arvojournals.org/article.aspx?articleid=2124532']

        links25043027 = ['https://www.nature.com/articles/nature13577', 
                        'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4339042/']
                        # --> PMC should win

        links26717410 = ['https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002330',
                        'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696735']
                        # --> plos wins unless PMC sections added

        links19286929 = ['https://www.physiology.org/doi/full/10.1152/ajplung.90369.2008?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed',
                        'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2692799/']
                        # --> PMC should with win 'h2', other with 'h1', 'n/a' with 'h6'

        links21996730 = ['https://www.nature.com/articles/onc2011433',
                        'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3356601/']
                        # --> nature has unique sections, should win with ['Ethics', 'Accession codes']

        link1, text1, div1 = util.maxSectionMatch(links15452052)
        link2, text2, div2 = util.maxSectionMatch(links25043027)
        link3, text3, div3 = util.maxSectionMatch(links26717410)
        link4, text4, div4 = util.maxSectionMatch(links19286929)
        self.assertIn('iovs', link1)
        self.assertEqual('h6', div1)
        self.assertIn('pmc', link2)
        self.assertEqual('h2', div2)
        self.assertIn('plos', link3)
        self.assertEqual('h2', div3)
        self.assertIn('pmc', link4)
        self.assertEqual('h2', div4)

        secText1 = util.findSectionText(text1, sectionTags = ['methods', 
                                'procedures'], soupDiv = div1)
        self.assertIn('according to protocols', secText1.lower())
        self.assertIn('cells plated on coverslips', secText1.lower())

        secText1 = util.findSectionText(text1, sectionTags = ['methods', 
                                'results'], soupDiv = div1)
        self.assertIn('according to protocols', secText1.lower())
        self.assertNotIn('gene expression patterns were', secText1.lower())

        secText1 = util.findSectionText(text1, sectionTags = ['results'], 
                                        soupDiv = div1)
        self.assertIn('gene expression patterns were', secText1.lower())
        self.assertIn('several corneal maintenance', secText1.lower())

        secText1 = util.findSectionText(text1, sectionTags = ['methods', 
                                'procedures'], soupDiv = 'h2')
        self.assertEqual('n/a', secText1.lower())

        secText1 = util.findSectionText(text1, sectionTags = ['xyz'], 
                                        soupDiv = 'h2')
        self.assertEqual('n/a', secText1.lower())


        link3, _, _ = util.maxSectionMatch(links26717410, 
                        possibleSections = ['references', 'data availability', 
                            'funding statement', 'abbreviations', 
                            'materials and methods', 'results', 'introduction', 
                            'discussion'])

        self.assertIn('pmc', link3)
        link4, _, div4 = util.maxSectionMatch(links19286929, soupAttempts = ['h1'])
        self.assertIn('physiology', link4)
        self.assertEqual('h1', div4)
        link4, _, div4 = util.maxSectionMatch(links19286929, soupAttempts = ['h6'])
        self.assertEqual('n/a', link4)
        self.assertEqual('n/a', div4)

        link5, text5, div5 = util.maxSectionMatch(links21996730, 
                            possibleSections = ['ethics', 'accession codes'])
                            
        self.assertIn('nature', link5)
        self.assertIn('ethics', text5)
        self.assertEqual('h2', div5)
        

    def text_pmidSectionExtraction(self):
        """ Check: Valid text returns from same paper with multiple sectionIDs.
            pmcPreference argument. Null return from paper wih no links on pubmed
            page or without access to full paper. """
    
        pmid1 = '15452052'
        pmid2 = '26717410'
        
        url1 = 'https://www.ncbi.nlm.nih.gov/pubmed/{0}'.format(pmid1)
        pubText1 = requests.get(url1).text
        url2 = 'https://www.ncbi.nlm.nih.gov/pubmed/{0}'.format(pmid2)
        pubText2 = requests.get(url2).text

        meth1 = util.pmidSectionExtraction(pmidText = pubText1, sectionID = 'methods')
        res1 = util.pmidSectionExtraction(pmidText = pubText1, sectionID = 'results')
        self.assertIn('according to protocols', meth1.lower())
        self.assertIn('gene expression patterns were', res1.lower())

        meth2_a = util.pmidSectionExtraction(pmidText = pubText2, 
                                    sectionID = 'methods', pmcPreference = True)
        meth2_b = util.pmidSectionExtraction(pmidText = pubText2, 
                                    sectionID = 'methods', pmcPreference = False)

        self.assertNotEqual(meth2_a, meth2_b)
        self.assertIn('mice were mantained following', meth2_a.lower())
        self.assertIn('mice were mantained following', meth2_b.lower())
        
        null1 = '31341269'
        nullURL1 = 'https://www.ncbi.nlm.nih.gov/pubmed/{0}'.format(null1)
        nullText1 = requests.get(nullURL1).text

        self.assertEqual('n/a', util.pmidSectionExtraction(pmidText = nullText1, 
                        sectionID = 'methods'))

        null2 = '13054692'
        nullURL2 = 'https://www.ncbi.nlm.nih.gov/pubmed/{0}'.format(null2)
        nullText2 = requests.get(nullURL2).text

        self.assertEqual('n/a', util.pmidSectionExtraction(pmidText = nullText2, 
                        sectionID = 'methods'))
        

if __name__ == '__main__':

    unittest.main()
