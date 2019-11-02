""" Unit-testing for age-related functions. Keep < 30sec if possible """

import unittest, sys, os, re, warnings, requests

sys.path.append('./src/')
import setup_metadataExtract as util

gsmURL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='

class TestDataProcess(unittest.TestCase):

    def test_addAgeStrings(self):
        """ Check: Convert numbers, use hand-written test cases. """
        
        case1 = '10 weeks old'
        case2 = '10 weeks for 5 weeks'
        case3 = '6-10 weeks for 7 days'
        case4 = '6-8 days for 7 days'
        nullCase = '6-10 eons for 7 iotas'

        self.assertEqual(10, util.addAgeStrings(0, case1, convertFrom = 'week', 
            convertTo = 'week'))
        self.assertEqual(15, util.addAgeStrings(0, case2, convertFrom = 'week', 
            convertTo = 'week'))
        self.assertEqual(8, util.addAgeStrings(0, case3, convertFrom = 'week', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case4, convertFrom = 'week', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case1, convertFrom = 'day', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case2, convertFrom = 'day', 
            convertTo = 'week'))
        self.assertEqual(1, util.addAgeStrings(0, case3, convertFrom = 'day', 
            convertTo = 'week'))
        self.assertEqual(14, util.addAgeStrings(0, case4, convertFrom = 'day', 
            convertTo = 'day'))
        self.assertEqual(0, util.addAgeStrings(0, case1, convertFrom = 'month', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case2, convertFrom = 'month', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case3, convertFrom = 'month', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, case4, convertFrom = 'month', 
            convertTo = 'week'))
        self.assertEqual(0, util.addAgeStrings(0, nullCase, convertFrom = 'week', 
            convertTo = 'week'))


    def test_numericTimeConvert(self):
        """ Check: Convert numbers, use hand-written test cases. Catch exceptions
            on bad entries and name conversions where possible. """
        
        case1a = '10 weeks old'
        case1b = '10-week-old'
        case1c = 'ten-week old'
        case2 = '10 weeks for 5 weeks'
        case3 = '6-10 weeks for 7 days'
        case4 = '6-8 days for 7 days'
        case5 = '4day-10mos for 1 week'
        case6a = '6.5-10.5 weeks'
        case6b = '6.5-10weeks'
        nullCase1 = '6-10 eons for 7 iotas'
        
        self.assertEqual(10, util.numericTimeConvert(case1a, convertTo = 'week'))
        self.assertEqual(70, util.numericTimeConvert(case1a, convertTo = 'Days'))
        self.assertEqual('n/a', util.numericTimeConvert(case1a, convertTo = 'Days', 
            checkConverts = ['day', 'month']))
        self.assertEqual(10, util.numericTimeConvert(case1b, convertTo = 'week'))
        self.assertEqual(70, util.numericTimeConvert(case1b, convertTo = 'Days'))
        self.assertEqual('n/a', util.numericTimeConvert(case1b, convertTo = 'Days', 
            checkConverts = ['day', 'month']))
        self.assertEqual(10, util.numericTimeConvert(case1c, convertTo = 'week'))
        self.assertEqual(70, util.numericTimeConvert(case1c, convertTo = 'Days'))
        self.assertEqual('n/a', util.numericTimeConvert(case1c, convertTo = 'Days', 
            checkConverts = ['day', 'month']))
        self.assertEqual(41, int(util.numericTimeConvert(case5)))
        self.assertEqual(8.5, util.numericTimeConvert(case6a))
        self.assertEqual(8.25, util.numericTimeConvert(case6b))
        
        with self.assertRaises(ValueError):
            util.numericTimeConvert(case1a, convertTo = 'hrs')
        with self.assertRaises(ValueError):
            util.numericTimeConvert(case1a, convertTo = 'week', 
                checkConverts = ['hrs', 'weeks'])
    
        self.assertEqual(15, util.numericTimeConvert(case2, convertTo = 'week'))
        self.assertEqual(9, util.numericTimeConvert(case3, convertTo = 'week'))
        self.assertEqual(2, util.numericTimeConvert(case4, convertTo = 'week'))
        self.assertEqual('n/a', util.numericTimeConvert(nullCase1, convertTo = 'week'))
        
        for text in ['7-8week', '7-8 weeks', '7wks-8wks', '7-8wk', '7weeks-8weeks', 
                    '6.5wks-8.5wks', '7wks-8 weeks']:
            self.assertEqual(7.5, util.numericTimeConvert(text, convertTo = 'week'))

        for text in ['10-12week', '10-12 weeks', '10wks-12wks', '10-12wk', 
                    '10.5weeks-11.5weeks', '10weeks-12weeks', '10weeks-12 wks']:
            self.assertEqual(11, util.numericTimeConvert(text, convertTo = 'week'))
        
        for text in ['4-10day', '4-10 days', '4days-10days', '4-10days', 
                       '4days-10days', '4days-10 days', '4d-10days']:
            self.assertEqual(1, util.numericTimeConvert(text, convertTo = 'week'))
        
        
    def test_geoAgeExtract(self):
        """ Check: Proper age on hand-picked test samples. Failed pickups on 
            alternative parse IDs. Check age extraction from GSE if age cannot
            be detected from sample. 
            
            TODO: Full text checks"""

        samp1, samp2, samp3 = 'GSM400641', 'GSM937915', 'GSM1402452'
        samp4, samp5, samp6 = 'GSM1418737', 'GSM1282831', 'GSM172972'
        samp1Text = requests.get(gsmURL + samp1).text
        samp2Text = requests.get(gsmURL + samp2).text
        samp3Text = requests.get(gsmURL + samp3).text
        samp4Text = requests.get(gsmURL + samp4).text
        samp5Text = requests.get(gsmURL + samp5).text
        samp6Text = requests.get(gsmURL + samp6).text
        """
        v1, s1 = util.geoAgeExtract(samp1Text)
        self.assertAlmostEqual(4.28571428, v1)
        self.assertEqual('Sample', s1)

        v2, s2 = util.geoAgeExtract(samp2Text)
        self.assertEqual('n/a', v2)
        self.assertEqual('n/a', s2)

        v3, s3 = util.geoAgeExtract(samp3Text, tryAgePMID = False)
        self.assertEqual('n/a', v3)
        self.assertEqual('Study', s3)

        v4, s4 = util.geoAgeExtract(samp4Text)
        self.assertEqual(13.5, v4)
        self.assertEqual('Sample', s4)

        v5, s5 = util.geoAgeExtract(samp5Text)
        self.assertEqual(16, v5)
        self.assertEqual('Sample', s5)

        v6a, s6a = util.geoAgeExtract(samp6Text, tryAgeStudy = False)
        self.assertEqual('n/a', v6a)
        self.assertEqual('Sample', s6a)
        
        v6b, s6b = util.geoAgeExtract(samp6Text, tryAgeStudy = True)
        self.assertEqual(13.5, v6b)
        self.assertEqual('Study', s6b)

        v1, s1 = util.geoAgeExtract(samp1Text, parseAgeIDs = ['Characteristics'])
        self.assertAlmostEqual(4.28571428, v1)
        
        v2, s2 = util.geoAgeExtract(samp2Text, parseAgeIDs = ['Characteristics'])
        self.assertEqual('n/a', v2)
        
        v4, s4 = util.geoAgeExtract(samp4Text, parseAgeIDs = ['Characteristics'])
        self.assertEqual(1, v4)
        
        v5a, s5a = util.geoAgeExtract(samp5Text, parseAgeIDs = ['Characteristics'],
            tryAgeStudy = False)
        self.assertEqual('n/a', v5a)
        
        v5b, s5b = util.geoAgeExtract(samp5Text, parseAgeIDs = ['Characteristics'],
            tryAgeStudy = True)
        self.assertEqual(16, v5b)
        """

if __name__ == '__main__':

    unittest.main()
