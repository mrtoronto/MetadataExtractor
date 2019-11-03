""" Unit-testing for age-related functions. Keep < 30sec if possible """

import unittest, sys, os, re, warnings, requests

sys.path.append('./src/')
import setup_metadataExtract as util

gsmURL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='

class TestDataProcess(unittest.TestCase):

    def test_enumAgeStrings(self):
        """ Check: Convert numbers, use hand-written test cases. """
        
        case1 = '10 weeks old'
        case2 = '10 weeks for 5 weeks'
        case3 = '6-10 weeks for 7 days'
        case4 = '6-8 days for 7 days'
        nullCase = '6-10 eons for 7 iotas'

        val, durs = util.enumAgeStrings([], [], case1, convertFrom = 'week', 
            convertTo = 'week')
        self.assertEqual(10, val[0])
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], case2, convertFrom = 'week', 
            convertTo = 'week')
        self.assertEqual(10, val[0])
        self.assertEqual(5, durs[0])

        val, durs = util.enumAgeStrings([], [], case3, convertFrom = 'week', 
            convertTo = 'week')
        self.assertEqual(8, val[0])
        self.assertEqual([], durs)
        
        val, durs = util.enumAgeStrings([], [], case4, convertFrom = 'week', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], case1, convertFrom = 'day', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)
        
        val, durs = util.enumAgeStrings([], [], case2, convertFrom = 'day', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)
        
        val, durs = util.enumAgeStrings([], [], case3, convertFrom = 'day', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual(1, durs[0])

        val, durs = util.enumAgeStrings([], [], case4, convertFrom = 'day', 
            convertTo = 'week')
        self.assertEqual(1, val[0])
        self.assertEqual(1, durs[0])
        
        val, durs = util.enumAgeStrings([], [], case1, convertFrom = 'month', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], case2, convertFrom = 'month', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], case3, convertFrom = 'month', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], case4, convertFrom = 'month', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)

        val, durs = util.enumAgeStrings([], [], nullCase, convertFrom = 'week', 
            convertTo = 'week')
        self.assertEqual([], val)
        self.assertEqual([], durs)


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
        case7 = '20wks-1.5yrs'
        nullCase1 = '6-10 eons for 7 iotas'
        
        val, flag = util.numericTimeConvert(case1a, convertTo = 'week')
        self.assertEqual(10, val)
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case1a, convertTo = 'Days')
        self.assertEqual(70, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case1a, convertTo = 'Days',
            checkConverts = ['day', 'month'])
        self.assertEqual('n/a', val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case1b, convertTo = 'week')
        self.assertEqual(10, val)
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case1b, convertTo = 'Days')
        self.assertEqual(70, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case1b, convertTo = 'Days',
            checkConverts = ['day', 'month'])
        self.assertEqual('n/a', val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case1c, convertTo = 'week')
        self.assertEqual(10, val)
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case1c, convertTo = 'Days')
        self.assertEqual(70, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case1c, convertTo = 'Days',
            checkConverts = ['day', 'month'])
        self.assertEqual('n/a', val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case5, convertTo = 'week')
        self.assertEqual(21, int(val))
        self.assertEqual(True, flag)

        val, flag = util.numericTimeConvert(case5, convertTo = 'week',
            flagRange = False)
        self.assertEqual(21, int(val))
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case6a, convertTo = 'week')
        self.assertEqual(8.5, val)
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case6b, convertTo = 'week')
        self.assertEqual(8.25, val)
        self.assertEqual(False, flag)

        val, flag = util.numericTimeConvert(case7, convertTo = 'week')
        self.assertEqual(49, val)
        self.assertEqual(True, flag)
        
        with self.assertRaises(ValueError):
            util.numericTimeConvert(case1a, convertTo = 'hrs')
        with self.assertRaises(ValueError):
            util.numericTimeConvert(case1a, convertTo = 'week', 
                checkConverts = ['hrs', 'weeks'])
    
        val, flag = util.numericTimeConvert(case2, convertTo = 'week')
        self.assertEqual(15, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case3, convertTo = 'week')
        self.assertEqual(9, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(case4, convertTo = 'week')
        self.assertEqual(2, val)
        self.assertEqual(False, flag)
        
        val, flag = util.numericTimeConvert(nullCase1, convertTo = 'week')
        self.assertEqual('n/a', val)
        self.assertEqual(False, flag)

        
        for text in ['7-8week', '7-8 weeks', '7wks-8wks', '7-8wk', '7weeks-8weeks', 
                    '6.5wks-8.5wks', '7wks-8 weeks']:
            val, flag = util.numericTimeConvert(text, convertTo = 'week')
            self.assertEqual(7.5, val)

        for text in ['10-12week', '10-12 weeks', '10wks-12wks', '10-12wk', 
                    '10.5weeks-11.5weeks', '10weeks-12weeks', '10weeks-12 wks']:
            val, flag = util.numericTimeConvert(text, convertTo = 'week')
            self.assertEqual(11, val)
        
        for text in ['4-10day', '4-10 days', '4days-10days', '4-10days', 
                       '4days-10days', '4days-10 days', '4d-10days']:
            val, flag = util.numericTimeConvert(text, convertTo = 'week')
            self.assertEqual(1, val)
        

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
        
        v1, s1, f1 = util.geoAgeExtract(samp1Text)
        self.assertAlmostEqual(4.28571428, v1)
        self.assertEqual('Sample', s1)
        self.assertEqual(False, f1)

        v2, s2, f2 = util.geoAgeExtract(samp2Text)
        self.assertEqual('n/a', v2)
        self.assertEqual('n/a', s2)
        self.assertEqual(False, f2)

        v3, s3, f3 = util.geoAgeExtract(samp3Text, tryAgePMID = False)
        self.assertEqual('n/a', v3)
        self.assertEqual('Study', s3)
        self.assertEqual(False, f3)

        v4, s4, f4 = util.geoAgeExtract(samp4Text)
        self.assertEqual(13.5, v4)
        self.assertEqual('Sample', s4)
        self.assertEqual(False, f4)

        v5, s5, f5 = util.geoAgeExtract(samp5Text)
        self.assertEqual(16, v5)
        self.assertEqual('Sample', s5)
        self.assertEqual(False, f5)

        v6a, s6a, f6 = util.geoAgeExtract(samp6Text, tryAgeStudy = False)
        self.assertEqual('n/a', v6a)
        self.assertEqual('Sample', s6a)
        self.assertEqual(False, f6)
        
        v6b, s6b, _ = util.geoAgeExtract(samp6Text, tryAgeStudy = True)
        self.assertEqual(6.75, v6b)
        self.assertEqual('Study', s6b)

        v1, s1, _ = util.geoAgeExtract(samp1Text, parseAgeIDs = ['Characteristics'])
        self.assertAlmostEqual(4.28571428, v1)
        
        v2, s2, _ = util.geoAgeExtract(samp2Text, parseAgeIDs = ['Characteristics'])
        self.assertEqual('n/a', v2)
        
        v4, s4, _ = util.geoAgeExtract(samp4Text, parseAgeIDs = ['Characteristics'])
        self.assertEqual(1, v4)
        
        v5a, s5a, _ = util.geoAgeExtract(samp5Text, parseAgeIDs = ['Characteristics'],
            tryAgeStudy = False)
        self.assertEqual('n/a', v5a)
        
        v5b, s5b, _ = util.geoAgeExtract(samp5Text, parseAgeIDs = ['Characteristics'],
            tryAgeStudy = True)
        self.assertEqual(16, v5b)
        

if __name__ == '__main__':

    unittest.main()
