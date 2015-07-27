#!/usr/bin/env python

# write a correct test!
import unittest

import pygimli as pg

class TestIterBug(unittest.TestCase):
    
    def test_MissingRefCounter(self):
        """
        """
        a = pg.RVector(10, 1)

        # das geht schief wegen fehlendem referenzcounter. der Iter nutzt das
        # temporäre Object a(0,9) das nicht weiter gezählt wird
        # ich glaub unsere pygimli objecte brauchen einen refcounter wenn sie von py aus generiert werden
        # print((a(0, 9).beginPyIter()[0], "!=", a[0]))
        self.assertNotEqual(a(0, 9).beginPyIter()[0], a[0])

        # das geht weil wir einen eigenen iter definieren der die referenz hällt
        # see: class VectorIter: in pygimli.__init__.py
        i = 0
        for ai in a(0, 9):
            self.assertEqual(ai, a[i])
            i = i +1
            
if __name__ == '__main__':
    
    unittest.main()


