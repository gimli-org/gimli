#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestFrameworks(unittest.TestCase):
    def test_Fit(self):
        """
        """
        func = lambda x, a, b, c: a + b * x + c * x**2
        x = np.linspace(0, 1, 10)
        data = func(x, 1.5, 2, 2.5)
        model, response = pg.frameworks.fit(func, data, x=x)
        np.testing.assert_allclose(model, [1.5, 2, 2.5])
        np.testing.assert_allclose(data, response)

        func = lambda t, a, b: a*np.exp(b*t)
        t = np.linspace(1, 2, 100)
        data = func(t, 1.1, 2.2)
        model, response = pg.frameworks.fit(func, data, t=t)

        # pg.plt.plot(t, data, '.', label='data')
        # pg.plt.plot(t, response, label='response')
        # pg.wait()

        np.testing.assert_allclose(model, [1.1, 2.2])
        np.testing.assert_allclose(data, response)


if __name__ == '__main__':

    #test = TestFrameworks()
    #test.test_Fit()

    unittest.main()