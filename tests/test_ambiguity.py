""" Tests for ambiguity.py
"""
import unittest

from ambiguity import ambiguity
import numpy as np
import scipy.signal
import os


class TestAmbiguity(unittest.TestCase):
    """ Simple unit tests for ambiguity
    """

    SIMPLE_ARGS = {
        "u_basic": np.ones((1, 51)),
        "fcode": True,
        "f_basic": np.dot(0.0031,
                          np.array(np.arange(-25, 26))),
        "F": 6,
        "K": 360,  # 60,
        "T": 1.1,
        "N": 360,  # 60,
        "sr": 10,
        "plot1_file": "fig_1.png",
        "plot2_file": "fig_2.png",
        "plot_format": "png",
    }

    _OUTPUT_DIR = 'test_plots'

    _OUTPUT_FORMAT = "png"

    _TEST_CASES = {

        "Pulse": {
            'u_basic': np.array([[1]]),
            'fcode': False,
            'f_basic': None,
            'F': 4,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'LFM': {
            'u_basic': np.ones((1, 51)),
            'fcode': True,
            'f_basic': np.multiply(0.0031, np.array(np.arange(-25, 26))),
            'F': 6,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'Weighted LFM': {
            'u_basic': np.array([np.sqrt(scipy.signal.chebwin(51, 50))]),
            'fcode': True,
            'f_basic': np.multiply(0.0031, np.array(np.arange(-25, 26))),
            'F': 6,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'Costas 7': {
            'u_basic': np.ones((1, 7)),
            'fcode': True,
            'f_basic': np.array([[4, 7, 1, 6, 5, 2, 3]]),
            'F': 12,
            'K': 60,
            'T': 1.1,
            'N': 80,
            'sr': 20,
        },

        'Barker 13': {
            'u_basic': np.array([[1, 1, 1, 1,
                                  1, -1, -1, 1,
                                  1, -1, 1, -1,
                                  1]]),
            'fcode': False,
            'f_basic': None,
            'F': 10,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'Frank 16': {
            'u_basic': np.array([[1, 1, 1, 1, 1,
                                  1j, -1, -1j, 1,
                                  -1, 1, -1, 1,
                                  -1j, -1, 1j]]),
            'fcode': False,
            'f_basic': None,
            'F': 10,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'P4 25': {
            'u_basic': np.exp(1j * np.pi * (1./25. *
                              np.power(np.array([np.arange(0, 25)]), 2) -
                                       np.array([np.arange(0, 25)]))),
            'fcode': False,
            'f_basic': None,
            'F': 15,
            'K': 80,
            'T': 1.1,
            'N': 80,
            'sr': 20,
        },

        'Complementary Pair': {
            'u_basic': np.array([[1, 1, -1, 0,
                                  0, 0, 0, 0,
                                  0, 1, 1j, 1]]),
            'fcode': False,
            'f_basic': None,
            'F': 10,
            'K': 60,
            'T': 1.1,
            'N': 60,
            'sr': 10,
        },

        'Pulse Train 1': {
            'u_basic': np.array([[1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1]]),
            'fcode': False,
            'f_basic': None,
            'F': 16,
            'K': 80,
            'T': 1.05,
            'N': 100,
            'sr': 10,
        },

        'Pulse Train 2': {
            'u_basic': np.array([[1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1]]),
            'fcode': False,
            'f_basic': None,
            'F': 12,
            'K': 80,
            'T': 0.042,
            'N': 60,
            'sr': 10,
        },

        'Stepped Freq. Pulse Train': {
            'u_basic': np.array([[1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1, 0, 0, 0, 0,
                                  1]]),
            'fcode': True,
            'f_basic': np.multiply(0.78, np.array([[0, 0, 0, 0, 0,
                                                    1, 0, 0, 0, 0,
                                                    2, 0, 0, 0, 0,
                                                    3, 0, 0, 0, 0,
                                                    4, 0, 0, 0, 0,
                                                    5]])),
            'F': 12,
            'K': 80,
            'T': 0.042,
            'N': 60,
            'sr': 10,
        },

        'Weighted Stepped Freq. Pulse Train': {
            'u_basic': np.multiply(np.conj(np.sqrt(scipy.signal.chebwin(36, 50))),
                                   np.array([[1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1, 0, 0, 0, 0,
                                              1]])),
            'fcode': True,
            'f_basic': np.multiply(0.7, np.array([[0, 0, 0, 0, 0,
                                                   1, 0, 0, 0, 0,
                                                   2, 0, 0, 0, 0,
                                                   3, 0, 0, 0, 0,
                                                   4, 0, 0, 0, 0,
                                                   5, 0, 0, 0, 0,
                                                   6, 0, 0, 0, 0,
                                                   7]])),
            'F': 16,
            'K': 70,
            'T': 0.03,
            'N': 60,
            'sr': 5,
        },

    }

    def setUp(self):
        if not os.path.exists(self._OUTPUT_DIR):
            os.makedirs(self._OUTPUT_DIR)

    def test_ambiguity_fcode1(self):
        """ Test basic functionality
        """
        args = self.SIMPLE_ARGS.copy()
        args['fcode'] = True
        args['plot1_file'] = os.path.join(self._OUTPUT_DIR,
                                          "fig_1_fcode1.png")
        args['plot2_file'] = os.path.join(self._OUTPUT_DIR,
                                          "fig_2_fcode1.png")
        print(args)
        self.assertTrue(ambiguity(**args))

    def test_ambiguity_fcode0(self):
        """ Test code path when no frequency coding.
        """
        args = self.SIMPLE_ARGS.copy()
        args['fcode'] = False
        args['plot1_file'] = os.path.join(self._OUTPUT_DIR,
                                          "fig_1_fcode0.png")
        args['plot2_file'] = os.path.join(self._OUTPUT_DIR,
                                          "fig_2_fcode0.png")
        print(args)
        self.assertTrue(ambiguity(**args))

    def test_ambiguity_signals(self):
        """ Test all the given sample signals
        """
        xtn = self._OUTPUT_FORMAT
        for signal_name in self._TEST_CASES.keys():
            print(signal_name)
            plot1_file = os.path.join(self._OUTPUT_DIR,
                                      signal_name + "_fig_1." + xtn)
            plot2_file = os.path.join(self._OUTPUT_DIR,
                                      signal_name + "_fig_2." + xtn)
            args = self._TEST_CASES[signal_name]
            args['plot1_file'] = plot1_file
            args['plot2_file'] = plot2_file
            args['plot_format'] = xtn
            print(args)
            self.assertTrue(ambiguity(**args))
