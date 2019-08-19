import numpy as np


def assert_array_equals_with_index_error_message(test_instance, expected, actual):
    """

    :param test_instance: unittest.TestCase
    :param expected: numpy.ndarray
    :param actual: numpy.ndarray
    :return:
    """
    self = test_instance
    ex_sh = expected.shape
    ac_sh = actual.shape
    msg = '\nArrays have not the same size.\n' + \
          'expected: ' + str(ex_sh) + '\n' + \
          'actual  : ' + str(ac_sh)
    self.assertEqual(ex_sh, ac_sh, msg)

    ex_it = np.nditer(expected, flags=['multi_index'])
    while not ex_it.finished:
        idx = ex_it.multi_index
        exp = expected.item(idx)
        act = actual.item(idx)
        msg = '\nValue at index ' + str(idx) + " is not equal to expected.\n" + \
              'expected: ' + str(exp) + '\n' + \
              'actual  : ' + str(act)
        self.assertAlmostEqual(exp, act, msg=msg)

        ex_it.iternext()
