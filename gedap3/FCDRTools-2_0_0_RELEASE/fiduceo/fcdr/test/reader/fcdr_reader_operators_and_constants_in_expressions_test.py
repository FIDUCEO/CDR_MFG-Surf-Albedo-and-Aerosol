import unittest as ut
import xarray as xr
import numpy as np
import fiduceo.fcdr.test.test_utils as tu
from fiduceo.fcdr.reader.fcdr_reader import FCDRReader as R
import math

t, f = True, False
pi = math.pi


class FCDRReaderExpressionOperatorsAndConstantsTest(ut.TestCase):

    def test_and(self):
        expression = 'a & b'  # the meaning is "a and b"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_bool_variable()
        ds['b'] = get_two_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray([[t, f, f, f], [f, f, f, f], [t, f, f, f]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_or(self):
        expression = 'a | b'  # the meaning is "a or b"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_bool_variable()
        ds['b'] = get_two_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray([[t, f, t, t], [t, f, t, t], [t, f, t, f]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_not(self):
        expression = '~ a'  # the meaning is "not a"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], np.asarray([f, t, f, t]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_and_with_vertical_1D(self):
        expression = 'a & b'  # the meaning is "a and b"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_vertical_one_dim_bool_variable()
        ds['b'] = get_two_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray([[f, f, f, f], [f, f, f, t], [f, f, f, f]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_or_with_vertical_1D(self):
        expression = 'a | b'  # the meaning is "a or b"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_vertical_one_dim_bool_variable()
        ds['b'] = get_two_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray([[t, f, f, t], [t, t, t, t], [t, f, f, f]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_not_with_vertical_1D(self):
        expression = '~ a'  # the meaning is "not a"
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_vertical_one_dim_bool_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y'], np.asarray([t, f, t]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_less_than(self):
        expression = 'a < 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[t, t, t, t], [t, t, t, f], [t, t, f, f]],
             [[t, f, f, f], [f, f, f, f], [f, f, f, f]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_less_than_or_equal(self):
        expression = 'a <= 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[t, t, t, t], [t, t, t, t], [t, t, t, f]],
             [[t, t, f, f], [t, f, f, f], [f, f, f, f]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_equal(self):
        expression = 'a == 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[f, f, f, f], [f, f, f, t], [f, f, t, f]],
             [[f, t, f, f], [t, f, f, f], [f, f, f, f]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_not_equal(self):
        expression = 'a != 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[t, t, t, t], [t, t, t, f], [t, t, f, t]],
             [[t, f, t, t], [f, t, t, t], [t, t, t, t]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_greater_than_or_equal(self):
        expression = 'a >= 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[f, f, f, f], [f, f, f, t], [f, f, t, t]],
             [[f, t, t, t], [t, t, t, t], [t, t, t, t]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_greater_than(self):
        expression = 'a > 5'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_three_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['z', 'y', 'x'], np.asarray(
            [[[f, f, f, f], [f, f, f, f], [f, f, f, t]],
             [[f, f, t, t], [f, t, t, t], [t, t, t, t]]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_unary_negation(self):
        expression = '- a'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[-10, -20, -30, -40], [-13, -23, -33, -43], [-16, -26, -36, -46]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_addition(self):
        expression = 'a + b'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        ds['b'] = get_vertical_one_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[210, 220, 230, 240], [313, 323, 333, 343], [416, 426, 436, 446]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_subtraction(self):
        expression = 'b - a'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        ds['b'] = get_vertical_one_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[190, 180, 170, 160], [287, 277, 267, 257], [384, 374, 364, 354]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_multiplication(self):
        expression = 'a * b'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        ds['b'] = get_vertical_one_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[2000, 4000, 6000, 8000], [3900, 6900, 9900, 12900], [6400, 10400, 14400, 18400]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_division(self):
        expression = 'b / a'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        ds['b'] = get_vertical_one_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        loaded_ = ds['v']
        # if two arrays of type integer are divided, the result data type can be integer or float
        # the expected preparation depends on returned datatype
        if np.issubdtype(loaded_.dtype, np.integer):
            expected = xr.Variable(['y', 'x'], np.asarray(
                [[20, 10, 6, 5],
                 [23, 13, 9, 6],
                 [25, 15, 11, 8]]))
        else:
            expected = xr.Variable(['y', 'x'], np.asarray(
                [[20, 10, 6.666666666666667, 5],
                 [23.076923076923077, 13.043478260869565, 9.090909090909092, 6.976744186046512],
                 [25, 15.384615384615385, 11.11111111111111, 8.695652173913043]]))
        tu.assert_array_equals_with_index_error_message(self, expected, loaded_)

    def test_power(self):
        expression = 'a ** 2'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[100, 400, 900, 1600], [169, 529, 1089, 1849], [256, 676, 1296, 2116]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_modulo(self):
        expression = 'a % 7'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_two_dim_int_variable()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['y', 'x'], np.asarray(
            [[3, 6, 2, 5], [6, 2, 5, 1], [2, 5, 1, 4]]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_sinus(self):
        expression = 'sin(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_radians()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], np.asarray(
            [0,
             0.70710678,
             0.989821441881,
             0.70710678,
             -0.70710678,
             -0.989821441881]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_cosine(self):
        expression = 'cos(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_radians()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], np.asarray(
            [1,
             0.70710678,
             0.142314838273,
             -0.70710678,
             -0.70710678,
             -0.142314838273]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_tangent(self):
        expression = 'tan(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = get_one_dim_radians()
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], np.asarray(
            [0,
             1,
             6.95515277177,
             -1,
             1,
             6.95515277177]))
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_sinus(self):
        expression = 'arcsin(sin(a))'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0, 0.5, 1])
        R._load_virtual_variable(ds, 'v')
        expected = ds['a']
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_cosine(self):
        expression = 'arccos(cos(a))'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [1, 0.5, 0])
        R._load_virtual_variable(ds, 'v')
        expected = ds['a']
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_tangent(self):
        expression = 'arctan(tan(a))'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [1, 0.5, 0])
        R._load_virtual_variable(ds, 'v')
        expected = ds['a']
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_arctan2(self):
        expression = 'arctan2(a, b)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        ds['b'] = xr.Variable(('x',), [0.8, 0.7, 0.6])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.atan2(0.2, 0.8),
                                        math.atan2(0.4, 0.7),
                                        math.atan2(0.5, 0.6), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_hyperbolic_sine(self):
        expression = 'sinh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.sinh(0.2),
                                        math.sinh(0.4),
                                        math.sinh(0.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_hyperbolic_cosine(self):
        expression = 'cosh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.cosh(0.2),
                                        math.cosh(0.4),
                                        math.cosh(0.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_hyperbolic_tangent(self):
        expression = 'tanh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.tanh(0.2),
                                        math.tanh(0.4),
                                        math.tanh(0.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_hyperbolic_sine(self):
        expression = 'arcsinh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.asinh(0.2),
                                        math.asinh(0.4),
                                        math.asinh(0.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_hyperbolic_cosine(self):
        expression = 'arccosh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2, 4, 5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.acosh(2),
                                        math.acosh(4),
                                        math.acosh(5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_inverse_hyperbolic_tangent(self):
        expression = 'arctanh(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [0.2, 0.4, 0.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.atanh(0.2),
                                        math.atanh(0.4),
                                        math.atanh(0.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_natural_logarithm(self):
        expression = 'log(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2.2, 3.4, 4.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.log(2.2),
                                        math.log(3.4),
                                        math.log(4.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_base_10_logarithm(self):
        expression = 'log10(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2.2, 3.4, 4.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.log10(2.2),
                                        math.log10(3.4),
                                        math.log10(4.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_natural_logarithm_1_plus_x(self):
        expression = 'log1p(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2.2, 3.4, 4.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.log1p(2.2),
                                        math.log1p(3.4),
                                        math.log1p(4.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_exponential(self):
        expression = 'exp(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2.2, 3.4, 4.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.exp(2.2),
                                        math.exp(3.4),
                                        math.exp(4.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_exponential_minus_one(self):
        expression = 'expm1(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [2.2, 3.4, 4.5])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [math.expm1(2.2),
                                        math.expm1(3.4),
                                        math.expm1(4.5), ])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_square_root(self):
        expression = 'sqrt(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [4, 9, math.pow(12, 2)])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [2, 3, 12])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_absolute_value(self):
        expression = 'abs(a)'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(('x',), [-2.2, -3.4, 4.5 * -1])
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(('x',), [2.2, 3.4, 4.5])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_replace_pi(self):
        expression = 'a*pi'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(['x'], np.asarray([1, 2, 3, 4, 0.01]))
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], [3.14159265359, 6.28318530718, 9.42477796077, 12.5663706144, 0.0314159265359])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_dont_replace_pi_if_it_part_of_a_word(self):
        expression = 'lpit + pit + lpi'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['lpit'] = xr.Variable(['x'], np.asarray([1, 2, 3, 4, 0.01]))
        ds['pit'] = xr.Variable(['x'], np.asarray([1, 2, 3, 4, 0.01]))
        ds['lpi'] = xr.Variable(['x'], np.asarray([1, 2, 3, 4, 0.01]))
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], [3, 6, 9, 12, 0.03])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])

    def test_number_constants_with_scientific_notation(self):
        expression = 'a + 257e-3'
        ds = xr.Dataset()
        ds['v'] = create_virtual_variable(expression)
        ds['a'] = xr.Variable(['x'], np.asarray([1, 2, 3, 4, 0.01]))
        R._load_virtual_variable(ds, 'v')
        expected = xr.Variable(['x'], [1.257, 2.257, 3.257, 4.257, 0.267])
        tu.assert_array_equals_with_index_error_message(self, expected, ds['v'])


def get_one_dim_radians():
    pi_4 = pi / 4
    pi_2 = pi / 2.2
    return xr.Variable(['x'], np.asarray([0,
                                          pi_4,
                                          pi_2,
                                          pi - pi_4,
                                          pi + pi_4,
                                          pi + pi_2]))


def get_two_dim_bool_variable():
    return xr.Variable(['y', 'x'], np.asarray(
        [[t, f, f, t], [f, f, f, t], [t, f, f, f]]))


def get_one_dim_bool_variable():
    return xr.Variable(['x'], np.asarray(
        [t, f, t, f]))


def get_vertical_one_dim_bool_variable():
    return xr.Variable(['y'], np.asarray(
        [f, t, f]))


def get_three_dim_int_variable():
    return xr.Variable(['z', 'y', 'x'], np.asarray(
        [[[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]],
         [[4, 5, 6, 7], [5, 6, 7, 8], [6, 7, 8, 9]]]))


def get_two_dim_int_variable():
    return xr.Variable(['y', 'x'], np.asarray(
        [[10, 20, 30, 40], [13, 23, 33, 43], [16, 26, 36, 46]]))


def get_vertical_one_dim_int_variable():
    return xr.Variable(['y'], np.asarray(
        [200, 300, 400]))


def create_virtual_variable(expression):
    v_var = xr.Variable([], np.NaN)
    v_var.attrs["virtual"] = "true"
    v_var.attrs["expression"] = expression
    return v_var
