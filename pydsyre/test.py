# Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
#
# This file is part of the heyoka.py library.
#
# This Source Code Form is subject to the terms of the Mozilla
# Public License v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest as _ut

class expression_class_test(_ut.TestCase):
    def runTest(self):
        self.test_expression_construction()
        self.test_random_constants()

    def test_expression_construction(self):
        from pydsyre import expression, get_reverse_kernel_map
        kernels = ["sum","mul","sin", "inv", "diff"]
        ex = expression(nvars=1, ncons=2, kernels=kernels)
        ker_idx = ex.get_kernels_idx()
        map = get_reverse_kernel_map()
        for idx, ker_name in zip(ker_idx, kernels):
            self.assertEqual(map[idx], ker_name)

    def test_expression_random_constants(self):
       from pydsyre import expression
       N_test=123
       ex = expression(nvars=1, ncons=N_test, kernels=["sum","mul","sin"])
       cons = ex.random_constants(-1.234,1.234)
       self.assertEqual(len(cons), N_test)
       for con in cons:
           self.assertGreater(con, -1.234)
           self.assertLess(con, 1.234)

def run_test_suite():
    import numpy as np

    retval = 0

    suite = _ut.TestLoader().loadTestsFromTestCase(expression_class_test)
    #suite.addTest(real128_test_case())


    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError("One or more tests failed.")