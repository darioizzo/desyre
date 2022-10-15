# Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
#
# This file is part of the heyoka.py library.
#
# This Source Code Form is subject to the terms of the Mozilla
# Public License v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import unittest as _ut

class expression_class_test(_ut.TestCase):

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

    def test_expression_random_genotype(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "inv", "diff"]
        ex = expression(nvars=1, ncons=2, kernels=kernels)
        geno = ex.random_genotype(length = 20)
        self.assertEqual(len(geno), 20 * 3)
        ex.check_genotype(geno)

    def test_expression_remove_nesting(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin"]
        ex = expression(nvars=1, ncons=0, kernels=kernels)
        geno = [0,1,1,2,1,0,2,2,0]
        geno2 = ex.remove_nesting(geno)
        self.assertNotEqual(geno, geno2)

    def test_expression_phenotype(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "diff"]
        ex = expression(nvars=1, ncons=0, kernels=kernels)
        geno = [2, 0, 0, 3, 1, 0, 1, 2, 0] # ['x0','sin(x0)','(sin(x0)-x0)','((sin(x0)-x0)*x0)']
        vars = [1.]
        phen = ex.phenotype(geno=geno, vars=vars, cons=[])
        gt = [1.0, 0.8414709848078965, -0.1585290151921035, -0.1585290151921035]
        self.assertEqual(len(phen), 4)
        for it1,it2 in zip(gt, phen):
            self.assertAlmostEqual(it1,it2)

    def test_expression_sphenotype(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "diff"]
        ex = expression(nvars=1, ncons=0, kernels=kernels)
        geno = [2, 0, 0, 3, 1, 0, 1, 2, 0] # ['x0','sin(x0)','(sin(x0)-x0)','((sin(x0)-x0)*x0)']
        gt = ['x0','sin(x0)','(sin(x0)-x0)','((sin(x0)-x0)*x0)']
        phen = ex.sphenotype(geno=geno)
        self.assertEqual(len(phen), 4)
        for it1,it2 in zip(gt, phen):
            self.assertEqual(it1,it2)

    def test_expression_complexity(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "diff"]
        ex = expression(nvars=1, ncons=0, kernels=kernels)
        geno = [0, 0, 0, 0, 0, 1, 3, 1, 2, 3, 0, 2, 3, 2, 0, 0, 0, 0, 3, 0, 6, 3, 3, 1, 1, 2, 4, 1, 7, 6, 2, 0, 1, 0, 4, 11, 3, 12, 10, 1, 1, 3, 2, 13, 4, 0, 13, 4, 1, 15, 7, 3, 3, 10, 3, 13, 17, 1, 14, 8]
        complexity = ex.complexity(geno = geno)
        self.assertEqual(len(complexity), 21)
        gt = [1, 3, 5, 9, 7, 7, 3, 5, 13, 13, 9, 2, 10, 20, 13, 21, 28, 27, 19, 48, 27]
        for it1,it2 in zip(gt, complexity):
            self.assertEqual(it1,it2)

    def test_expression_mse_ddmse(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "diff"]
        ex = expression(nvars=3, ncons=3, kernels=kernels)
        geno = ex.random_genotype(length = 200)
        cons = ex.random_constants(lb =-1.,ub = 1.)
        xs = [[1.,2.,3.], [0.1,0.2,-0.3]]
        ys = [0.,1]
        mse = ex.mse(geno=geno, cons=cons, xs = xs, ys = ys)
        mse2, grad, hess = ex.ddmse(geno=geno, cons=cons, xs = xs, ys = ys)
        for it1,it2 in zip(mse, mse2):
            self.assertEqual(it1,it2)

    def test_expression_check_genotype(self):
        from pydsyre import expression
        kernels = ["sum", "mul", "sin", "diff"]
        ex = expression(nvars=3, ncons=3, kernels=kernels)
        geno = ex.random_genotype(length = 200)
        ex.check_genotype(geno)

    

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