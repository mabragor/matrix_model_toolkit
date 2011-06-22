load 'loop_eqn_gen.sage'

import re

class GaussianSolver(object):
    def __init__(self):
        self._not_interesting_list = var('Lambda, beta, z0, y(z0)')
        self._loop_eqn_gen = LoopEquation()
    
    def debug(self, message):
        '''
        Prints out a message surrounded by the box of asterisks, for it to be
        noticeable.
        
        @type message: string
        @param message: The message we want to print, usually to know, where we
        are.
        '''
        
        m_len = len(message)
        
        tmp_str = reduce(lambda x, y: x + y, ['*' for i in xrange(m_len + 4)])
        print tmp_str
        print '* ' + message + ' *'
        print tmp_str
    
    def _set_initial_conditions(self):
        '''
        Set initial conditions for iterative procedure, like rho_{1,0,0} and
        its derivatives.
        '''
        # Define basic variables to work with
        var ('z0 beta T2')
        
        # Define y(x) 'special' function.
        def deriv(self, *args, **kwds):
            if kwds['diff_param'] == 0:
                # print args, kwds
                return args[0] / self(args[0])
            else:
                print 'SOME CRAPTASTIC ERROR'
                
        y = function('y', nargs=1, derivative_func=deriv)
        
        # Define zero order 1-point resolvent. Ensure that both n_s_a and
        # n_d_a versions are defined.
        function('rho_1_0_0', nargs=1)
        function('rho_0_1_0', nargs=1)
        
        ans0 = [rho_1_0_0(z0) == (z0 - y(z0)) / 2 / beta,
                rho_0_1_0(z1) == (z1 - y(z1)) / 2 / beta]
                
        self._precompute_derivatives(ans0)
        
    
    def solve(self, max_order, max_diff_points=None):
        '''
        @type  max_order: int
        @param max_order: Maximum order of rho_1_0, that we want to find.
        
        @type  max_diff_points: int
        @param max_diff_points: At the highest order, we evaluate resolvents
        with only max_diff_points number of non-coincident points.
        Defaults to None.
        '''
        self._set_initial_conditions()
        
        
        # We should translate apparent order to the real order of equation
        # by subtracting -2.
        for cur_ord in xrange(1, max_order + 1):
            print '*************************************'
            print '* Starting calculation of new order *'
            print '*************************************'   
            
            ans = self._get_and_solve_eqns(cur_ord, 0)
            num_iters = max([cur['n_s_a'] + cur['n_d_a'] - 1 for cur in ans[1]])
            if cur_ord == max_order and max_diff_points != None:
                num_iters = min(num_iters, max_diff_points - 1)
            self._precompute_derivatives(ans[0])
            
            for n_d_a in xrange(1, num_iters + 1):
                ans = self._get_and_solve_eqns(cur_ord, n_d_a)
                self._precompute_derivatives(ans[0])
    
    def _z_to_y(self, foo):
        '''
        Rewrites as much as possible in terms of y(z)'s instead of z's.
        
        Accepts arguments of different types and behavior is different.
        '''
        
        if type(foo) == str:
            return self._z_to_y_func(foo)
        elif type(foo) == sage.symbolic.expression.Expression:
            return self._z_to_y_expr(foo)
        else:
            print "THIS JUST SHOULD NOT HAPPEN!"
            raise Exception
            
    def _z_to_y_str(self, expr_str):
        '''
        Rewrites as much as possible in terms of y(z)'s instead of z's.
        
        Accepts formula in string format. A low-level function in a sense that
        it does not take into account necessity to hold form of denominators.
        
        @type  expr_str: string
        @param expr_str: EXAMPLE: "z^2 - 4 * beta * Lambda"
        
        @rtype : string
        @return: Simplified expression, also in form of a string.
        '''
        
        pattern = re.compile('z(?P<var_num>\d+)\^(?P<power>\d+)')
        
        def repl_func(match):
            var_num = int(match.group('var_num'))
            var = eval('z%d' % var_num)
            power = int(match.group('power'))
            
            # print var, power
            
            if power % 2 == 0:
                return ('(y(z%d)^2 + 4 * beta * Lambda)^%d' %
                        (var_num, power / 2))
            else:
                return ('(y(z%d)^2 + 4 * beta * Lambda)^%d * z%d' %
                        (var_num, (power - 1) / 2, var_num))
                        
        new_expr_str = pattern.sub(repl_func, expr_str)
        
        return new_expr_str 
   
    def _z_to_y_expr(self, expr):
        '''
        @type  expr: symbolic expression
        @param expr: symbolic expression, we want to have as much y's and as
        little z's, as possible.
        
        @rtype: symbolic expression
        @return: Simplified version of the expression
        '''
        self.debug('Starting simplification procedure...')
        
        # First we isolate the denominator
        expr_den = expr.denominator()
        expr_num = expr.numerator().expand()
        
        print 'Denomiantor: ', expr_den
        print 'Numberator: ', expr_num
        
        expr_str = str(expr_num)
        new_expr_str = self._z_to_y_str(expr_str)
        
        new_expr = symbolic_expression(new_expr_str)
        new_expr = new_expr.simplify_full()
        
        # We plug the denominator back.
        new_expr = new_expr / expr_den
        
        return new_expr
    
    def _y_to_rads_str(self, expr_str):
        '''
        Converts all y's in the string form of expression into their
        exact radical values.
        
        @attention: This procedure is expresimental, just to test if bruteforce
        simplification will work well.
        '''
        
        pattern = re.compile('y\(z(?P<var_num>\d+)\)')
        
        def repl_func(match):
            var_num = int(match.group('var_num'))
            var = eval('z%d' % var_num)
            power = int(match.group('power'))
            
            return 'sqrt(z%d^2 - 4 * beta * Lambda)' %  var_num
                        
        new_expr_str = pattern.sub(repl_func, expr_str)
        
        return new_expr_str
    
    def _rads_to_y_str(self, expr_str):
        '''
        Rewrites explicit radicals in the compact form of y's.
        
        Replacement proceeds 
        '''
        
        pattern = re.compile('y\(z(?P<var_num>\d+)\)')
        
        def repl_func(match):
            var_num = int(match.group('var_num'))
            var = eval('z%d' % var_num)
            power = int(match.group('power'))
            
            return 'sqrt(z%d^2 - 4 * beta * Lambda)' %  var_num
                        
        new_expr_str = pattern.sub(repl_func, expr_str)
        
        return new_expr_str
    
    def _z_to_y_func(self, func_name):
        '''
        @type  func_name: string
        @param func_name: Name of the function, expression we want to simplify.
        
        @rtype : callable symbolic expression
        @return: Simplified version of the function.
        '''
        func = globals().get(func_name, None)
        
        expr = func(*func.arguments())
        
        new_expr = self._z_to_y_expr(expr)
        
        new_func = new_expr.function(*func.arguments())
        
        globals().update({func_name: new_func})
        
        return eval(func_name)
    
    def _precompute_derivatives(self, eqn_list):
        '''
        Once any new resolvent is computed, its derivatives in coincident
        points should be computed.
        Assignments of newly found function are made on the fly.
        
        This procedure computes drhodiff, drho and d2rho for any rho
        
        @type  eqn_list: list
        @param eqn_list: List of equation, defining newly computed resolvents.
        
        @rtype : list
        @return: List of newly computed resolvents, populated with their
        derivatives in coincident points.
        '''
        # This will contain instances of functions assigned
        result = []
        
        # We will use this extensively,so we proxy it.
        g = globals()
        
        for eqn in eqn_list:
            # Instances of functions, assigned in one round.
            local_result = []
            
            print 'eqn: ', eqn
            op_str = str(eqn.lhs())
            res_params = self._extract_resolvent_params(op_str)
            # These will be central objects, so we proxy them.
            n_s_a, n_d_a, order = (res_params['n_s_a'], res_params['n_d_a'],
                                   res_params['order'])
            ops = res_params['var_list']
            print n_s_a, n_d_a, order, ops
            print 
            
            g.update({'rho_%d_%d_%d' % (n_s_a, n_d_a, order):
                      eqn.rhs().function(*ops)})
            self._z_to_y('rho_%d_%d_%d' % (n_s_a, n_d_a, order))
            
            # Curing the special case, where we have only one z0
            if n_s_a == 1:
                g.update({'rho_%d_%d_%d' % (n_s_a - 1, n_d_a + 1, order):
                          eval('rho_%d_%d_%d' % (n_s_a, n_d_a, order))})
            
            f = eval('rho_%d_%d_%d' % (n_s_a, n_d_a, order))

            # print f(z0, z1)
            local_result.append(f)
                        
            # If n_s_a == 1, regardless of the circumstances, we can
            # calculate drho_1_k_ and d2rho_1_k_ and drhodiff_0_1+k_
            # This IF block is written in non-object-oriented fashion.
            # Be careful!!!
            print 'Performing differentiations w.r.t z0'
            if n_s_a == 1:
                if self._not_exists('drho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                    print 'Calculate drho'
                    # Calculate derivative of the function w.r.t z0
                    tmp_func = eval('rho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    new_func = self._smart_diff(tmp_func, ops[0])
                    
                    # Define new function
                    g.update({'drho_%d_%d_%d' % (n_s_a, n_d_a, order):
                              new_func})           
                    self._z_to_y('drho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    
                    # Add a link to the instance of newly created function
                    # to our array.
                    tmp_func = eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    local_result.append(tmp_func)
                               
                if self._not_exists('d2rho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                    print 'Calculate d2rho'
                    # Calculate second derivative.
                    tmp_func = eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    new_func = self._smart_diff(tmp_func, ops[0])
                    
                    #Define new function.
                    g.update({'d2rho_%d_%d_%d' % (n_s_a, n_d_a, order):
                              new_func})           
                    self._z_to_y('d2rho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    
                    # Add a link to the instance of newly created function
                    # to our array.
                    tmp_func = eval('d2rho_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order))
                    local_result.append(tmp_func)
                        
                if self._not_exists('drhodiff_%d_%d_%d' %
                                    (n_s_a - 1, n_d_a + 1, order)):
                    
                    print 'Calculate drhodiff'
                    g.update({'drhodiff_%d_%d_%d' %
                              (n_s_a - 1, n_d_a + 1, order):
                              eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))})
                    self._z_to_y('drhodiff_%d_%d_%d' %
                                 (n_s_a - 1, n_d_a + 1, order))
                    
                    tmp_func = eval('drhodiff_%d_%d_%d' %
                                    (n_s_a - 1, n_d_a + 1, order))
                    
                    local_result.append(tmp_func)
            
            # If n_d_a > 0, we can calculate
            # drho_i+1_k-1_, d2rho_i+1_k-1_ and drhodiff_i_k_
            
            # First some preparations are necessary:
            
            # Returns first argument, not equal to z0.
            def correct_op():
                if n_s_a == 0:
                    return ops[0]
                else:
                    return ops[1]
            
            # We may proceed with calculation.
            print 'Performing differentiations w.r.t second argument'
            
            if n_d_a > 0:
                # print 'Im here'
                if self._not_exists('drhodiff_%d_%d_%d' %
                                    (n_s_a, n_d_a, order)):
                    print 'Calculate drhodiff'
                    
                    tmp_func = eval('rho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    new_func = self._smart_diff(tmp_func, correct_op())
                    
                    g.update({'drhodiff_%d_%d_%d' % (n_s_a, n_d_a, order):
                              new_func})
                    self._z_to_y('drhodiff_%d_%d_%d' % (n_s_a, n_d_a, order))
                    
                    local_result.append(new_func)
                    
                if self._not_exists('drho_%d_%d_%d' % 
                                    (n_s_a + 1, n_d_a - 1, order)):
                    print 'Calculate drho'
                    
                    tmp_func = eval('drhodiff_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order))
                    # Variable name is hardcoded. Do something.
                    # print tmp_func(*ops), {str(correct_op()): z0}
                    # First we substitute desired argument by a 'wildcard'
                    # variable 'x'.
                    var('_x')
                    pre_limit_expr = tmp_func(*ops).subs({correct_op():
                                                          z0 + _x})
                    # Then we take the limit, taking into account, that the
                    # result could be a number.
                    print pre_limit_expr
                    tmp_rhs = symbolic_expression(pre_limit_expr.series(_x, 1).coefficient(_x, 0))
                    
                    self.debug('After taking the limit this equals: ')
                    print tmp_rhs
                    
                    # print 'tmp_rhs', tmp_rhs
                    # We generate correct argument list for resulting
                    # expression.
                    # print ops
                    if ops[0] == z0:
                        arg_list = [ops[0]] + ops[2:]
                    else:
                        arg_list = [z0] + ops[1:]
                    
                    # print arg_list
                    
                    new_func = tmp_rhs.function(*arg_list)
                    
                    # print new_func
                    
                    g.update({'drho_%d_%d_%d' % (n_s_a + 1, n_d_a - 1, order):
                              new_func})    
                    self._z_to_y('drho_%d_%d_%d' %
                                 (n_s_a + 1, n_d_a - 1, order))
                    
                    self.debug('Passed calculation of drho')
                              
                    local_result.append(new_func)
                    
                if self._not_exists('d2rhodiff_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order)):
                    print 'Calculate d2rhodiff'
                    
                    tmp_func = eval('drhodiff_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order))
                    new_func = self._smart_diff(tmp_func, correct_op())
                    
                    g.update({'d2rhodiff_%d_%d_%d' % (n_s_a, n_d_a, order):
                              new_func})
                    self._z_to_y('d2rhodiff_%d_%d_%d' % (n_s_a, n_d_a, order))
                    
                    # We do not append this function, since it does not appear
                    # in the loop equations.
                    # local_result.append(tmp_func)
                
                if self._not_exists('d2rho_%d_%d_%d' % 
                                    (n_s_a + 1, n_d_a - 1, order)):
                    print 'Calculate d2rho'
                    
                    tmp_func = eval('d2rhodiff_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order))
                    # print tmp_func(*ops)
                    pre_limit_expr = tmp_func(*ops).subs({correct_op(): z0 + _x})
                    tmp_rhs = symbolic_expression(pre_limit_expr.series(_x, 1).coefficient(_x, 0))
                    # print 'tmp_rhs', tmp_rhs
                    
                    if ops[0] == z0:
                        arg_list = [ops[0]] + ops[2:]
                    else:
                        arg_list = [z0] + ops[1:]
                    
                    new_func = tmp_rhs.function(*arg_list)
                    
                    g.update({'d2rho_%d_%d_%d' %
                              (n_s_a + 1, n_d_a - 1, order):
                              new_func})
                    self._z_to_y('d2rho_%d_%d_%d' %
                                 (n_s_a + 1, n_d_a - 1, order))
                    
                    tmp_func = eval('d2rho_%d_%d_%d' %
                                    (n_s_a + 1, n_d_a - 1, order))
                    local_result.append(tmp_func)
            
            print local_result
            result.extend(local_result)
        
        return result
    
    def _smart_diff(self, func, arg):
        '''
        Fixed differentiation procedure, which differentiates y(x) completely.
        
        @type  func: callable symbolic expression
        @param func: Function, we want to differentiate.
        
        @type  arg: variable
        @param arg: Variable, w.r.t which differentiation should be performed.
        
        @rtype : callable symbolic expression
        @return: The result of differentiation.
        '''
        new_expr = str(func.derivative(arg)(*func.arguments()))
        
        # We substitute D[0](y)(x) by x/y(x).
        # Note that this is so only for Gaussian model!
        find_pattern = re.compile('D\[0\]\(y\)\((?P<y_arg>[a-zA-Z0-9]+)\)')
        
        subs_pattern = '(\g<y_arg> / y(\g<y_arg>))'
        
        new_corr_expr = find_pattern.sub(subs_pattern, new_expr)
        
        result = symbolic_expression(new_corr_expr).function(*func.arguments())
        return result
    
    def _get_and_solve_eqns(self, order, n_d_a):
        '''
        Generates closed system of equations and then solves it.
        
        @type  order: int
        @param order: Common order of generated equations.
        
        @type  n_d_a: int
        @param n_d_a: Number of arguments other than z0 in equations. Number
        of z0's differs from equation to equation, and this is the key of the
        construction.
        
        @rtype : tuple(list of answers,
        list of names of newly defined functions)
        @return: 
        '''
        
        # First we generate list of variables and functions, we do not want to
        # solve for.
        add_args = [var('z%d' % i) for i in range(1, n_d_a + 1)]
        add_y_args = [y(var('z%d' % i)) for i in range(1, n_d_a + 1)]
        local_not_interesting_list = (list(self._not_interesting_list) + 
                                      add_args +
                                      add_y_args)
        print local_not_interesting_list
        
        # Here we collect all equations.
        all_eqns = set([])
        
        # First we generate 'first' equation and estimate number of
        # indeterminates.
        l_e_g = self._loop_eqn_gen
        l_e_g.order = -2 - n_d_a + order
        l_e_g.loop_equation(1, n_d_a)
        first_eqn = l_e_g.to_connected()
        all_eqns.add(first_eqn)
        indet_list = self._get_indets_list(first_eqn)
        
        print indet_list
        
        # Now we add required number of equations
        n_s_a = 2
        while n_s_a < len(indet_list) + 1:
            l_e_g.order = -1 - n_s_a - n_d_a + order
            l_e_g.loop_equation(n_s_a, n_d_a)
            cur_eqn = l_e_g.to_connected()
            all_eqns.add(cur_eqn)
            
            cur_indet_list = self._get_indets_list(cur_eqn)
            indet_list.update(cur_indet_list)
            
            n_s_a += 1
            print indet_list
        
        # solve currently cannot solve for functions like f(x), so we should
        # make a substitution and write f_x instead.
        
        print 'Making subs...'
        
        s_all_eqns, s_indet_list = self._make_substitution(all_eqns, indet_list)
        
        # All equations should be linear, and we make use of it here
        s_all_eqns = self._simplify_eqn_form(s_all_eqns, s_indet_list)
        
        # Now we solve
        print 'Solving equations...'
        
        print list(s_all_eqns), tuple(s_indet_list)
        
        all_ans = solve(list(s_all_eqns), *s_indet_list)
        
        print 'Got an answer:'
        if type(all_ans[0]) == list:
            all_ans = all_ans[0]
        print all_ans
        
        return (all_ans, [self._get_resolvent_params_from_name(str(cur_indet))
                          for cur_indet in indet_list])
    
    def _simplify_eqn_form(self, eqns, indets):
        '''
        Custom made simplification procedure for equations.
        
        Exploits the fact that system is linear at each order of perturbation
        theory.
        
        @type  eqns: set
        @param eqns: Set of equations.
        
        @type  indets: set
        @param indets: Set of indeterminates, we want to solve for.
        
        @rtype : set
        @return: Set of more-or-less simplified equations.
        '''
        
        simpl_eqns = set([])
        for eqn in eqns:
            simpl_eqn = 0
            for indet in indets:
                coeff = eqn.collect(indet).coefficient(indet, 1)
                print coeff
                coeff = coeff.simplify()
                
                simpl_eqn += coeff * indet
                eqn -= coeff * indet
            
            # First we simplify just in case
            eqn = eqn.full_simplify().factor()
            
            eqn = self._z_to_y(eqn)
            
            eqn = eqn.partial_fraction(y(z0))
            print '*******************************************'
            print '* Simplified form of free coefficient is: *'
            print '*******************************************'
            print eqn
            
            simpl_eqn += eqn
            
            simpl_eqns.add(simpl_eqn)
        
        return simpl_eqns
        
    
    def _make_substitution(self, eqns, indets):
        '''
        Substitutes all functions rho_i_j_k(z0,.,zk) by rho_i_j_k_z0_..._zk
        symbolic variables. Needed for solve method to work.
        
        @type  eqns: set
        @param eqns: Equations, that we want to solve.
        
        @type  indets: set
        @param indets: For what we want to solve?
        
        @rtype : tuple(new_eqns, new_indets)
        @return: A tuple, containing equations, where substitutions were made
        and new indeterminates to solve for.
        '''
        
        new_eqns = set([])
        new_indets = set([])
        
        for cur_eqn in eqns:
            cur_eqn_str = str(cur_eqn)
            for cur_indet in indets:
                ind_data = self._get_resolvent_params_from_name(str(cur_indet))
                print 'ind_data: ', ind_data
                
                func_pattern = 'rho_%d_%d_%d\(%s\)' % (ind_data['n_s_a'],
                                                       ind_data['n_d_a'],
                                                       ind_data['order'],
                                                       ind_data['var_list'])
                
                tmp_str = re.sub(',[ ]?', '_', ind_data['var_list'])
                repl_pattern = 'rho_%d_%d_%d_%s' % (ind_data['n_s_a'],
                                                    ind_data['n_d_a'],
                                                    ind_data['order'],
                                                    tmp_str)
                
                print func_pattern, repl_pattern
                
                # Define variable with new name
                var(repl_pattern)
                new_indets.add(eval(repl_pattern))
                
                # Substitute occurencies of a function by new varname in the
                # equation.
                cur_eqn_str = re.sub(func_pattern, repl_pattern, cur_eqn_str)
            print cur_eqn_str
            new_eqns.add(symbolic_expression(cur_eqn_str))
                
        return (new_eqns, new_indets)
    
    def _get_resolvent_params_from_name(self, res_name):
        '''
        Get parameters like number of arguments and order from string
        representation of resolvent name
        
        @type  res_name: string
        @param res_name: rho_i_j_k(z0,.,zk)
        
        @rtype : dict
        @return: A dictionary containing following information
            - 'n_s_a' : Number of z0 arguments
            - 'n_d_a' : Number of arguments other than z0
            - 'order' : Which genus contribution is it?
            - 'var_list' : String representing on that arguments function was
            called.
        '''
        
        pattern = re.compile('rho_' +
                             '(?P<n_s_a>\d+)_' +
                             '(?P<n_d_a>\d+)_' +
                             '(?P<order>\d+)' +
                             '\((?P<var_list>[a-z0-9, ]+)\)')
        
        m = pattern.match(res_name)
        
        result = {'n_s_a': int(m.group('n_s_a')),
                  'n_d_a': int(m.group('n_d_a')),
                  'order': int(m.group('order')),
                  'var_list': m.group('var_list')}
        
        return result
    
    def _get_indets_list(self, equation):
        '''
        Get list of indeterminate resolvents in equation.
        
        @type  equation: Expression
        @param equation: Symbolic expression of equation
        
        @rtype : set(function)
        @return: Set of symbolic functions.
        '''
        print 'loop equation: ', equation
        
                             
        equation_str = str(equation)
        # print equation_str
        
        pattern = re.compile('rho_' +
                             '(?P<n_s_a>\d+)_' +
                             '(?P<n_d_a>\d+)_' +
                             '(?P<order>\d+)' +
                             '\((?P<var_list>[a-z0-9, ]+)\)')
        
        match_list = []
        search_start = 0
        while True:
            cur_match = pattern.search(equation_str, search_start)
            if cur_match == None:
                break
            else:
                match_list.append(cur_match)
                search_start = cur_match.end()
                        
        # print match_list
        
        func_list = [globals()['rho_' + cur_p.group('n_s_a') +
                               '_' + cur_p.group('n_d_a') +
                               '_' + cur_p.group('order')]
                     for cur_p in match_list]
        
        def _l(x):
            if type(x) == tuple:
                return x
            else:
                return (x,)
        
        arg_list = [_l(var(cur_p.group('var_list'))) for cur_p in match_list]
        
        print arg_list
        
        indet_list = [cur_func(*cur_arg)
                      for cur_func, cur_arg in zip(func_list, arg_list)]
        
        return set(indet_list)
                
    def _not_exists(self, func_name):
        '''
        Determines, if the expression for the function was already found and
        assigned. If even symbolic function with this name does not exist, it
        defines it.
        
        @type  func_name: string
        @param func_name: Name of the function to check.
        
        @rtype : bool
        @return: True if only symbolic function exists, False if some
        expression is already assigned. 
        '''
        try:
            tmp = eval(func_name)
        except NameError:
            # print func_name
            function(func_name)
            return True
        
        if (tmp == func_name.__repr__()):
            return False
        else:
            return True
        
    def _extract_resolvent_params(self, func_name):
        '''
        Extracts information about resolvent from its name.
        
        @type  func_name: string
        @param func_name: Name of the resolvent.
        
        @rtype : dict
        @return: Dictoinary, containing following fields:
            - 'n_s_a': Number of z0 in resolvent.
            - 'n_d_a': Number of arguments, other than z0 in resolvent.
            - 'order': Which genus contribution it is?
        '''
        
        print 'func_name: ', func_name
        name_pattern = re.compile('rho_(?P<n_s_a>\d+)_' + 
                                      '(?P<n_d_a>\d+)_' +
                                      '(?P<order>\d+)[_(]'
                                      '(?P<vars>[a-zA-Z0-9_]+)')
        m = name_pattern.match(func_name)
        print m
        result = {'n_s_a': int(m.group('n_s_a')),
                  'n_d_a': int(m.group('n_d_a')),
                  'order': int(m.group('order'))}
                  
        tmp_var = var(re.sub('_', ',', m.group('vars')))
        if type(tmp_var) == tuple:
            pass
        else:
            tmp_var = (tmp_var,)
        result.update({'var_list': list(tmp_var)})
        
        return result

###############################################################################
# Testing part
###############################################################################

a = GaussianSolver()
var('z0 z1')
# function('rho_1_0_1', nargs=1)
# function('rho_1_1_0', nargs=2)
# function('drho_1_1_0', nargs=2)
# function('d2rho_1_1_0', nargs=2)
# function('drhodiff_0_2_0', nargs=2)

# b = [rho_1_1_0(z0, z1) == z0 + z1 ** 2]

# a._precompute_derivatives(b)
# a._set_initial_conditions()
# a._get_and_solve_eqns(-1,0)
a.solve(3,0)



