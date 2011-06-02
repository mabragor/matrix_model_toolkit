load 'loop_eqn_gen.sage'

import re

class GaussianSolver(object):
    def _set_initial_conditions(self):
        '''
        Set initial conditions for iterative procedure, like rho_{1,0,0} and
        its derivatives.
        '''
        # Define basic variables to work with
        var ('z0 beta')
        
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
        ans0 = [rho_1_0_0(x) == (x - y(x)) / 2 / beta,
                rho_0_1_0(x) == (x - y(x)) / 2 / beta]
        
    
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
        for eqn in eqn_list:
            # Instances of functions, assigned in one round.
            local_result = []
            
            print eqn.lhs(), eqn.lhs().operator()
            op = eqn.lhs().operator()
            res_params = self._extract_resolvent_params(str(op))
            # These will be central objects, so we proxy them.
            n_s_a, n_d_a, order = (res_params['n_s_a'], res_params['n_d_a'],
                                   res_params['order'])
            ops = eqn.lhs().operands()
            
            print 'I"m here'
            print 'I"m here'
            
            globals().update({str(eqn.lhs().operator()): eqn.rhs().function(*ops)})           
            
            print 'I"m here'
            
            f = eval('rho_%d_%d_%d' % (n_s_a, n_d_a, order))

            print f(z1, z0)
            local_result.append(f)
                        
            # If n_s_a == 1, regardless of the circumstances, we can
            # calculate drho_1_k_ and d2rho_1_k_ and drhodiff_0_1+k_
            # This IF block is written in non-object-oriented fashion.
            # Be careful!!!
            if n_s_a == 1:
                if self._not_exists('drho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                    print 'Calculate drho'
                    # Calculate derivative of the function w.r.t z0
                    tmp_rhs = diff(eqn.rhs(), ops[0])
                    print tmp_rhs, tmp_rhs.function(*ops)
                    
                    # Define new function
                    globals().update({'drho_%d_%d_%d' % (n_s_a, n_d_a, order):
                                      tmp_rhs.function(*ops)})           
            
                    # Add a link to the instance of newly created function
                    # to our array.
                    tmp_func = eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    local_result.append(tmp_func)
                               
                if self._not_exists('d2rho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                    print 'Calculate d2rho'
                    # Calculate second derivative.
                    tmp_func = eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))
                    tmp_rhs = diff(tmp_func(*ops), ops[0])
                    
                    #Define new function.
                    globals().update({'d2rho_%d_%d_%d' % (n_s_a, n_d_a, order): tmp_rhs.function(*ops)})           
                    
                    # Add a link to the instance of newly created function
                    # to our array.
                    tmp_func = eval('d2rho_%d_%d_%d' % 
                                    (n_s_a, n_d_a, order))
                    local_result.append(tmp_func)
                        
                if self._not_exists('drhodiff_%d_%d_%d' % (n_s_a - 1, n_d_a + 1,
                                                           order)):
                    
                    print 'Calculate drhodiff'
                    globals().update({'drhodiff_%d_%d_%d' % (n_s_a - 1, n_d_a + 1, order):
                                        eval('drho_%d_%d_%d' % (n_s_a, n_d_a, order))})
                    
                    tmp_func = eval('drhodiff_%d_%d_%d' % (n_s_a - 1, n_d_a + 1, order))
                    
                    local_result.append(tmp_func)  
            '''
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
            
            if n_d_a > 1:
                if self._not_exists('drhodiff_%d_%d_%d' %
                                    (n_s_a, n_d_a, order)):
                    tmp_func = eval('drhodiff_%d_%d_%d' %
                                    (n_s_a, n_d_a, order))
                    tmp_rhs = diff(eqn.rhs(), correct_op())
                    
                    tmp_eqn = tmp_func(*ops) == tmp_rhs
                    
                    local_result.append(tmp_eqn)
                else:
                    # Added for combatibility with the following code.
                
                    if self._not_exists('drho_%d_%d_%d' % 
                                        (n_s_a + 1, n_d_a - 1, order)):
                        tmp_func = eval('drho_%d_%d_%d' % 
                                        (n_s_a + 1, n_d_a - 1, order))
                        # Variable name is hardcoded. Do something.
                        tmp_rhs = limit(local_result[-1].rhs(),
                                        correct_op() == z0)
                        if ops[0] == z0:
                            arg_list = ops
                        else:
                            arg_list = [z0] + ops
                        
                        tmp_eqn = tmp_func(*arg_list) == tmp_rhs
                        
                        local_result.append                
                        
                    
                    if self._not_exists('d2rhodiff_%d_%d_%d' % 
                                        (n_s_a + 1, n_d_a - 1, order)):
                        if self._not_exists('d2rhodiff_%d_%d_%d' % 
                                            (n_s_a, n_d_a, order)):
               ''' 
            
                
             
    def _not_exists(self, func_name):
        try:
            tmp = eval(func_name)
        except NameError:
            return False
        
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
        
        name_pattern = re.compile('rho_(?P<n_s_a>\d+)_' + 
                                      '(?P<n_d_a>\d+)_' +
                                      '(?P<order>\d+)')
        m = name_pattern.match(func_name)
        result = {'n_s_a': int(m.group('n_s_a')),
                  'n_d_a': int(m.group('n_d_a')),
                  'order': int(m.group('order'))}
        
        return result

###############################################################################
# Testing part
###############################################################################

a = GaussianSolver()
var('z0 z1')
function('rho_1_1_0', nargs=2)
function('drho_1_1_0', nargs=2)
function('d2rho_1_1_0', nargs=2)
function('drhodiff_0_2_0', nargs=2)

b = [rho_1_1_0(z0, z1) == z0 + z1 ** 2]

a._precompute_derivatives(b)


