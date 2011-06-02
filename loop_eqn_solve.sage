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
        
        This procedure computes drhodiff, drho and d2rho for any rho
        
        @type  eqn_list: list
        @param eqn_list: List of equation, defining newly computed resolvents.
        
        @rtype : list
        @return: List of newly computed resolvents, populated with their
        derivatives in coincident points.
        '''
        result = []
        
        for eqn in eqn_list:
            local_result = []
            
            local_result.append(eqn)
            
            print eqn.lhs(), eqn.lhs().operator()
            res_params = self._extract_resolvent_params(eqn.lhs().operator())
            # These will be central objects, so we proxy them.
            n_s_a, n_d_a, order = (res_params['n_s_a'], res_params['n_d_a'],
                                   res_params['order'])
            
            _vars = eqn.lhs().operands()
            
            
            # If n_s_a == 1, regardless of the circumstances, we can
            # calculate drho_1_k_ and d2rho_1_k_ and drhodiff_0_1+k_
            if n_s_a == 1:
                if self._not_exists('drho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                    # Common logic for all such code pieces.
                    # Calculate derivative, then assemble string that
                    # represents symbolic equation, then eval it.
                    
                    tmp_rhs = diff(eqn.rhs(), z0)._repr_()
                    
                    tmp_eqn = 'drho_%d_%d_%d' % (n_s_a, n_d_a, order) + 
                              ' == ' + tmp_rhs
                               
                if self._not_exists('d2rho_%d_%d_%d' % (n_s_a, n_d_a, order)):
                if self._not_exists('drhodiff_%d_%d_%d' % (n_s_a - 1, n_d_a + 1,
                                                           order)):  
            
            # If n_d_a > 0, we can calculate
            # drho_i+1_k-1_, d2rho_i+1_k-1_ and drhodiff_i_k_
            if n_d_a > 1:
                if self._not_exists('drhodiff_%d_%d_%d' %
                                    (n_s_a, n_d_a, order)):
                    if self._not_exists('drho_%d_%d_%d' % 
                                        (n_s_a + 1, n_d_a - 1, order)):
                
                    if self._not_exists('d2rhodiff_%d_%d_%d' % 
                                        (n_s_a + 1, n_d_a - 1, order)):
                        if self._not_exists('d2rhodiff_%d_%d_%d' % 
                                            (n_s_a, n_d_a, order)):
                
            
            # General case
            if n_s_a > 0 and n_d_a > 0:
            elif n_s_a > 0 and n_d_a == 0:
                if n_s_a == 1:
                    # We still may do something
                else:
                    # We do nothing.
                    pass
            elif n_s_a == 0 and n_d_a > 0:
                # We can calculate derivatives, as usual.
            elif n_s_a == 0 and n_d_a == 0:
                # This situation just should not occur - something strange
                # happened.    
             
    def _not_exists(func_name):
        try:
            eval(func_name)
        except NameError:
            return True
        else:
            return False
        
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
        
a = LoopEquation()
a.loop_equation(1,0)
b = a.to_connected(); b