'''
In this module various methods concerning loop equations for matrix models are
contained.
'''

DEBUG = True

class LoopEquation (object):
    '''
    @type  beta: Symbol
    @ivar beta: the beta in the power of Wandermonde determinant.
    '''
    def __init__(self, order=0):
        self.beta = var('beta')
        self.order = order
    
    def _wand_part(self, num_same_points, num_diff_points):
        '''
        Part of loop equation, coming from differentiation of Wandermonde
        determinant.
        
        @type  num_same_points: int
        @param num_same_points: Number of z0's in general term of
            loop equation.
        @type  num_diff_points: int
        @param num_diff_points: Number of zi's for i > 0 in general term of
            loop equation.
        '''
        
        var_list = self._vars(range(num_diff_points + 1))
        # print var_list, var_list_en
        
        # num_points = num_same_points + num_diff_points
        
        result = beta * self.r(num_same_points + 1,
                                num_diff_points)(*var_list)
        result += (1 - beta) * self.rd(num_same_points,
                                        num_diff_points)(*var_list)
        
        return result
    
    def _tail_part(self, num_same_points, num_diff_points):
        '''
        Part of loop equation, coming from differentiation of the 'tail' -
        product of 1/(z_i - lambda).
        
        For optimization, rdd and rd functions (partial derivatives w.r.t
        certain z are precomputed)
        
        @type  num_points: int
        @param num_points: Number of z's in loop equation.
        '''
        
        if num_same_points == 1:
            result = 0
        else:
            var_list = self._vars(range(num_diff_points + 1))
            result = (1/2 * (num_same_points - 1) *
                      self.rdd(num_same_points, num_diff_points)(*var_list))
                      
        r = self.r(num_same_points, num_diff_points)
        for i in range(1, num_diff_points + 1):
            result += (r(var_list)    
        
        
    def _potential_part(self, num_points):
        '''
        Part of loop equation, coming from differentiation of e^{-V/g}.
        
        @type  num_points: int
        @param num_points: Number of z's in loop equation.
        '''
    
    def debug(self, message):
        if DEBUG:
            print message
    
    def r(self, num_same_points, num_diff_points):
        '''
        Disconnected n-point resolvent of order 'order' in h.
        
        For now returns symbolic function.
        
        Repeated argument at the beginning (z0) should be provided into
        resulting function only once. Multiplicity is encoded in the function
        name.
        
        EXAMPLE:
        r_i_j_k (z0, z1) =(in old notation)= r_(i+j)_k (z0,z0,..,z0, z1)
        '''
        
        func_name = 'r_%d_%d_%d' % (num_same_points, num_diff_points,
                                    self.order)
        
        self.debug(func_name)
        
        num_points = num_same_points + num_diff_points
        
        tmp_func = function(func_name, nargs=num_diff_points + 1,
                            latex_name='r_{%d,%d,%d}' % (num_same_points,
                                                         num_diff_points,
                                                         self.order))
        
        return tmp_func
    
    def rd(self, num_same_points, num_diff_points):
        '''
        First derivative w.r.t z0 of disconnected n-point resolvent of order
        'order' in h.
        
        For now returns symbolic function.
        
        Repeated argument at the beginning (z0) should be provided into
        resulting function only once. Multiplicity is encoded in the function
        name.
        
        EXAMPLE:
        rd_i_j_k (z0, z1) =(in old notation)= rd_(i+j)_k (z0,z0,..,z0, z1)
        '''
        
        func_name = 'rd_%d_%d_%d' % (num_same_points, num_diff_points,
                                    self.order)
        
        self.debug(func_name)
        
        num_points = num_same_points + num_diff_points
        
        tmp_func = function(func_name, nargs=num_diff_points + 1,
                            latex_name="r'_{%d,%d,%d}" % (num_same_points,
                                                         num_diff_points,
                                                         self.order))
        
        return tmp_func
    
    def rdd(self, num_same_points, num_diff_points):
        '''
        Second derivative w.r.t z0 of disconnected n-point resolvent of order
        'order' in h.
        
        For now returns symbolic function.
        
        Repeated argument at the beginning (z0) should be provided into
        resulting function only once. Multiplicity is encoded in the function
        name.
        
        EXAMPLE:
        rdd_i_j_k (z0, z1) =(in old notation)= rdd_(i+j)_k (z0,z0,..,z0, z1)
        '''
        
        func_name = 'rdd_%d_%d_%d' % (num_same_points, num_diff_points,
                                    self.order)
        
        self.debug(func_name)
        
        num_points = num_same_points + num_diff_points
        
        tmp_func = function(func_name, nargs=num_diff_points + 1,
                            latex_name="rdd'_{%d,%d,%d}" % (num_same_points,
                                                         num_diff_points,
                                                         self.order))
        
        return tmp_func
    
    
    # def rd(self, 
    
    def _vars(self, num_list):
        '''
        Generates list of symbolic variables from list of integers by following
        rule:
        
            [i, j, ...] -> [zi, zj, ...].
            
        @type  num_list: [int, int, ...]
        @param num_list: List of integers
        @rtype : [Expression, Expression, ...]
        @return: List of (newly declared, if necessary) symbolic variables.
        '''
        
        name_list = map(lambda i: 'z{0}'.format(i), num_list)
        var_list = map(var, name_list)
        return var_list
    
    