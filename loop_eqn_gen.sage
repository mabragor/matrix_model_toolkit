'''
In this module various methods concerning loop equations for matrix models are
contained.
'''

import re
from itertools import product

DEBUG = True

class LoopEquation (object):
    '''
    @type  beta: Symbol
    @ivar beta: the beta in the power of Wandermonde determinant.
    '''
    def __init__(self, order=-2, potential=x**2):
        self.beta = var('beta')
        self.order = order
        self.V = potential
    
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
        
        r = self._symb_func('r', self.order)
        rd = self._symb_func('rd', self.order)
        
        result = beta * r(num_same_points + 1,
                          num_diff_points)(*var_list)
        result += (beta - 1) * rd(num_same_points,
                                  num_diff_points)(*var_list)
        
        return result
    
    def _tail_part(self, num_same_points, num_diff_points):
        '''
        Part of loop equation, coming from differentiation of the 'tail' -
        product of 1/(z_i - lambda).
        
        For optimization, rdd and rd functions (partial derivatives w.r.t
        certain z are precomputed)
        
        @type  num_same_points: int
        @param num_same_points: Number of z0's in loop equation.
        
        @type  num_diff_points: int
        @param num_diff_points: Number of points in loop equation, other than
        z0.
        '''
        
        # Declare all variables.
        var_list = self._vars(range(num_diff_points + 1))
        
        # First, we add second derivatives in the coincident points.
        # Note, that here number of SAME arguments is smaller than in
        # _wand_part by one.
        if num_same_points == 1:
            result = 0
        else:
            rdd = self._symb_func('rdd', self.order)
            result = (1/2 * (num_same_points - 1) *
                      rdd(num_same_points - 1, num_diff_points)(*var_list))
    
        # And now we add derivatives in different points.
        # Peculiarity is: wether we should add 'z0' to the list of arguments.
        # Only add it, if it is in fact present (num_same_args > 1)              
        
        # First of all, if there are no diff_points, do nothing
        if num_diff_points == 0:
            return result
        
        # And else execute the following code.
        
        def zero_arg():
            if num_same_points > 1:
                return self._vars([0])
            else:
                return []
	
        for i in range(1, num_diff_points + 1):
            var_diff = zero_arg() + self._vars(range(1, num_diff_points + 1))
        
            var_same = (self._vars(range(0, i)) +
                        self._vars(range(i+1, num_diff_points +  1)))
        
            var_ord_changed = (zero_arg() + self._vars([i]) +
                               self._vars(range(1, i)) +
                               self._vars(range(i+1, num_diff_points + 1)))
            # print var_diff, var_same, var_ord_changed
        
            r = self._symb_func('r', self.order)
            rddiff = self._symb_func('rddiff', self.order)
            
            result -= (1 / (var_list[i] - var_list[0]) ** 2 * 
                       (r(num_same_points - 1, num_diff_points)(*var_diff) - 
                        r(num_same_points, num_diff_points - 1)(*var_same))
                       )
        
            result += (1 / (var_list[i] - var_list[0]) *
                       rddiff(num_same_points - 1, num_diff_points)
                             (*var_ord_changed)
                       )
    
        return result
        
        
    def _potential_part(self, num_same_points, num_diff_points):
        '''
        Part of loop equation, coming from differentiation of e^{-V/g}.
        
        @type  num_same_points: int
        @param num_same_points: Number of z0's in loop equation.
        
        @type  num_diff_points: int
        @param num_diff_points: Number of points in loop equation, other than
        z0.
        '''
        
        # Define all the variables resolvents can depend on.
        var_list = self._vars(range(num_diff_points + 1))
        
        # First, we add a dirty hack to handle only Gaussian potential
        # Will implement generic potential later.
        if bool(self.V == x^2):
            var('T2')
            
            z = var_list[0]
            r = self._symb_func('r', self.order + 1)
            result = - z * r(num_same_points, num_diff_points)(*var_list)
            
            
            # Add a contribution of particular (zero) correlator
            def zero_arg():
                if num_same_points > 1:
                    return self._vars([0])
                else:
                    return []
            
            # Correctly handle 1/g^2 term.
            if num_diff_points == 0 and num_same_points == 1:
                if self.order == -2:
                    multiplier = 1
                else:
                    multiplier = 0
            else:
                r = self._symb_func('r', self.order + 2)
                var_diff = zero_arg() + self._vars(range(1, num_diff_points + 1)
                                                   )
                multiplier = r(num_same_points - 1, num_diff_points)(*var_diff)
            
            var ('Lambda')
            result += Lambda * multiplier
            
            return result
        else:
            return 0 
    
    def debug(self, message):
        if DEBUG:
            print message
    
    def _symb_func(self, function_name, order):
        '''
        Returns symbol of the function, defined by funcname, if you wish
        to obtain the
        expression of disconnected resolvent in terms of connected ones, use
        '_r' method instead.
        
        Repeated argument at the beginning (z0) should be provided into
        resulting function only once. Multiplicity is encoded in the function
        name.
        
        @type  order: int
        @param order: The least power of string constant g in this function.
        If order is negative, then function name contains 'm' instead of '-'.
        
        EXAMPLE:
        r_i_j_k (z0, z1) =(in old notation)= r_(i+j)_k (z0,z0,..,z0, z1)
        '''
        def wrapper(num_same_points, num_diff_points):
        
            # Generate name and latex representation for the function.
            # Note the substitution of minus sign by 'm' in function name.
            func_name = '%s_%d_%d_' % (function_name, num_same_points,
                                         num_diff_points)
            def ord_name():
                if order < 0:
                    return 'm' + str(abs(order))
                else:
                    return str(order)
            
            func_name += ord_name()
        
            latex_repr = '%s_{%d,%d,%d}' % (function_name, num_same_points,
                                            num_diff_points, order)
            
            # self.debug(func_name)
        
            num_points = num_same_points + num_diff_points
        
            def _lambda ():
                if num_same_points == 0:
                    return 0
                else:
                    return 1
                    
            tmp_func = function(func_name, nargs=num_diff_points + _lambda(),
                                latex_name=latex_repr)
            
            return tmp_func
        return wrapper

    def _vars(self, num_list):
        '''
        Generates list of symbolic variables from list of integers by following
        rule:
            
        [i, j, ...] -> [zi, zj, ...].
            
        @type  num_list: [int, int, ...]
        @param num_list: List of integers
        @rtype: [Expression, Expression, ...]
        @return: List of (newly declared, if necessary) symbolic variables.
        '''
        
        name_list = map(lambda i: 'z{0}'.format(i), num_list)
        var_list = map(var, name_list)
        return var_list
    
    def loop_equation(self, num_same_points, num_diff_points):
        '''
        Generates loop equation for given number of points and specified order.
        Order is specified in 'order' instance attribute.
        
        'equation' instance attribute is updated on successs, containing newly
        obtained loop equation.
        
        @type  num_same_points: int
        @param num_same_points: Number of equal points in the equation. Should
        be greater of equal to 1.
        
        @type  num_diff_points: int
        @param num_diff_points: Number of mutually distint points in the
        equation. They are also different from 'z0', quantitiy of which is
        specified by 'num_same_points' parameter. 
        
        @rtype: Expression
        @return: Loop equation in terms of disconnected resolvents. To get to
        the view in terms of connected ones, use 'to_connected' method.
        '''
        self.equation = (self._wand_part(num_same_points, num_diff_points) +
                         self._tail_part(num_same_points, num_diff_points) +
                         self._potential_part(num_same_points, num_diff_points))
        
        return self.equation
                
    def to_connected(self, equation=None):
        '''
        Rewrites specified loop equation in terms of connected resolvents,
        instead of disconnected ones. On success, equation_conn attribute is
        updated.
        
        @type  equation: Expression
        @param equation: Loop equation to rewrite. If not specified,
        self.equation attribute is used.
        
        @rtype: Expression
        @return: Loop equation in terms of connected resolvents.
        '''

        # If no equation is provided 'equation' instance attribute is used.
        # Note that we convert to str.
        if equation == None:
            equation = self.equation._repr_()
        else:
            equation = equation._repr_()

        # First we define regular expression, that matches all names of
        # disconnected resolvents and their derivatives.
        funcname_re = re.compile('(?P<func_name>r[a-zA-Z0-9]*)_' +
                                 '(?P<num_same_args>\d+)_' + 
                                 '(?P<num_diff_args>\d+)_' +
                                 '(?P<order>m?\d+)')
        
        # Then we define, with what to replace every function name
        # task is to write '_' before this name.
        replace_pattern = ('self._\g<func_name>(' + 
                           '\g<num_same_args>,' + 
                           '\g<num_diff_args>,' +
                           '\'\g<order>\')')
                           
        # And finally perform the replacement and evaluate resulting sage
        # command.
        sage_command = funcname_re.sub(replace_pattern, equation)
        # print sage_command
        self.equation_conn = eval(preparse(sage_command))

        return self.equation_conn

    @staticmethod
    def beads_into_glasses(num_args, num_groups):
        '''
        Place num_args distinct beads into num_groups equivalent glasses.
        
        @type  num_args: int
        @param num_args: Number of integers (distinct beads) to divide into
        groups.
        
        @type  num_groups: int
        @param num_groups: Number of groups (equivalent glasses)
        to divide into.
        
        @rtype: set
        @return: Set of frozensets of tuples of numbers, where:
            - tuple describes ID numbers of beads in the given glass;
            - frozenset describes any given partition of beads;
            - set enumerates all possible partitions.
        '''
        
        # Generate masks, that describe, which integer goes where.
        pre_mask_list = list(product(range(num_groups), repeat=num_args))
        
        mask_list = []
        for i in xrange(len(pre_mask_list)):
            # Some groups have zero beads in them, this is not a mask.
            if set(pre_mask_list[i]) != set(range(num_groups)):
                pass
            else:
                mask_list.append(pre_mask_list[i])
        
        # Then for each generated mask we partition our integer list.
        # We do not add duplicates.
        int_partition_set = set([])
        for cur_mask in mask_list:
            cur_partition = [[] for i in xrange(num_groups)]
            
            for i in xrange(num_args):
                cur_partition[cur_mask[i]].append(i+1)
        
            # Only hashable types (such as tuple) can be members of sets,
            # so we move to them.
            cur_partition = map(lambda x: tuple(x), cur_partition)
            cur_partition = frozenset(cur_partition)
            
            int_partition_set.add(cur_partition)
            
        return int_partition_set
        

    def _r(self, num_same_args, num_diff_args, str_order):
        '''
        This is a wrapper for abstract L{_conn_resolv} method, that serves to
        express disconnected resolvent in terms of connected ones.
        '''
        return self._conn_resolv('rho', None, None,
                                 num_same_args, num_diff_args, str_order)
                                 
    def _rd(self, num_same_args, num_diff_args, str_order):
        '''
        This is a wrapper for abstract L{_conn_resolv} method, that serves to
        express first derivative of disconnected resolvent in coincident point
        in terms of connected ones.
        '''
        return self._conn_resolv('rho', 'drho', num_same_args, 
                                 num_same_args, num_diff_args, str_order)
    
    def _rdd(self, num_same_args, num_diff_args, str_order):
        '''
        This is a wrapper for abstract L{_conn_resolv} method, that serves to
        express second derivative of disconnected resolvent in coincident point
        in terms of connected ones.
        '''
        return self._conn_resolv('rho', 'd2rho', num_same_args,
                                 num_same_args, num_diff_args, str_order)
    
    def _rddiff(self, num_same_args, num_diff_args, str_order):
        '''
        This is a wrapper for abstract L{_conn_resolv} method, that serves to
        express first derivative of disconnected resolvent in first
        non-coincident point in terms of connected ones.
        '''
        return self._conn_resolv('rho', 'drhodiff', num_same_args + 1,
                                 num_same_args, num_diff_args, str_order)
    
    
    def _conn_resolv(self, base_func_name, mod_func_name, magic_number,
                     num_same_args, num_diff_args, str_order):
        '''
        Expression for disconnected resolvent is terms of connected ones.
        Already known connected resolvents are substituted by their value,
        while not known are left unevaluated.
        
        @type  base_func_name: string
        @param base_func_name: Name to use for all functions, which do not
        happen to be called upon argument number 'magic_number'.
        
        @type  mod_func_name: string
        @param mod_func_name: Name to use for the function, that is called
        upon argument number 'magic_number'.
        
        @type  magic_number: int
        @param magic_number: Number of argument, that alteres function name.
        Usually, this means, that differentiation is performed with respect to
        it. This should be manifested in change of function name.
        
        @type  num_same_args: int
        @param num_same_args: Number of z0 arguments in the resolvent.
        
        @type  num_diff_args: int
        @param num_diff_args: Number of arguments, different from z0 and
        mutually distinct.
        
        @type  str_order: string
        @param str_order: String that describes order of the resolvent.
        Minus sign is encoded by 'm' in the beginning of the string. Apart
        from that, str_order is the string representation of the absolute
        value of order.
        
        @rtype: Expression
        @return: Symbolic expression, that expresses disconnected resolvent in
        terms of connected ones.
        '''
        
        # First we recover order of disconnected resolvent.
        if str_order[0] == 'm':
            order = -int(str_order[1:])
        else:
            order = int(str_order) 
        
        # Total number  of arguments of the resolvent is equal to.
        num_args = num_same_args + num_diff_args
        
        # Then we generate list of all possible partitions of arguments of
        # the resolvent.
        partition_list = [self.beads_into_glasses(num_args, num_groups) 
                          for num_groups in xrange(1, num_args + 1)]
        
        # print partition_list
        
        def wrapper(*points):
            '''
            Splits 'points' into tuples according to just generated
            partition_list.
            '''
            
            result = 0
            for i in xrange(num_args):
                # Total order of connected resolvents is equal to.
                num_groups = i + 1
                total_order = order - (2 * (num_args - num_groups) - num_args)
                
                # If total order is negative, this partitioning does not
                # contribute
                if total_order < 0:
                    continue
                
                # All possible distributions of orders among resolvents.
                order_mask = self._order_mask(total_order, num_groups)
                
                for cur_partition, cur_ord_mask in product(partition_list[i],
                                                           order_mask):
                    # Split arguments into groups, corresponding to
                    # particular connected resolvents.
                    split_args = self._split_args(points, num_same_args,
                                                  num_diff_args,
                                                  cur_partition,
                                                  magic_number)
                    
                    # Generate list of all funcnames.
                    f = self._generate_func_names
                    func_names = f(base_func_name, split_args, cur_ord_mask, 
                                   mod_func_name)
                    
                    
                    # Finally, calculate functions on their arguments.
                    # And add result to overall result.         
                    only_args = [elem['points'] for elem in split_args]
                    ready_to_apply = zip(func_names, only_args)
                    
                    tmp_result = mul([tmp[0](*(tmp[1]))
                                     for tmp in ready_to_apply])
                                     
                    result += tmp_result
            
            return result
                    
        return wrapper

    def _split_args(self, points, n_s_a, n_d_a, partition, magic_number):
        '''
        @type  points: tuple
        @param points: Disconnected resolvent argument list to distirbute over
        connected resolvents.
        
        @type  n_s_a: int
        @param n_s_a: Multiplicity of z0 argument.
        
        @type  n_d_a: int
        @param n_d_a: Number of arguments other than z0.
        
        @type  partition: frozenset(tuples())
        @param partition: How arguments are partitioned.
        '''
        def get_arg(arg_num):
            '''
            In 'points' tuple first argument actually represents multiple
            first arguments, where multiplicity is stated by num_same_args
            variable. This function takes this peculiarity into account.
            
            @type  arg_num: int
            @param arg_num: Argument number, if we write arguments in 'flat 
            form' - 'z0, z0, z0, z1,., zn'
            
            @rtype: Variable
            @return: Correct argument, taking into account multiplicity of
            the first argument.
            '''
            # print points, arg_num, n_s_a
                
            if arg_num < n_s_a + 1:
                return points[0]
            else:
                # n_s_a equal to zero is the special case.
                return points[arg_num - n_s_a + (-1) * int(n_s_a == 0)]
            
        # First divide all arguments flatly.
        list_vars = []
        # We must mark, what partition has 'magic' argument, in which
        # differentiation is performed.
        has_magic_number_array = []
        for little_arg_list in partition:
            flag = (magic_number != None and magic_number in little_arg_list)
            if flag:
                has_magic_number_array.append(True)
            else:
                has_magic_number_array.append(False)
            
            # Generate flat list of arguments for particular resolvent.
            tmp_list = [get_arg(i) for i in little_arg_list]
            list_vars.append(tmp_list)
            
        # print list_vars, has_magic_number_array
        
        # Then in each group assemble first arguments, which are equal to each
        # other and calculate corresponding n_s_a and n_d_a.
        # Write the result as a dictionary.
        list_vars = [self._flat_list_to_compressed_dict(x) for x in list_vars]
        for i in xrange(len(has_magic_number_array)):
            list_vars[i].update({'has_magic_number': has_magic_number_array[i]})
        
        # print list_vars
        
        return list_vars
        
        
    def _flat_list_to_compressed_dict(self, flat_list):
        '''
        Converts 'flat' representation of the argument list to the compressed
        form.
        
        @type  flat_list: list
        @param flat_list: List of arguments in the 'flat' form, that is
        [z0, z0, z0,., z1, z2, ...]
        
        @rtype: dict
        @return: A dictionary, containing following information
            - 'n_s_a': Number of z0's in arglist;
            - 'n_d_a': Number of all other arguments;
            - 'points': List of all arguments, where z0, if present at all, is
            placed at the beginning of the list 1 time.
        '''
        zero_arg = self._vars([0])[0]
        n_s_a = flat_list.count(zero_arg)
        n_d_a = len(flat_list) - n_s_a
        
        compressed_arg_list = [zero_arg] * int(n_s_a != 0) + flat_list[n_s_a:]
        
        result = {'n_s_a': n_s_a, 'n_d_a': n_d_a, 'points': compressed_arg_list}
        
        return result

    @staticmethod
    def _generate_func_names(base_func_name, split_args, ord_mask, 
                             mod_func_name=None):
        '''
        Ensures that all function names to appear in equation, do exist.
        
        For purposes of incorporating derivatives in coincident points and
        derivatives in general, three optional arguments can be specified.
        If given, they determine, which resolvent name should be altered, and
        how.
        
        @type  base_func_name: string
        @param base_func_name: Template to use for start of almost all
        function names.
        
        @type  split_args: list[dict]
        @param split_args: Argument lists for individual resolvents in
        compressed form.
        
        @type  mod_func_name: string
        @param mod_func_name: Optional. What name to use, if function is called
        on an argument with number 'magic_number'
        
        @type  ord_mask: list(int)
        @param ord_mask: What are the orders of connected resolvents?
        
        @rtype: list(functions)
        @return: List of functions, ready to be called on corresponding
        argument lists. If functions were not evaluated before, they are only
        symbolic names, otherwise, expressions would be obtained upon calling.
        
        @attention: Any correctness checks are OMITTED!!! Use at your own risk.
        '''
        
        funcs = []
        
        for i in xrange(len(split_args)):
            if split_args[i]['has_magic_number']:
                pre_func_name = mod_func_name
            else:
                pre_func_name = base_func_name
                
            func_name = '%s_%d_%d_%d' % (pre_func_name,
                                             split_args[i]['n_s_a'],
                                             split_args[i]['n_d_a'],
                                             ord_mask[i])
                                             
            try:
                # Check, if function exists.
                cur_f = eval(func_name)
            
            except NameError:
                # If not defined, reserve symbolic name.
                cur_f = function(func_name, nargs = len(split_args[i]['points']))
            
            finally:
                funcs.append(cur_f)
        
        return funcs
        
        
    @classmethod
    def _order_mask(cls, total_order, num_groups):
        '''
        Generates all possible distributions of orders among
        connected resolvents.
        
        @type  total_order: int
        @param total_order: Sum of all orders of the resolvents.
        
        @type  num_groups: int
        @param num_groups: How many resolvents there are?
        
        @rtype: list
        @return: All possible distributions (each order can be greater or
        equal than zero) as list of lists.
        '''
        
        # If division is not required, just return in the specified required 
        # format.
        if num_groups == 1:
            return [[total_order]]
        
        result  = []
        for i in xrange(total_order + 1):
            # Generate all distributions on smaller number of groups.
            tmp_list = cls._order_mask(total_order - i, num_groups - 1)
            
            # Append to each group kept number of units.
            tmp_list = map(lambda x: x + [i], tmp_list)
            
            # Add everything to the resuling list.
            result.extend(tmp_list)
            
        return result

    @staticmethod
    def _if_exists_if_not_define(func_name, num_same_args,
                                 num_diff_args, order):
        
        # Get the string representation of order
        str_order = 'm' * ((1 + sgn(order)) / 2) + str(abs(order))

        compiled_func_name = '%s_%d_%d_%s' % (func_name,
                                              num_same_args,
                                              num_diff_args,
                                              str_order)

        try:
            type(eval(compiled_func_name))
        except NameError:
            function(compiled_func_name, nargs=num_sage_args + num_diff_args)



