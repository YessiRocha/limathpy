from sympy import symbols, Function, Eq, Derivative, dsolve, solve, init_printing, exp, sin, cos, tan, plot_parametric
t = symbols('t')
y = Function('y')
x = Function('x')
C1 = symbols('C1')
C2 = symbols('C2')

#Solving ordinary differential equations.


#Solving first order ordinary differential equations of the form p(t)y'(t) + q(t)y(t) = g(t).

def first_ode(vector): 
    """A function that returns the general solution of a first order differential equation of the form p(t)y'(t) + q(t)y(t) = g(t);
       p, q and g are functions which depend on t.
    Args: 
        vector (list): a list of the form [p(t), q(t), g(t)].
    Returns:
        Eq: the general solution of the equation, C1 is a constant which depends on some initial condition.
    Example: 
        >>> from limathpy import first_ode
        >>> first_ode([t, 2, 2 + t])
        Eq(y(t), C1/t**2 + t/3 + 1)"""
    eq = Eq(vector[0]*Derivative(y(t), t) + vector[1]*y(t), vector[2])
    sol = dsolve(eq, y(t))
    return sol


def solve_first_ode(vector, init_cond = [1, 0]): 
    """A function that returns the solution of a first order differential equation of the form p(t)y'(t) + q(t)y(t) = g(t) 
       for some given initial conditions; p, q and g are functions which depend on t.
    Args: 
        vector (list): a list of the form [p(t), q(t), g(t)].
        init_cond (list): a list of the form [t1, y(t1)], for some t1. Defaults to y(t1 = 1) = 0.
    Returns:
        Eq: the solution of the equation depending on the initial condition.
    Example: 
        >>> from limathpy import solve_first_ode
        >>> solve_first_ode([t, 2, 2 + t], [1, 1])
        Eq(y(t), t/3 + 1 - 4/(3*t**2))"""
    equation = first_ode(vector).rhs
    evalu = equation.subs({t: init_cond[0]})
    lin = Eq(evalu, 0)
    solu = solve(lin, C1)
    final = equation.subs({C1: solu[0]})
    return Eq(y(t), final)


#Solving first order ordinary differential equations of the form r(t)y''(t) + p(t)y'(t) + q(t)y(t) = g(t).

def second_ode_const(vector): 
    """A function that returns the solution of a second order differential equation of the form r(t)y''(t) + p(t)y'(t) + q(t)y(t) = g(t); 
       r, p, q and g are functions which depend on t.
    Args: 
        vector (list): a list of the form [r(t), p(t), q(t), g(t)].
    Returns:
        Eq: the general solution of the equation, C1 and C2 are constants which depend on some initial conditions.
    Example: 
        >>> from limathpy import second_ode_const
        >>> second_ode_const([t**2, 2*t, 0, 1])
        Eq(y(t), C1 + C2/t + log(t))"""
    eq = Eq(vector[0]*Derivative(y(t), t, 2) + vector[1]*Derivative(y(t), t) + vector[2]*y(t), vector[3])
    sol = dsolve(eq, y(t))
    return sol


def solve_2nd_ode(vector, init_cond = [[1, 2], [0, 1]]): #y(1)=0, y(0)=1 vector=(y'', y', y, g(t)
    """A function that returns the solution of a second order differential equation of the form r(t)y''(t) + p(t)y'(t) + q(t)y(t) = g(t) 
       for some given initial conditions; r, p, q and g are functions which depend on t.
    Args: 
        vector (list): a list of the form [r(t), p(t), q(t), g(t)].
        init_cond (list of two lists): a matrix of the form [[t1, t2], [y(t1), y(t2)]], for some t1, t2. Defaults to y(t1 = 1) = 0 and y(t2 = 2) = 1.
    Returns:
        Eq: the general solution of the equation, C1 and C2 are constants which depend on some initial conditions.
    Example: 
        >>> from limathpy import solve_2nd_ode
        >>> solve_2nd_ode([t**2, 2*t, 0, 1], [[1, 0], [2, 0]])
        Eq(y(t), log(t) - 2*log(2) + 2*log(2)/t)"""
    equation = second_ode_const(vector).rhs
    evalu1 = equation.subs({t: init_cond[0][0]})
    evalu2 = equation.subs({t: init_cond[1][0]})
    eq1 = Eq(evalu1, init_cond[0][1])
    eq2 = Eq(evalu2, init_cond[1][1])
    const = solve((eq1, eq2))
    final = equation.subs(const)
    return Eq(y(t), final)


#Solving ordinary differential equations systems.

def system_ode(matrix):
    """A function that, given a 2x2 matrix (list of two lists), returns the general sotutions of the associated ordinary differential equations system.
    Args: 
        matrix (list of two lists): a list of two lists of the form [[t1, t2], [t3, t3]], 
                                    where you obtain the following system x'(t) = t1*x(t) + t2*y(t); y'(t) = t3*x(t) + t4*y(t). 
    Returns:
        tuple: a tuple of the form (x(t), y(t)), with x(t) and y(t) the general solutions of the system. C1 and C2 are constants that depend on some initial conditions.
    Example:
    >>> from limathpy import system_ode
    >>> system_ode([[1, 0], [0, -3]])
    (C1*exp(t), C2*exp(-3*t))
    """
    t = symbols('t')
    x, y = symbols('x y', cls = Function)
    eq1 = Eq(Derivative(x(t), t), matrix[0][0]*x(t) + matrix[0][1]*y(t))
    eq2 = Eq(Derivative(y(t), t), matrix[1][0]*x(t) + matrix[1][1]*y(t))
    sols = dsolve((eq1, eq2))
    return sols[0].rhs, sols[1].rhs


def lin_system(matrix, init_cond = [[1, 1], [0, 1]]): 
    """A function that, given an ordinary differential equations system, returns the linear system for some given initial conditions. 
    Args: 
        matrix (list of two lists): a list of two lists of the form [[t1, t2], [t3, t3]], 
                                    where you obtain the following system x'(t) = t1*x(t) + t2*y(t); y'(t) = t3*x(t) + t4*y(t). 
        init_cond (list of two lists): a matrix of the form [[t1, t2], [x(t1), y(t2)]], for some t1, t2. Defaults to x(t1 = 1) = 0 and x(t2 = 1) = 1.
    Returns:
        tuple: a tuple of the form (Eq_1, Eq_2), with Eq_1 and Eq_2 equations (Eq) with C1 and C2 as variables that are to be found using the initial conditions. 
    Example:
    >>> from limathpy import lin_system
    >>> lin_system([[1, 1], [0, -3]], [[0, 0], [0, 1]])
    (Eq(-C1/4 + C2, 0), Eq(C1, 1))"""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    sols = system_ode(matrix)
    lin1 = Eq(sols[0].subs({t: init_cond[0][0]}), init_cond[1][0])
    lin2 = Eq(sols[1].subs({t: init_cond[0][1]}), init_cond[1][1])
    return lin1, lin2


def solve_system_ed(matrix, init_cond = [[1, 1], [0, 1]]):
    """A function that, given a 2x2 matrix (list of two lists), returns the sotutions of the associated ordinary differential equations system 
    for some given initial conditions.
    Args: 
        matrix (list of two lists): a list of two lists of the form [[t1, t2], [t3, t3]], 
                                    where you obtain the following system x'(t) = t1*x(t) + t2*y(t); y'(t) = t3*x(t) + t4*y(t). 
        init_cond (list of two lists): a matrix of the form [[t1, t2], [y(t1), y(t2)]], for some t1, t2. Defaults to y(t1 = 1) = 0 and y(t2 = 2) = 1.                         
    Returns:
        tuple: a tuple of the form (x(t), y(t)), with x(t) and y(t) the general solutions of the system. C1 and C2 are constants that depend on some initial conditions.
    Example:
    >>> from limathpy import solve_system_ed
    >>> solve_system_ed([[1, 1], [0, -3]], [[0, 0], [0, 1]])
    (exp(t)/4 - exp(-3*t)/4, exp(-3*t))    
    """
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    C1, C2 = symbols('C1 C2')
    sys_ed = system_ode(matrix)
    sys_lin = lin_system(matrix, init_cond)
    dict_sols = solve(sys_lin)
    expr1 = sys_ed[0].subs(dict_sols)
    expr2 = sys_ed[1].subs(dict_sols)
    return expr1, expr2

#Retrato de fase.

def phase_portrait(matriz, lim_initialconditions=2):
    system = system_ode(matriz)
    p = plot_parametric((0,0), (t, 0, 0), show=False, title = 'Phase portrait')
    for i in range(0, lim_initialconditions):
        for j in range(0, lim_initialconditions):
            const = lin_system(matriz, [i, j])
            expr1 = system[0].subs({C1: const[1].rhs, C2: const[0].rhs})#los pone al revés al c1 y c2 el sis de arriba 
            expr2 = system[1].subs(({C1: const[1].rhs, C2: const[0].rhs}))
            p1 = plot_parametric((expr1, expr2), (t, 0, 10), show = False)
            p.append(p1[0])
    return p.show()
