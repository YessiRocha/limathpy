# my file of the calculus of several variables
import sympy as sp


def partial_derivate(expr, var):
    """Function that returns a sympy expression that represents the
       partial derivative of the given function.
       
    Args:
        expr: Any sympy function.
        var: Variable with respect to whict it is derived.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> z = symbols('z')
        >>> partial_derivate(5*x*y - x*sp.cos(z) + z**8*y, y)
        5*x + z**8"""
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    partial = sp.diff(expr, var)
    return partial


def gradient(expr, var):
    """Function that returns a list representing the function and
       the gradient of the given function.
       
    Args:
        expr: Any sympy function.
        var: List of variables in order.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> z = symbols('z')
        >>> gradient(x*sp.cos(3*z)**2 + sp.sin(y), [x, y, z])
        [x*cos(3*z)**2 + sin(y), cos(3*z)**2, cos(y), -6*x*sin(3*z)*cos(3*z)]
        """
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    gradient = [expr]
    for i in range(len(var)):
        gradient.append(sp.diff(gradient [0], var[i]))
    return gradient


def jacobian(function , var):
    """Function that returns a matrix that represents the Jacobian
       matrix of the given coordinate functions.
       
    Args:
        function: Any list of sympy functions.
        var: List of variables in order.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> z = symbols('z')
        >>> jacobian([x*y*z, x**2 + y**2, sp.sin(x*y*z)], [x, y, z])
        Matrix([
        [           y*z,            x*z,            x*y],
        [           2*x,            2*y,              0],
        [y*z*cos(x*y*z), x*z*cos(x*y*z), x*y*cos(x*y*z)]])
    """
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    m = len(function)
    n = len(var)
    mat = sp.zeros(m, n)
    for i in range(m):
        for j in range(n):
            mat[i, j] = sp.diff(function[i], var[j])
    return mat


def hessian(function, var):
    """Function that returns a matrix that represents the Hessian
       matrix of the given coordinate functions.
       
    Args:
        function: Any list of sympy functions.
        var: List of variables in order.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> z = symbols('z')
        >>> hessian([x*y*z, x*sp.cos(y) -z, x*y*sp.sin(z)], [x, y, z])
        Matrix([
        [0,         0,           0],
        [0, -x*cos(y),           0],
        [0,         0, -x*y*sin(z)]])
    """
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    m = len(function)
    n = len(var)
    mat = sp.zeros(m, n)
    for i in range(m):
        for j in range(n):
            jacb = sp.diff(function[i], var[j])
            mat[i, j] = sp.diff(jacb, var[j])
    return mat


def divergence(expr, var):
    """Function that returns the Divergence of a given vector field.
    
    Args:
        expr: Any list of sympy functions.
        var: List of variables in order.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> divergence([sp.cos(x*2) + y, x + 2*y**4], [x, y])
        8*y**3 - 2*sin(2*x)"""
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    n = len(expr)
    partials = [0]*n
    divergence = 0
    for i in range(n):
        partials[i] = sp.diff(expr[i], var[i])
        divergence = divergence + partials[i]
    return divergence


def laplacian(expr, var):
    """Function that returns the Laplacian of a function of several
       variables.
    
    Args:
        expr: Any sympy function.
        var: List of variables in order.
        
    Returns:
        list: with the reverse of the function.
        
    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> y = symbols('y')
        >>> z = symbols('z')
        >>> laplacian(2*x*y + y**3 - sp.cos(z), [x, y, z])
        6*y + cos(z)"""
    x = sp.symbols('x')
    y = sp.symbols('y')
    z = sp.symbols('z')
    n = len(var)
    partial1 = [0]*n
    partial2 = [0]*n
    laplacian = 0
    for i in range(n):
        partial1[i] = sp.diff(expr, var[i])
        partial2[i] = sp.diff(partial1[i], var[i])
        laplacian = laplacian + partial2[i]
    return laplacian
