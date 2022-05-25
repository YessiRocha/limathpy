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
