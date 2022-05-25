"""This module contains functions to solve some calclus problems in one variable. """

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


def n_derivatives(expr, n=1):
    """Function that returns a list with the n derivates of an expression,
    with n given.

    Args: 
        expr: Any sympy function
        n (int, optional): The number of derivatives required.
        Defaults to one.

    Returns:
        list: the funtion and the indicated derivatives.

    Example:
        >>> from sympy import symbols
        >>> x = symbols('x')
        >>> n_derivatives(x**4, 4)
        [x**4, 4*x**3, 12*x**2, 24*x, 24]"""
    derivatives = [expr]
    x = sp.symbols('x')
    for i in range(n):
        derivatives.append(sp.Derivative(derivatives[-1], x).doit())
    return derivatives


def graph_fyd(expression):
    """Function that graphs an expression given as a string and its derivative
    on the same plane.

    Args: 
        expression (str): Expression of a function.

    Example:
        >>> import matplotlib.pyplot as plt
        >>> graph_fyd("x**2")

    .. image:: graph_fyd.png
      :align: center"""
    x = sp.symbols('x')
    expr = sp.sympify(expression)
    deriv = sp.diff(expr, x)
    f, f_prime = sp.lambdify(x, expr, 'numpy'), sp.lambdify(x, deriv, 'numpy')
    domain = np.linspace(-10, 10)
    f_eval = f(domain)
    f_prime_eval = f_prime(domain)
    if type(f_prime_eval) == float or type(f_prime_eval) == int:
        for j in range(len(domain) - 1):
            f_prime_eval = np.append(f_prime_eval, [f_prime(domain)])

    fig, ax = plt.subplots()
    ax.set_title("Function and derivative")
    ax.plot(domain, f_eval, label=expression)
    ax.plot(domain, f_prime_eval, label='Derivative')
    ax.set_xlabel("$x$")
    ax.legend(loc='center',
              bbox_to_anchor=(0.78, -0.13),
              shadow=True,
              ncol=2)
    plt.draw_if_interactive()


def tangent_line(expression, x_0):
    """Function that gives the equation of a tangent line to a function about
    a given point.

    Args:
        expression (str): Expression of a function.
        x_0: Value of the x coordinate for the point of tangency of the line.

    Returns:
        A sympy equation of the tangent line to the function through
        the given point.

    Example:
        >>> tangent_line("x**2", 1)
        Eq(y, 2*x - 1)"""
    x, y = sp.symbols('x y')
    expr = sp.sympify(expression)
    y_0 = expr.subs({x: x_0})
    deriv = sp.diff(expr, x)
    slope = deriv.subs({x: x_0})
    line = sp.Eq(y, slope * (x - x_0) + y_0)
    return line


def root_f(expression, number=0):
    """Function that returns if an expression evaluated to a given number n,
    is zero.

    Args:
        expression (str): Expression of a function.
        number (optional): The real number at which the function is evaluated.
        Defaults to zero.

    Returns:
        bool: True if the given number is the root of the function,
        False otherwise.

    Example:
        >>> root_f("x**2", 0)
        True"""
    x = sp.symbols('x')
    expr = sp.sympify(expression)
    return expr.subs({x: number}) == 0


def revolution_area(expression, lower_bound, upper_bound):
    """Function that calculates the area, over an interval, of a surface of
    revolution whose axis of rotation is the x or y-axis.

    Args:
        expression: A sympy function that generate the surface of revolution.
        lower_bound: The lower bound of the interval over which the surface
        of revolution is defined.
        upper_bound: The upper bound of the interval over which the surface
        of revolution is defined.

    Returns:
        A numerical sympy expression of the area of the surface of revolution
        on the given interval.

    Example:
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> revolution_area(x**2, 0, 2)
    pi*(-asinh(4) + 132*sqrt(17))/32"""
    x = sp.symbols('x')
    expr = expression * sp.sqrt(1 + (expression.diff(x)) ** 2)
    surface_area = 2*sp.pi*sp.integrate(expr, (x, lower_bound, upper_bound))
    return sp.simplify(surface_area)


def reverse_func(expression):
    """Function that returns the inverse of a given expression.

    Args:
        expression: A sympy function.

    Returns:
        list: with the reverse of the function.
        it can have length greater than 1 due to the domain in which the
        function is injective.

    Example:
    >>> from sympy import symbols
    >>> x = symbols('x')
    >>> reverse_func(x**2)
    [-sqrt(y), sqrt(y)]"""
    x, y = sp.symbols('x y')
    rev = sp.solve(sp.Eq(y, expression), x)
    return rev


class TestLimitDiverges(Exception):
    def __init__(self, message="The quotient limit for this test does "
                               "not converge."):
        super().__init__(message)


def seq_converg(expression):
    """Function that determines whether a sequence converges to zero or diverges
    using the quotient limit test.

    The limit of :math:`a_{n+1}/a_n` is considered, for a sequence
    :math:`(a_n)_{n \in \mathbb{N}}` of positive real numbers.

    Args:
        expression: A sympy function in terms of n.

    Returns:
        Message: Indicates if the sequence converges to zero, diverges
        or nothing can be said about it.

    Example:
    >>> from sympy import symbols
    >>> n = symbols('n')
    >>> seq_converg(1/2**n)
    'The sequence 2**(-n) converges to zero.'"""
    n = sp.symbols('n')
    expr_2 = expression.subs({n: n + 1})
    r = sp.limit_seq(expr_2 / expression, n)
    if 0 < r < 1:
        return f"The sequence {expression} converges to zero."
    elif r > 1:
        return f"The sequence {expression} diverges."
    elif r == 1:
        return f"Nothing can be said about the sequence {expression}, " \
               f"try another method."
    else:
        raise TestLimitDiverges


def seri_converg(expression):
    """Function that determines whether a serie from 1 to infinity converges
    or diverges using the quotient limit test.

    The limit of :math:`a_{n+1}/a_n` is considered, for a sequence
    :math:`(a_n)_{n \in \mathbb{N}}`.

    Args:
        expression: A sympy function in terms of n.

    Returns:
        Message: Indicates if the serie from 1 to infinity converges,
        diverges or nothing can be said about it.

    Example:
    >>> from sympy import symbols
    >>> n = symbols('n')
    >>> seri_converg(1/2**n)
    'The series from 1 to infinity of the sequence 2**(-n) converges.'"""
    n = sp.symbols('n')
    expr_2 = expression.subs({n: n + 1})
    r = sp.limit_seq(expr_2 / expression, n)
    if 0 < r < 1:
        return f"The infinite series of the sequence {expression} " \
               f"converges."
    elif r > 1:
        return f"The infinite series of the sequence {expression} " \
               f"diverges."
    elif r == 1:
        return f"Nothing can be said about the infinite series " \
               f"{expression}, try another method."
    else:
        raise TestLimitDiverges
