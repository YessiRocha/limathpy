# my file of calculus in one variable
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
        >>> from sympy import diff, symbols, Eq, sympify
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
    >>> import sympy as sp
    >>> x = sp.symbols('x')
    >>> revolution_area(x**2, 0, 2)
    """
    x = sp.symbols('x')
    expr = expression * sp.sqrt(1 + (expression.diff(x)) ** 2)
    surface_area = 2*sp.pi*sp.integrate(expr, (x, lower_bound, upper_bound))
    return sp.simplify(surface_area)
