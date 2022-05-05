#my file of calculus in one variable
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt 


x = sp.symbols('x')
def n_derivates(expr, n=1):
    """Function that returns a list with the n derivates of an expression, with n given
    Args: 
        param1 = a sympy expression of a function
        param2 = the number of derivatives required (this parameter is optional, by default returns one derivative)
    Returns: 
        A list with the funtion and the indicated derivatives"""
    derivates = [expr]
    for i in range(n):
        derivates.append(sp.Derivative(derivates[-1], x).doit())
    return derivates
  
  
  def graph_fyd(expresion):
    """Function that graphs an expression given as a string and its derivative in the same plane
    Args:
        param1 = expression of a function as a string
    Returns:
        The graph of the function and its derivative in the same plane """
    x = sp.symbols('x')
    expr = sp.sympify(expresion)
    deriv = sp.diff(expr, x)
    f, f_prime = sp.lambdify(x, expr, 'numpy'), sp.lambdify(x, deriv, 'numpy')
    domain = np.linspace(-10, 10)
    f_eval = f(domain)
    f_prime_eval = f_prime(domain)
    if type(f_prime_eval) == float or type(f_prime_eval) == int:
        for j in range(len(domain)-1):
            f_prime_eval=np.append(f_prime_eval, [f_prime(domain)])
            
    
    fig, ax = plt.subplots()
    ax.set_title("Función y derivada")
    ax.plot(domain, f_eval)
    ax.plot(domain, f_prime_eval)
    ax.set_xlabel("$x$")
    plt.show()
