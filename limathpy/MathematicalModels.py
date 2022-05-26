#File to MathematicalModels
from sympy import plot, symbols, Function, Eq, Derivative, dsolve, solve
from matplotlib import pyplot as plt
from celluloid import Camera
from IPython.display import HTML
import numpy as np
import scipy.integrate as spi
import pylab as pl

def sistem(matriz):
    """Given a list of lists, a system of differential equations returns.
    Args:
        matriz: is a simple structured text parser project

    Example:
    >>> from matplotlib import pyplot as plt
    >>> import numpy as np
    >>> from limathpy import diagram
    >>> diagrama(f, 0.1, 100)""""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    eq1 = Eq(Derivative(x(t), t), matriz[0][0]*x(t) + matriz[0][1]*y(t))
    eq2 = Eq(Derivative(y(t), t), matriz[1][0]*x(t) + matriz[1][1]*y(t))
    sols = dsolve((eq1, eq2))
    return sols[0].rhs, sols[1].rhs
    

def linear_system(matriz, cond_inic):
    """Given a system of differential equations, the linear system returns to :math:`t=0"""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    sols = sistem(matriz)

    lineal1 = Eq(sols[0].subs({t:0}), cond_inic[0])
    lineal2 = Eq(sols[1].subs({t:0}), cond_inic[1])
    return lineal1, lineal2

def sistema_ed(matriz, cond_inic):
    """Given a matrix and initial conditions, the solution of the differential equation system returns."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    C1, C2 = symbols('C1 C2')
    sis_ed = sistem(matriz)
    sis_lin = linear_system(matriz, cond_inic)
    dict_sols = solve(sis_lin)
    expr1 = sis_ed[0].subs(dict_sols)
    expr2 = sis_ed[1].subs(dict_sols)
    return expr1, expr2
    
    
def diagram(par, x0, it):
    """ A function that, returns a spiderweb diagram of some function.
    Args:
        par: is a simple structured text parser project
        x0: initial condition
        it: number of steps

    Example:
    >>> from matplotlib import pyplot as plt
    >>> import numpy as np
    >>> from limathpy import diagram
    >>> diagrama(f, 0.1, 100)"""
    def f(x):
        return par*x*(1-x)
    fig, ax = plt.subplots()
    camera = Camera(fig)
    x = [x0]
    y = [x0]
    s = np.arange(0, 1, 0.01)

    for i in range(it):
        ax.plot(s, f(s), color='blue')
        ax.plot(s, s, color='black')
        x.append(x[2*i])
        x.append(f(x[2*i]))
        y.append(f(y[2*i]))
        y.append(f(y[2*i]))
        ax.plot(x, y, color='red')
        camera.snap()
    return camera.animate()    


def fibonacci(n):
    """A function that, returns the n-th Fibonacci number
    Args:
        n (int): the integer number 
    Example:
        >>> from limathpy import fibonacci
        >>> [fibonacci(n) for n in range(1, 20)]
        [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584]"""
    if n == 1:
        return 0
    elif n == 2:
        return 1
    else:
        return fibonacci(n-1) + fibonacci(n-2)



# population size
#N=1
#beta=1.4247
#gamma=0.14286
# time step
#TS=1.0 
#ND=70.0
#S0=1-1e-6
#I0=1e-6
#INPUT = (S0, I0, 0.0)


def diff_eqs(INP,t):  
    Y=np.zeros((3))
    V = INP
    """Diferential equations"""
    Y[0] = - beta * V[0] * V[1]
    Y[1] = beta * V[0] * V[1] - gamma * V[1]
    Y[2] = gamma * V[1]
    return Y   # For odeint

#t_start = 0.0; t_end = ND; t_inc = TS
#t_range = np.arange(t_start, t_end+t_inc, t_inc)
#RES = spi.odeint(diff_eqs,INPUT,t_range)

#Grafic
#pl.plot(RES[:,0]*N, '-g', label='Susceptible')
#pl.plot(RES[:,2]*N, '-k', label='Rerecovered')
#pl.plot(RES[:,1]*N, '-r', label='Infected')
#pl.legend(loc=0)
#pl.title('Basic Model SIR')
#pl.xlabel('Time')
#pl.savefig('sirpy')
