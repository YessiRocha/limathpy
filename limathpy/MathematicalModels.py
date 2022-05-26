#File to MathematicalModels
from matplotlib import pyplot as plt
from celluloid import Camera
import numpy as np

#Lotka-Volterra predator-prey model


def rungekutta2_fg(f, g, t0, x0, y0, h, samples):
    """ A function that, returns a system of ordinary differential equations with 2nd Order Runge Kutta.
    Args:
        f: first differential equation
        g: second differential equation
        t0: initial condition of the observation time
        x0: initial condition of the number of dams
        y0: initial condition of the number of predators
        h: algorithm parameter
        samples (int):total number of samples

    Example:
    >>> import numpy as np
    >>> from limathpy import rungekutta2_fg
    >>> print(' [ ti, xi, yi]')
    >>> print(table)"""
    size = samples + 1
    table = np.zeros(shape=(size,3), dtype=float)
    table[0] = [t0, x0, y0]
    ti = t0
    xi = x0
    yi = y0
    for i in range(1, size, 1):
        K1x = h * f(ti, xi, yi)
        K1y = h * g(ti, xi, yi)
        K2x = h * f(ti+h, xi + K1x, yi+K1y)
        K2y = h * g(ti+h, xi + K1x, yi+K1y)
        xi = xi + (1/2)*(K1x+K2x)
        yi = yi + (1/2)*(K1y+K2y)
        ti = ti + h
        table = np.array(table)
        table[i] = [ti, xi, yi]
        #Parameters of the equations
        a = 0.5
        b = 0.7
        c = 0.35
        d = 0.35
        f = lambda t, x, y : a*x - b*x*y
        g = lambda t, x, y : - c*y + d*x*y
        t0 = 0
        x0 = 2
        y0 = 1
        h = 0.5
        samples = 101
        table = rungekutta2_fg(f, g, t0, x0, y0, h, samples)
        ti = table[:,0]
        xi = table[:,1]
        yi = table[:,2]
    np.set_printoptions(precision=6)
    return(table)
    
    
def diagram(par, x0, it):
    """A function that, returns a spiderweb diagram of some function.

    Args:
        par: is a simple structured text parser project
        x0: initial condition
        it: number of steps

    Example:
    >>> from matplotlib import pyplot as plt
    >>> from celluloid import Camera
    >>> from IPython.display import HTML
    >>> import numpy as np
    >>> from limathpy import diagram
    >>> anim = diagram(3.8, 0.1, 200)
    >>> HTML(anim.to_html5_video())"""
    def f(x):
        return par*x*(1-x)
    fig, ax = plt.subplots()
    camera = Camera(fig)
    x = [x0]
    y = [x0]
    s = np.arange(0, 1, 0.01)

    for i in range(it):
        ax.plot(s, f(s), color = 'blue')
        ax.plot(s, s, color = 'black')
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
