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
    >>> a = 0.5
    >>> b = 0.7
    >>> c = 0.35
    >>> d = 0.35
    >>> f = lambda t, x, y: a*x - b*x*y
    >>> g = lambda t, x, y: - c*y + d*x*y
    >>> t0 = 0
    >>> x0 = 2
    >>> y0 = 1
    >>> h = 0.5
    >>> samples = 101
    >>> table = rungekutta2_fg(f, g, t0, x0, y0, h, samples)
[[ 0.        2.        1.      ]
 [ 0.5       1.754875  1.16975 ]
 [ 1.        1.457533  1.302069]
 [ 1.5       1.167405  1.373599]
 [ 2.        0.922773  1.381103]
 [ 2.5       0.734853  1.33689 ]
 [ 3.        0.598406  1.258434]
 [ 3.5       0.502789  1.161433]
 [ 4.        0.43776   1.05747 ]
 [ 4.5       0.39535   0.954156]
 [ 5.        0.36995   0.856056]
 [ 5.5       0.357857  0.765649]
 [ 6.        0.356772  0.684061]
 [ 6.5       0.365425  0.611594]
 [ 7.        0.383295  0.548072]
 [ 7.5       0.41044   0.493068]
 [ 8.        0.447386  0.446051]
 [ 8.5       0.495053  0.406481]
 [ 9.        0.554706  0.373878]
 [ 9.5       0.627897  0.347872]
 [10.        0.716394  0.328248]
 [10.5       0.822045  0.315003]
 [11.        0.946561  0.308414]
 [11.5       1.091121  0.309138]
 [12.        1.255739  0.318359]
 [12.5       1.438245  0.33799 ]
 [13.        1.632757  0.370937]
 [13.5       1.827612  0.421359]
 [14.        2.003072  0.494714]
 [14.5       2.130143  0.59701 ]
 [15.        2.173551  0.73223 ]
 [15.5       2.103102  0.896997]
 [16.        1.913257  1.074276]
 [16.5       1.637137  1.233405]
 [17.        1.334464  1.342914]
 [17.5       1.058838  1.387297]
 [18.        0.836776  1.371067]
 [18.5       0.670949  1.310046]
 [19.        0.552642  1.221667]
 [19.5       0.470823  1.1201  ]
 [20.        0.416031  1.015227]
 [20.5       0.381277  0.913283]
 [21.        0.361785  0.817857]
 [21.5       0.354481  0.73078 ]
 [22.        0.357532  0.652774]
 [22.5       0.370003  0.583897]
 [23.        0.391633  0.523839]
 [23.5       0.422684  0.472107]
 [24.        0.463858  0.42815 ]
 [24.5       0.516226  0.39144 ]
 [25.        0.581185  0.361532]
 [25.5       0.660398  0.338109]
 [26.        0.755699  0.321032]
 [26.5       0.868935  0.310397]
 [27.        1.001686  0.306619]
 [27.5       1.154789  0.310545]
 [28.        1.32756   0.323625]
 [28.5       1.51659   0.348143]
 [29.        1.71397   0.387494]
 [29.5       1.905037  0.446412]
 [30.        2.066237  0.530796]
 [30.5       2.165082  0.646391]
 [31.        2.166022  0.795122]
 [31.5       2.045845  0.96879 ]
 [32.        1.813899  1.144051]
 [32.5       1.517631  1.287402]
 [33.        1.218981  1.371905]
 [33.5       0.962067  1.390516]
 [34.        0.762497  1.353795]
 [34.5       0.616756  1.27932 ]
 [35.        0.514275  1.183658]"""
    size = samples + 1
    table = np.zeros(shape=(size, 3), dtype=float)
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
        table[i] = [ti, xi, yi]
    table = np.array(table)
    ti = table[:,0]
    xi = table[:,1]
    yi = table[:,2]
    np.set_printoptions(precision=6)
    print(' [ ti, xi, yi]')
    return table
        
    
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
