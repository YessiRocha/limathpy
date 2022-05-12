from sympy import symbols, Function, Eq, Derivative, dsolve, solve, plot_parametric

#Resolución de ecuaciones diferenciales ordinarias




#Resolución de sistemasde ecuaciones diferenciales ordinarias

def sistema(matriz):
    """Given a 2x2 Matrix it returns an ordinary differential equations system.

    Examples
    --------
    >>> from limathpy import sistema
    >>> sistema([[1, -1], [1, 0]])
    ((C1/2 - sqrt(3)*C2/2)*exp(t/2)*cos(sqrt(3)*t/2) - (sqrt(3)*C1/2 + C2/2)*exp(t/2)*sin(sqrt(3)*t/2), C1*exp(t/2)*cos(sqrt(3)*t/2) - C2*exp(t/2)*sin(sqrt(3)*t/2))

    """
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    eq1 = Eq(Derivative(x(t), t), matriz[0][0]*x(t) + matriz[0][1]*y(t))
    eq2 = Eq(Derivative(y(t), t), matriz[1][0]*x(t) + matriz[1][1]*y(t))
    sols = dsolve((eq1, eq2))
    return sols[0].rhs, sols[1].rhs


def sistema_lineal(matriz, cond_inic):
    """Given an ordinary differential equations system, it returns the linear system when :math:`t=0`.
    
    Examples
    --------
    >>> from limathpy import sistema_lineal
    >>> sistema_lineal([[1,4], [1, 0]], [0,1])
    (Eq(C1*(1 - sqrt(17))/2 + C2*(1 + sqrt(17))/2, 0), Eq(C1 + C2, 1))
    
    """
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    sols = sistema(matriz)

    lineal1 = Eq(sols[0].subs({t:0}), cond_inic[0])
    lineal2 = Eq(sols[1].subs({t:0}), cond_inic[1])
    return lineal1, lineal2


def sistema_ed(matriz, cond_inic):
    """Given a 2x2 Matrix and some initial conditions, it returns the solution of the associated system of differential equations.
    
    Examples
    --------
    >>> from limathpy import sistema_ed
    >>> sistema_ed([[1,-2], [1, 0]], [1,1])
    (-3*sqrt(7)*exp(t/2)*sin(sqrt(7)*t/2)/7 + exp(t/2)*cos(sqrt(7)*t/2), sqrt(7)*exp(t/2)*sin(sqrt(7)*t/2)/7 + exp(t/2)*cos(sqrt(7)*t/2))
    
    """
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    C1, C2 = symbols('C1 C2')
    sis_ed = sistema(matriz)
    sis_lin = sistema_lineal(matriz, cond_inic)
    dict_sols = solve(sis_lin)
    expr1 = sis_ed[0].subs(dict_sols)
    expr2 = sis_ed[1].subs(dict_sols)
    return expr1, expr2

def phase_portrait(matriz, ci = [2, 1]):
    """Given a 2x2 Matrix and some initial conditions, it returns the phase portrait of the associated system of differential equations.
    
    Examples
    --------
    >>> from limathpy import phase_portrait
    >>> phase_portrait([[1.5, 5], [-1, 0]], ci = [2, 1])
    aquí va una imagen del retrato de fase 
    
    """
    t, a = symbols('t a')
    sol = sistema_ed(matriz, ci) 
    lista_de_soluciones0 = [sol[0].rhs.subs({a:i}) for i in range(-10, 16)]
    lista_de_soluciones1 = [sol[1].rhs.subs({a:i}) for i in range(-10, 16)]
    p1 = plot_parametric(sol[0].rhs.subs({a:1}), sol[1].rhs.subs({a:1}), (t, 0, 10), show = False)
    lista_de_gráficas = [plot_parametric(lista_de_soluciones0[k], lista_de_soluciones1[k], (t, 0, 10), show = False) 
                     for k in range(26)]
    for gráfica in lista_de_gráficas:
        p1.append(gráfica[0])

    return p1.show()
