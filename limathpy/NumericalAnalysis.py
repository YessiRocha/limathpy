#File to Numerical analysis 
from sympy import plot, symbols, Function, Eq, Derivative, dsolve, solve

def sistema(matriz):
    """Dada un lista de listas, regresa un sistema de ecuaciones diferenciales."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    eq1 = Eq(Derivative(x(t), t), matriz[0][0]*x(t) + matriz[0][1]*y(t))
    eq2 = Eq(Derivative(y(t), t), matriz[1][0]*x(t) + matriz[1][1]*y(t))
    sols = dsolve((eq1, eq2))
    return sols[0].rhs, sols[1].rhs

def sistema_lineal(matriz, cond_inic):
    """Dado un sistema de ecuaciones diferenciales, regresa el sistema lineal en t=0."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    sols = sistema(matriz)

    lineal1 = Eq(sols[0].subs({t:0}), cond_inic[0])
    lineal2 = Eq(sols[1].subs({t:0}), cond_inic[1])
    return lineal1, lineal2

def sistema_ed(matriz, cond_inic):
    """Dada una matriz y condiciones iniciales, regresa la solución del sistema de ed."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    C1, C2 = symbols('C1 C2')
    sis_ed = sistema(matriz)
    sis_lin = sistema_lineal(matriz, cond_inic)
    dict_sols = solve(sis_lin)
    expr1 = sis_ed[0].subs(dict_sols)
    expr2 = sis_ed[1].subs(dict_sols)
    return expr1, expr2
    
sols = sistema_ed(([[0, 1], [-1, 0]]), [1, -2])
sols


x = symbols('x')
p1 = plot(x*x)
