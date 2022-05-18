from sympy import symbols, Function, Eq, Derivative, dsolve, solve, init_printing, exp, sin, cos, tan, plot_parametric
t = symbols('t')
y = Function('y')
x = Function('x')
C1 = symbols('C1')
C2 = symbols('C2')

#Resolución de ecuaciones diferenciales ordinarias.

def first_ode(vector): #vector=(función que acompaña a la y'(t), función que acompaña a la y(t), función que depende de t)
    """Dado un vector, regresa solución general de una ecuación diferencial lineal ordinaria de primer orden."""
    eq = Eq(vector[0]*Derivative(y(t), t) + vector[1]*y(t), vector[2])
    sol = dsolve(eq, y(t))
    return sol


def solve_first_ode(vector, cond_inic=[1,0]): #cond_inic es y(t=1)=0 ; vector=(función que acompaña a la y'(t), función que acompaña a la y(t), función que depende de t)
    equation = first_ode(vector).rhs
    evalu = equation.subs({t: cond_inic[0]})
    lin = Eq(evalu, 0)
    solu = solve(lin, C1)
    final = equation.subs({C1: solu[0]})
    return Eq(y(t), final)

#Para resolver EDO de segundo orden de la forma y'' + p(t)y' + q(t)y + g(t)=0.

def second_ode_const(vector): #vector=(función que acompaña a la y''(t),función que acompaña a la y'(t), función que acompaña a la y(t), función que depende de t)
    """Dado un vector, resuelve una ecuación diferencial ordinaria de segundo orden con constantes C1 y C2."""
    eq = Eq(vector[0]*Derivative(y(t), t, 2) + vector[1]*Derivative(y(t), t) + vector[2]*y(t), vector[3])
    sol = dsolve(eq, y(t))
    return sol


def solve_2nd_edo(vector, cond_inic=[[1, 0], [0, 1]]): #y(1)=0, y(0)=1 vector=(y'', y', y, g(t)
    """Dado un vector y condiciones iniciales resuelve una ecuación diferencial ordinaria de segundo orden con constantes."""
    equation = second_ode_const(vector).rhs
    evalu1 = equation.subs({t: cond_inic[0][0]})
    evalu2 = equation.subs({t: cond_inic[1][0]})
    eq1 = Eq(evalu1, cond_inic[0][1])
    eq2 = Eq(evalu2, cond_inic[1][1])
    const = solve((eq1, eq2))
    final = equation.subs(const)
    return Eq(y(t), final)

#Para resolver sistemas de ecuaciones diferenciales ordinarias.

def sistema(matriz):
    """Dada un lista de listas, regresa un sistema de ecuaciones diferenciales."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    eq1 = Eq(Derivative(x(t), t), matriz[0][0]*x(t) + matriz[0][1]*y(t))
    eq2 = Eq(Derivative(y(t), t), matriz[1][0]*x(t) + matriz[1][1]*y(t))
    sols = dsolve((eq1, eq2))
    return sols[0].rhs, sols[1].rhs


def sistema_lineal(matriz, cond_inic=[[1, 2], [3, 4]]): #x(1)=3, y(2)=4 vector=(y'', y', y, g(t)
    """Dado un sistema de ecuaciones diferenciales, regresa el sistema lineal en t=0."""
    t = symbols('t')
    x, y = symbols('x y', cls=Function)
    sols = sistema(matriz)
    lineal1 = Eq(sols[0].subs({t:cond_inic[0][0]}), cond_inic[1][0])
    lineal2 = Eq(sols[1].subs({t:cond_inic[0][1]}), cond_inic[1][1])
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

#Retrato de fase.

def phase_portrait(matriz, lim_initialconditions=2):
    system = sistema(matriz)
    p = plot_parametric((0,0), (t, 0, 0), show=False, title = 'Phase portrait')
    for i in range(0, lim_initialconditions):
        for j in range(0, lim_initialconditions):
            const = sistema_lineal(matriz, [i, j])
            expr1 = system[0].subs({C1: const[1].rhs, C2: const[0].rhs})#los pone al revés al c1 y c2 el sis de arriba 
            expr2 = system[1].subs(({C1: const[1].rhs, C2: const[0].rhs}))
            p1 = plot_parametric((expr1, expr2), (t, 0, 10), show = False)
            p.append(p1[0])
    return p.show()
