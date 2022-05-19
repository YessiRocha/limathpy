import sympy as sp
import numpy as np
import matplotlib.pyplot as plt 
from scipy import linalg


def matrix_n(matriz, n):
    """Powers of a matrix.

    Args:
        matriz: matrix of sympy
        n (int): the power.
        
    Returns:
        list: List of n powers of a given matrix.
        
    Examples:
        >>> import sympy as sp
        >>> matrix_n(sp.Matrix([[1,2],[3,4]]),2)
        [Matrix([
        [1, 2],
        [3, 4]]), Matrix([
        [ 7, 10],
        [15, 22]])]
    """
    pot_matriz = [matriz**i for i in range(1,n+1)]
    return pot_matriz


def integers_list(lista):
    """Integers in a list.
    
    Args:
        lista (list): any list
    
    Returns:
        string: 'True' if all values in the given list are integers or 
        'False' in other case.
    
    Examples:
        >>> integers_list([1,3,0.5])
        False
    """
    floats = [float(num) for num in lista]
    list1 = [flotante.is_integer() for flotante in floats]
    return all(list1)


def int_eigvals_n(n):
    """Values of n for which [[1,n],[1,1]] has integer eigenvalues.
    
    Args:
        n (int): the element (1,2) of the matrix.
    
    Returns:
        list: A list with the values between 0 and n for which [[1,n],[1,1]] has 
        positive eigenvalues.
    
    Examples: 
        >>> import sympy as sp
        >>> int_eigvals_n(100)
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
    """
    enes=[]
    for i in range(0,n+1):
        matriz = sp.Matrix([[1,i],[1,1]])
        eigenvalores = list(matriz.eigenvals())
        if integers_list(eigenvalores) == True:
            enes.append(i)           
    return enes


def int_eigenvalues(matriz):
    """Examine if the eigenvalues are integers
    
    Args:
        matriz (matrix) = square sympy Matrix.
        
    Returns:
        string:
        'The matrix has all its integer eigenvalues' or 
        'The matrix has not all its integer eigenvalues.'
        
    Examples:
        >>> import sympy as sp
        >>> int_eigenvalues(sp.Matrix([[1,2],[1,3]]))
        The Matrix([[1, 2], [1, 3]]) has not all its integer eigenvalues
        
    """
    eigenvalues=list(matriz.eigenvals())
    if integers_list(eigenvalues) == True:
        print(f"The {matriz} has all its integer eigenvalues")
    else:
        print(f"The {matriz} has not all its integer eigenvalues")


def inner_product(vector1, vector2):
    """Inner product. 

    Args:
        vector1 (matrix): Sympy matrix 2x1
        vector1 (matrix): Sympy matrix 2x1.

    Returns:
        (float): usual inner product in :math:`\mathbb{R}^n`
    
    Example:
        >>> import sympy as sp
        >>> v1=sp.Matrix([0,1])
        >>> v2=sp.Matrix([1,2])
        >>> inner_product(v1,v2)
        2.0
    """
    return float(list(vector1.T*vector2)[0])


def change_basis(base1, base2):
    """Change of basis matrix, a matrix that translates vector
    representations from one basis, such as the standard coordinate 
    system, to another basis.
    
    Args:
        base1 (matriz): Sympy matrix as a representation of a base
        base2 (matriz): Sympy matrix as a representation of a base.

    Returns:
        matrix: Change of basis matrix from base1 to base2
    
    Example:
        >>> # If we want the base change matrix from B1={(3,1),(2,-1)} to B2={(2,4),(-5,3)}
        >>> # we use the matrix representations of each base as in this example 
        >>> import sympy as sp
        >>> B1=sp.Matrix([[3,2],[1,-1]])
        >>> B2=sp.Matrix([[2,-5],[4,3]])
        >>> B1, B2, change_basis(B1,B2)
        (Matrix([
        [3,  2],
        [1, -1]]),
        Matrix([
        [2, -5],
        [4,  3]]),
        Matrix([
        [ 7/13,  1/26],
        [-5/13, -5/13]])) 
    """
    eb2=base2.inv()*base1
    return eb2


def graph_solution(expr1, expr2):
    """ Graphical Solution to a 2 x 2 System of Equations.
    
    Args: 
        expr1 (expression): equation expressed in terms of :math:`x`
        expr2 (expression): equation expressed in terms of :math:`x`.
    
    Example:
        >>> import sympy as sp
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt 
        >>> expr1 = x - 2
        >>> expr2 = 2*x  +3
        >>> graph_solution(expr1,expr2)
        .. image:: graph_solution.png
          :align: center

    """
    x, y = sp.symbols('x y')
    a=sp.solve([y-expr1,y-expr2], [x, y])
    rect1 = sp.sympify(expr1)
    rect2 = sp.sympify(expr2, x)
    f1, f2 = sp.lambdify(x, rect1, 'numpy'), sp.lambdify(x, rect2, 'numpy')
    if a == [0]:
        domain = np.linspace(-5, 5)
        image = np.linspace(-5, 5)
    else:
        domain = np.linspace(-5+float(a[x]), 5+float(a[x]))
        image = np.linspace(-5+float(a[y]), 5+float(a[y]))
    inter = [a[x],a[y]]
    F_eval = f1(domain)
    G_eval = f2(domain)
    if type(G_eval) == float or type(G_eval) == int:
        for j in range(len(domain) - 1):
            G_eval = np.append(G_eval, [f2(domain)])
    fig, ax = plt.subplots()
    ax.set_title("Graphic solution")
    ax.plot(domain, image)
    ax.plot(domain, F_eval, color='blue', label=expr1)
    ax.plot(domain, G_eval, color='green', label=expr2)
    ax.plot(a[x], a[y], '.', color='black', linewidth = 0.25)
    plt.axhline(a[y], color='gray', linewidth = 0.5, linestyle='dashed')
    plt.axvline(a[x], color='gray', linewidth = 0.5, linestyle='dashed')
    plt.annotate((a[x],a[y]), (a[x],a[y]), fontsize=12, color='black')
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    plt.draw_if_interactive()
