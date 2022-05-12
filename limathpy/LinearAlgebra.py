import sympy as sp


def matrix_n(matriz, n):
    """Powers of a matrix.

    Args:
        parameter1 = matrix of sympy
        parameter2 = n (the power)
        
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
        param1 = list
    
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
        param1 = integer, the element (1,2) of the matrix
    
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
    """Integer eigenvalues
    
    Args:
        param1 = square sympy Matrix
        
    Returns:
        string:
        'The matrix has all its eigenvalues positive' or 
        'The matrix has not all its eigenvalues positive.'
        
    Examples:
        >>> import sympy as sp
        >>> int_eigenvalues(sp.Matrix([[1,2],[1,3]]))
        The Matrix([[1, 2], [1, 3]]) has not all its eigenvalues positive
        
    """
    eigenvalues=list(matriz.eigenvals())
    if integers_list(eigenvalues) == True:
        print(f"The {matriz} has all its eigenvalues positive")
    else:
        print(f"The {matriz} has not all its eigenvalues positive")

        
def inner_product(vector1,vector2):
    """Documentation"""
    return list(vector1.T*vector2)[0]

def change_basis(base1,base2):
    eb2=base2.inv()*base1
    return eb2
