import sympy as sp

def matrix_n(matriz, n):
    """Powers of a matrix.
    
    Args:
        parameter1 = matrix of sympy
        parameter2 = n (the power)
        
    Returns:
        List of n powers of a given matrix.
        
    Examples:
        >>> import sympy as sp
        >>> matrix_n(sp.Matrix([[1,2],[3,4]]),2)
        [Matrix([
        [1, 2],
        [3, 4]]),
        Matrix([
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
        'True' if all values in the given list are integers or 
        'False' in other case.
    
    Examples:
        >>> integers_list([1,3,0.5])
        False
    """
    floats = [float(num) for num in lista]
    list1 = [flotante.is_integer() for flotante in floats]
    return all(list1)


def matrices_eigenvalores_ent(n):
    """Función que regresa los valores n para los cuales la matriz [[1,n],[1,1]] tiene eigenvalores enteros"""
    enes=[]
    for i in range(0,n+1):
        matriz = Matrix([[1,i],[1,1]])
        eigenvalores = list(matriz.eigenvals())
        if val_ent_list(eigenvalores) == True:
            enes.append(i)      
    return enes


def mat_eigen_ent(matriz):
    """Integer eigenvalues
    Args:
        param1=square sympy Matrix
    Returns:
        'Matrix tiene todos sus valores propios enteros or 
        Matrix no tiene todos sus valores propios enteros' """
    eigenvalues=list(matriz.eigenvals())
    if val_ent_list(eigenvalues) == True:
        print(f"{matriz} tiene todos sus valores propios enteros")
    else:
        print(f"{matriz} no tiene todos sus valores propios enteros")

        
def camb_base(matriz):
    pass


def inner_product(vector1,vector2):
    """Documentation"""
    return list(vector1.T*vector2)[0]
