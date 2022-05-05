import sympy as sp

def potencias_matriz(matriz, n):
    """Powers of a matrix 
       Args:
            parameter1 = matrix of sympy
            parameter2 = n 
       Returns:
            The n powers of a given matrix in a list."""
    pot_matriz = [matriz**i for i in range(1,n+1)]
    return pot_matriz


def val_ent_list(lista):
    """Función que regresa True si la lista tiene todos sus valores enteros y False si no"""
    flotantes = [float(num) for num in lista]
    lista1 = [flotante.is_integer() for flotante in flotantes]
    return all(lista1)


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
