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
