#File to Numerical analysis 
def suma(x, y, *extras):
  if extras:
    return f"{extras[0]} {x+y}"
  else:
    return x+y

suma(4, 5, "La suma es"), suma(4, 5)
