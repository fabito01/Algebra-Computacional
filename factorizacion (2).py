# -*- coding: utf-8 -*-
"""
@author: ALBERTO PEÑA PEÑALVER Y FABIO TORRES MARTINEZ
"""


# para esta actividad necesitamos la anterior (cuerpos_finitos.py)
# completa, pero quitando las funciones de factorizaciÃ³n en anillo_fp_x
# y en anillo_fq_x, ya que las implementaremos aquÃ­ (no se olviden de
# quitarlas)
import cuerpos_finitos as cf
import random

'''
FACTORIZACION Fpx
'''

# square-free factorization
# input: fpx --> anillo_fp_x
# input: f --> polinomio fabricado por fpx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mÃ³nicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mÃ³nicos
# distintos entre si, ei >= 1, c en fp, entonces g = f1 * f2 * ... * fr

def sqfree_fact_fpx(fpx, f): #FASE 1
    f = fpx.reducir(f)
    f = fpx.monico(f)
    l = []
    sol = fpx.uno()
    p = fpx.fp.p
    while not fpx.es_uno(f):  
        g = fpx.div(f, fpx.gcd(f, derivar_fp(fpx, f))) #g = f/ gcd(f,f')
        while not fpx.es_uno(g): #siguiendo el pseudocódigo visto en clase
            f = fpx.div(f, g) 
            h = fpx.gcd(f, g)
            m = fpx.div(g, h)
            if not fpx.es_uno(m):
                l.append(m)
            g = h
        if not fpx.es_uno(f):
            f = entre_p(fpx,p,f)
    for x in l: #multiplicamos todos los factores con multiplicidad 1 obtenidos
        sol = fpx.mult(sol, x)
    return sol


def derivar_fp(fpx, f): #hacemos la derivada
    l = list(f)
    l.pop(0)
    for i in range(len(l)):
        l[i] = (l[i] * (i+1)) % fpx.fp.p
    return fpx.reducir(tuple(l))

def entre_p(fpx,p, f): #para cuando se anule la derivada, dividimos todos los exponentes entre p
    a = []
    for i in range(((len(f)-1)//p)+1):
        a.append(f[p*i]) #solo tendrá posiciones potencia de p (pues se anuló la derivada)
    return fpx.reducir(fpx.elem_de_tuple(tuple(a)))

# distinct-degree factorization
# input: fpx --> anillo_fp_x
# input: g --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mÃ³nicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mÃ³nicos de h de grado = i, el Ãºltimo hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr

def didegr_fact_fpx(fpx, g): #FASE 2
    l = []
    g = fpx.monico(g)
    h = fpx.mod((0, 1), g) # x
    while not fpx.es_uno(g): #en la iteracion i, obtenemos el producto de factores irreducibles de grado i  
        h = fpx.pot_mod(h, fpx.fp.p, g) 
        f = fpx.gcd(fpx.resta(h, (0, 1)), g) #x^p -x (mod g)
        if not fpx.es_uno(f):
            l.append(f)
            g = fpx.div(g, f) #le quitamos los factores de grado i de nuestro polinomio
            h = fpx.mod(h, g) #redefinimos h
        else:
            l.append(fpx.uno()) #si no tiene, añadimos el 1 para indicarlo y que trabaje eficientemente en la fase 3 
    return l
didegr_fact_fpx

# equal-degree factorization
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: r --> int
# input: h --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mÃ³nicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mÃ³nicos de grado = r


def eqdegr_fact_fpx(fpx,r,h): #FASE 3
    return fase3(fpx, r, h, []) #hacemos una llamada a una funcion auxilar que trabaje con una lista acumuladora 

def fase3(fpx, d, gd,lista):
    if fpx.grado(gd) < d: #no debemos continuar
        return None
    if fpx.grado(gd) == d: #hemos acabado, añadimos nuestro factor a la lista de factores
        lista.append(gd)
        return lista
    #construimos nuestro polinomio aleatorio, de grado menor que nuestro polinomio de entrada
    grado_a = random.randint(1, fpx.grado(gd)-1)
    i = 0
    a = []
    while i < grado_a:
        a.append(fpx.fp.aleatorio())
        i += 1
    a.append(1)
    a = fpx.reducir(tuple(a))
    aux = fpx.pot_mod(a, (fpx.fp.p**d - 1)//2, gd) #construimos los polinomios A(x), B(x), C(x) tal y como hemos visto en clase
    A = fpx.gcd(a, gd)
    B = fpx.gcd(fpx.resta(aux, fpx.uno()), gd)
    C = fpx.gcd(fpx.suma(aux, fpx.uno()), gd)
    fase3(fpx,d, A, lista), fase3(fpx,d, B, lista), fase3(fpx, d, C, lista) #llamamos recursivamente a fase3
    return lista

# multiplicidad de factor irreducible mÃ³nico
# input: fpx --> anillo_fp_x
# input: f --> polinomio de fpx (objeto opaco) no nulo
# input: u --> polinomio irreducible mÃ³nico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fpx(fpx, f, u): #dividimos f entre cada factor. Lo hará tantas veces como sea su multiplicidad.
    i = 0
    while fpx.es_igual(fpx.mod(f, u), fpx.cero()):
        f = fpx.div(f, u)
        i += 1
    return i

# factorizaciÃ³n de Cantor-Zassenhaus
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: f --> polinomio de fpx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorizaciÃ³n completa de f en irreducibles mÃ³nicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fpx(fpx, f):                     # mantener esta implementacion
    g = sqfree_fact_fpx(fpx, f)
    h = didegr_fact_fpx(fpx, g)
    irreducibles = []
    for r in range(len(h)):
        if fpx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fpx(fpx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fpx(fpx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para aÃ±adir la funciÃ³n de factorizaciÃ³n de Cantor-Zassenhaus
# como un mÃ©todo de la clase anillo_fp_x
cf.anillo_fp_x.factorizar = fact_fpx

'''
FACTORIZACION Fqx
'''
# square-free factorization
# input: fqx --> anillo_fq_x
# input: f --> polinomio fabricado por fqx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mÃ³nicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mÃ³nicos
# distintos entre si, ei >= 1, c en fq, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fqx(fqx, f): #FASE 1
    f = fqx.reducir(f)
    f = fqx.monico(f)
    l = []
    sol = fqx.uno()
    p = fqx.fq.fpx.fp.p
    while not fqx.es_uno(f):
        g = fqx.div(f, fqx.gcd(f, derivar_fq(fqx, f))) #g = f/ gcd(f,f')
        while not fqx.es_uno(g): #siguiendo el pseudocódigo visto en clase
            f = fqx.div(f, g)
            h = fqx.gcd(f, g)
            m = fqx.div(g, h)
            if not fqx.es_uno(m):
                l.append(m)
            g = h
        if not fqx.es_uno(f):
            f = entre_p(fqx,p,f) #funciona el entre_p() de antes
    for x in l: #multiplicamos todos los factores con multiplicidad 1 obtenidos
        sol = fqx.mult(sol, x)
    return sol

def derivar_fq(fqx, f): #hacemos la derivada
    l = list(f)
    l.pop(0)
    for i in range(len(l)):
        l[i] = (fqx.fq.mult(l[i], (i+1,))) #ligero cambio respecto a fpx a la hora de multiplicar
    return fqx.reducir(tuple(l))


# distinct-degree factorization
# input: fqx --> anillo_fq_x
# input: g --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mÃ³nicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mÃ³nicos de h de grado = i, el Ãºltimo hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fqx(fqx, g): #FASE 2
    g = fqx.monico(g)
    l = []
    h = fqx.mod(((), (1,)), g) # x
    while not fqx.es_uno(g): #en la iteracion i, obtenemos el producto de factores irreducibles de grado i  
        h = fqx.pot_mod(h, fqx.fq.q, g) #x^q -x (mod g)
        f = fqx.gcd(fqx.resta(h, ((), (1,))), g)
        if not fqx.es_uno(f):
            l.append(f)
            g = fqx.div(g, f)  #le quitamos los factores de grado i de nuestro polinomio
            h = fqx.mod(h, g)  #redefinimos h
        else:
            l.append(fqx.uno())  #si no tiene, añadimos el 1 para indicarlo y que trabaje eficientemente en la fase 3 
    return l

# equal-degree factorization
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: r --> int
# input: h --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mÃ³nicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mÃ³nicos de grado = r
def eqdegr_fact_fqx(fqx, r, h): #FASE 3
    return tuple(fase33(fqx,r,h,[])) #hacemos una llamada a una funcion auxilar que trabaje con una lista acumuladora 

def fase33(fqx, d: int, gd, lista): 
    if fqx.grado(gd) < d: #no debemos continuar
        return None
    if fqx.grado(gd) == d: #hemos acabado, añadimos nuestro factor a la lista de factores
        lista.append(gd)
        return lista
    #construimos nuestro polinomio aleatorio, de grado menor que nuestro polinomio de entrada
    grado_a = random.randint(1, fqx.grado(gd)-1)
    i = 0
    a = []
    for i in range(grado_a+1):
        a.append(fqx.fq.aleatorio()) 
    a = fqx.reducir(tuple(a))
    aux = fqx.pot_mod(a, (fqx.fq.q**d - 1)//2, gd) #construimos los polinomios A(x), B(x), C(x) tal y como hemos visto en clase
    A = fqx.gcd(a, gd)
    B = fqx.gcd(fqx.resta(aux, fqx.uno()), gd)
    C = fqx.gcd(fqx.suma(aux, fqx.uno()), gd)
    fase33(fqx, d, A, lista), fase33(fqx, d, B, lista), fase33(fqx, d, C, lista) #llamamos recursivamente a fase3
    return lista

# multiplicidad de factor irreducible mÃ³nico
#input: fqx --> anillo_fq_x
# input: f --> polinomio de fqx (objeto opaco) no nulo
# input: u --> polinomio irreducible mÃ³nico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fqx(fqx, f, u):  #dividimos f entre cada factor. Lo hará tantas veces como sea su multiplicidad.
    i = 0
    while fqx.es_igual(fqx.mod(f, u), fqx.cero()):
        f = fqx.div(f, u)
        i += 1
    return i

# factorizaciÃ³n de Cantor-Zassenhaus
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: f --> polinomio de fqx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorizaciÃ³n completa de f en irreducibles mÃ³nicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fqx(fqx, f):                     # mantener esta implementaciÃ³n
    g = sqfree_fact_fqx(fqx, f)
    h = didegr_fact_fqx(fqx, g)
    irreducibles = []
    for r in range(len(h)):
        if fqx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fqx(fqx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fqx(fqx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para aÃ±adir la funciÃ³n de factorizaciÃ³n de Cantor-Zassenhaus
# como un mÃ©todo de la clase anillo_fq_x
cf.anillo_fq_x.factorizar = fact_fqx
