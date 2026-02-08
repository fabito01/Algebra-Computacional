# -*- coding: utf-8 -*-
"""
@author: ALBERTO PEÑA PEÑALVER Y FABIO TORRES MARTINEZ
"""

# Se deben implementar 4 clases: cuerpo_fp, anillo_fp_x, cuerpo_fq, anillo_fq_x
# respetando el esqueleto de abajo. Se pueden aÃ±adir funciones auxiliares, pero
# las indicadas deben estar. Los elementos de cada uno de estos cuerpos/anillos
# son objetos opacos.
import re
import random

'''
Definimos 3 funciones auxiliares externas a la clase, que nos permitirán saber si un número es primo,
ordenar 2 tuplas de mayor a menor longitud y por último factorizar un número entero.
'''


def es_primo(p):  # True si es primo, False si no lo es
    i = 2
    primo = True
    while i*i <= p and primo:
        if p % i == 0:
            primo = False
        else:
            i += 1
    return primo


def mayor(a, b):  # toma 2 tuplas y devuelve primero el mayor longitud y luego el menor longitud
    if len(a) < len(b):
        a, b = b, a
    return a, b


# nos devuelve los primos de la factorización del entero, sin sus exponentes
def factorizacion_int(a):
    i = 2
    l = []
    while a != 1:
        if a % i == 0 and es_primo(i):
            l.append(i)
            a //= i
        else:
            i += 1
    return list(set(l))


class cuerpo_fp:
    def __init__(self, p):  # construye el cuerpo de p elementos Fp = Z/pZ
        if p <= 1 or not es_primo(p):
            raise ValueError("p tiene que ser un número primo")
        self.p = p

    def cero(self):        # devuelve el elemento 0
        return 0

    def uno(self):         # devuelve el elemento 1
        return 1

    def elem_de_int(self, n):  # fabrica el elemento dado por la clase de n
        return n % self.p

    def elem_de_str(self, s):  # fabrica el elemento a partir de un string (parser)
        return int(s) % self.p

    def conv_a_int(self, a):  # devuelve un entero entre 0 y p-1
        return a % self.p

    def conv_a_str(self, a):  # pretty-printer
        return str(a % self.p)

    def suma(self, a, b):     # a+b
        return (a+b) % self.p

    def inv_adit(self, a):    # -a
        return (-a) % self.p

    def mult(self, a, b):     # a*b
        return (a*b) % self.p

    def pot(self, a, k):      # a^k (k entero)
        if k == 0:
            return 1
        # si el k fuese negativo, entonces calculamos (a^-1)^k= (a^k)^-1
        if k < 0:
            a = self.inv_mult(a)
            k = -k
        aux = self.pot(a, k//2)  # potencicacion rapida
        aux = self.mult(aux, aux)
        if k % 2 == 1:
            aux = self.mult(aux, a)
        return aux

    def inv_mult(self, a):    # a^(-1)
        if a == 0:
            raise ValueError("No tiene inverso multiplicativo")
        return self.pot(a, self.p - 2)

    def es_cero(self, a):     # a == 0
        return a == self.cero()

    def es_uno(self, a):      # a == 1
        return a == self.uno()

    def es_igual(self, a, b):  # a == b
        return a == b

    def aleatorio(self):      # fabrica un elemento aleatorio con prob uniforme
        return random.randint(0, self.p - 1)

    def tabla_suma(self):     # devuelve la matriz de pxp (lista de listas) de la suma
        A = [[self.suma(i, j) for j in range(0, self.p)]
             for i in range(0, self.p)]
        "para que lo ponga bonito"
        # for fila in A:
        #   print(" ".join(f"{x:6.0f}" for x in fila))
        return A

    def tabla_mult(self):     # devuelve la matriz de pxp (lista de listas) de la mult
        A = [[self.mult(i, j) for j in range(0, self.p)]
             for i in range(0, self.p)]
        "para que lo ponga bonito"
        # for fila in A:
        #    print(" ".join(f"{x:6.0f}" for x in fila))
        return A

    def tabla_inv_adit(self):  # devuelve una lista de long p con los inv_adit
        return [self.inv_adit(i) for i in range(0, self.p)]

    # devuelve una lista de long p con los inv_mult (en el i­ndice 0 pone un '*')
    def tabla_inv_mult(self):
        return ["*"] + [self.inv_mult(i) for i in range(1, self.p)]

    def cuadrado_latino(self, a):  # cuadrado latino a*i+j (con a != 0)
        if a % self.p == 0:
            raise ValueError("a no puede ser 0")
        A = [[self.suma(self.mult(a, i), j)for j in range(0, self.p)] for i in range(0, self.p)] 
        return A


class anillo_fp_x:
    def __init__(self, fp, var='x'):  # construye el anillo Fp[var]
        self.fp = fp
        self.var = var

    def cero(self):                 # ()
        return tuple()

    def uno(self):                  # (1,)
        return (self.fp.uno(),)

    def reducir(self, a):
        b = list(a)
        # la primera condición por si fuese el pol nulo
        while len(b) > 0 and (b[-1]) == self.fp.cero():
            b.pop(-1)
        return tuple(b)

    # fabrica un polinomio a partir de la tupla (a0, a1, ...) #lista de objetos opacos de p
    def elem_de_tuple(self, a):
        return self.reducir(tuple([self.fp.conv_a_int(x) for x in a]))

    def elem_de_int(self, a):       # fabrica un polinomio a partir de los di­gitos de a en base p
        pol = []
        if a < 0:
            raise ValueError('Tiene que ser un entero positivo')
        if a == 0:
            return self.cero()
        while a > 0:
            pol.append(a % self.fp.p)
            a = a//self.fp.p
        return tuple(pol)

    def elem_de_str(self, s):  # Hemos usado ia    # fabrica un polinomio a partir de un string (parser) (coherente con el pretty-printer)
        s = s.replace(" ", "")  # quitamos espacios
        s = s.replace("*", "")  # quitamos multiplicaciones explícitas
        var = self.var
        # Diccionario para guardar coeficientes por exponente
        coefs = {}
        terminos = re.findall(r'[+-]?[^+-]+', s)
        for t in terminos:
            if t == '':
                continue
            signo = 1
            if t[0] == '+':
                t = t[1:]
            elif t[0] == '-':
                t = t[1:]
                signo = -1
            if var in t:
                if '^' in t:
                    coef_str, exp_str = t.split(var + '^')
                    exp = int(exp_str)
                elif t.endswith(var):
                    coef_str = t[:-1]
                    exp = 1
                else:
                    coef_str = '1'
                    exp = 1
                coef = self.fp.elem_de_str(
                    coef_str) if coef_str else self.fp.elem_de_str('1')
            else:
                coef = self.fp.elem_de_str(t)
                exp = 0
            # aplicar signo
            # asegurar tipo del cuerpo
            coef = self.fp.suma(coef, self.fp.elem_de_str('0'))
            if signo == -1:
                coef = self.fp.resta(self.fp.elem_de_str('0'), coef)
            coefs[exp] = coef
        # Construir tupla de coeficientes ordenada por grado
        grado_max = max(coefs.keys()) if coefs else 0
        resultado = tuple(coefs.get(i, self.fp.elem_de_str('0'))
                          for i in range(grado_max + 1))
        return self.reducir(resultado)

    # devuelve la tupla de coeficientes #esos coefs luego los interpretaré con la clase fp
    def conv_a_tuple(self, a):
        return self.reducir(a)

    def conv_a_int(self, a):        # devuelve el entero correspondiente al polinomio
        s = 0
        m = len(a)
        for i in range(0, m):
            s += (a[i]*self.fp.p**i)
        return s

    def conv_a_str(self, a):  # hemos usado ia       # pretty-printer
        partes = []
        for i, a in enumerate(a):
            if a == 0:
                continue
            if i == 0:
                parte = f"{a}"
            elif i == 1:
                parte = f"{'' if a == 1 else '-' if a == -1 else a}x"
            else:
                parte = f"{'' if a == 1 else '-' if a == -1 else a}x^{i}"
            partes.append(parte)
        if not partes:
            return "0"
        s = " + ".join(partes)
        s = s.replace("+ -", "- ")
        return s

    def suma(self, a, b):  # a+b
        m, n = mayor(a, b)
        ln = len(n)
        g = [self.fp.suma(m[i], n[i]) for i in range(ln)]
        for j in range(ln, len(m)):
            g.append(m[j] % self.fp.p)
        return self.reducir(tuple(g))

    def resta(self, a, b):
        return self.suma(a, self.inv_adit(b))

    def inv_adit(self, a):          # -a
        return tuple([(-a[i]) % self.fp.p for i in range(0, len(a))])

    def mult(self, a, b):           # a*b
        if a == self.cero() or b == self.cero():
            return self.cero()
        d1, d2 = len(a), len(b)
        c = (d1+d2-1)*[0]  # calculamos de antemano el grado del producto
        for i in range(d1):
            for j in range(d2):
                c[i+j] = (c[i+j] + (a[i]*b[j])) % self.fp.p
        # no hace falta reducir porque sabiamos el grado del producto
        return tuple(c)

    def mult_por_escalar(self, a, e):  # a*e (con e en Z/pZ)
        if self.fp.es_cero(e):  # descatamos casos triviales
            return self.cero()
        elif self.fp.es_uno(e):
            return a
        else:
            # multiplicamos por el pol de grado 0 con coef e
            return self.reducir(self.mult(a, (e,)))

    def divmod(self, a, b):        # devuelve q,r tales que a=bq+r y deg(r)<deg(b)
        a = self.reducir(a)
        b = self.reducir(b)
        la = len(a)
        lb = len(b)
        if la < lb:  # si el grado de a es menor que el de b
            return self.cero(), a
        q = [0]*(la-lb+1)  # inicializamos cociente
        r = a
        lamda = self.fp.inv_mult(b[-1])
        while len(r) >= len(b):  # mientras el grado de r sea mayor o igual que el de b
            # coeficiente por el que multiplicaremos b (ya hemos invertido b[-1])
            jan = lamda*(r[-1])
            # lo añadimos al cociente en su posición correspondiente
            q[len(r)-lb] = jan % self.fp.p
            r = self.reducir(self.resta(r, self.mult(self.mult_por_escalar(
                b, jan), tuple([0]*((len(r)-lb))+[1]))))  # iteramos la division
        return tuple(q), r

    def div(self, a, b):           # q
        return self.divmod(a, b)[0]

    def mod(self, a, b):           # r
        return self.divmod(a, b)[1]

    def grado(self, a):            # deg(a)
        return len(a)-1

    def monico(self, a):
        if a[-1] == self.fp.uno():
            return a
        else:
            return self.mult_por_escalar(a, self.fp.inv_mult(a[-1]))

    def gcd(self, a, b):           # devuelve g = gcd(a,b) mÃ³nico
        a, b = mayor(a, b)
        if a == b:
            return a
        if (a or b) == self.uno():
            return self.uno()
        if a == self.cero():
            return self.monico(b)
        if b == self.cero():
            return self.monico(a)
        r = self.mod(a, b)  # tomamos el resto
        while (r != self.cero()):  # si necesitamos volver a iterar
            # calculamos el siguiente resto que usaremos para hacer en la siguiente el gcd(b,r1)
            r1 = self.mod(b, r)
            b = r
            r = r1
        return self.monico(b)

    def gcd_ext(self, a, b):       # devuelve g,x,y tales que g=ax+by, g=gcd(a,b) monico
        if a == self.cero():
            return self.mult_por_escalar(b, self.fp.inv_mult(b[-1])), self.cero(), (self.fp.inv_mult(b[-1]),)
        if b == self.cero():
            return self.mult_por_escalar(a,self.fp.inv_mult(a[-1])), (self.fp.inv_mult(a[-1]),), self.cero()
        q, r = self.divmod(a, b)  # r = a-bq
        # g = bx + ry= bx + (a -bq)y = ay + b(x -qy)
        g, x, y = self.gcd_ext(b, r)
        x = self.resta(x, self.mult(q, y))
        return self.monico(g), y, x

    def inv_mod(self, a, b):       # devuelve x tal que ax = 1 mod b -> ax + bk = 1
        g, x, y = self.gcd_ext(a, b)
        if not self.es_uno(g):
            raise Exception(a, b, g, 'No es invertible')
        return x

    def pot_mod(self, a, k, b):    # a^k mod b -> a^k = x (mod b) -> a^k = x + bm -> a^k +bn = x
        if k == 0:
            return self.uno()
        if k < 0:
            k = -k
            a = self.inv_mod(a, b)
        aux = self.pot_mod(a, k//2, b)  # potenciacion rapida
        aux = self.mult(aux, aux)
        if k % 2 == 1:
            aux = self.mult(aux, a)
        return self.mod(aux, b)

    def es_cero(self, a):          # a == ()
        return a == self.cero()

    def es_uno(self, a):           # a == (1,)
        return a == self.uno()

    def es_igual(self, a, b):      # a == b
        return a == b

    def es_irreducible(self, f):   # test de irreducibilidad de Rabin
        n = self.grado(f)
        if n <= 1:  # los polinomios de grado 1 son irreducibles
            return True
        pi = factorizacion_int(n)
        p = self.fp.p
        g = self.pot_mod((0, 1), p**n, f)
        if not self.es_igual(g, (0, 1)):  # no hace falta hacer modulo a x pues deg(f)>=2
            # (x^p^n = x -> x^p^n - x = 0 (mod f)) si no es cero, entonces es porque no le divide f, luego no cumple la primera premisa de Rabin
            return False
        rabin = True
        i = 0
        while i < len(pi) and rabin:
            di = n//pi[i]
            h = self.resta(self.pot_mod((0, 1), p**di, f),
                           (0, 1))  # (x^p^di - x)
            if not self.es_uno(self.gcd(f, h)):
                rabin = False  # si gcd no es uno ya sabemos que no cumple la segunda premisa de Rabin
            i += 1
        return rabin
    

class cuerpo_fq:
    def __init__(self, fp, g, var='ç'):  # construye el cuerpo Fp[var]/<g(var)>
        # g es objeto fabricado por fp
        self.fpx = anillo_fp_x(fp, 'x')
        if self.fpx.grado(g) <= 0 or not self.fpx.es_irreducible(g):
            raise ValueError(
                'g debe ser polinomio de grado mayor o igual q 1 e irreducible')
        self.g = g
        self.n = self.fpx.grado(g)
        self.q = self.fpx.fp.p**self.n
        self.var = var

    def cero(self):                # 0
        return self.fpx.cero()

    def uno(self):                 # 1
        return self.fpx.uno()

    def elem_de_tuple(self, a):    # fabrica elemento a partir de tupla de coeficientes
        return self.fpx.mod(self.fpx.elem_de_tuple(a), self.g)

    def elem_de_int(self, a):      # fabrica elemento a partir de entero
        return self.fpx.mod(self.fpx.elem_de_int(a), self.g)

    def elem_de_str(self, s):  # hemos usado ia  # fabrica elemento parseando string
        # quitar espacios y tratar negativos mod p
        s = s.replace(" ", "").replace("-", "+2*")
        base = self.fpx.fp.p
        var = self.var[0]  # ej: 'ç'
        n = self.n        # grado del polinomio irreducible
        # Si el usuario mete directamente un número
        if s.isdigit():
            return (int(s) % base, ) + (0,) * (n - 1)
        # Inicializamos el vector de coeficientes
        coef = [0] * n
        # Separamos por '+'
        terminos = s.split('+')
        for t in terminos:
            if t == '':
                continue
            if var not in t:
                # término constante
                coef[0] = int(t) % base
            else:
                # término con variable
                if '^' in t:
                    a, exp = t.split('^')
                    exp = int(exp)
                    a = a.replace(var, '')
                    if a == '':
                        c = 1
                    elif a == '*':
                        c = 1
                    else:
                        c = int(a.replace('*', '')) % base
                    if exp < n:
                        coef[exp] = c
                else:
                    # término como '2ç' o 'ç'
                    a = t.replace(var, '')
                    if a == '':
                        c = 1
                    elif a == '*':
                        c = 1
                    else:
                        c = int(a.replace('*', '')) % base
                    coef[1] = c
        return tuple(coef)

    def conv_a_tuple(self, a):     # devuelve tupla de coeficientes sin ceros "extra"
        return self.fpx.mod(self.fpx.conv_a_tuple(a), self.g)

    def conv_a_int(self, a):       # devuelve el entero correspondiente
        return self.fpx.conv_a_int(self.fpx.mod(a, self.g))

    def conv_a_str(self, a):  # hemos usado ia     # pretty-printer
        var = self.var[0]
        base = self.fpx.fp.p
        partes = []
        for i, c in enumerate(a):
            if c % base == 0:
                continue
            if i == 0:
                partes.append(f"{c}")
            elif i == 1:
                if c == 1:
                    partes.append(f"{var}")
                else:
                    partes.append(f"{c}{var}")
            else:
                if c == 1:
                    partes.append(f"{var}^{i}")
                else:
                    partes.append(f"{c}{var}^{i}")
        if not partes:
            return "0"
        return " + ".join(partes)

    def suma(self, a, b):          # a+b
        return self.fpx.mod(self.fpx.suma(a, b), self.g)

    def inv_adit(self, a):         # -a
        return self.fpx.mod(self.fpx.inv_adit(a), self.g)

    def mult(self, a, b):          # a*b
        return self.fpx.mod(self.fpx.mult(a, b), self.g)

    def pot(self, a, k):           # a^k (k entero)
        return self.fpx.pot_mod(a, k, self.g)

    def inv_mult(self, a):         # a^(-1)
        return self.fpx.pot_mod(a, -1, self.g)

    def es_cero(self, a):          # a == 0
        return a == self.cero()

    def es_uno(self, a):           # a == 1
        return a == self.uno()

    def es_igual(self, a, b):      # a == b
        return a == b

    def aleatorio(self):           # devuelve un elemento aleatorio con prob uniforme
        l = random.randint(1, self.n - 1)
        a = [random.randint(0, self.fpx.fp.p - 1) for _ in range(l)]
        return self.fpx.reducir(tuple(a))

    def elementos(self):  # nos devuelve una lista con todos los elementos de Fq
        v = [[]]
        for i in range(self.n):
            w = []
            for i in v:
                for j in range(self.fpx.fp.p):
                    w.append(tuple([j]+list(i)))
            v = w
        return v

    def tabla_suma(self):          # matriz de qxq correspondiente a la suma (con la notacion int)
        v = self.elementos()
        A = [[0]*self.q for _ in range(self.q)]
        for i in range(self.q):
            for j in range(self.q):
                A[i][j] = (self.fpx.conv_a_int(self.suma(v[i], v[j])))
        return A

    def tabla_mult(self):          # matriz de qxq correspondiente a la mult (con la notacion int)
        v = self.elementos()
        A = [[0]*self.q for _ in range(self.q)]
        for i in range(self.q):
            for j in range(self.q):
                A[i][j] = (self.fpx.conv_a_int(self.mult(v[i], v[j])))
        return A

    def tabla_inv_adit(self):      # lista de inv_adit (con la notacion int)
        return [self.fpx.conv_a_int(self.fpx.inv_adit(pol)) for pol in self.elementos()]

    def tabla_inv_mult(self):      # lista de inv_mult (con la notacion int)
        u = self.elementos()
        _ = u.pop(0)
        return ["*"] + [self.conv_a_int(self.fpx.inv_mod(self.fpx.reducir(pol), self.g)) for pol in u]

    def cuadrado_latino(self, a):  # cuadrado latino para a != 0 (con notacion int)
        if self.fpx.mod(a, self.g) == self.cero():
            raise ValueError("a no puede ser 0")
        w = self.elementos()
        A = [[self.fpx.conv_a_int(self.suma(self.mult(a, i), j))for j in w] for i in w]
        return A


class anillo_fq_x:
    # Fq[var], var debe ser distinta que la de fq
    def __init__(self, fq, var='x'):
        self.fq = fq
        self.var = var

    def cero(self):                # ()
        return self.fq.cero()

    def uno(self):                 # ((1,),)
        return (self.fq.uno(),)

    def elem_de_tuple(self, a):  # fabrica elemento a partir de una tupla
        return self.reducir(tuple([self.fq.elem_de_tuple(k) for k in a]))

    def elem_de_int(self, a):  # fabrica elemento a partir de entero
        pol = []
        if a < 0:
            raise ValueError('Tiene que ser un entero positivo')
        if a == 0:
            return self.cero()
        while a > 0:
            pol.append(a % self.fq.q)
            # a0 + a1*q + a2*q^2+...an*g^n -> (a0,a1,a2,...,an)
            a = a//self.fq.q
        # (a0,a1,a2,...an) -> (elem_de_int(a0), elem_de_int(a1),...,elem_de_int(an))
        c = [self.fq.fpx.elem_de_int(x) for x in pol]
        return self.reducir(tuple(c))

    def elem_de_str(self, s):  # hemos usado ia #parser
        s = s.replace(" ", "")
        var_x = self.var  # normalmente 'x'
        coeficientes = {}
        # --- Separar términos respetando paréntesis ---
        partes = []
        buffer = ""
        nivel = 0
        for ch in s:
            if ch == '(':
                nivel += 1
            elif ch == ')':
                nivel -= 1
            if ch == '+' and nivel == 0:
                partes.append(buffer)
                buffer = ""
            else:
                buffer += ch
        if buffer:
            partes.append(buffer)
        # --- Procesar cada término ---
        for t in partes:
            if not t:
                continue
            if var_x not in t:
                # término constante
                coef = self.fq.elem_de_str(t.strip('()'))
                exp = 0
            else:
                if '^' in t:
                    a, exp = t.split(f'{var_x}^')
                    exp = int(exp)
                elif t.endswith(var_x):
                    a = t[:-len(var_x)]
                    exp = 1
                else:
                    a = t
                    exp = 0
                a = a.strip('*')
                if a in ('', '('):
                    coef = self.fq.elem_de_str('1')
                else:
                    a = a.strip('()')
                    coef = self.fq.elem_de_str(a)
            coeficientes[exp] = coef
        # --- Crear tupla de coeficientes ---
        if not coeficientes:
            return (self.fq.cero(),)
        grado = max(coeficientes.keys())
        res = []
        for i in range(grado + 1):
            coef = coeficientes.get(i, self.fq.cero())
            res.append(coef)
        return tuple(res)

    def conv_a_tuple(self, a):  # devuelve la tupla de coeficientes
        return self.reducir(a)

    def conv_a_int(self, a):    # devuelve el entero correspondiente al polinomio
        s = 0
        for i in range(len(a)):
            s += (self.fq.fpx.conv_a_int(a[i])*self.fq.q**i)
        return s

    def conv_a_str(self, a):  # hemos usado ia #pretty printer
        var_x = self.var
        partes = []
        for i, coef in enumerate(a):
            if coef == tuple(0 for _ in coef):  # elemento 0 del cuerpo
                continue
            coef_str = self.fq.conv_a_str(coef)
            if i == 0:
                partes.append(f"{coef_str}")
            elif i == 1:
                partes.append(f"({coef_str}){var_x}")
            else:
                partes.append(f"({coef_str}){var_x}^{i}")
        if not partes:
            return "0"
        return " + ".join(partes)

    def reducir(self, a):
        b = list(a)
        for i in range(len(b)):
            b[i] = self.fq.fpx.reducir(b[i])
        # la primera condición por si fuese el pol nulo
        while len(b) > 0 and (b[-1]) == self.fq.cero():
            b.pop(-1)
        return tuple(b)

    def suma(self, a, b):
        a, b = mayor(a, b)
        lb = len(b)
        g = [self.fq.suma(a[i], b[i]) for i in range(lb)]
        for j in range(lb, len(a)):
            g.append(self.fq.fpx.mod(a[j], self.fq.g))
        return self.reducir(tuple(g))

    def resta(self, a, b):
        return self.suma(a, self.inv_adit(b))

    def inv_adit(self, a):
        return self.reducir(tuple([self.fq.inv_adit(k) for k in a]))

    def mult(self, a, b):
        if a == self.cero() or b == self.cero():
            return self.cero()
        d1, d2 = len(a), len(b)
        # calculamos la longitud del producto previamente
        c = (d1+d2-1)*[self.cero()]
        for i in range(d1):
            for j in range(d2):
                c[i+j] = self.fq.fpx.reducir(self.fq.suma(c[i+j],
                                             self.fq.mult(a[i], b[j])))
        return tuple(c)

    def mult_por_escalar(self, a, e):
        return self.reducir(tuple([self.fq.mult(k, e) for k in a]))

    def divmod(self, a, b):  # analogo al de fpx
        a = self.reducir(a)
        b = self.reducir(b)
        la = len(a)
        lb = len(b)
        if la < lb:
            return self.cero(), a
        q = [tuple() for _ in range(la-lb+1)]
        r = a
        lamda = self.fq.inv_mult(b[-1])
        while len(r) >= len(b):
            jan = self.fq.fpx.mult(lamda, (r[-1]))
            q[len(r)-lb] = self.fq.elem_de_tuple(jan)
            r = self.reducir(self.resta(r, self.mult(self.mult_por_escalar(
                b, jan), (tuple([() for _ in range(len(r)-lb)]+[(1,)])))))
        return tuple(q), r

    def div(self, a, b):  # q
        return self.divmod(a, b)[0]

    def mod(self, a, b):  # r
        return self.divmod(a, b)[1]

    def grado(self, a):
        return len(a) - 1

    def monico(self, a):
        if a[-1] == self.fq.uno():
            return a
        else:
            return self.mult_por_escalar(a, self.fq.inv_mult(a[-1]))

    def gcd(self, a, b):  # identico al de fpx
        a, b = mayor(a, b)
        if a == b:
            return a
        if (a or b) == self.uno():
            return self.uno()
        if a == self.cero():
            return self.monico(b)
        if b == self.cero():
            return self.monico(a)
        r = self.mod(a, b)
        while (r != self.cero()):
            r1 = self.mod(b, r)
            b = r
            r = r1
        return self.monico(b)

    def gcd_ext(self, a, b):     # devuelve g,x,y tales que g=ax+by, g=gcd(a,b) monico
        if a == self.cero():
            return self.mult_por_escalar(b, self.fq.inv_mult(b[-1])), self.cero(), (self.fq.inv_mult(b[-1]),)
        if b == self.cero():
            return self.mult_por_escalar(a, self.fq.inv_mult(a[-1])), (self.fq.inv_mult(a[-1]),), self.cero()
        q, r = self.divmod(a, b)  # r = a-bq
        # g = bx + ry= bx + (a -bq)y = ay + b(x -qy)
        g, x, y = self.gcd_ext(b, r)
        x = self.resta(x, self.mult(q, y))
        return self.monico(g), y, x

    def inv_mod(self, a, b):  # devuelve x tal que ax = 1 mod b -> ax + bk = 1
        g, x, y = self.gcd_ext(a, b)
        if not self.es_uno(g):
            raise Exception(a, b, g, 'No es invertible')
        return x

    def pot_mod(self, a, k, b):  # identico al de fpx
        if k == 0:
            return self.uno()
        if k < 0:
            k = -k
            a = self.inv_mod(a, b)
        aux = self.pot_mod(a, k//2, b)
        aux = self.mult(aux, aux)
        if k % 2 == 1:
            aux = self.mult(aux, a)
        return self.mod(aux, b)

    def es_cero(self, a):  # a == ()
        return a == self.cero()

    def es_uno(self, a):  # a == ((1,),)
        return a == self.uno()

    def es_igual(self, a, b):  # a == b
        return a == b

    def es_irreducible(self, f):  # analogo a fpx
        n = self.grado(f)
        if n <= 1:  # los polinomios de grado 1 son irreducibles
            return True
        pi = factorizacion_int(n)
        q = self.fq.q
        g = self.pot_mod(((), self.fq.uno()), q**n, f)
        # no se necesita hacer modulo a x porque sabemos ya que deg(f)>=2
        if not self.es_igual(g, ((), self.fq.uno())):
            return False
        rabin = True
        i = 0
        while i < len(pi) and rabin:
            di = n//pi[i]
            h = self.resta(self.pot_mod(((), self.fq.uno()),
                           q**di, f), ((), self.fq.uno()))
            if not self.es_uno(self.gcd(f, h)):
                rabin = False
            i += 1
        return rabin
    
    
    

