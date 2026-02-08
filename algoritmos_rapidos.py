# para esta actividad necesitamos la primera (cuerpos_finitos.py) completa
import cuerpos_finitos as cf

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio (objeto opaco creado por fpx)
# output: f*g calculado usando el método de Karatsuba



def fp_x_mult_karatsuba(fpx, f, g):
    n = len(f)
    m = len(g)
    if n < m:
        f, g, n, m = g, f, m, n
    #añadimos mult escuela?
    mitad = n//2
    if n <= 20:
        c = fpx.mult(f, g)
    elif m >= mitad:
        #dividimos f
        flow = f[:mitad]
        fpx.reducir(flow)
        fhi = f[mitad:]
        #dividimos g
        glow = g[:mitad]
        fpx.reducir(glow)
        ghi = g[mitad:]
       
       
        c0 = fp_x_mult_karatsuba(fpx, flow, glow)
        c2 = fp_x_mult_karatsuba(fpx, fhi, ghi)
        fsum = fpx.suma(flow, fhi)
        gsum = fpx.suma(glow, ghi)
       
        c1 = fpx.resta(fp_x_mult_karatsuba(fpx, fsum, gsum),fpx.suma(c0, c2))
        c = [fpx.fp.cero()]*(n+m-1)
        for i in range (len(c0)):
            c[i] = c0[i]
        for i in range (len(c1)):
            c[i+mitad] = fpx.fp.suma(c[i+mitad],c1[i]) #hace módulo????
        for i in range (len(c2)):
            c[i+2*mitad] = fpx.fp.suma(c[i+2*mitad],c2[i]) #hace módulo????
    else:
        flow = f[:mitad]
        fpx.reducir(flow)
        fhi = f[mitad:]
        c0 = fpx.elem_de_tuple(fp_x_mult_karatsuba(fpx, flow, g))
        c1 = fp_x_mult_karatsuba(fpx, fhi, g)
        c = [fpx.fp.cero()]* (n+m-1)
        for i in range(len(c0)):
            c[i]= c0[i]
        for i in range(len(c1)):
            c[i+mitad]= fpx.fp.suma(c[i+mitad],c1[i])
    return fpx.reducir(tuple(c)) # reducir o no reducir, esa es la cuestion


# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.mult_fast = fp_x_mult_karatsuba


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_inf_vec(fp, n, a, b):
    sol = fp_x_mult_karatsuba(cf.anillo_fp_x(fp,'x'), a, b)
    return sol[:n]



# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_sup_vec(fp, n, a, b):
    a = list(a)
    a.reverse()
    a = tuple(a)
    sol = fp_x_mult_karatsuba(cf.anillo_fp_x(fp,'x'), a, b)
    return sol[n-1:]


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fp (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_vec(fp, n, a, b):
    row_U = a[:n]
    res_U = fp_toep_sup_vec(fp, n, row_U, b)
    col_L = (fp.cero(),) + a[n:]
    res_L = fp_toep_inf_vec(fp, n, col_L, b)
    sol = []
    for i in range(n):
        v = fp.suma(res_L[i],res_U[i])
        sol.append(v)
       
    return tuple(sol)

#print('a',fp_toep_vec(cf.cuerpo_fp(5), 3, (1,2,3,4,5), (3,2,3)))
#print('b',fp_toep_vec(cf.cuerpo_fp(5), 5, (1,4,2,1,3,2,3,1,4), (2,2,1,3,4)))

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fp_toep_inf_inv(fp, n, a):
    if n == 1:
        return (fp.inv_mult(a[0]),)
    m = (n + 1) // 2


    inv_top = fp_toep_inf_inv(fp, m, a[:m])
    
    inv_A1_lleno = list(inv_top) + [fp.cero()] * (n - m)

    prod = fp_toep_inf_vec(fp, n, a, tuple(inv_A1_lleno))
    residuo = prod[m:] 


    correction = fp_toep_inf_vec(fp, n - m, inv_top[:n-m], residuo)

    inv_bottom = []
    for val in correction:
        inv_bottom.append(fp.inv_adit(val))
    return inv_top + tuple(inv_bottom)   
    
    # if n == 1:
    #     return (fp.inv_mult(a[0]),)
    # m = (n + 1) // 2
    # inv_A1 = fp_toep_inf_inv(fp, m, a[:m])

    # #inv_A1_lleno = list(inv_A1) + [fp.cero()] * (n - m)
    # if n%2 != 0 :
    #     prod = fp_toep_inf_vec(fp,n-m,a[m:] + (fp.cero(),),inv_A1)
    #     resul = fp_toep_inf_vec(fp, n - m, inv_A1, prod)
    #     for val in resul:
    #         inv_A1 = inv_A1 + (fp.inv_adit(val),)
    #     return inv_A1

    # else:
    #     prod = fp_toep_inf_vec(fp, n-m, a[m:],inv_A1)
    #     resul = list(fp_toep_inf_vec(fp, n - m, inv_A1, prod))
    #     for val in resul[:(n-m-1)]:
    #         inv_A1 = inv_A1 + (fp.inv_adit(val),)
    #     return inv_A1
    
    
    # if n == 1:
    #     return (fp.inv_mult(a[0]),)
    # m = (n + 1) // 2
    # inv_A1 = fp_toep_inf_inv(fp, m, a[:m])
    # #inv_A1_lleno = list(inv_A1) + [fp.cero()] * (n - m)
    # if n%2 != 0 :
    #     prod = fp_toep_inf_vec(fp,m,a[m:] + (fp.cero(),),inv_A1)
    #     #print(prod)
    #     resul = fp_toep_inf_vec(fp, m, inv_A1, prod)
    #     #print(resul)
    #     for val in resul[:n]:
    #         inv_A1 = inv_A1 + (fp.inv_adit(val),)
    #     return inv_A1


    # else:
    #     prod = fp_toep_inf_vec(fp, n-m, a[m:],inv_A1)
    #     resul = list(fp_toep_inf_vec(fp, m, inv_A1, prod))
    #     for val in resul:
    #         inv_A1 = inv_A1 + (fp.inv_adit(val),)
    #     return inv_A1
        

    #residuo = prod[m:] # Tomamos solo la parte inferior
    
    

# fp = cf.cuerpo_fp(5)
# a = (1, 2, 0)
# inv_a = fp_toep_inf_inv(fp, 3, a)
# print(f"Columna inversa calculada: {inv_a}")




# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fp_toep_sup_inv(fp, n, a):
    if n == 1:
        return (fp.inv_mult(a[0]),)

    m = (n + 1) // 2

    inv_A1 = fp_toep_sup_inv(fp, m, a[:m])

    inv_A1_lleno = list(inv_A1) + [fp.cero()] * (n - m)
    
    prod = fp_toep_inf_vec(fp, n, a, tuple(inv_A1_lleno))
    
    residuo = prod[m:]

    correction = fp_toep_inf_vec(fp, n - m, inv_A1[:n-m], residuo)

    inv_right = []
    for val in correction:
        inv_right.append(fp.inv_adit(val))

    return inv_A1 + tuple(inv_right)
    
# fp = cf.cuerpo_fp(5)
# a = (1,2,3,8,5)
# inv_a = fp_toep_sup_inv(fp, 5, a)
# print(f"Columna inversa calculada: {inv_a}")   
   
   
# fp = cf.cuerpo_fp(5)
# a = (3,2,3)
# inv_a = fp_toep_sup_inv(fp, 3, a)
# print(f"Columna inversa calculada: {inv_a}")   
   
      
    
  
    
    
    

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio no nulo (objeto opaco creado por fpx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fp_x_divmod(fpx, f, g):
    f = fpx.reducir(f)
    g = fpx.reducir(g)
    n = len(f) - 1 # Grado de f
    m = len(g) - 1 # Grado de g

    # Caso grado f < grado g
    if n < m:
        return (fpx.fp.cero(),), f

    # tamaño del cociente q
    k = n - m + 1


    f_rev = list(f)
    f_rev.reverse()  #(a0,a1,...,an) --> (an,an-1,..,a1,a0)
    
    g_rev = list(g) 
    g_rev.reverse() #(b0,b1,...,bm) --> (bm,bm-1,...,b0)


    col_g = g_rev[:k]
    if len(col_g) < k:
        col_g = col_g + [fpx.fp.cero()] * (k - len(col_g))
    
    vec_f = f_rev[:k]

    inv_col_g = fp_toep_inf_inv(fpx.fp, k, tuple(col_g))

    q_rev = fp_toep_inf_vec(fpx.fp, k, inv_col_g, tuple(vec_f))

    q_list = list(q_rev)
    q_list.reverse()
    q = tuple(q_list)

    gq = fp_x_mult_karatsuba(fpx, g, q)
    
    r = fpx.resta(f, gq)

    q = fpx.reducir(q)
    r = fpx.reducir(r)

    return q, r

# fpx= cf.anillo_fp_x(cf.cuerpo_fp(7),'x')
# f = (5, 1, 0, 2, 3)
# g = (4, 0, 1)
# print('a',fp_x_divmod(fpx,f,g))

# fpx= cf.anillo_fp_x(cf.cuerpo_fp(3),'x')
# f = (1,1,0,2,1)
# g = (1, 0, 1)
# print('b',fp_x_divmod(fpx,f,g))


# fpx= cf.anillo_fp_x(cf.cuerpo_fp(5),'x')
# f = (2,1,4,3)
# g = (1,2)
# print('c',fp_x_divmod(fpx,f,g))


# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.divmod_fast = fp_x_divmod

# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# utilizar el algoritmo de Cooley-Tuckey
def fp_fft(fp, g, k, a):
    n = 2**k
    if n == 1:
        return a
    
    half_n = n // 2
    p_par = a[::2]
    p_imp = a[1::2]
    
    g2 = fp.mult(g, g)
    
    t_par = fp_fft(fp, g2, k - 1, p_par)
    t_imp = fp_fft(fp, g2, k - 1, p_imp)
    
    vec = [0] * n
    pot = fp.uno()
    for i in range(half_n):
        aux = fp.mult(pot, t_imp[i])
        vec[i] = fp.suma(t_par[i], aux)
        vec[i + half_n] = fp.suma(t_par[i], fp.inv_adit(aux))
        pot = fp.mult(pot, g)
        
    return tuple(vec)

# fp = cf.cuerpo_fp(97)
# fpx = cf.anillo_fp_x(fp,'x')

# pol1 = tuple([1, 2, 3] + [0]*13) 

# tra1 = fp_fft(fp, 8, 4, pol1)
# print(tra1)
# print(f"Polinomio entrada (len 16): {pol1}")
# print(f"Resultado Transformada: {tra1}")







# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fp_ifft(fp, g, k, a):
    n = 2**k
    pol =fp_fft(fp,fp.inv_mult(g),k,a)
    pol=list(pol)
    for i in range(n):
        pol[i] = fp.mult(pol[i],fp.inv_mult(n)) 
    return tuple(pol)


# pol2 = tuple([1,1] + [0]*14)

# pol = fp_x_mult_karatsuba(fpx,pol1,pol2)
# print(pol)
# tra2 = fp_fft(fp,8,4,pol2)
# print(tra2)
# # tra1 = list(tra1)
# # tra2 = list(tra2)
# tra = [fp.mult(tra1[i],tra2[i]) for i in range(16)]
# tra = tuple(tra)
# # print(tra)
# xxx = fp_ifft(fp,8,4,tra)
# print(xxx)




# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio (objeto opaco creado por fqx)
# output: f*g calculado usando el método de Karatsuba
def fq_x_mult_karatsuba(fqx, f, g):
    n = len(f)
    m = len(g)
    if n < m:
        f, g, n, m = g, f, m, n
    mitad = n//2
    if n <= 20:
        c = fqx.mult(f, g)
    elif m >= mitad:
        #dividimos f
        flow = f[:mitad]
        fqx.reducir(flow)
        fhi = f[mitad:]
        #dividimos g
        glow = g[:mitad]
        fqx.reducir(glow)
        ghi = g[mitad:]
       
       
        c0 = fq_x_mult_karatsuba(fqx, flow, glow)
        c2 = fq_x_mult_karatsuba(fqx, fhi, ghi)
        fsum = fqx.suma(flow, fhi)
        gsum = fqx.suma(glow, ghi)
       
        c1 = fqx.resta(fq_x_mult_karatsuba(fqx, fsum, gsum),fqx.suma(c0, c2))
        c = [fqx.fq.cero()]*(n+m-1)
        for i in range (len(c0)):
            c[i] = c0[i]
        for i in range (len(c1)):
            c[i+mitad] = fqx.fq.suma(c[i+mitad],c1[i]) #hace módulo????
        for i in range (len(c2)):
            c[i+2*mitad] = fqx.fq.suma(c[i+2*mitad],c2[i]) #hace módulo????
    else:
        flow = f[:mitad]
        fqx.reducir(flow)
        fhi = f[mitad:]
        c0 = fqx.elem_de_tuple(fq_x_mult_karatsuba(fqx, flow, g))
        c1 = fq_x_mult_karatsuba(fqx, fhi, g)
        c = [fqx.fq.cero()]* (n+m-1)
        for i in range(len(c0)):
            c[i]= c0[i]
        for i in range(len(c1)):
            c[i+mitad]= fqx.fq.suma(c[i+mitad],c1[i])
    return fqx.reducir(tuple(c)) #reducir o no reducir

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.mult_fast = fq_x_mult_karatsuba



# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_inf_vec(fq, n, a, b):
    sol = fq_x_mult_karatsuba(cf.anillo_fq_x(fq,'x'), a, b)
    return sol[:n]

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_sup_vec(fq, n, a, b):
    a = list(a)
    a.reverse()
    a = tuple(a)
    sol = fq_x_mult_karatsuba(cf.anillo_fq_x(fq,'x'), a, b)
    return sol[n-1:]

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fq (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_vec(fq, n, a, b):
    row_U = a[:n]
    res_U = fq_toep_sup_vec(fq, n, row_U, b)
   
    col_L = (fq.cero(),) + a[n:]
    res_L = fq_toep_inf_vec(fq, n, col_L, b)
    sol = []
    for i in range(n):
        v = fq.suma(res_L[i],res_U[i])
        sol.append(v)
       
    return tuple(sol)





# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fq_toep_inf_inv(fq, n, a):
    if n == 1:
        return (fq.inv_mult(a[0]),)
    m = (n + 1) // 2


    inv_top = fq_toep_inf_inv(fq, m, a[:m])
    
    inv_A1_lleno = list(inv_top) + [fq.cero()] * (n - m)

    prod = fq_toep_inf_vec(fq, n, a, tuple(inv_A1_lleno))
    residuo = prod[m:] 


    correction = fq_toep_inf_vec(fq, n - m, inv_top[:n-m], residuo)

    inv_bottom = []
    for val in correction:
        inv_bottom.append(fq.inv_adit(val))
    return inv_top + tuple(inv_bottom)   

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fq_toep_sup_inv(fq, n, a):
    if n == 1:
        return (fq.inv_mult(a[0]),)

    m = (n + 1) // 2

    inv_A1 = fq_toep_sup_inv(fq, m, a[:m])

    inv_A1_lleno = list(inv_A1) + [fq.cero()] * (n - m)
    
    prod = fq_toep_inf_vec(fq, n, a, tuple(inv_A1_lleno))
    
    residuo = prod[m:]

    correction = fq_toep_inf_vec(fq, n - m, inv_A1[:n-m], residuo)

    inv_right = []
    for val in correction:
        inv_right.append(fq.inv_adit(val))

    return inv_A1 + tuple(inv_right)



# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio no nulo (objeto opaco creado por fqx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fq_x_divmod(fqx, f, g):
    f = fqx.reducir(f)
    g = fqx.reducir(g)
    n = len(f) - 1 # Grado de f
    m = len(g) - 1 # Grado de g

    # Caso grado f < grado g
    if n < m:
        return (fqx.fq.cero(),), f

    # tamaño del cociente q
    k = n - m + 1


    f_rev = list(f)
    f_rev.reverse()  #(a0,a1,...,an) --> (an,an-1,..,a1,a0)
    
    g_rev = list(g) 
    g_rev.reverse() #(b0,b1,...,bm) --> (bm,bm-1,...,b0)


    col_g = g_rev[:k]
    if len(col_g) < k:
        col_g = col_g + [fqx.fq.cero()] * (k - len(col_g))
    
    vec_f = f_rev[:k]

    inv_col_g = fq_toep_inf_inv(fqx.fq, k, tuple(col_g))

    q_rev = fq_toep_inf_vec(fqx.fq, k, inv_col_g, tuple(vec_f))

    q_list = list(q_rev)
    q_list.reverse()
    q = tuple(q_list)

    gq = fq_x_mult_karatsuba(fqx, g, q)
    
    r = fqx.resta(f, gq)

    q = fqx.reducir(q)
    r = fqx.reducir(r)

    return q, r

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.divmod_fast = fq_x_divmod



# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a q-1
# input: a -> tupla de longitud n de elementos de fq
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# utilizar el algoritmo de Cooley-Tuckey
def fq_fft(fq, g, k, a):
    n = 2**k
    if n == 1:
        return a
    
    half_n = n // 2
    p_par = a[::2]
    p_imp = a[1::2]
    
    g2 = fq.mult(g, g)
    
    t_par = fq_fft(fq, g2, k - 1, p_par)
    t_imp = fq_fft(fq, g2, k - 1, p_imp)
    
    vec = [0] * n
    pot = fq.uno()
    for i in range(half_n):
        aux = fq.mult(pot, t_imp[i])
        vec[i] = fq.suma(t_par[i], aux)
        vec[i + half_n] = fq.suma(t_par[i], fq.inv_adit(aux))
        pot = fq.mult(pot, g)
        
    return tuple(vec)

# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fq
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fq_ifft(fq, g, k, a):
    n = 2**k
    pol =fq_fft(fq,fq.inv_mult(g),k,a)
    pol=list(pol)
    val_n = fq.fpx.fp.elem_de_int(n)
    elemento_n =(val_n,)
    for i in range(n):
        pol[i] = fq.mult(pol[i],fq.inv_mult(elemento_n)) 
    return tuple(pol)


