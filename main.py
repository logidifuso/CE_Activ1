import parametros as params
from individuo import Individuo

import numpy as np
from bitarray import bitarray
import random
import time
import pandas as pd


#####################################################################################
#                               INICIALIZACIÓN
#####################################################################################
# Generación de la población inicial y la ordena de mayor
# a menor fitness (es decir de valores de más bajos al evaluar la función
# a valores más altos)
# -------------------------------------------------------------------------------------
def genera_genes(_long_gen, _num_genes):
    """
    Genera el material genético inicial
    :param _long_gen: longitud de los genes
    :param _num_genes: número de dimensiones, define el número de genes a generar
    :return: devuelve una lista (de genes) donde cada elemento (gen) está
    representado por un bitarray
    """
    individuo = []

    for _i in range(_num_genes):
        dim_j = bitarray()
        for j in range(_long_gen):
            dim_j.append(random.getrandbits(1))
        individuo.append(dim_j)
    return individuo


def asigna_probabilidades(_pobl):
    """
    Asignación de probabilidades (lineales) de selección para cruze (atributo de
    los objetos de la clase Individuo).
    :param _pobl: lista que contiene los individuos de una generación dada
    :return:
    """
    _pobl.sort()  # todo: Es necesario ordenarla? o ya viene siempre ordenada?? en Inicio viene
    pos = 0
    acum = 0
    for el in _pobl:
        el.set_prob_lin(pos)
        acum += el.get_prob_padre()
        el.set_prob_padre_acumulada(acum)
        pos += 1
    return


def genera_poblacion_inicial(_long_gen, _num_genes, _mu):
    """
    Devuelve una lista de instancias de la clase Individuo con el número de
    genes y la longitud de estos especificados.
    :param _long_gen: longitud de los genes
    :param _num_genes: número de genes de cada individuo
    :param _mu:
    :return:
    """
    _poblacion = []
    for _i in range(_mu):
        g = genera_genes(_long_gen, _num_genes)  # Genera genes de un nuevo individuo
        indiv = Individuo(g)  # Creacción de individuo
        indiv.set_fitness()  # Cálculo del fitness
        _poblacion.append(indiv)  # Annade nuevo individuo a la población
    _poblacion.sort()

    asigna_probabilidades(_poblacion)  # Asigna las probabilidades a cada individuo
    return _poblacion


def evalua_poblacion(pobl, _fitness_optimo):
    """
    Retorna las estadísticas principales de la población dada
    La tupla devuelta se compone de:
        1. exito = True o False
        2. media_fitness = media de los fitness de esa generación
        3. mejor = mejor fitness en la generación
        4. peor = peor fitness en la generación
        5. desviación = desviación típica del fitness en la generación
    :param pobl:
    :param _fitness_optimo: mejor fitness posible con la codificación usada
    :return: devuelve (exito, mejor_indiv, media_fitness, mejor, peor, desviacion)
    """
    fitnesses = []
    for el in pobl:
        fitnesses.append(el.get_fitness())
    media_fitness = np.average(fitnesses)
    peor = max(fitnesses)
    mejor = min(fitnesses)
    # varianza = np.var(fitnesses)
    desviacion = np.std(fitnesses)
    mejor_indiv = min(pobl)
    if mejor <= _fitness_optimo:
        exito = True
    else:
        exito = False
    return exito, mejor_indiv, media_fitness, mejor, peor, desviacion


def nueva_generacion(_poblacion, _mu, _lambda):
    """
    Función que realiza las operaciones necesarias para crear una nueva generación
    Partiendo de una población inicial con los valores de fitness, probabilidad de
    selección individual para cruze y probabilidad acumulada para cada individuo reliza
    el resto del proceso:
    1. Selección de padres para reproducción usando el algoritmo estocástico universal (SUS)
    2. Cruze
    3. Mutación
    4. Selección hijos (relevo generacional) y Elitismo
    5. Igualmente completa el resto de información para cada individuo:
        - Valor de fitness
        - Probabilidad selección padre y probabilidad acumulada

    :param _poblacion:
    :param _mu:
    :param _lambda:
    :return:
    """
    # ############################## Selección de los padres  ##########################
    seleccion_padres = 0
    indice = 0
    r = random.uniform(0, 1 / _lambda)
    lista_padres = []

    while seleccion_padres < _lambda:
        while r <= _poblacion[indice].get_prob_padre_acumulada():
            lista_padres.append(_poblacion[indice])
            r = r + 1 / _lambda
            seleccion_padres += 1
        indice += 1

    # ############################# CRUZE Y MUTACIÓN       ############################
    random.shuffle(lista_padres)  # Barajamos los padres --> cruze aleatorio
    elite = max(_poblacion)  # Se reserva el mejor para aplicar elitismo
    _poblacion = []  # Reseteo de la población - relevo generacional

    k = 0
    while k < (_mu - 1):
        hijos = lista_padres[k].cruze(lista_padres[k + 1])  # Generamos los hijos (por parejas)
        _poblacion.append(Individuo(hijos[0]))
        _poblacion.append(Individuo(hijos[1]))
        _poblacion[k].mutacion()  # Se muta cada uno
        _poblacion[k + 1].mutacion()  # de los hijos
        _poblacion[k].set_fitness()
        _poblacion[k + 1].set_fitness()
        k += 2

    mejor_hijo = max(_poblacion)

    if elite > mejor_hijo:
        peor_hijo = min(_poblacion)
        _poblacion.remove(peor_hijo)
        _poblacion.append(elite)

    _poblacion.sort()
    pos = 0
    acum = 0
    for _elem in _poblacion:
        _elem.set_prob_lin(pos)
        acum += _elem.get_prob_padre()
        _elem.set_prob_padre_acumulada(acum)
        pos += 1

    return _poblacion


def experimento(_params):
    # Fijación de parámetros #todo: más tarde habría que cambiar el paso a la clase individuo, para hacerlo más legible
    num_genes = int(_params['dimen'])
    long_gen = int(_params['num_bits'])
    xi_inic = float(_params['Interval_min'])
    xi_fin = float(_params['Interval_max'])
    mu = int(_params['tamano_poblacion'])
    flag_funcion = _params['seleccion_func']
    numero_de_generaciones = int(_params['num_de_generaciones'])
    numero_de_ejecuciones = int(_params['num_de_ejecuciones'])

    # Pasa los parámetros del experimento actual a la clase Individuo
    Individuo.num_genes = int(_params['dimen'])
    Individuo.long_gen = int(_params['num_bits'])
    Individuo.xi_inic = float(_params['Interval_min'])
    Individuo.xi_fin = float(_params['Interval_max'])
    Individuo.mu = int(_params['tamano_poblacion'])
    Individuo.pm = float(_params['prob_mutacion'])
    Individuo.flag_funcion = _params['seleccion_func']
    Individuo.gray = _params['gray']
    Individuo.s = float(_params['s'])

    _fitness_optimo = params.fitness_optimo(flag_funcion, xi_inic, xi_fin, long_gen, num_genes)

    _experimento = []

    _AES = 0
    _SR = 0
    _MBF = 0
    start = time.time()

    for _ejec in range(numero_de_ejecuciones):
        # Generación de la población inicial
        poblacion = genera_poblacion_inicial(long_gen, num_genes, mu)
        ejecucion = list([])  # Inicia la lista de valores para registrar el desempenno del algoritmo
        ejecucion.append(evalua_poblacion(poblacion, _fitness_optimo))  # Anexamos los datos de la población inicial
        _primer_hit = False
        for _gener in range(numero_de_generaciones):

            poblacion = nueva_generacion(poblacion, mu, mu)  # (exito, media, mejor, peor, desviacion)
            datos_generacion = evalua_poblacion(poblacion, _fitness_optimo)
            ejecucion.append(datos_generacion)

            if (datos_generacion[0] is True) and (_primer_hit is False):
                _primer_hit = _gener

        if _primer_hit is not False:
            _AES += _primer_hit * numero_de_generaciones
            _SR += 1
        _experimento.append(ejecucion)

    end = time.time()

    _AES /= numero_de_ejecuciones
    _SR = (100 * _SR) / numero_de_ejecuciones
    _experimento = np.asarray(_experimento)
    for _i in range(numero_de_ejecuciones):
        _MBF += _experimento[_i, numero_de_generaciones - 1, 3]
    _MBF /= numero_de_ejecuciones
    tiempo = end - start

    # A continuación imprimimos los resultados #todo: esto se puede eliminar
    """
    print("Tiempo requerido: %s" % tiempo)
    print("\nNumero de ""runs"": ", numero_de_ejecuciones)
    print("Usando como criterio de exito un valor de fitness máximo =", _fitness_optimo)
    print("El AES obtenido es:", _AES)
    print("El SR (Success Rate) obtenido es: " + str(_SR) + "%")
    print("El MBF obtenido es: ", _MBF)
    input("Pulsa Enter para continuar...")
    """

    return tiempo, _fitness_optimo, _AES, _SR, _MBF, _experimento


l_params = params.lee_parametros("parametros.csv")


for exper in l_params:
    result_experimento = experimento(exper)

    nombre_exp = exper['nombre_experimento']
    tiempo_exp = result_experimento[0]
    optimo = result_experimento[1]
    AES = result_experimento[2]
    SR = result_experimento[3]
    MBF = result_experimento[4]
    datos_evolucion = result_experimento[5]

    dic_result = dict(exper)
    dic_result["tiempo"] = result_experimento[0]
    dic_result["optimo"] = result_experimento[1]
    dic_result["AES"] = result_experimento[2]
    dic_result["SR"] = result_experimento[3]
    dic_result["MBF"] = result_experimento[4]
    # Encapsulo la lista [result_experimento[5]]  para que pandas no intente "desgranarla"
    #temp = datos_evolucion[:, :, 2:]
    #dic_result["datos_evolucion"] = temp.tolist()
    #print("linea 273-borrar, los datos de evolucion son:", dic_result["datos_evolucion"])
    print(dic_result)

    # Creacción de un dataframe de pandas con toda la información del experimento
    # y guardado en un archivo .csv
    df = pd.DataFrame(dic_result, index=[dic_result['nombre_experimento']])
    df.to_csv(str(dic_result['nombre_experimento'])+'.csv')


    print("El AES es:", AES)
    print("El SR es:", SR)
    print("El mejor fitness posible es:", optimo)
    print("El MBF es:", MBF)

    matriz_mejores = datos_evolucion[:, :, 3]
    vector_media_mejores = np.average(matriz_mejores, axis=0)

    import matplotlib.pyplot as plt

    t = np.arange(0, int(exper["num_de_generaciones"])+1)
    s = vector_media_mejores

    fig, ax = plt.subplots()
    ax1 = ax.plot(t, s)

    plt.ylim([0.00307, max(vector_media_mejores)])                      # todo: el fitness optimo en lugar del numero!!

    ax.set(xlabel='Generacion', ylabel='Valor de la función de "fitness"',
           title='Evolución del mejor individuo de la población')
    ax.grid()

    plt.show()



df_leido = pd.read_csv("Exper1.csv")
df_leido = df_leido["MBF"]
df_leido = df_leido[0]          # Importante desencapsular !!!
print(df_leido)
input()

#datos_evolucion = np.asarray(df_leido)
#print(datos_evolucion)
datos_evolucion = datos_evolucion[0]


matriz_mejores = datos_evolucion[1]
vector_media_mejores = np.average(matriz_mejores, axis=0)

# import matplotlib          #todo: quitar
#import matplotlib.pyplot as plt

t = np.arange(0, int(exper["num_de_generaciones"]) + 1)
s = vector_media_mejores

fig, ax = plt.subplots()
ax1 = ax.plot(t, s)

plt.ylim([0.00307, max(vector_media_mejores)])  # todo: el fitness optimo en lugar del numero!!

ax.set(xlabel='Generacion', ylabel='Valor de la función de "fitness"',
       title='Evolución del mejor individuo de la población')
ax.grid()

plt.show()

print("Y hasta aquí llego ahora")

"""
dimen,
num_bits,
Interval_min,
Interval_max,
tamano_poblacion,
prob_mutacion,
seleccion_func,
s,
num_de_generaciones,
num_de_ejecuciones


2,8,-10,10,10,0.01,"Esfera",1.5,50,100
"""
