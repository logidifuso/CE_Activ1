import parametros as params
import graficos_progreso as gp
from individuo import Individuo

import numpy as np
from bitarray import bitarray
import random
import time
import pandas as pd
import json
import argparse

#####################################################################################
#                        Definiciones Funciones principales
#####################################################################################

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
    # Fijación de parámetros
    num_genes = _params['dimen']
    long_gen = _params['num_bits']
    xi_inic = _params['Interval_min']
    xi_fin = _params['Interval_max']
    mu = _params['tamano_poblacion']
    flag_funcion = _params['seleccion_func']
    numero_de_generaciones = _params['num_de_generaciones']
    numero_de_ejecuciones = _params['num_de_ejecuciones']

    # Pasa los parámetros del experimento actual a la clase Individuo
    Individuo.num_genes = num_genes
    Individuo.long_gen = long_gen
    Individuo.xi_inic = xi_inic
    Individuo.xi_fin = xi_fin
    Individuo.mu = mu
    Individuo.pm = _params['prob_mutacion']
    Individuo.flag_funcion = _params['seleccion_func']
    Individuo.gray = _params['gray']
    Individuo.s = _params['s']

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
        _MBF += _experimento[_i, numero_de_generaciones, 3]
    _MBF /= numero_de_ejecuciones
    tiempo = end - start

    return tiempo, _fitness_optimo, _AES, _SR, _MBF, _experimento


#####################################################################################
#                               Rutina principal
#####################################################################################


# ______________ Parseado del fichero de parámetros y lectura ______________________
# Parseado
parser = argparse.ArgumentParser()
parser.add_argument("-csv_params", "--csv_params",
                    type=str, help="Ruta y fichero de parámetros. Formato .csv")
args = parser.parse_args()
if args.csv_params is not None:
    direccion = args.csv_params
    print(direccion)
else:
    direccion = "parametros.csv"

# Lectura
l_params = params.lee_parametros(direccion)

# _______________________  Ejecucion del experimento  ______________________________
for exper in l_params:

    # Typecasting: El módulo csv devuelve strings, por tanto hay que convertir los tipos manualmente
    exper['dimen'] = int(exper['dimen'])
    exper['num_bits'] = int(exper['num_bits'])
    exper['Interval_min'] = float(exper['Interval_min'])
    exper['Interval_max'] = float(exper['Interval_max'])
    exper['tamano_poblacion'] = int(exper['tamano_poblacion'])
    exper['prob_mutacion'] = float(exper['prob_mutacion'])
    exper['s'] = float(exper['s'])
    exper['num_de_generaciones'] = int(exper['num_de_generaciones'])
    exper['num_de_ejecuciones'] = int(exper['num_de_ejecuciones'])

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
    # Selección de los valores: media_fitness, mejor, peor y desviación para grabar
    dic_result["datos_evolucion"] = result_experimento[5][:, :, 2:].tolist()

    #####################################################################################
    #                       Guardado de los resultados y gráficas
    #####################################################################################

    # Creacción de un dataframe de pandas con toda la información del experimento
    # y guardado en un archivo .csv
    nombre_fichero_resultados = "".join(["./resultados/",
                                         str(dic_result['nombre_experimento'])
                                         ])

    df = pd.DataFrame(dic_result)
    df.to_csv(nombre_fichero_resultados + '.csv')

    # Guardado de los resultados del experimento en un archivo .json
    with open(nombre_fichero_resultados + '.json', 'w') as outfile:
        json.dump(dic_result, outfile)

    # #############################   Gráficas ####################################

    # matriz_mejores = datos_evolucion[:, :, 3]
    # v_media_mejores = np.average(matriz_mejores, axis=0)
    v_media_medias = np.average(datos_evolucion[:, :, 2], axis=0)
    v_media_mejores = np.average(datos_evolucion[:, :, 3], axis=0)
    v_media_peores = np.average(datos_evolucion[:, :, 4], axis=0)

    gp.graf_medias_fitness_por_generacion(exper['num_de_generaciones']+1,
                                          v_media_medias,
                                          v_media_mejores,
                                          v_media_peores)

    gp.graf_mejor_fitness_por_generacion(exper['num_de_generaciones']+1,
                                         exper['num_de_ejecuciones'],
                                         exper['prob_mutacion'],
                                         datos_evolucion)

    print("El AES es:", AES)
    print("El SR es:", SR)
    print("El mejor fitness posible es:", optimo)
    print("El MBF es:", MBF)

    """
    
    
    # Gráfica con la media del mejor fitness (MBF por generación)
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
    """


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
