#import parametros as params
from individuo import Individuo
import time
#import numpy as np
from bitarray import bitarray
#import random

"""
def nueva_generacion(_poblacion, _mu, _lambda):
    #####################################
    # - Selección de los padres
    #####################################
    # Se usa algoritmo estocástico universal (SUS)
    seleccion_padres = 0
    # En nuestro caso lamdda = al tamanno de la población,
    # pero dejo la variable para posibles futuros experimentos
    indice = 0
    r = random.uniform(0, 1 / _lambda)
    lista_padres = []

    while seleccion_padres < _lambda:
        while r <= _poblacion[indice].get_prob_padre_acumulada():
            lista_padres.append(_poblacion[indice])
            r = r + 1 / _lambda
            seleccion_padres += 1
        indice += 1

    ##############################################
    # Parte 5 - CRUZE Y MUTACIÓN
    ##############################################
    random.shuffle(lista_padres)  # Barajamos los padres --> cruze aleatorio
    elite = max(_poblacion)  # Reservo el mejor de la población por si debemos aplicar elitismo
    _poblacion = []  # Reseteo de la población - aplicamos relevo generacional

    k = 0
    while k < (_mu - 1):
        hijos = lista_padres[k].cruze(lista_padres[k + 1])  # Generamos los hijos (por parejas)
        _poblacion.append(Individuo(hijos[0]))
        _poblacion.append(Individuo(hijos[1]))
        _poblacion[k].mutacion()        # Se muta cada uno
        _poblacion[k + 1].mutacion()    # de los hijos
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
    for elem in _poblacion:
        elem.set_prob_lin(pos)
        acum += elem.get_prob_padre()
        elem.set_prob_padre_acumulada(acum)
        pos += 1

    return _poblacion

"""

_params = dict

def experimento(_params, _fitness_optimo):
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
    Individuo.s = float(_params['s'])

    _experimento = []
    AES = 0
    SR = 0
    MBF = 0
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
            AES += _primer_hit * numero_de_generaciones
            SR += 1
        _experimento.append(ejecucion)

    end = time.time()

    AES /= numero_de_ejecuciones
    SR = (100*SR) / numero_de_ejecuciones
    _experimento = np.asarray(_experimento)
    for _i in range(numero_de_ejecuciones):
        MBF += _experimento[_i, numero_de_generaciones - 1, 3]
    MBF /= numero_de_ejecuciones
    tiempo = end - start
    # A continuación imprimimos los resultados #todo: esto se puede eliminar
    print("Tiempo requerido: %s" % tiempo)
    print("\nNumero de ""runs"": ", numero_de_ejecuciones)
    print("Usando como criterio de exito un valor de fitness máximo =", _fitness_optimo)
    print("El AES obtenido es:", AES)
    print("El SR (Success Rate) obtenido es: " + str(SR) + "%")
    print("El MBF obtenido es: ", MBF)
    input("Pulsa Enter para continuar...")

    return tiempo, AES, SR, MBF, _experimento
