
import csv
from bitarray import bitarray

from individuo import Individuo


def lee_parametros(file_datos):
    """
    Lee un archivo .csv con los valores de los parámetros a usar por el
    algoritmo y devuelve una lista de diccionarios donde cada diccionario
    contiene los valores de los parámetros para cada uno de los experimentos
    a realizar.

    Los keys de cada uno de los diccionarios devueltos son:
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

    :param file_datos: archivo .csv con los valores de los parámetros para cada experimento
    :return: lista de diccionarios ordenados
    """
    lista_experimentos = []     # TODO: en Any.. es posible limitarlo a un tipo único (e.g.: int ?)
    with open(file_datos) as datos:
        diccionarios = csv.DictReader(datos, delimiter=',')
        for linea in diccionarios:
            lista_experimentos.append(linea)
    return lista_experimentos


def fitness_optimo(funcion, x_min, x_max, num_bits, dimen):
    """
    Devuelve el mejor fitness posible que se puede obtener con la amplitud
    del intervalo (x_max - x_min) dada y el número de bits usados

    :param funcion: Determina la función que se pretende maximizar
    :param x_min: Inicio del intervalo de maximización para cada dimensión
    :param x_max: Final del intervalo de maximización para cada dimensión
    :param num_bits: número de bits usados para codificar el valor real de la
                     variable xi
    :param dimen: Número de dimensiones de la función a evaluar
    :return: devuelve el mejor fitness posible que se puede obtener con la
            amplitud del intervalo (x_max - x_min) dada y el número de bits
            usados
    """

    delta = (x_max - x_min)/(2**num_bits - 1)

    if funcion == "Esfera":
        optimo = bin(round((-x_min)/delta))
    elif funcion == "Schwefel":
        optimo = bin(round((420.9687 - x_min) / delta))
    else:
        return 1

    gen_optimo = bitarray(optimo[2:])
    genes_optimos = []
    for j in range(dimen):
        genes_optimos.append(gen_optimo)

    indiv_optimo = Individuo(genes_optimos)
    indiv_optimo.set_fitness()
    mejor_posible = indiv_optimo.get_fitness()

    return mejor_posible
