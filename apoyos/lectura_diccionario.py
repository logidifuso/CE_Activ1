import csv
from collections import OrderedDict
from typing import List, Any


def lee_parametros(file_datos):
    """
    Lee un archivo .csv con los valores de los parámetros a usar por el
    algoritmo y devuelve una lista de diccionarios donde cada diccionario
    contiene los valores de los parámetros para cada uno de los experimentos
    a realizar
    :param file_datos: archivo .csv con los valores de los parámetros para cada experimento
    :return: lista de diccionarios ordenados
    """
    lista_experimentos= [] #TODO: en Any.. es posible limitarlo a un tipo único (e.g.: int ?)
    with open(file_datos) as datos:
        diccionarios = csv.DictReader(datos, delimiter=',')
        for linea in diccionarios:
            lista_experimentos.append(linea)
    return lista_experimentos

