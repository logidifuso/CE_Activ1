import parametros as params
from individuo import Individuo

import numpy as np
from bitarray import bitarray
import random


l_params = params.lee_parametros("parametros.csv")

for exper in l_params:
    dicc_params = exper

    # Fijación de parámetros #todo: más tarde habría que cambiar el paso a la clase individuo, para hacerlo más legible
    num_genes = int(exper['dimen'])
    long_gen = int(exper['num_bits'])
    xi_inic = float(exper['Interval_min'])
    xi_fin = float(exper['Interval_max'])
    mu = int(exper['tamano_poblacion'])
    flag_funcion = exper['seleccion_func']
    numero_de_generaciones = int(exper['num_de_generaciones'])
    numero_de_ejecuciones = int(exper['num_de_ejecuciones'])

    # Pasa los parámetros del experimento actual a la clase Individuo
    Individuo.num_genes = int(exper['dimen'])
    Individuo.long_gen = int(exper['num_bits'])
    Individuo.xi_inic = float(exper['Interval_min'])
    Individuo.xi_fin = float(exper['Interval_max'])
    Individuo.mu = int(exper['tamano_poblacion'])
    Individuo.pm = float(exper['prob_mutacion'])
    Individuo.flag_funcion = exper['seleccion_func']
    Individuo.s = float(exper['s'])

    _experimento = []
    fitness_optimo = params.fitness_optimo(flag_funcion, xi_inic, xi_fin, long_gen, num_genes)