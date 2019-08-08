
import numpy as np
import random
import math
from bitarray import bitarray


class Individuo(object):
    """
    Implementación de la estructura de datos correspondiente a cada individuo y los métodos
    utilizados.

    Atributos:
        * genes: genotipo del individuo
        * fitness: El valor de fitness representa la variación con respecto a la mejor
        * aproximación posible. A menor valor de fitness, mejor individuo
        * mating_prob: Probabilidad de cruze
        * prob_acumulada: Probabilidad de cruze acumulada (para algoritmo estocástico uni-
          versal (SUS)

    Métodos:
        * cruze: Implementa el cruze entre 2 individuos
    """
    tag = 0
    num_genes = None
    long_gen = None
    xi_inic = None
    xi_fin = None
    pm = None
    flag_funcion = None
    s = None
    mu = None

    def __init__(self, genes):
        self.genes = genes
        self.fitness = None
        self.prob_lin = None
        self.prob_acumulada = None
        self.id = Individuo.tag
        Individuo.tag += 1

    def __lt__(self, other):
        """
        "Less than" method, para ordenar los individuos en las poblaciones según su valor de
        fitness. A mayor valor de fitness --> "menor" es el individuo, pues tal como está
        definido queremos minimizar el valor de fitness
        :param other:
        :return:
        """
        return self.fitness > other.fitness

    def __str__(self):
        return str(self.genes) + " Fitness:" + str(self.fitness) + "  Prob indivdual: " + str(
            self.prob_lin) + " Prob Acum: " + str(self.prob_acumulada)

    #######################################################################################
    # Método que implementa el cruze entre 2 individuos
    #######################################################################################
    def cruze(self, otro):
        """
        Implementación de cruze de 2 puntos entre los individuos "self" y "otro"
        :param otro:
        :return: retorna los genes de 2 hijos (hijo1, hijo2)
        """
        punto1 = 0
        punto2 = 0

        hijo1 = []
        hijo2 = []

        while punto1 >= punto2:
            punto1 = random.randint(0, Individuo.long_gen)
            punto2 = random.randint(0, Individuo.long_gen)

        for i in range(Individuo.num_genes):
            gen_hijo1 = self.get_genes()[i][:]
            gen_hijo1[punto1:punto2] = otro.get_genes()[i][punto1:punto2]

            gen_hijo2 = otro.get_genes()[i][:]
            gen_hijo2[punto1:punto2] = self.get_genes()[i][punto1:punto2]

            hijo1.append(gen_hijo1)
            hijo2.append(gen_hijo2)
        return hijo1, hijo2

    ############################################################################################
    # Método que implementa las mutaciones
    ############################################################################################
    def mutacion(self):
        """
        Implementación de la mutación. El procedimiento consiste en:
        1. Crear un array de la longitud correspondiente con probabilidad de 1's = pm
        2. Convertir a bitarray
        3. Operación XOR entre los genes y el bitarray de mutación
        :return:
        """
        for i in range(Individuo.num_genes):
            array_mutacion = np.random.choice([0, 1], size=Individuo.long_gen, p=[1 - Individuo.pm, Individuo.pm])
            array_mutacion = bitarray(array_mutacion.tolist())
            self.get_genes()[i] = self.get_genes()[i] ^ array_mutacion

    ###########################################################################################
    # Evaluación y escritura del fitness correspondiente a la función escogida (esfera/Schwefel) 
    ###########################################################################################
    def set_fitness(self):
        """
        Método que evalúa el fitness del individuo y lo asigna al correspondiente
        atributo "fitness"
        Hace uso de las variables de clase:
        - flag.funcion: Función esfera o Schwefel
        - xi_inin, xi_fin : Inicio y fin del intervalo de evaluación de la función
        - long_gen: longitud de los genes del individuo
        - num_genes: número de genes del individuo. Igual al número de dimensiones
        que se están evaluando en la función.
        :return:
        """
        def bitarray_a_entero(cadena_bits):
            n = cadena_bits.to01()  # Convierte el bitarray a una cadena 0 y 1's
            return int(n, 2)  # Coge la cadena como entero en base 2 y la convierte en entero

        def f_esfera(genotipo, a1, a2, long_bits):
            tot = 0
            factor = (a2 - a1) / (2 ** long_bits - 1)
            for i in genotipo:
                tot += (a1 + factor * bitarray_a_entero(i)) ** 2
            return tot

        def f_schwefel(genotipo, a1, a2, long_bits, dim_n):
            tot = 0
            factor = (a2 - a1) / (2 ** long_bits - 1)
            for i in genotipo:
                aux = (a1 + factor * bitarray_a_entero(i))
                tot += aux * (math.sin(math.sqrt(abs(aux))))
            return 418.9829 * dim_n - tot
            # Introduzco un control de excepciones?? --> Opcional..

        if Individuo.flag_funcion == "Esfera":
            self.fitness = f_esfera(self.genes, Individuo.xi_inic, Individuo.xi_fin, Individuo.long_gen)
        elif Individuo.flag_funcion == "Schwefel":
            self.fitness = f_schwefel(self.genes, Individuo.xi_inic, Individuo.xi_fin, Individuo.long_gen,
                                      Individuo.num_genes)
        return

    ##########################################################################################
    # Funciones "setter" y "getter" para encapsulado
    ##########################################################################################
    def set_prob_lin(self, rank):
        """
        Fija la probabilidad lineal de selección para cruze en función del rango (posición
        en la lista de población previamente ordenada
        :param rank:
        :return:
        """
        self.prob_lin = (2 - Individuo.s) / Individuo.mu + (2 * rank * (Individuo.s - 1)) / \
                        (Individuo.mu * (Individuo.mu - 1))

    def set_prob_padre_acumulada(self, prob):
        self.prob_acumulada = prob
        return

    def get_genes(self):
        return self.genes

    def get_fitness(self):
        return self.fitness

    def get_prob_padre(self):
        return self.prob_lin

    def get_prob_padre_acumulada(self):
        return self.prob_acumulada

    def get_id(self):
        return self.id
