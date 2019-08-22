import numpy as np
import matplotlib.pyplot as plt


def graf_mejor_aproximacion(x_ref, y_ref, problema_tipo):
    fig, ax = plt.subplots()

    line1, = ax.plot(x_ref, y_ref, label=problema_tipo)
    # print("\n\nEl fenotipo del mejor individuo que vamos a evaluar es:\n", mejor_individuo.get_fenotipo())
    # input("Pulsa enter para continuar")
    line2, = ax.plot(x_ref, evaluar_fenotipo(mejor_individuo, x_ref),
                     "o", label='Mejor individuo')

    ax.set_xlim([min(x_ref), max(x_ref)])

    ax.set_title("Función objetivo y aproximación")
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.grid()

    ax.legend()
    plt.show()
    return


# #############################################################################
# ############### PLOTS DE PROGRESO O CONVERGENCIA   ##########################
###############################################################################

# Medias del fitness por generación
def graf_medias_fitness_por_generacion(max_generaciones, _media_medias_fitness,
                                       _media_mejor_fitness, _media_peor_fitness):
    """

    :param max_generaciones:
    :param _media_medias_fitness:
    :param _media_mejor_fitness:
    :param _media_peor_fitness:
    :return:
    """
    x = np.arange(max_generaciones)
    fig, ax = plt.subplots()

    ax.plot(x, _media_medias_fitness, label='Media del fitness')
    ax.plot(x, _media_mejor_fitness, label='Media del mejor')
    # line3 = ax.plot(x, mejor_mejores_fitness, label='Mejores individuo')
    # line4 = ax.plot(x, peor_peores_fitness, label='Peores individuos')
    ax.plot(x, _media_peor_fitness, label='Media del peor')

    ax.set_ylim([0.001, 10])
    ax.set_xlim([0, max_generaciones])

    ax.set_title("Medias de los fitness por generacion")
    ax.set_xlabel("Generacion")
    ax.set_ylabel("Valor del error")

    ax.set_yscale('log')

    ax.grid(True)
    ax.legend()
    plt.show(block=True)

    return


# Variación máxima del fitness por generación
def graf_delta_fitess_por_generacion(max_generaciones, _media_mejor_fitness,
                                     _mejor_mejores_fitness, _peor_peores_fitness):
    """

    :param max_generaciones:
    :param _media_mejor_fitness:
    :param _mejor_mejores_fitness:
    :param _peor_peores_fitness:
    :return:
    """
    x = np.arange(max_generaciones)
    fig, ax = plt.subplots()

    # line1 = ax.plot(x, media_medias_fitness, label='Media del fitness')
    ax.plot(x, _media_mejor_fitness, label='Media del mejor')
    ax.plot(x, _mejor_mejores_fitness, label='Mejores individuo')
    ax.plot(x, _peor_peores_fitness, label='Peores individuos')
    # line5 = ax.plot(x, media_peor_fitness, label='Media del peor')

    ax.set_ylim([0.001, 10])
    ax.set_xlim([0, max_generaciones])

    ax.set_title("Variación máxima del fitness por generación")
    ax.set_xlabel("Generacion")
    ax.set_ylabel("Valor del error")

    ax.set_yscale('log')

    ax.grid(True)
    ax.legend()
    plt.show(block=False)
    return


# Media de la desviación típica del fitness por generación
def graf_media_desviacion_por_generacion(max_generaciones, _media_desviacion):
    x = np.arange(max_generaciones)
    fig, ax = plt.subplots()

    ax.plot(x, _media_desviacion, label='Media de la desviacion típica')

    ax.set_ylim([0.0001, 0.1])
    ax.set_xlim([0, max_generaciones])

    ax.set_title("Media de la desviación típica de la función de error")
    ax.set_xlabel("Generacion")
    ax.set_ylabel("Valor del error")

    ax.set_yscale('log')

    ax.grid(True)
    ax.legend()
    plt.show(block=False)
    return


# Mejor fitness por generación de cada ejecución
def graf_mejor_fitness_por_generacion(max_generaciones, _num_ejecuciones,
                                      p_mutacion, _ejecuciones):

    """

    :param max_generaciones:
    :param _num_ejecuciones:
    :param p_mutacion:
    :param _ejecuciones:
    :return:
    """

    x = np.arange(max_generaciones)
    fig, ax = plt.subplots()

    for i in range(_num_ejecuciones):
        ej_mejor_fitness = _ejecuciones[i, :, 3]
        etiqueta = str("Run %i" % i)
        ax.plot(x, ej_mejor_fitness, label=etiqueta)

    ax.set_ylim([0.001, 10])
    ax.set_xlim([0, max_generaciones])

    ax.set_title("Mejor fitness por generacion")
    ax.set_xlabel("Generacion")
    ax.set_ylabel("Valor del error")

    ax.set_yscale('log')

    ax.grid(True)
    ax.legend()

    #texto = '\n'.join(r'p_mutacion=%.3f' % p_mutacion)
    #texto = str(p_mutacion)
    texto = ''.join(r'p_mutacion=%.3f' % p_mutacion)
    ax.text(2, 0.002, texto, fontsize=14, verticalalignment="top") #la posición del texto se da en coordenadas de datos!

    plt.show(block=False)
    nombre_archivo = "".join(["./resultados/mejor_fitness_generacion",
                              "pmut", str(p_mutacion).replace(".", "")
                              ])
    plt.savefig(nombre_archivo)

    return
