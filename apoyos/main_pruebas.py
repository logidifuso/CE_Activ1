#mport parametros_pruebas
#import individuo

import lectura_diccionario
import individuo_prueba


lista_params = lectura_diccionario.lee_parametros("parametros.csv")


diccionario = lista_params[0]
print(lista_params[0])
print(diccionario['prob_mutacion'])
prob_mutacion = diccionario['prob_mutacion']




print("El valor de la variable en la clase al principio es:", individuo_pruebas.IndivPruebas.variable)

individuo_pruebas.IndivPruebas.variable = "Pepito"

print("El valor de la variable en la clase tras cambiar es:", individuo_pruebas.IndivPruebas.variable)

"""
for elem is lista_params:
    dicc_params = el
"""

"""

params = {
    'TAMANO_POBLACION': 2000,
    'LONG_MAX_GENOTIPO': 240,
    'MAX_WRAPS': 2,
    'MAX_VAL_CODON': 256,

    'MIN_LONG_FENOTIPO_INICIAL': 1,
    'MAX_LONG_FENOTIPO_INICIAL': 8,

    'ARCHIVO_GRAMATICA': 'gramatica_nucleos.bnf',
    'PROBLEMA_TIPO': 'Problema1',

    'OPCION_SELECCION': 'Torneo',
    'TAMANO_TORNEO': 2,

    'NUM_EJECUCIONES': 8,
    'MAX_GENERACIONES': 200,

    'U': 0.1,
    'K0': 1,
    'K1': 10,
    'S': 1.8,

    'p_mutacion': 0.1,  # todo: decidir si es una constante o se usa en algo memético
    'p_cruze': 0.901

}
"""