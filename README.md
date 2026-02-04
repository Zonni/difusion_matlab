# 1D Neutron Diffusion Solver (FEM & FDM)

Este repositorio contiene una implementación en **MATLAB** orientada a objetos para resolver la ecuación de difusión de neutrones en 1D. El código permite calcular el factor de multiplicación efectivo ($k_{eff}$) y los modos del flujo neutrónico ($\phi$) para reactores homogéneos y heterogéneos.

## Características

* **Modelos de Energía:**
    * **1 Grupo (1g):** Difusión simple monoenergética.
    * **2 Grupos (2g):** Modelo acoplado Rápido-Térmico.
* **Métodos Numéricos:**
    * **FEM:** Método de Elementos Finitos con polinomios de Lagrange de orden variable.
    * **FDM:** Método de Diferencias Finitas (centradas).
* **Geometría:** Soporte para mallas no uniformes y configuraciones heterogéneas.

## Estructura del Código

El proyecto se organiza en clases modulares:

* **`Malla1D`**: Define la discretización espacial y tamaño de celdas.
* **`Materiales1D1g` / `Materiales1D2g`**: Gestiona las secciones eficaces ($D, \Sigma_a, \nu\Sigma_f$) y matrices de scattering ($\Sigma_{s12}$).
* **`ElementoFinito`**: Genera funciones de forma simbólicas, derivadas y cuadratura de Gauss.
* **`ProblemaDifusion1D1g` / `ProblemaDifusion1D2g`**:
    * Ensamblaje de matrices globales $A$ (Rigidez/Pérdidas) y $B$ (Producción).
    * Aplicación de condiciones de contorno (Dirichlet).
    * Solución del problema de autovalores generalizado.
    * Visualización de flujos.

## Uso

Los scripts principales configuran la geometría, asignan materiales y ejecutan el solver. Ejemplos incluidos:

* `FEM_1G_HOM_N14.m`: Caso 1 Grupo Homogéneo.
* `FEM_1G_HET_N14.m`: Caso 1 Grupo Heterogéneo.
* `FEM_2G_HOM_N14.m`: Caso 2 Grupos Homogéneo.
* `FEM_2G_HET_N14.m`: Caso 2 Grupos Heterogéneo.

## Requisitos

* MATLAB (con Symbolic Math Toolbox para la generación de elementos finitos).

## Autores

Proyecto realizado por 

- Jaime Gómez (Universitat de València)

Bajo la supervisión de:

- Antoni Vidal-Ferràndiz (Universitat Politècnica de València)
- Amanda Carreño  (Universitat Politècnica de València)
