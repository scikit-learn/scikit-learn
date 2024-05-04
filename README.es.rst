.. -*- mode: rst -*-

|Azure| |CirrusCI| |Codecov| |CircleCI| |Nightly wheels| |Black| |PythonVersion| |PyPi| |DOI| |Benchmark|

.. |Azure| image:: https://dev.azure.com/scikit-learn/scikit-learn/_apis/build/status/scikit-learn.scikit-learn?branchName=main
   :target: https://dev.azure.com/scikit-learn/scikit-learn/_build/latest?definitionId=1&branchName=main

.. |CircleCI| image:: https://circleci.com/gh/scikit-learn/scikit-learn/tree/main.svg?style=shield
   :target: https://circleci.com/gh/scikit-learn/scikit-learn

.. |CirrusCI| image:: https://img.shields.io/cirrus/github/scikit-learn/scikit-learn/main?label=Cirrus%20CI
   :target: https://cirrus-ci.com/github/scikit-learn/scikit-learn/main

.. |Codecov| image:: https://codecov.io/gh/scikit-learn/scikit-learn/branch/main/graph/badge.svg?token=Pk8G9gg3y9
   :target: https://codecov.io/gh/scikit-learn/scikit-learn

.. |Nightly wheels| image:: https://github.com/scikit-learn/scikit-learn/workflows/Wheel%20builder/badge.svg?event=schedule
   :target: https://github.com/scikit-learn/scikit-learn/actions?query=workflow%3A%22Wheel+builder%22+event%3Aschedule

.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/scikit-learn.svg
   :target: https://pypi.org/project/scikit-learn/

.. |PyPi| image:: https://img.shields.io/pypi/v/scikit-learn
   :target: https://pypi.org/project/scikit-learn

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black

.. |DOI| image:: https://zenodo.org/badge/21369/scikit-learn/scikit-learn.svg
   :target: https://zenodo.org/badge/latestdoi/21369/scikit-learn/scikit-learn

.. |Benchmark| image:: https://img.shields.io/badge/Benchmarked%20by-asv-blue
   :target: https://scikit-learn.org/scikit-learn-benchmarks

.. |PythonMinVersion| replace:: 3.9
.. |NumPyMinVersion| replace:: 1.19.5
.. |SciPyMinVersion| replace:: 1.6.0
.. |JoblibMinVersion| replace:: 1.2.0
.. |ThreadpoolctlMinVersion| replace:: 3.1.0
.. |MatplotlibMinVersion| replace:: 3.3.4
.. |Scikit-ImageMinVersion| replace:: 0.17.2
.. |PandasMinVersion| replace:: 1.1.5
.. |SeabornMinVersion| replace:: 0.9.0
.. |PytestMinVersion| replace:: 7.1.2
.. |PlotlyMinVersion| replace:: 5.14.0

.. image:: https://raw.githubusercontent.com/scikit-learn/scikit-learn/main/doc/logos/scikit-learn-logo.png
  :target: https://scikit-learn.org/

**scikit-learn** es un módulo de Python para el aprendizaje automático
construido sobre SciPy y se distribuye bajo la licencia BSD de 3 cláusulas.

David Cournapeau comenzó el proyecto en 2007 como un proyecto de Google
Summer of Code, y desde entonces muchos voluntarios han contribuido.
Consulte la página
`About us <https://scikit-learn.org/dev/about.html#authors>`__
para obtener una lista de los colaboradores principales.

Actualmente es mantenido por un equipo de voluntarios.

Sitio web: https://scikit-learn.org

Instalación
-----------

Dependencias
~~~~~~~~~~~~

scikit-learn requiere:

- Python (>= |PythonMinVersion|)
- NumPy (>= |NumPyMinVersion|)
- SciPy (>= |SciPyMinVersion|)
- joblib (>= |JoblibMinVersion|)
- threadpoolctl (>= |ThreadpoolctlMinVersion|)

=======

**Scikit-learn 0.20 fue la última versión que soportó Python 2.7 y Python 3.4.**
Scikit-learn 1.0 y posteriores requieren Python 3.7 o más reciente.
Scikit-learn 1.1 y posteriores requieren Python 3.8 o más reciente.

Las capacidades de trazar de scikit-learn (es decir, funciones que comienzan
con ``plot_`` y clases que terminan con ``Display``) requieren Matplotlib
(>= |MatplotlibMinVersion|).
Para ejecutar los ejemplos, se requiere Matplotlib >= |MatplotlibMinVersion|.
Algunos ejemplos requieren scikit-image >= |Scikit-ImageMinVersion|,
algunos ejemplos requieren pandas >= |PandasMinVersion|,
lgunos ejemplos necesitan seaborn >= |SeabornMinVersion| y plotly >= |PlotlyMinVersion|.

Instalación por el Usuario
~~~~~~~~~~~~~~~~~~~~~~~~~~

Si ya tiene una instalación funcional de NumPy y SciPy, la forma más fácil de instalar scikit-learn es usando ``pip``::

    pip install -U scikit-learn

o ``conda``::

    conda install -c conda-forge scikit-learn

La documentación incluye instrucciones de instalación más detalladas:  `installation instructions <https://scikit-learn.org/stable/install.html>`_.


Registro de Cambios
-------------------

Consulte el `changelog <https://scikit-learn.org/dev/whats_new.html>`__ para ver
un historial de cambios notables en scikit-learn.

Desarrollo
----------

Damos la bienvenida a nuevos colaboradores de todo nivel de experiencia.
Los objetivos de la comunidad de scikit-learn son ser útiles, acogedores y
eficaces. La guía `Development Guide <https://scikit-learn.org/stable/developers/index.html>`_
tiene información detallada sobre cómo contribuir con código, documentación,
pruebas y más. Hemos incluido información básica en este README.

Enlaces Importantes
~~~~~~~~~~~~~~~~~~~

- Repositorio oficial del código fuente: https://github.com/scikit-learn/scikit-learn
- Descargar versiones: https://pypi.org/project/scikit-learn/
- Sistema de seguimiento de tareas: https://github.com/scikit-learn/scikit-learn/issues

Código Fuente
~~~~~~~~~~~~~

Puede verificar las últimas fuentes con el comando::

    git clone https://github.com/scikit-learn/scikit-learn.git

Contribución
~~~~~~~~~~~~

Para aprender más sobre cómo hacer una contribución a scikit-learn,
consulte nuestra guía de contribución
`Contributing guide
<https://scikit-learn.org/dev/developers/contributing.html>`_.

Pruebas
~~~~~~~

Después de la instalación, puede lanzar la suite de pruebas desde fuera del
directorio fuente (necesitará tener ``pytest`` >= |PyTestMinVersion| instalado)::

    pytest sklearn

Consulte la página web https://scikit-learn.org/dev/developers/contributing.html#testing-and-improving-test-coverage
para obtener más información.

    La generación de números aleatorios puede controlarse durante las pruebas
    al configurarse la variable de entorno ``SKLEARN_SEED``.

Envío de Pull Request
~~~~~~~~~~~~~~~~~~~~~

Antes de abrir un Pull Request, revise la página completa de
contribución (Contributing) para asegurarse de que su código cumple con nuestras
reglas: https://scikit-learn.org/stable/developers/index.html

Historia del Proyecto
---------------------

El proyecto fue iniciado en 2007 por David Cournapeau como un proyecto de Google
Summer of Code, y desde entonces muchos voluntarios han contribuido. Vea la página
`About us <https://scikit-learn.org/dev/about.html#authors>`__  para obtener una
lista de los colaboradores principales.

El proyecto actualmente es mantenido por un equipo de voluntarios.

**Nota**: `scikit-learn` anteriormente se conocía como `scikits.learn`.

Ayuda y Soporte
----------------

Documentación
~~~~~~~~~~~~~

- Documentación HTML (versión estable): https://scikit-learn.org
- Documentación HTML (versión de desarrollo): https://scikit-learn.org/dev/
- Preguntas frecuentes: https://scikit-learn.org/stable/faq.html

Comunicación
~~~~~~~~~~~~

Lista de correo: https://mail.python.org/mailman/listinfo/scikit-learn
Logos & Branding: https://github.com/scikit-learn/scikit-learn/tree/main/doc/logos
Blog: https://blog.scikit-learn.org
Calendario: https://blog.scikit-learn.org/calendar/
Twitter: https://twitter.com/scikit_learn
Stack Overflow: https://stackoverflow.com/questions/tagged/scikit-learn
Discusiones en GitHub: https://github.com/scikit-learn/scikit-learn/discussions
Sitio web: https://scikit-learn.org
LinkedIn: https://www.linkedin.com/company/scikit-learn
YouTube: https://www.youtube.com/channel/UCJosFjYm0ZYVUARxuOZqnnw/playlists
Facebook: https://www.facebook.com/scikitlearnofficial/
Instagram: https://www.instagram.com/scikitlearnofficial/
TikTok: https://www.tiktok.com/@scikit.learn
Mastodon: https://mastodon.social/@sklearn@fosstodon.org
Discord: https://discord.gg/h9qyrK8Jc8

Citación
~~~~~~~~

Si utiliza scikit-learn en una publicación científica, agradeceríamos las citas: https://scikit-learn.org/stable/about.html#citing-scikit-learn
