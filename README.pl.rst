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

**scikit-learn** to moduł Pythona przeznaczony do uczenia maszynowego,
zbudowany na bazie SciPy, rozpowszechniany na licencji 3-Clause BSD.

Projekt został rozpoczęty w 2007 roku przez Davida Cournapeau jako projekt
w ramach Google Summer of Code, a od tego czasu swój wkład wniosło wielu
wolontariuszy. Listę głównych współtwórców można znaleźć na stronie
`About us <https://scikit-learn.org/dev/about.html#authors>`__ .

Obecnie projektem zarządza zespół wolontariuszy.

Strona internetowa: https://scikit-learn.org


Instalacja
----------

Zależności
~~~~~~~~~~

scikit-learn wymaga:

- Python (>= |PythonMinVersion|)
- NumPy (>= |NumPyMinVersion|)
- SciPy (>= |SciPyMinVersion|)
- joblib (>= |JoblibMinVersion|)
- threadpoolctl (>= |ThreadpoolctlMinVersion|)

=======

**Scikit-learn 0.20 była ostatnią wersją obsługującą Pythona 2.7 i Pythona 3.4.**
scikit-learn 1.0 i kolejne wymagają Pythona 3.7 lub nowszego.
scikit-learn 1.1 i kolejne wymagają Pythona 3.8 lub nowszego.

Możliwości tworzenia wykresów scikit-learn (tj. funkcje zaczynające się od ``plot_``
i klasy kończące się na ``Display``) wymagają Matplotlib (>= |MatplotlibMinVersion|).
Do uruchomienia przykładów wymagany jest Matplotlib >= |MatplotlibMinVersion|.
Niektóre przykłady wymagają scikit-image >= |Scikit-ImageMinVersion|,
niektóre przykłady wymagają pandas >= |PandasMinVersion|,
niektóre przykłady wymagają seaborn >= |SeabornMinVersion|
oraz plotly >= |PlotlyMinVersion|.

Instalacja dla Użytkowników
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Jeśli masz już działającą instalację NumPy i SciPy,
najłatwiejszym sposobem na instalację scikit-learn jest użycie ``pip``::

  pip install -U scikit-learn

lub ``conda``::

    conda install -c conda-forge scikit-learn

Dokumentacja zawiera bardziej szczegółowe instrukcje instalacji `installation instructions <https://scikit-learn.org/stable/install.html>`_.

Dziennik Zmian
--------------

Zapoznaj się z dziennikiem zmian `changelog <https://scikit-learn.org/dev/whats_new.html>`__,
aby zobaczyć historię znaczących zmian w scikit-learn.

Rozwój
------

Zapraszamy do współpracy nowych współtwórców na każdym poziomie doświadczenia.
Cele społeczności scikit-learn to bycie pomocnym, przyjaznym i skutecznym.
Przewodnik dla programistów `Development Guide <https://scikit-learn.org/stable/developers/index.html>`_
zawiera szczegółowe informacje o współtworzeniu kodu, dokumentacji, testów
i więcej. Podstawowe informacje zostały zawarte w tym README.

Ważne Linki
~~~~~~~~~~~

- Oficjalne repozytorium kodu źródłowego: https://github.com/scikit-learn/scikit-learn
- Pobieranie wydań: https://pypi.org/project/scikit-learn/
- Śledzenie zgłoszeń: https://github.com/scikit-learn/scikit-learn/issues

Kod źródłowy
~~~~~~~~~~~~

Możesz sprawdzić najnowsze źródła za pomocą polecenia::

    git clone https://github.com/scikit-learn/scikit-learn.git

Współtworzenie
~~~~~~~~~~~~~~

Aby dowiedzieć się więcej o współtworzeniu scikit-learn, zapoznaj się z naszym
przewodnikiem dla współtwórców `Contributing guide
<https://scikit-learn.org/dev/developers/contributing.html>`_.

Testowanie
~~~~~~~~~~

Po instalacji możesz uruchomić zestaw testów spoza katalogu źródłowego
(musisz mieć zainstalowane ``pytest`` >= |PyTestMinVersion|)::

    pytest sklearn

Sprawdź stronę internetową https://scikit-learn.org/dev/developers/contributing.html#testing-and-improving-test-coverage
aby uzyskać więcej informacji.

    Generację liczb losowych podczas testowania można kontrolować ustawiając
    zmienną środowiskową ``SKLEARN_SEED``.

Składanie Wniosku o Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Przed otwarciem Pull Requesta, zapoznaj się ze stroną o wkładzie, aby upewnić
się, że Twój kod spełnia nasze wytyczne: https://scikit-learn.org/stable/developers/index.html

Historia Projektu
-----------------

Projekt został rozpoczęty w 2007 roku przez Davida Cournapeau jako projekt
w ramach Google Summer of Code, a od tego czasu swój wkład wniosło wielu
wolontariuszy. Listę głównych współtwórców można znaleźć na stronie
`About us <https://scikit-learn.org/dev/about.html#authors>`__ .

Projekt obecnie jest utrzymywany przez zespół wolontariuszy.

**Uwaga**: scikit-learn był wcześniej znany jako `scikits.learn`.

Pomoc i Wsparcie
----------------

Dokumentacja
~~~~~~~~~~~~

- Dokumentacja HTML (stabilna wersja): https://scikit-learn.org
- Dokumentacja HTML (wersja deweloperska): https://scikit-learn.org/dev/
- FAQ: https://scikit-learn.org/stable/faq.html

Komunikacja
~~~~~~~~~~~

- Lista Mailingowa: https://mail.python.org/mailman/listinfo/scikit-learn
- Logotypy i Branding: https://github.com/scikit-learn/scikit-learn/tree/main/doc/logos
- Blog: https://blog.scikit-learn.org
- Kalendarz: https://blog.scikit-learn.org/calendar/
- Twitter: https://twitter.com/scikit_learn
- Stack Overflow: https://stackoverflow.com/questions/tagged/scikit-learn
- Dyskusje na GitHubie: https://github.com/scikit-learn/scikit-learn/discussions
- Strona Internetowa: https://scikit-learn.org
- LinkedIn: https://www.linkedin.com/company/scikit-learn
- YouTube: https://www.youtube.com/channel/UCJosFjYm0ZYVUARxuOZqnnw/playlists
- Facebook: https://www.facebook.com/scikitlearnofficial/
- Instagram: https://www.instagram.com/scikitlearnofficial/
- TikTok: https://www.tiktok.com/@scikit.learn
- Mastodon: https://mastodon.social/@sklearn@fosstodon.org
- Discord: https://discord.gg/h9qyrK8Jc8

Cytowanie
~~~~~~~~~

W przypadku wykorzystania scikit-learn w pracach naukowych, prosimy o cytowanie: https://scikit-learn.org/stable/about.html#citing-scikit-learn
