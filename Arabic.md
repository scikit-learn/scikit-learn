.. -*- mode: rst -*-

|Azure|_ |Travis|_ |Codecov|_ |CircleCI|_ |Nightly wheels|_ |Black|_ |PythonVersion|_ |PyPi|_ |DOI|_

.. |Azure| image:: https://dev.azure.com/scikit-learn/scikit-learn/_apis/build/status/scikit-learn.scikit-learn?branchName=main
.. _Azure: https://dev.azure.com/scikit-learn/scikit-learn/_build/latest?definitionId=1&branchName=main

.. |Travis| image:: https://api.travis-ci.com/scikit-learn/scikit-learn.svg?branch=main
.. _Travis: https://travis-ci.com/scikit-learn/scikit-learn

.. |Codecov| image:: https://codecov.io/gh/scikit-learn/scikit-learn/branch/main/graph/badge.svg?token=Pk8G9gg3y9
.. _Codecov: https://codecov.io/gh/scikit-learn/scikit-learn

.. |CircleCI| image:: https://circleci.com/gh/scikit-learn/scikit-learn/tree/main.svg?style=shield&circle-token=:circle-token
.. _CircleCI: https://circleci.com/gh/scikit-learn/scikit-learn

.. |Nightly wheels| image:: https://github.com/scikit-learn/scikit-learn/workflows/Wheel%20builder/badge.svg?event=schedule
.. _`Nightly wheels`: https://github.com/scikit-learn/scikit-learn/actions?query=workflow%3A%22Wheel+builder%22+event%3Aschedule

.. |PythonVersion| image:: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue
.. _PythonVersion: https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue

.. |PyPi| image:: https://badge.fury.io/py/scikit-learn.svg
.. _PyPi: https://badge.fury.io/py/scikit-learn

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
.. _Black: https://github.com/psf/black

.. |DOI| image:: https://zenodo.org/badge/21369/scikit-learn/scikit-learn.svg
.. _DOI: https://zenodo.org/badge/latestdoi/21369/scikit-learn/scikit-learn


.. |PythonMinVersion| replace:: 3.7
.. |NumPyMinVersion| replace:: 1.14.6
.. |SciPyMinVersion| replace:: 1.1.0
.. |JoblibMinVersion| replace:: 0.11
.. |ThreadpoolctlMinVersion| replace:: 2.0.0
.. |MatplotlibMinVersion| replace:: 2.2.2
.. |Scikit-ImageMinVersion| replace:: 0.14.5
.. |PandasMinVersion| replace:: 0.25.0
.. |SeabornMinVersion| replace:: 0.9.0
.. |PytestMinVersion| replace:: 5.0.1

.. image:: doc/logos/scikit-learn-logo.png
  :target: https://scikit-learn.org/

**scikit-learn** وحدة بايثون للتعلم الآلي المبنية على سايباي و يتم توزيعها بموجب ترخيص بي اس دي المكون من 3 فقرات  .

بدا هذا المشروع في عام 2007 بواسطة ديفيد كورنابيو في صيف جوجل من مشروع كود, و منذ ذلك الحين ساهم العديد من المتطوعين
انظر الى صفحة  
 `About us <https://scikit-learn.org/dev/about.html#authors>`__ 
 للحصول على قائمة المساهمين الاساسيين 

الموقع تحت الصيانة من قبل المتطوعين.

Website: https://scikit-learn.org

التثبييت
------------

المعتمدات
~~~~~~~~~~~~

scikit-learn requires:

- Python (>= |PythonMinVersion|)
- NumPy (>= |NumPyMinVersion|)
- SciPy (>= |SciPyMinVersion|)
- joblib (>= |JoblibMinVersion|)
- threadpoolctl (>= |ThreadpoolctlMinVersion|)

=======

**Scikit-learn 0.20 كان الاصدار الاخير الذي يدعم بايثون 2.7 و بايثون 3.4**
scikit-learn 0.23 هذا الاصدار و ما بعده يتطلب بايثون 3.6 و ما بعده.
scikit-learn 1.0 هذا الاصدار و ما بعده يتطلب بايثون 3.7 و ما بعده.

Scikit-learn plotting capabilities (i.e., functions start with ``plot_`` and
classes end with "Display") require Matplotlib (>= |MatplotlibMinVersion|).
For running the examples Matplotlib >= |MatplotlibMinVersion| is required.
A few examples require scikit-image >= |Scikit-ImageMinVersion|, a few examples
require pandas >= |PandasMinVersion|, some examples require seaborn >=
|SeabornMinVersion|.

تعليمات التثبييت للمستخدم
~~~~~~~~~~~~~~~~~

اذا كان لديك تثبييت فعال لل نمباي و سايباي فإن اسهل طريقة لتثبييت سايكيت ليرن هي استخدام  ``pip``   ::

    pip install -U scikit-learn

او ``conda``::

    conda install -c conda-forge scikit-learn

هذا الملف يحتوي على تفاصيل اكثر `installation instructions <https://scikit-learn.org/stable/install.html>`_.


سجل التغييرات 
---------

انظر الى `changelog <https://scikit-learn.org/dev/whats_new.html>`__
للحصول على تاريخ من التغييرات الملحوظة في سايكيت ليرن

التطوير
-----------

نرحب بالمساهمين الجدد من جميع مستويات الخبرة.  أهداف مجتمع سايكيت ليرن مفيدة و مرحبة و فعالة

`Development Guide <https://scikit-learn.org/stable/developers/index.html>`_

هذا الرابط يحتوي معلومات مفصلة حول اكواد المساهمة و التوثيقات و الاختبارات و المزيد. هناك معلومات اساسية مضافة في الملف ريدمي

روابط مهمة
~~~~~~~~~~~~~~~

- موقع المشروع الرسمي: https://github.com/scikit-learn/scikit-learn
- لتحميل الاصدارات: https://pypi.org/project/scikit-learn/
- لتتبع المشكلات: https://github.com/scikit-learn/scikit-learn/issues

الكود الاصلي
~~~~~~~~~~~

يمكنك الاطلاع على الكود الاصلي من هنا::

    git clone https://github.com/scikit-learn/scikit-learn.git

المساهمة
~~~~~~~~~~~~

لمعرفة المزيد حول تقديم مساهمة في سايكيت ليرن يرجى الاطلاع على
`Contributing guide
<https://scikit-learn.org/dev/developers/contributing.html>`_.

الاختبارات
~~~~~~~

بعد التثبييت, يمكنك تشغييل مجموعة الاختبار من خارج المصدر (تتطلب الحصول على ``pytest`` >= |PyTestMinVersion| مثبت)::

    pytest sklearn

اطلع على https://scikit-learn.org/dev/developers/advanced_installation.html#testing
للحصول على المزيد من المعلومات.

    يمكن التحكم في توليد الارقام العشوائية اثناء الاختبار عن طريق ضبط ``SKLEARN_SEED`` environment variable.

لتسليم Pull request
~~~~~~~~~~~~~~~~~~~~~~~~~

قبل ارسال طلب pull request الق نظرة على ملف صفحة المساهمة الكاملة للتاكد من توافق التعليمات البرمجية الخاصة بك مع ارشاداتنا: https://scikit-learn.org/stable/developers/index.html

تاريخ المشروع
---------------

بدا هذا المشروع في عام 2007 بواسطة ديفيد كورنابيو في صيف جوجل من مشروع كود, و منذ ذلك الحين ساهم العديد من المتطوعين
انظر الى صفحة  
 `About us <https://scikit-learn.org/dev/about.html#authors>`__ 
 للحصول على قائمة المساهمين الاساسيين 

الموقع تحت الصيانة من قبل المتطوعين.

**ملاحظة**: `scikit-learn` كان يدعى سابقا ب  `scikits.learn`.

لتقديم الدعم و المساعدة
----------------

التوثيقات
~~~~~~~~~~~~~

- HTML documentation (stable release): https://scikit-learn.org
- HTML documentation (development version): https://scikit-learn.org/dev/
- FAQ: https://scikit-learn.org/stable/faq.html

للتواصل
~~~~~~~~~~~~~

- Mailing list: https://mail.python.org/mailman/listinfo/scikit-learn
- Gitter: https://gitter.im/scikit-learn/scikit-learn
- Twitter: https://twitter.com/scikit_learn
- Stack Overflow: https://stackoverflow.com/questions/tagged/scikit-learn
- Github Discussions: https://github.com/scikit-learn/scikit-learn/discussions
- Website: https://scikit-learn.org

الاقتباس
~~~~~~~~

اذا كنت تستخدم سايكيت ليرن في منشور علمي يرجى الاشارة عبر: https://scikit-learn.org/stable/about.html#citing-scikit-learn
