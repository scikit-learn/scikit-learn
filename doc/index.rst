
..
    We are putting the title as a raw HTML so that it doesn't appear in
    the contents
    
.. raw:: html

    <h1>scikit-learn: machine learning in Python</h1>
    <style type="text/css">
    p {
        margin: 7px 0 7px 0 ;
    }
    span.linkdescr a {
        color:  #3E4349 ;
    }
    </style>

..  
   Here we are building a banner: a javascript selects randomly 4 images in 
   the list

.. only:: html

    .. |banner1| image:: auto_examples/svm/images/plot_oneclass_1.png
       :height: 140
       :target: auto_examples/svm/plot_oneclass.html

    .. |banner2| image:: auto_examples/cluster/images/plot_ward_structured_vs_unstructured_2.png
       :height: 140
       :target: auto_examples/cluster/plot_ward_structured_vs_unstructured.html

    .. |banner3| image:: auto_examples/gaussian_process/images/plot_gp_regression_1.png
       :height: 140
       :target: auto_examples/gaussian_process/plot_gp_regression.html

    .. |banner4| image:: auto_examples/cluster/images/plot_lena_ward_segmentation_1.png
       :height: 140
       :target: auto_examples/cluster/plot_lena_ward_segmentation.html

    .. |banner5| image:: auto_examples/svm/images/plot_svm_nonlinear_1.png
       :height: 140
       :target: auto_examples/svm/plot_svm_nonlinear.html

    .. |banner6| image:: auto_examples/applications/images/plot_species_distribution_modeling_1.png
       :height: 140
       :target: auto_examples/applications/plot_species_distribution_modeling.html

    .. |banner7| image:: auto_examples/gaussian_process/images/plot_gp_probabilistic_classification_after_regression_1.png
       :height: 140
       :target: auto_examples/gaussian_process/plot_gp_probabilistic_classification_after_regression.html

    .. |banner8| image:: auto_examples/ensemble/images/plot_forest_importances_faces_1.png
       :height: 140
       :target: auto_examples/ensemble/plot_forest_importances_faces.html

    .. |banner9| image:: auto_examples/svm/images/plot_weighted_samples_1.png
       :height: 140
       :target: auto_examples/svm/plot_weighted_samples.html

    .. |banner10| image:: auto_examples/linear_model/images/plot_sgd_weighted_samples_1.png
       :height: 140
       :target: auto_examples/linear_model/plot_sgd_weighted_samples.html

    .. |banner11| image:: auto_examples/cluster/images/plot_kmeans_digits_1.png
       :height: 140
       :target: auto_examples/cluster/plot_kmeans_digits.html

    .. |banner12| image:: auto_examples/decomposition/images/plot_faces_decomposition_2.png
       :height: 140
       :target: auto_examples/decomposition/plot_faces_decomposition.html

    .. |banner13| image:: auto_examples/decomposition/images/plot_faces_decomposition_3.png
       :height: 140
       :target: auto_examples/decomposition/plot_faces_decomposition.html

    .. |banner14| image:: auto_examples/images/plot_lda_vs_qda_1.png
       :height: 140
       :target: auto_examples/plot_lda_vs_qda.html

    .. |center-div| raw:: html

        <div style="text-align: center; margin: -7px 0 -10px 0;" id="banner">

    .. |end-div| raw:: html

        </div>

        <SCRIPT>
        // Function to select 4 imgs in random order from a div
        function shuffle(e) {       // pass the divs to the function
          var replace = $('<div>');
          var size = 4;
          var num_choices = e.size();

          while (size >= 1) {
            var rand = Math.floor(Math.random() * num_choices);
            var temp = e.get(rand);      // grab a random div from our set
            replace.append(temp);        // add the selected div to our new set
            e = e.not(temp); // remove our selected div from the main set
            size--;
            num_choices--;
          }
          $('#banner').html(replace.html() ); // update our container div 
                                              // with the new, randomized divs
        }
        shuffle ($('#banner a.external'));
        </SCRIPT>

    |center-div| |banner1| |banner2| |banner3| |banner4| |banner5| |banner6| |banner7| |banner8| |banner9| |banner10| |banner11| |banner12| |banner13| |banner14| |end-div|


.. topic:: Easy-to-use and general-purpose machine learning in Python

    ``scikit-learn`` is a Python module integrating classic machine
    learning algorithms in the tightly-knit scientific Python
    world (`numpy <http://numpy.scipy.org>`_, `scipy
    <http://www.scipy.org>`_, `matplotlib
    <http://matplotlib.sourceforge.net/>`_).
    It aims to provide simple and efficient solutions to learning
    problems, accessible to everybody and reusable in various
    contexts: **machine-learning as a versatile tool for science and
    engineering**.


.. raw:: html

  <table class="contentstable" style="width: 100% ; margin-top: -8px">
    <tr valign="top"><td width="28%">
      <p class="biglink"><a class="biglink" href="supervised_learning.html">
                Supervised learning</a><br/>
         <span class="linkdescr">
                <a href="modules/svm.html">Support vector machines</a>,
                <a href="modules/linear_model.html">linear models</a>,
                <a href="modules/naive_bayes.html">naives Bayes</a>,
                <a href="modules/gaussian_process.html">Gaussian process</a>...
         </span></p>
    </td><td align="center" width="32%">
      <p class="biglink"><a class="biglink" href="unsupervised_learning.html">
        Unsupervised learning</a><br/>
         <span class="linkdescr">
                <a href="modules/clustering.html">Clustering</a>,
                <a href="modules/mixture.html">Gaussian mixture models</a>,
                <a href="modules/manifold.html">manifold learning</a>,
                <a href="modules/decomposition.html">matrix factorization</a>,
                <a href="modules/covariance.html">covariance</a>...
         </span></p>
    </td><td align="right" width="30%">
      <p class="biglink"><a class="biglink" href="index.html#user-guide">
        And much more</a><br/>
         <span class="linkdescr">
                <a href="model_selection.html">Model selection</a>,
                <a href="datasets/index.html">datasets</a>,
                <a href="modules/feature_extraction.html">feature extraction...</a>
                <strong>See below</strong>.</span></p>
    </td></tr>
  </table>

**License:** Open source, commercially usable: **BSD license** (3 clause)

.. include:: includes/big_toc_css.rst

Documentation for scikit-learn **version** |release|. For other versions and
printable format, see :ref:`documentation_resources`.

User Guide
==========

.. toctree::
   :maxdepth: 2

   user_guide.rst

Example Gallery
===============

.. toctree::
   :maxdepth: 2

   auto_examples/index


Development
===========
.. toctree::
   :maxdepth: 2

   developers/index
   developers/performance
   developers/utilities
   developers/debugging
   about

.. toctree::
   :hidden:

   support
   whats_new
