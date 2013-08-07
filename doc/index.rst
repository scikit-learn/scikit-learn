.. We want the previous/next button to work on the user guide pages and on the
   API. We have to put the doctree so that sphinx populates the
   rellinks. Yet, we don't want it to be displayed on the main page, hence we
   don't display through the css.

.. raw:: html

   <div class="no-display">

.. toctree::

    tutorial/index
    user_guide
    auto_examples/index
    support
    whats_new
    presentations
    about
    documentation
    datasets/index
    datasets/covtype
    datasets/labeled_faces
    datasets/mldata
    datasets/olivetti_faces
    datasets/twenty_newsgroups
    modules/classes
    testimonials/testimonials
    developers/index
    developers/debugging
    developers/maintainer
    developers/performance
    developers/utilities
    install
    tutorial/basic/tutorial
    tutorial/machine_learning_map/index



.. raw:: html

   </div>


.. This is were the real work starts.


.. raw:: html

    <!-- Block section -->
    <div class="container-index">
    <div class="container index-upper">
    <div class="row-fluid">

    <!-- Classification -->
    <div class="span4 box">
    <h2 >

:ref:`Classification <supervised-learning>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Identifying to which set of categories a new observation belong
    to.</p>
    <div class="box-links">
    <strong>Applications</strong>: Spam detection, Image recognition.</br>
    <strong>Algorithms</strong>:&nbsp;

:ref:`SVM<svm>`, :ref:`nearest neighbors<classification>`, :ref:`random forest<forest>`, ...

.. raw:: html

    <small class="float-right box-example-links">

:ref:`Examples<general_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>

    <!-- Regression -->
    <div class="span4 box">
    <h2>

:ref:`Regression <supervised-learning>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Predicting a continuous value for a new example.</p>
    <div class="box-links">
    <strong>Applications</strong>: Drug response, Stock prices.</br>
    <strong>Algorithms</strong>:&nbsp;

:ref:`SVR<svm>`, :ref:`ridge regression<ridge_regression>`, :ref:`Lasso<lasso>`, ...

.. raw:: html

    <small class="float-right box-example-links">

:ref:`Examples<general_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>

    <!-- Clustering -->
    <div class="span4 box">
    <h2>

:ref:`Clustering<clustering>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Automatic grouping of similar objects into sets.</p>
    <div class="box-links">
    <strong>Applications</strong>: Customer segmentation, Grouping experiment outcomes</br>
    <strong>Algorithms</strong>:&nbsp;

:ref:`k-Means<k_means>`, :ref:`spectral clustering<spectral_clustering>`, :ref:`mean-shift<mean_shift>`, ...

.. raw:: html

    <small class="float-right example-links">

:ref:`Examples<cluster_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>

    <!-- row -->
    </div>
    <div class="row-fluid">

    <!-- Dimension reduction -->
    <div class="span4 box">
    <h2>

:ref:`Dimensionality reduction<decompositions>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Reducing the number of random variables to consider.</p>
    <div class="box-links">
    <strong>Applications</strong>: Visualization, Increased efficiency</br>
    <strong>Algorithms</strong>:&nbsp;

:ref:`PCA<PCA>`, :ref:`Isomap<isomap>`, :ref:`non-negative matrix factorization<NMF>`.

.. raw:: html

    <small class="float-right example-links">

:ref:`Examples<decomposition_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>

    <!-- Model selection -->
    <div class="span4 box">
    <h2>

:ref:`Model selection<model_selection>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Comparing, validating and choosing parameters and models.</p>
    <div class="box-links">
    <strong>Goal</strong>: Improved accuracy via parameter tuning</br>
    <strong>Modules</strong>:&nbsp;

:ref:`grid search<grid_search>`, :ref:`cross validation<cross_validation>`, :ref:`metrics<model_evaluation>`.

.. raw:: html

    <small class="float-right example-links">

:ref:`Examples<general_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>


    <!-- Preprocessing -->
    <div class="span4 box">
    <h2>

:ref:`Preprocessing<preprocessing>`

.. raw:: html

    </h2>
    <blockquote>
    <p>Feature extraction and normalization.</p>
    <div class="box-links">
    <strong>Application</strong>: Transforming input data such as text for use with machine learning algorithms.</br>
    <strong>Modules</strong>:&nbsp;

:ref:`preprocessing<preprocessing>`, :ref:`feature extraction<feature_extraction>`.

.. raw:: html

    <span class="example-links">
    <small class="float-right example-links">

:ref:`Examples<general_examples>`

.. raw:: html

    </small>
    </div>
    </blockquote>
    </div>

    <!-- row -->
    </div>
    </div> <!-- container -->


    <div class="container index-lower">
        <div class="row-fluid">
            <!-- News -->
            <div class="span4">
                <h4>News</h4>
                <ul>
                <li><em>On-going development:</em>
                <a href="whats_new.html"><em>What's new</em> (changelog)</a>
                </li>
                <li><em>August 2013.</em> scikit-learn 0.14 is available for download (<a href="whats_new.html">Changelog</a>).
                </li>
                <li><em>July 22-28th, 2013: international sprint.</em>
                During this week-long sprint, we gathered most of the core
                developers in Paris.
                We want to thank our sponsors, our
                hosts <a href="http://www.telecom-paristech.fr/">Télécom ParisTech</a>
                and <a href="http://www.tinyclues.com/">tinyclues</a>, and
                donations that helped fund this event.

                </ul>
            </div>

            <!-- Community -->
            <div class="span4">
                <h4>Community</h4>
                <ul>
                <li><em>Questions?</em> See <a href="http://stackoverflow.com/questions/tagged/scikit-learn">stackoverflow</a> # scikit-learn</li>
                <li><em>Mailing list:</em> <a href="https://lists.sourceforge.net/lists/listinfo/scikit-learn-general">scikit-learn-general@lists.sourceforge.net</a></li>
                <li><em>IRC:</em> #scikit-learn @ <a href="http://webchat.freenode.net/">freenode</a></li>
                </ul>

                <form target="_top" id="paypal-form" method="post" action="https://www.paypal.com/cgi-bin/webscr">
                    <input type="hidden" value="_s-xclick" name="cmd">
                    <input type="hidden" value="74EYUMF3FTSW8" name="hosted_button_id">
                </form>

                <a class="btn btn-warning btn-big" onclick="document.getElementById('paypal-form').submit(); return false;">Help us, <strong>donate!</strong></a>
                <a class="btn btn-warning btn-big cite-us" href="./about.html#citing-scikit-learn"><strong>Cite us!</strong></a>

                <small style="display: block; margin-top: 10px"><a href="about.html#funding">Read more about donations</a></small>
            </div>

            <!-- who using -->
            <div class="span4">
                <h4>Who uses scikit-learn?</h4>

                <div id="testimonials_carousel" class="carousel slide">
                    <div class="carousel-inner">
                        <div class="active item">
                          <img src="_images/inria.jpg" class="thumbnail" />
                          <p>
                          <em>"We use scikit-learn to support leading-edge basic research [...]"</em>
                          </p>
                        </div>
                        <div class="item">
                          <img src="_images/evernote.png" class="thumbnail" />
                          <p>
                          <em>"For these tasks, we relied on the excellent scikit-learn package for Python."</em>
                          </p>
                        </div>
                        <div class="item">
                          <img src="_images/telecomparistech.jpg"
                               class="thumbnail" />
                          <p>
                          <em>"The great benefit of scikit-learn is its fast learning curve [...]"</em>
                          </p>
                        </div>
                        <div class="item">
                          <img src="_images/aweber.png" class="thumbnail" />
                          <p>
                          <em>"It allows us to do AWesome stuff we would not otherwise accomplish"</em>
                          </p>
                        </div>
                    </div>
                </div>
                <p align="right">
                <small class="example-link">
                <a href="testimonials/testimonials.html">More testimonials</a>
                </small>
                </p>
            </div>

        </div>
    </div>
    </div>


    <script>
      $('#testimonials_carousel').carousel()
    </script>


