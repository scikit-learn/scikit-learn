.. We want the previous/next button to work on the user guide pages and on the
   API. We have to put the doctree so that sphinx populates the
   rellinks. Yet, we don't want it to be displayed on the main page, hence we
   don't display through the css.

.. raw:: html

   <div class="no-display">

.. toctree::

    tutorials
    user_guide
    auto_examples/index
    support
    whats_new
    presentations


.. raw:: html

   </div>


.. This is were the real work starts.


.. raw:: html

    <!-- Block section -->
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

:ref:`Model Selection<model_selection>`

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
    <p>Creating and normalizing features.</p>
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
            <div class="span6">
                <h4>News</h4>
                <ul>
                <li><em>July 2013.</em>
                The sprint is drawing to a close and we will push a 0.14
                release candidate soon. Or new shiny website is up, too.
                During this week-long sprint, we gathered most of the core
                developers in Paris for a blend of fun and productivity.
                <!--
                Here are some of the biggest changes in the upcoming version:
                <ul>
                <li>Python 3 support</li>
                <li>Ensembles of Randomized Trees speed improvements</li>
                <li>Restricted Boltzman Machines</li>
                <li>Missing data imputation</li>
                <li>Bi-clustering</li>
                </ul>
                -->

                We want to thank again our sponsors, in particular Télécom ParisTech 
                and Tinyclues, who hosted our sprint during the week and the weekend respectively.
                Private donations also helped fund this event.

                </li>
                <li><em>July 2013.</em>
                The sprint has started! It is great to meet so many core developers in person. Also,
                some new contributors managed to make it, which makes for a great and productive crowd.
                We kick-started the sprint with merges of some very important contributions and hope
                we can continue this pace through the week.
                </li>
                <!--
                    <li><em>July 2013.</em> The scikit-learn international code sprint is around the corner! Please, sponsor us.</li>
                    <li><em>February 2013.</em> scikit-learn 0.13.1 is available for download. <a href="whats_new.html">Changelog</a> .</li>
                    -->
                </ul>
            </div>

            <!-- Community -->
            <div class="span6">
                <h4>Community</h4>
                <ul>
                <li><em>Questions?</em> See <a href="http://stackoverflow.com/questions/tagged/scikit-learn">stackoverflow</a> # scikit-learn for usage questions</li>
                <li><em>Mailing list:</em> scikit-learn-general@lists.sourceforge.net</li>
                <li><em>IRC:</em> #scikit-learn @ <a href="http://webchat.freenode.net/">freenode</a></li>
                <li><em>Help us:</em> <button class="btn btn-small" onclick="document.getElementById('paypal-form').submit(); return false;">Donate!</button> (<a href="about.html#funding">read more</a>)</li>
                    <form target="_top" id="paypal-form" method="post" action="https://www.paypal.com/cgi-bin/webscr">
                    <input type="hidden" value="_s-xclick" name="cmd">
                    <input type="hidden" value="74EYUMF3FTSW8" name="hosted_button_id">
                    </form>
                </ul>
            </div>

            <!-- who using -->
            <!--
            <div class="span4">
                <h4>Who is using scikit-learn?</h4>

                </h4>
                <div id="myCarousel" class="carousel slide">
                    <div class="carousel-inner">
                        <div class="active item"><img src="_images/inria.jpg" class="thumbnail" /><br /> <em>-- Great stuff!</em></div>
                        <div class="item"><img src="_static/img/google.png" class="thumbnail" /><br /> <em>-- So good!</em></div>
                    </div>
                    <div style="margin-top: 5px"><a href="#">More testimonials</a></div>
                </div>
                <script>$('#myCarousel').carousel()</script>
            </div>
            -->

        </div>
    </div>
