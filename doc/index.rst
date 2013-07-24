.. raw:: html

    <!-- Block section -->
    <div class="container", style="width:100%;">
    <div class="row-fluid">

    <!-- Classification -->
    <div class="span4 box">
    <h2 >

:ref:`Classification <supervised-learning>`

.. raw:: html

    </h2>
    <blockquote>
    <p class="box-tagline">Identifying to which set of categories a new observation belong
    to.</p>
    <strong>Applications</strong>: Spam detection, image recognition.</br>
    <strong>Algorithms</strong>:

:ref:`SVM<svm>`, :ref:`nearest neighbors<classification>`, :ref:`random forest<forest>`, ...

.. raw:: html

    <span class="example-links">
      <small class="float-right"><a href="#">Examples</a></small>
    </span>
    </blockquote>
    </div>

    <!-- Regression -->
    <div class="span4 box">
    <h2>

:ref:`Regression <supervised-learning>`

.. raw:: html

    </h2>
    <blockquote>
    <p class="box-tagline">Predicting a continuous value for a new example.</p>
    <strong>Applications</strong>: drug response, stock prices.</br>
    <strong>Algorithms</strong>:

:ref:`SVR<svm>`, :ref:`ridge regression<ridge_regression>`, :ref:`Lasso<lasso>`, ...

.. raw:: html

    <span class="example-links">
      <small class="float-right"><a href="#">Examples</a></small>
    </span>
    </blockquote>
    </div>

    <!-- Clustering -->
    <div class="span4 box">
    <h2>

:ref:`Clustering<clustering>`

.. raw:: html

    </h2>
    <blockquote>
    <p class="box-tagline">Automatic grouping of similar objects into sets.</p>
    <strong>Applications</strong>: customer segmentation, grouping experiment outcomes</br>
    <strong>Algorithms</strong>:

:ref:`k-Means<k_means>`, :ref:`spectral clustering<spectral_clustering>`, :ref:`mean-shift<mean_shift>`, ...

.. raw:: html

    <span class="example-links">
      <small class="float-right"><a href="#">Examples</a></small>
    </span>
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
    <p class="box-tagline">Reducing the number of random variables to consider.</p>
    <strong>Applications</strong>: visualization, increased efficiency</br>
    <strong>Algorithms</strong>:

:ref:`PCA<PCA>`, :ref:`Isomap<isomap>`, :ref:`non-negative matrix factorization<NMF>`, ...

.. raw:: html

    <span class="example-links">
      <small class="float-right"><a href="#">Examples</a></small>
    </span>
    </blockquote>
    </div>

    <!-- Model selection -->
    <div class="span4 box">
    <h2>

:ref:`Model Selection<model_selection>`

.. raw:: html

    </h2>
    <blockquote>
    <p class="box-tagline">Comparing, validating and choosing parameters and models.</p>
    <strong>Goal</strong>: Improved accuracy via parameter tuning</br>
    <strong>Relevant modules</strong>:
    
:ref:`grid search<grid_search>`, :ref:`cross validation<cross_validation>`, :ref:`evaluation metrics<model_evaluation>`

.. raw:: html

    <span class="example-links">
      <small class="float-right"><a href="#">Examples</a></small>
    </span>
    </blockquote>
    </div>


    <!-- Preprocessing -->
    <div class="span4 box">
    <h2>
    
:ref:`Preprocessing<preprocessing>`

.. raw:: html

    </h2>
    <blockquote>
    <p class="box-tagline">Creating and normalizing features.</p>
    <strong>Application</strong>: transforming input data such as text for use with machine learning algorithms.</br>
    <strong>Relevant modules</strong>:

:ref:`preprocessing<preprocessing>`, :ref:`feature extraction<feature_extraction>`

.. raw:: html

    <small class="float-right"><a href="#">Examples</a></small>
    </blockquote>
    </div>
    <!-- row -->
    </div>
    </div> <!-- container -->


    <div class="container" style="padding-top: 40px; width:100%">
        <div class="row-fluid">
            <!-- News -->
            <div class="span3" style="border-right: 1px solid #CCC; padding-right:5px">
                <h4 class="no-bg">News</h4>
                <ul>
                <li>The scikit-learn international code sprint is around the corner! Please, sponsor us.</li>
                <li>scikit-learn 0.13.1 is available for download.</li>
                </ul>
            </div>

            <!-- Sponsors -->
            <div class="span3" style="border-right: 1px solid #CCC; padding-right:5px">
                <h4 class="no-bg">Sponsors/Donations</h4>
                <p>Any donations are very welcome!</p>
                <form target="_top" method="post" action="https://www.paypal.com/cgi-bin/webscr">
                <input type="hidden" value="_s-xclick" name="cmd">
                <input type="hidden" value="74EYUMF3FTSW8" name="hosted_button_id">
                <input border="0" type="image" style="margin: 0 auto; position: relative; left: 6%;" alt="PayPal - The safer, easier way to pay online!" name="submit" src="https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif">
                <img border="0" width="1" height="1" src="https://www.paypalobjects.com/en_US/i/scr/pixel.gif" alt="">
                </form>
                <a href="#">Read more here ...</a>
            </div>

            <!-- Community -->
            <div class="span3" style="border-right: 1px solid #CCC; padding-right:5px">
                <h4 class="no-bg">Community</h4>
                <ul>
                <li>Appication and usage questions are best posted on <a href="#">stackoverflow.com</a> with tag sklearn.</li>
                <li>The mailing list for general discussions is scikit-learn-general@lists.sourceforge.net</li>
                <li>There is a #scikit-learn IRC channel on freenode that is frequented by devs and user.</li>
                </ul>
            </div>

            <!-- who using -->
            <div class="span3">
                <h4 class="no-bg">Who is using it</h4>
                <div id="myCarousel" class="carousel slide">
                    <ol class="carousel-indicators">
                    <li data-target="#myCarousel" data-slide-to="0" class="active"></li>
                    <li data-target="#myCarousel" data-slide-to="1"></li>
                    <li data-target="#myCarousel" data-slide-to="2"></li>
                    </ol>
                    <!-- Carousel items -->
                    <div class="carousel-inner">
                        <div class="active item"><img style="height:70px" src="img/inria.jpg"/></div>
                        <div class="item"><img style="height:70px" src="img/google.png"/></div>
                        <div class="item"><img style="height:70px" src="img/telecom.jpg"/></div>
                    </div>
                </div>
                <script>$('.carousel').carousel()</script>
            </div>

        </div>
    </div>
