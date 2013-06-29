

<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>scikit-learn/examples/applications/plot_stock_market.py at master · scikit-learn/scikit-learn</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://a248.e.akamai.net/assets.github.com/images/modules/logos_page/Octocat.png">
    <link rel="assets" href="https://a248.e.akamai.net/assets.github.com/">
    <link rel="xhr-socket" href="/_sockets" />
    
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="4869318" name="octolytics-actor-id" /><meta content="ccoovrey" name="octolytics-actor-login" /><meta content="e7ec2a90e7a4c651ace2744342b6573f111e37aa1cae3182d43b226b258c4945" name="octolytics-actor-hash" />

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="U8307+Ij9dtPLUiTYtXY8elEEDggpAoJuK1AB3BThXA=" name="csrf-token" />

    <link href="https://a248.e.akamai.net/assets.github.com/assets/github-be6069435cb0250bd316958375c0de713801bdf5.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://a248.e.akamai.net/assets.github.com/assets/github2-91a602bfb8b2a53f90aa7aea127a69e4f9ef1f11.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://a248.e.akamai.net/assets.github.com/assets/frameworks-1f72571b966545f4e27481a3b0ebbeeed4f2f139.js" type="text/javascript"></script>
      <script src="https://a248.e.akamai.net/assets.github.com/assets/github-3b51dd74a94c713c22309e373955e5fa02a3bb65.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="f5ee7511175547406ed413cec6ae4c9b">

        <link data-pjax-transient rel='permalink' href='/scikit-learn/scikit-learn/blob/b9877a736517355b6ab4bc1599ef7b26696bb00a/examples/applications/plot_stock_market.py'>

  <meta property="og:title" content="scikit-learn"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/scikit-learn/scikit-learn"/>
  <meta property="og:image" content="https://a248.e.akamai.net/assets.github.com/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="scikit-learn: machine learning in Python"/>

  <meta name="description" content="scikit-learn: machine learning in Python" />

  <meta content="365630" name="octolytics-dimension-user_id" /><meta content="scikit-learn" name="octolytics-dimension-user_login" /><meta content="843222" name="octolytics-dimension-repository_id" /><meta content="scikit-learn/scikit-learn" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="843222" name="octolytics-dimension-repository_network_root_id" /><meta content="scikit-learn/scikit-learn" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/scikit-learn/scikit-learn/commits/master.atom" rel="alternate" title="Recent Commits to scikit-learn:master" type="application/atom+xml" />

  </head>


  <body class="logged_in page-blob linux vis-public env-production  kill-the-chrome">

    <div class="wrapper">
      
      
      

      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/">
  <span class="mega-octicon octicon-mark-github"></span>
</a>

    <div class="divider-vertical"></div>

      <a href="/scikit-learn/scikit-learn/notifications" class="notification-indicator tooltipped downwards contextually-unread" title="You have unread notifications in this repository">
    <span class="mail-status unread"></span>
  </a>
  <div class="divider-vertical"></div>


      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="ccoovrey"
      data-repo="scikit-learn/scikit-learn"
      data-branch="master"
      data-sha="e627a0f22161ecce710a8d0a3aa01a2d870076c1"
  >

    <input type="hidden" name="nwo" value="scikit-learn/scikit-learn" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
            <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    

  

    <ul id="user-links">
      <li>
        <a href="/ccoovrey" class="name">
          <img height="20" src="https://secure.gravatar.com/avatar/2956a24d7299201437a3e7ea54b37c55?s=140&amp;d=https://a248.e.akamai.net/assets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png" width="20" /> ccoovrey
        </a>
      </li>

        <li>
          <a href="/new" id="new_repo" class="tooltipped downwards" title="Create a new repo">
            <span class="octicon octicon-repo-create"></span>
          </a>
        </li>

        <li>
          <a href="/settings/profile" id="account_settings"
            class="tooltipped downwards"
            title="Account settings (You have no verified emails)">
            <span class="octicon octicon-tools"></span>
          </a>
            <span class="settings-warning">!</span>
        </li>
        <li>
          <a class="tooltipped downwards" href="/logout" data-method="post" id="logout" title="Sign out">
            <span class="octicon octicon-log-out"></span>
          </a>
        </li>

    </ul>


<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo-create"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-list-unordered"></span> New organization</a>
  </li>



    <li class="section-title">
      <span title="scikit-learn/scikit-learn">This repository</span>
    </li>
    <li>
      <a href="/scikit-learn/scikit-learn/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
    </li>
</ul>

</div>


    
  </div>
</div>

      

      

<div class="flash-global flash-warn">
<div class="container">

    <h2>
      You don't have any verified emails.  We recommend <a href="https://github.com/settings/emails">verifying</a> at least one email.
    </h2>
    <p>
      Email verification helps our support team help you in case you have any email issues or lose your password.
    </p>












</div>
</div>



          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">


    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="U8307+Ij9dtPLUiTYtXY8elEEDggpAoJuK1AB3BThXA=" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="843222" />

    <div class="select-menu js-menu-container js-select-menu">
      <span class="minibutton select-menu-button  js-menu-target">
        <span class="js-select-button">
          <span class="octicon octicon-eye-watch"></span>
          Watch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-remove-close js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container">

            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for discussions in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-watch"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-unwatch"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for discussions in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

    <li class="js-toggler-container js-social-container starring-container ">
      <a href="/scikit-learn/scikit-learn/unstar" class="minibutton with-count js-toggler-target star-button starred upwards" title="Unstar this repo" data-remote="true" data-method="post" rel="nofollow">
        <span class="octicon octicon-star-delete"></span>
        <span class="text">Unstar</span>
      </a>
      <a href="/scikit-learn/scikit-learn/star" class="minibutton with-count js-toggler-target star-button unstarred upwards" title="Star this repo" data-remote="true" data-method="post" rel="nofollow">
        <span class="octicon octicon-star"></span>
        <span class="text">Star</span>
      </a>
      <a class="social-count js-social-count" href="/scikit-learn/scikit-learn/stargazers">1,384</a>
    </li>

        <li>
          <a href="/scikit-learn/scikit-learn/fork" class="minibutton with-count js-toggler-target fork-button lighter upwards" title="Fork this repo" rel="nofollow" data-method="post">
            <span class="octicon octicon-git-branch-create"></span>
            <span class="text">Fork</span>
          </a>
          <a href="/scikit-learn/scikit-learn/network" class="social-count">700</a>
        </li>


</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author">
            <a href="/scikit-learn" class="url fn" itemprop="url" rel="author"><span itemprop="title">scikit-learn</span></a></span
          ><span class="repohead-name-divider">/</span><strong
          ><a href="/scikit-learn/scikit-learn" class="js-current-repository js-repo-home-link">scikit-learn</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container
            ">

          <div class="repository-sidebar">

              

<div class="repo-nav repo-nav-full js-repository-container-pjax js-octicon-loaders">
  <div class="repo-nav-contents">
    <ul class="repo-menu">
      <li class="tooltipped leftwards" title="Code">
        <a href="/scikit-learn/scikit-learn" class="js-selected-navigation-item selected" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /scikit-learn/scikit-learn">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped leftwards" title="Issues">
          <a href="/scikit-learn/scikit-learn/issues" class="js-selected-navigation-item js-disable-pjax" data-gotokey="i" data-selected-links="repo_issues /scikit-learn/scikit-learn/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>389</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped leftwards" title="Pull Requests"><a href="/scikit-learn/scikit-learn/pulls" class="js-selected-navigation-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /scikit-learn/scikit-learn/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>110</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/scikit-learn/scikit-learn/wiki" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_wiki /scikit-learn/scikit-learn/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>


    </ul>
    <div class="repo-menu-separator"></div>
    <ul class="repo-menu">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/scikit-learn/scikit-learn/pulse" class="js-selected-navigation-item " data-pjax="true" data-selected-links="pulse /scikit-learn/scikit-learn/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/scikit-learn/scikit-learn/graphs" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_graphs repo_contributors /scikit-learn/scikit-learn/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/scikit-learn/scikit-learn/network" class="js-selected-navigation-item js-disable-pjax" data-selected-links="repo_network /scikit-learn/scikit-learn/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

    </ul>

  </div>
</div>


              <div class="only-with-full-nav">

                

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><strong>HTTPS</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/scikit-learn/scikit-learn.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/scikit-learn/scikit-learn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=clone">
  <h3><strong>SSH</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="git@github.com:scikit-learn/scikit-learn.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="git@github.com:scikit-learn/scikit-learn.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><strong>Subversion</strong> checkout URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/scikit-learn/scikit-learn" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/scikit-learn/scikit-learn" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>



<p class="clone-options">You can clone with
    <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
    <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
    <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>,
  and <a href="https://help.github.com/articles/which-remote-url-should-i-use">other methods.</a>
</p>




                  <a href="/scikit-learn/scikit-learn/archive/master.zip"
                     class="minibutton sidebar-button"
                     title="Download this repository as a zip file"
                     rel="nofollow">
                    <span class="octicon octicon-cloud-download"></span>
                    Download ZIP
                  </a>

              </div>
          </div>

          <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
            


<!-- blob contrib key: blob_contributors:v21:07366f5e266256412738cae9d8340609 -->
<!-- blob contrib frag key: views10/v8/blob_contributors:v21:07366f5e266256412738cae9d8340609 -->


      <p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

        <a href="/scikit-learn/scikit-learn/find/master" data-pjax data-hotkey="t" style="display:none">Show File Finder</a>

        <div class="file-navigation">
          


<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.6.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.6.X" rel="nofollow" title="0.6.X">0.6.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.7.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.7.X" rel="nofollow" title="0.7.X">0.7.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.8.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.8.X" rel="nofollow" title="0.8.X">0.8.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.9.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.9.X" rel="nofollow" title="0.9.X">0.9.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.10.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.10.X" rel="nofollow" title="0.10.X">0.10.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.11.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.11.X" rel="nofollow" title="0.11.X">0.11.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.12.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.12.X" rel="nofollow" title="0.12.X">0.12.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.13.X/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.13.X" rel="nofollow" title="0.13.X">0.13.X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian" rel="nofollow" title="debian">debian</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/invalid-n-folds/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="invalid-n-folds" rel="nofollow" title="invalid-n-folds">invalid-n-folds</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/master/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="master" rel="nofollow" title="master">master</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/py3k/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="py3k" rel="nofollow" title="py3k">py3k</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/sprint01/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="sprint01" rel="nofollow" title="sprint01">sprint01</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.12.0-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.12.0-1" rel="nofollow" title="debian/0.12.0-1">debian/0.12.0-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.11.0-2/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.11.0-2" rel="nofollow" title="debian/0.11.0-2">debian/0.11.0-2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.11.0-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.11.0-1" rel="nofollow" title="debian/0.11.0-1">debian/0.11.0-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.10.0-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.10.0-1" rel="nofollow" title="debian/0.10.0-1">debian/0.10.0-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.9.0.dfsg-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.9.0.dfsg-1" rel="nofollow" title="debian/0.9.0.dfsg-1">debian/0.9.0.dfsg-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.8.1.dfsg-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.8.1.dfsg-1" rel="nofollow" title="debian/0.8.1.dfsg-1">debian/0.8.1.dfsg-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.8.0.dfsg-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.8.0.dfsg-1" rel="nofollow" title="debian/0.8.0.dfsg-1">debian/0.8.0.dfsg-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.7.1.dfsg-3/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.7.1.dfsg-3" rel="nofollow" title="debian/0.7.1.dfsg-3">debian/0.7.1.dfsg-3</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.7.1.dfsg-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.7.1.dfsg-1" rel="nofollow" title="debian/0.7.1.dfsg-1">debian/0.7.1.dfsg-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.6.0.dfsg-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.6.0.dfsg-1" rel="nofollow" title="debian/0.6.0.dfsg-1">debian/0.6.0.dfsg-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.5-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.5-1" rel="nofollow" title="debian/0.5-1">debian/0.5-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.4-3/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.4-3" rel="nofollow" title="debian/0.4-3">debian/0.4-3</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.4-2/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.4-2" rel="nofollow" title="debian/0.4-2">debian/0.4-2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.4-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.4-1" rel="nofollow" title="debian/0.4-1">debian/0.4-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.3-4/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.3-4" rel="nofollow" title="debian/0.3-4">debian/0.3-4</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.3-3/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.3-3" rel="nofollow" title="debian/0.3-3">debian/0.3-3</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.3-2/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.3-2" rel="nofollow" title="debian/0.3-2">debian/0.3-2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.3-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.3-1" rel="nofollow" title="debian/0.3-1">debian/0.3-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/debian/0.2+svn625-1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="debian/0.2+svn625-1" rel="nofollow" title="debian/0.2+svn625-1">debian/0.2+svn625-1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.13-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.13-branching" rel="nofollow" title="0.13-branching">0.13-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.13.1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.13.1" rel="nofollow" title="0.13.1">0.13.1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.13/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.13" rel="nofollow" title="0.13">0.13</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.12-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.12-branching" rel="nofollow" title="0.12-branching">0.12-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.12.1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.12.1" rel="nofollow" title="0.12.1">0.12.1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.12/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.12" rel="nofollow" title="0.12">0.12</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.11-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.11-branching" rel="nofollow" title="0.11-branching">0.11-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.11-beta/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.11-beta" rel="nofollow" title="0.11-beta">0.11-beta</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.11/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.11" rel="nofollow" title="0.11">0.11</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.10-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.10-branching" rel="nofollow" title="0.10-branching">0.10-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.10/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.10" rel="nofollow" title="0.10">0.10</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.9-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.9-branching" rel="nofollow" title="0.9-branching">0.9-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.9/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.9" rel="nofollow" title="0.9">0.9</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.8-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.8-branching" rel="nofollow" title="0.8-branching">0.8-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.8.1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.8.1" rel="nofollow" title="0.8.1">0.8.1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.8/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.8" rel="nofollow" title="0.8">0.8</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.7-branching/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.7-branching" rel="nofollow" title="0.7-branching">0.7-branching</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.7.1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.7.1" rel="nofollow" title="0.7.1">0.7.1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.7/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.7" rel="nofollow" title="0.7">0.7</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.6-rc/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.6-rc" rel="nofollow" title="0.6-rc">0.6-rc</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.6.0/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.6.0" rel="nofollow" title="0.6.0">0.6.0</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.5.rc3/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.5.rc3" rel="nofollow" title="0.5.rc3">0.5.rc3</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.5.rc2/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.5.rc2" rel="nofollow" title="0.5.rc2">0.5.rc2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.5.rc/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.5.rc" rel="nofollow" title="0.5.rc">0.5.rc</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.5/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.5" rel="nofollow" title="0.5">0.5</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.4/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.4" rel="nofollow" title="0.4">0.4</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.3/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.3" rel="nofollow" title="0.3">0.3</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.2-beta/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.2-beta" rel="nofollow" title="0.2-beta">0.2-beta</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.2/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.2" rel="nofollow" title="0.2">0.2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.1-beta/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.1-beta" rel="nofollow" title="0.1-beta">0.1-beta</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/scikit-learn/scikit-learn/blob/0.1/examples/applications/plot_stock_market.py" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="0.1" rel="nofollow" title="0.1">0.1</a>
            </div> <!-- /.select-menu-item -->
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

          <div class="breadcrumb">
            <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/scikit-learn/scikit-learn" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">scikit-learn</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/scikit-learn/scikit-learn/tree/master/examples" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">examples</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/scikit-learn/scikit-learn/tree/master/examples/applications" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">applications</span></a></span><span class="separator"> / </span><strong class="final-path">plot_stock_market.py</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="examples/applications/plot_stock_market.py" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
          </div>
        </div>


        <div class="commit commit-loader file-history-tease js-deferred-content" data-url="/scikit-learn/scikit-learn/contributors/master/examples/applications/plot_stock_market.py">
          Fetching contributors…

          <div class="participation">
            <p class="loader-loading"><img alt="Octocat-spinner-32-eaf2f5" height="16" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" /></p>
            <p class="loader-error">Cannot retrieve contributors at this time</p>
          </div>
        </div>


        <div id="files" class="bubble">
          <div class="file">
            <div class="meta">
              <div class="info">
                <span class="icon"><b class="octicon octicon-file-text"></b></span>
                <span class="mode" title="File Mode">file</span>
                  <span>259 lines (214 sloc)</span>
                <span>8.216 kb</span>
              </div>
              <div class="actions">
                <div class="button-group">
                        <a class="minibutton tooltipped leftwards"
                           title="Clicking this button will automatically fork this project so you can edit the file"
                           href="/scikit-learn/scikit-learn/edit/master/examples/applications/plot_stock_market.py"
                           data-method="post" rel="nofollow">Edit</a>
                  <a href="/scikit-learn/scikit-learn/raw/master/examples/applications/plot_stock_market.py" class="button minibutton " id="raw-url">Raw</a>
                    <a href="/scikit-learn/scikit-learn/blame/master/examples/applications/plot_stock_market.py" class="button minibutton ">Blame</a>
                  <a href="/scikit-learn/scikit-learn/commits/master/examples/applications/plot_stock_market.py" class="button minibutton " rel="nofollow">History</a>
                </div><!-- /.button-group -->
              </div><!-- /.actions -->

            </div>
                <div class="blob-wrapper data type-python js-blob-data">
      <table class="file-code file-diff">
        <tr class="file-code-line">
          <td class="blob-line-nums">
            <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>

          </td>
          <td class="blob-line-code">
                  <div class="highlight"><pre><div class='line' id='LC1'><span class="sd">&quot;&quot;&quot;</span></div><div class='line' id='LC2'><br/></div><div class='line' id='LC3'><span class="sd">.. _stock_market:</span></div><div class='line' id='LC4'><br/></div><div class='line' id='LC5'><span class="sd">=======================================</span></div><div class='line' id='LC6'><span class="sd">Visualizing the stock market structure</span></div><div class='line' id='LC7'><span class="sd">=======================================</span></div><div class='line' id='LC8'><br/></div><div class='line' id='LC9'><span class="sd">This example employs several unsupervised learning techniques to extract</span></div><div class='line' id='LC10'><span class="sd">the stock market structure from variations in historical quotes.</span></div><div class='line' id='LC11'><br/></div><div class='line' id='LC12'><span class="sd">The quantity that we use is the daily variation in quote price: quotes</span></div><div class='line' id='LC13'><span class="sd">that are linked tend to cofluctuate during a day.</span></div><div class='line' id='LC14'><br/></div><div class='line' id='LC15'><br/></div><div class='line' id='LC16'><span class="sd">Learning a graph structure</span></div><div class='line' id='LC17'><span class="sd">--------------------------</span></div><div class='line' id='LC18'><br/></div><div class='line' id='LC19'><span class="sd">We use sparse inverse covariance estimation to find which quotes are</span></div><div class='line' id='LC20'><span class="sd">correlated conditionally on the others. Specifically, sparse inverse</span></div><div class='line' id='LC21'><span class="sd">covariance gives us a graph, that is a list of connection. For each</span></div><div class='line' id='LC22'><span class="sd">symbol, the symbols that it is connected too are those useful to explain</span></div><div class='line' id='LC23'><span class="sd">its fluctuations.</span></div><div class='line' id='LC24'><br/></div><div class='line' id='LC25'><span class="sd">Clustering</span></div><div class='line' id='LC26'><span class="sd">----------</span></div><div class='line' id='LC27'><br/></div><div class='line' id='LC28'><span class="sd">We use clustering to group together quotes that behave similarly. Here,</span></div><div class='line' id='LC29'><span class="sd">amongst the :ref:`various clustering techniques &lt;clustering&gt;` available</span></div><div class='line' id='LC30'><span class="sd">in the scikit-learn, we use :ref:`affinity_propagation` as it does</span></div><div class='line' id='LC31'><span class="sd">not enforce equal-size clusters, and it can choose automatically the</span></div><div class='line' id='LC32'><span class="sd">number of clusters from the data.</span></div><div class='line' id='LC33'><br/></div><div class='line' id='LC34'><span class="sd">Note that this gives us a different indication than the graph, as the</span></div><div class='line' id='LC35'><span class="sd">graph reflects conditional relations between variables, while the</span></div><div class='line' id='LC36'><span class="sd">clustering reflects marginal properties: variables clustered together can</span></div><div class='line' id='LC37'><span class="sd">be considered as having a similar impact at the level of the full stock</span></div><div class='line' id='LC38'><span class="sd">market.</span></div><div class='line' id='LC39'><br/></div><div class='line' id='LC40'><span class="sd">Embedding in 2D space</span></div><div class='line' id='LC41'><span class="sd">---------------------</span></div><div class='line' id='LC42'><br/></div><div class='line' id='LC43'><span class="sd">For visualization purposes, we need to lay out the different symbols on a</span></div><div class='line' id='LC44'><span class="sd">2D canvas. For this we use :ref:`manifold` techniques to retrieve 2D</span></div><div class='line' id='LC45'><span class="sd">embedding.</span></div><div class='line' id='LC46'><br/></div><div class='line' id='LC47'><br/></div><div class='line' id='LC48'><span class="sd">Visualization</span></div><div class='line' id='LC49'><span class="sd">-------------</span></div><div class='line' id='LC50'><br/></div><div class='line' id='LC51'><span class="sd">The output of the 3 models are combined in a 2D graph where nodes</span></div><div class='line' id='LC52'><span class="sd">represents the stocks and edges the:</span></div><div class='line' id='LC53'><br/></div><div class='line' id='LC54'><span class="sd">- cluster labels are used to define the color of the nodes</span></div><div class='line' id='LC55'><span class="sd">- the sparse covariance model is used to display the strength of the edges</span></div><div class='line' id='LC56'><span class="sd">- the 2D embedding is used to position the nodes in the plan</span></div><div class='line' id='LC57'><br/></div><div class='line' id='LC58'><span class="sd">This example has a fair amount of visualization-related code, as</span></div><div class='line' id='LC59'><span class="sd">visualization is crucial here to display the graph. One of the challenge</span></div><div class='line' id='LC60'><span class="sd">is to position the labels minimizing overlap. For this we use an</span></div><div class='line' id='LC61'><span class="sd">heuristic based on the direction of the nearest neighbor along each</span></div><div class='line' id='LC62'><span class="sd">axis.</span></div><div class='line' id='LC63'><span class="sd">&quot;&quot;&quot;</span></div><div class='line' id='LC64'><span class="k">print</span><span class="p">(</span><span class="n">__doc__</span><span class="p">)</span></div><div class='line' id='LC65'><br/></div><div class='line' id='LC66'><span class="c"># Author: Gael Varoquaux gael.varoquaux@normalesup.org</span></div><div class='line' id='LC67'><span class="c"># License: BSD 3 clause</span></div><div class='line' id='LC68'><br/></div><div class='line' id='LC69'><span class="kn">import</span> <span class="nn">datetime</span></div><div class='line' id='LC70'><br/></div><div class='line' id='LC71'><span class="kn">import</span> <span class="nn">numpy</span> <span class="kn">as</span> <span class="nn">np</span></div><div class='line' id='LC72'><span class="kn">import</span> <span class="nn">pylab</span> <span class="kn">as</span> <span class="nn">pl</span></div><div class='line' id='LC73'><span class="kn">from</span> <span class="nn">matplotlib</span> <span class="kn">import</span> <span class="n">finance</span></div><div class='line' id='LC74'><span class="kn">from</span> <span class="nn">matplotlib.collections</span> <span class="kn">import</span> <span class="n">LineCollection</span></div><div class='line' id='LC75'><br/></div><div class='line' id='LC76'><span class="kn">from</span> <span class="nn">sklearn</span> <span class="kn">import</span> <span class="n">cluster</span><span class="p">,</span> <span class="n">covariance</span><span class="p">,</span> <span class="n">manifold</span></div><div class='line' id='LC77'><br/></div><div class='line' id='LC78'><span class="c">###############################################################################</span></div><div class='line' id='LC79'><span class="c"># Retrieve the data from Internet</span></div><div class='line' id='LC80'><br/></div><div class='line' id='LC81'><span class="c"># Choose a time period reasonnably calm (not too long ago so that we get</span></div><div class='line' id='LC82'><span class="c"># high-tech firms, and before the 2008 crash)</span></div><div class='line' id='LC83'><span class="n">d1</span> <span class="o">=</span> <span class="n">datetime</span><span class="o">.</span><span class="n">datetime</span><span class="p">(</span><span class="mi">2003</span><span class="p">,</span> <span class="mo">01</span><span class="p">,</span> <span class="mo">01</span><span class="p">)</span></div><div class='line' id='LC84'><span class="n">d2</span> <span class="o">=</span> <span class="n">datetime</span><span class="o">.</span><span class="n">datetime</span><span class="p">(</span><span class="mi">2008</span><span class="p">,</span> <span class="mo">01</span><span class="p">,</span> <span class="mo">01</span><span class="p">)</span></div><div class='line' id='LC85'><br/></div><div class='line' id='LC86'><span class="n">symbol_dict</span> <span class="o">=</span> <span class="p">{</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;TOT&#39;</span><span class="p">:</span> <span class="s">&#39;Total&#39;</span><span class="p">,</span></div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;XOM&#39;</span><span class="p">:</span> <span class="s">&#39;Exxon&#39;</span><span class="p">,</span></div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CVX&#39;</span><span class="p">:</span> <span class="s">&#39;Chevron&#39;</span><span class="p">,</span></div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;COP&#39;</span><span class="p">:</span> <span class="s">&#39;ConocoPhillips&#39;</span><span class="p">,</span></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;VLO&#39;</span><span class="p">:</span> <span class="s">&#39;Valero Energy&#39;</span><span class="p">,</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;MSFT&#39;</span><span class="p">:</span> <span class="s">&#39;Microsoft&#39;</span><span class="p">,</span></div><div class='line' id='LC93'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;IBM&#39;</span><span class="p">:</span> <span class="s">&#39;IBM&#39;</span><span class="p">,</span></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;TWX&#39;</span><span class="p">:</span> <span class="s">&#39;Time Warner&#39;</span><span class="p">,</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CMCSA&#39;</span><span class="p">:</span> <span class="s">&#39;Comcast&#39;</span><span class="p">,</span></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CVC&#39;</span><span class="p">:</span> <span class="s">&#39;Cablevision&#39;</span><span class="p">,</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;YHOO&#39;</span><span class="p">:</span> <span class="s">&#39;Yahoo&#39;</span><span class="p">,</span></div><div class='line' id='LC98'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;DELL&#39;</span><span class="p">:</span> <span class="s">&#39;Dell&#39;</span><span class="p">,</span></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;HPQ&#39;</span><span class="p">:</span> <span class="s">&#39;HP&#39;</span><span class="p">,</span></div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;AMZN&#39;</span><span class="p">:</span> <span class="s">&#39;Amazon&#39;</span><span class="p">,</span></div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;TM&#39;</span><span class="p">:</span> <span class="s">&#39;Toyota&#39;</span><span class="p">,</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CAJ&#39;</span><span class="p">:</span> <span class="s">&#39;Canon&#39;</span><span class="p">,</span></div><div class='line' id='LC103'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;MTU&#39;</span><span class="p">:</span> <span class="s">&#39;Mitsubishi&#39;</span><span class="p">,</span></div><div class='line' id='LC104'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;SNE&#39;</span><span class="p">:</span> <span class="s">&#39;Sony&#39;</span><span class="p">,</span></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;F&#39;</span><span class="p">:</span> <span class="s">&#39;Ford&#39;</span><span class="p">,</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;HMC&#39;</span><span class="p">:</span> <span class="s">&#39;Honda&#39;</span><span class="p">,</span></div><div class='line' id='LC107'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;NAV&#39;</span><span class="p">:</span> <span class="s">&#39;Navistar&#39;</span><span class="p">,</span></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;NOC&#39;</span><span class="p">:</span> <span class="s">&#39;Northrop Grumman&#39;</span><span class="p">,</span></div><div class='line' id='LC109'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;BA&#39;</span><span class="p">:</span> <span class="s">&#39;Boeing&#39;</span><span class="p">,</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;KO&#39;</span><span class="p">:</span> <span class="s">&#39;Coca Cola&#39;</span><span class="p">,</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;MMM&#39;</span><span class="p">:</span> <span class="s">&#39;3M&#39;</span><span class="p">,</span></div><div class='line' id='LC112'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;MCD&#39;</span><span class="p">:</span> <span class="s">&#39;Mc Donalds&#39;</span><span class="p">,</span></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;PEP&#39;</span><span class="p">:</span> <span class="s">&#39;Pepsi&#39;</span><span class="p">,</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;KFT&#39;</span><span class="p">:</span> <span class="s">&#39;Kraft Foods&#39;</span><span class="p">,</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;K&#39;</span><span class="p">:</span> <span class="s">&#39;Kellogg&#39;</span><span class="p">,</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;UN&#39;</span><span class="p">:</span> <span class="s">&#39;Unilever&#39;</span><span class="p">,</span></div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;MAR&#39;</span><span class="p">:</span> <span class="s">&#39;Marriott&#39;</span><span class="p">,</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;PG&#39;</span><span class="p">:</span> <span class="s">&#39;Procter Gamble&#39;</span><span class="p">,</span></div><div class='line' id='LC119'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CL&#39;</span><span class="p">:</span> <span class="s">&#39;Colgate-Palmolive&#39;</span><span class="p">,</span></div><div class='line' id='LC120'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;NWS&#39;</span><span class="p">:</span> <span class="s">&#39;News Corp&#39;</span><span class="p">,</span></div><div class='line' id='LC121'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;GE&#39;</span><span class="p">:</span> <span class="s">&#39;General Electrics&#39;</span><span class="p">,</span></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;WFC&#39;</span><span class="p">:</span> <span class="s">&#39;Wells Fargo&#39;</span><span class="p">,</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;JPM&#39;</span><span class="p">:</span> <span class="s">&#39;JPMorgan Chase&#39;</span><span class="p">,</span></div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;AIG&#39;</span><span class="p">:</span> <span class="s">&#39;AIG&#39;</span><span class="p">,</span></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;AXP&#39;</span><span class="p">:</span> <span class="s">&#39;American express&#39;</span><span class="p">,</span></div><div class='line' id='LC126'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;BAC&#39;</span><span class="p">:</span> <span class="s">&#39;Bank of America&#39;</span><span class="p">,</span></div><div class='line' id='LC127'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;GS&#39;</span><span class="p">:</span> <span class="s">&#39;Goldman Sachs&#39;</span><span class="p">,</span></div><div class='line' id='LC128'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;AAPL&#39;</span><span class="p">:</span> <span class="s">&#39;Apple&#39;</span><span class="p">,</span></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;SAP&#39;</span><span class="p">:</span> <span class="s">&#39;SAP&#39;</span><span class="p">,</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CSCO&#39;</span><span class="p">:</span> <span class="s">&#39;Cisco&#39;</span><span class="p">,</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;TXN&#39;</span><span class="p">:</span> <span class="s">&#39;Texas instruments&#39;</span><span class="p">,</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;XRX&#39;</span><span class="p">:</span> <span class="s">&#39;Xerox&#39;</span><span class="p">,</span></div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;LMT&#39;</span><span class="p">:</span> <span class="s">&#39;Lookheed Martin&#39;</span><span class="p">,</span></div><div class='line' id='LC134'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;WMT&#39;</span><span class="p">:</span> <span class="s">&#39;Wal-Mart&#39;</span><span class="p">,</span></div><div class='line' id='LC135'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;WAG&#39;</span><span class="p">:</span> <span class="s">&#39;Walgreen&#39;</span><span class="p">,</span></div><div class='line' id='LC136'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;HD&#39;</span><span class="p">:</span> <span class="s">&#39;Home Depot&#39;</span><span class="p">,</span></div><div class='line' id='LC137'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;GSK&#39;</span><span class="p">:</span> <span class="s">&#39;GlaxoSmithKline&#39;</span><span class="p">,</span></div><div class='line' id='LC138'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;PFE&#39;</span><span class="p">:</span> <span class="s">&#39;Pfizer&#39;</span><span class="p">,</span></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;SNY&#39;</span><span class="p">:</span> <span class="s">&#39;Sanofi-Aventis&#39;</span><span class="p">,</span></div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;NVS&#39;</span><span class="p">:</span> <span class="s">&#39;Novartis&#39;</span><span class="p">,</span></div><div class='line' id='LC141'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;KMB&#39;</span><span class="p">:</span> <span class="s">&#39;Kimberly-Clark&#39;</span><span class="p">,</span></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;R&#39;</span><span class="p">:</span> <span class="s">&#39;Ryder&#39;</span><span class="p">,</span></div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;GD&#39;</span><span class="p">:</span> <span class="s">&#39;General Dynamics&#39;</span><span class="p">,</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;RTN&#39;</span><span class="p">:</span> <span class="s">&#39;Raytheon&#39;</span><span class="p">,</span></div><div class='line' id='LC145'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CVS&#39;</span><span class="p">:</span> <span class="s">&#39;CVS&#39;</span><span class="p">,</span></div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;CAT&#39;</span><span class="p">:</span> <span class="s">&#39;Caterpillar&#39;</span><span class="p">,</span></div><div class='line' id='LC147'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="s">&#39;DD&#39;</span><span class="p">:</span> <span class="s">&#39;DuPont de Nemours&#39;</span><span class="p">}</span></div><div class='line' id='LC148'><br/></div><div class='line' id='LC149'><span class="n">symbols</span><span class="p">,</span> <span class="n">names</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">(</span><span class="n">symbol_dict</span><span class="o">.</span><span class="n">items</span><span class="p">())</span><span class="o">.</span><span class="n">T</span></div><div class='line' id='LC150'><br/></div><div class='line' id='LC151'><span class="n">quotes</span> <span class="o">=</span> <span class="p">[</span><span class="n">finance</span><span class="o">.</span><span class="n">quotes_historical_yahoo</span><span class="p">(</span><span class="n">symbol</span><span class="p">,</span> <span class="n">d1</span><span class="p">,</span> <span class="n">d2</span><span class="p">,</span> <span class="n">asobject</span><span class="o">=</span><span class="bp">True</span><span class="p">)</span></div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="n">symbol</span> <span class="ow">in</span> <span class="n">symbols</span><span class="p">]</span></div><div class='line' id='LC153'><br/></div><div class='line' id='LC154'><span class="nb">open</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="n">q</span><span class="o">.</span><span class="n">open</span> <span class="k">for</span> <span class="n">q</span> <span class="ow">in</span> <span class="n">quotes</span><span class="p">])</span><span class="o">.</span><span class="n">astype</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">float</span><span class="p">)</span></div><div class='line' id='LC155'><span class="n">close</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">array</span><span class="p">([</span><span class="n">q</span><span class="o">.</span><span class="n">close</span> <span class="k">for</span> <span class="n">q</span> <span class="ow">in</span> <span class="n">quotes</span><span class="p">])</span><span class="o">.</span><span class="n">astype</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">float</span><span class="p">)</span></div><div class='line' id='LC156'><br/></div><div class='line' id='LC157'><span class="c"># The daily variations of the quotes are what carry most information</span></div><div class='line' id='LC158'><span class="n">variation</span> <span class="o">=</span> <span class="n">close</span> <span class="o">-</span> <span class="nb">open</span></div><div class='line' id='LC159'><br/></div><div class='line' id='LC160'><span class="c">###############################################################################</span></div><div class='line' id='LC161'><span class="c"># Learn a graphical structure from the correlations</span></div><div class='line' id='LC162'><span class="n">edge_model</span> <span class="o">=</span> <span class="n">covariance</span><span class="o">.</span><span class="n">GraphLassoCV</span><span class="p">()</span></div><div class='line' id='LC163'><br/></div><div class='line' id='LC164'><span class="c"># standardize the time series: using correlations rather than covariance</span></div><div class='line' id='LC165'><span class="c"># is more efficient for structure recovery</span></div><div class='line' id='LC166'><span class="n">X</span> <span class="o">=</span> <span class="n">variation</span><span class="o">.</span><span class="n">copy</span><span class="p">()</span><span class="o">.</span><span class="n">T</span></div><div class='line' id='LC167'><span class="n">X</span> <span class="o">/=</span> <span class="n">X</span><span class="o">.</span><span class="n">std</span><span class="p">(</span><span class="n">axis</span><span class="o">=</span><span class="mi">0</span><span class="p">)</span></div><div class='line' id='LC168'><span class="n">edge_model</span><span class="o">.</span><span class="n">fit</span><span class="p">(</span><span class="n">X</span><span class="p">)</span></div><div class='line' id='LC169'><br/></div><div class='line' id='LC170'><span class="c">###############################################################################</span></div><div class='line' id='LC171'><span class="c"># Cluster using affinity propagation</span></div><div class='line' id='LC172'><br/></div><div class='line' id='LC173'><span class="n">_</span><span class="p">,</span> <span class="n">labels</span> <span class="o">=</span> <span class="n">cluster</span><span class="o">.</span><span class="n">affinity_propagation</span><span class="p">(</span><span class="n">edge_model</span><span class="o">.</span><span class="n">covariance_</span><span class="p">)</span></div><div class='line' id='LC174'><span class="n">n_labels</span> <span class="o">=</span> <span class="n">labels</span><span class="o">.</span><span class="n">max</span><span class="p">()</span></div><div class='line' id='LC175'><br/></div><div class='line' id='LC176'><span class="k">for</span> <span class="n">i</span> <span class="ow">in</span> <span class="nb">range</span><span class="p">(</span><span class="n">n_labels</span> <span class="o">+</span> <span class="mi">1</span><span class="p">):</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">print</span><span class="p">(</span><span class="s">&#39;Cluster </span><span class="si">%i</span><span class="s">: </span><span class="si">%s</span><span class="s">&#39;</span> <span class="o">%</span> <span class="p">((</span><span class="n">i</span> <span class="o">+</span> <span class="mi">1</span><span class="p">),</span> <span class="s">&#39;, &#39;</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">names</span><span class="p">[</span><span class="n">labels</span> <span class="o">==</span> <span class="n">i</span><span class="p">])))</span></div><div class='line' id='LC178'><br/></div><div class='line' id='LC179'><span class="c">###############################################################################</span></div><div class='line' id='LC180'><span class="c"># Find a low-dimension embedding for visualization: find the best position of</span></div><div class='line' id='LC181'><span class="c"># the nodes (the stocks) on a 2D plane</span></div><div class='line' id='LC182'><br/></div><div class='line' id='LC183'><span class="c"># We use a dense eigen_solver to achieve reproducibility (arpack is</span></div><div class='line' id='LC184'><span class="c"># initiated with random vectors that we don&#39;t control). In addition, we</span></div><div class='line' id='LC185'><span class="c"># use a large number of neighbors to capture the large-scale structure.</span></div><div class='line' id='LC186'><span class="n">node_position_model</span> <span class="o">=</span> <span class="n">manifold</span><span class="o">.</span><span class="n">LocallyLinearEmbedding</span><span class="p">(</span></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">n_components</span><span class="o">=</span><span class="mi">2</span><span class="p">,</span> <span class="n">eigen_solver</span><span class="o">=</span><span class="s">&#39;dense&#39;</span><span class="p">,</span> <span class="n">n_neighbors</span><span class="o">=</span><span class="mi">6</span><span class="p">)</span></div><div class='line' id='LC188'><br/></div><div class='line' id='LC189'><span class="n">embedding</span> <span class="o">=</span> <span class="n">node_position_model</span><span class="o">.</span><span class="n">fit_transform</span><span class="p">(</span><span class="n">X</span><span class="o">.</span><span class="n">T</span><span class="p">)</span><span class="o">.</span><span class="n">T</span></div><div class='line' id='LC190'><br/></div><div class='line' id='LC191'><span class="c">###############################################################################</span></div><div class='line' id='LC192'><span class="c"># Visualization</span></div><div class='line' id='LC193'><span class="n">pl</span><span class="o">.</span><span class="n">figure</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="n">facecolor</span><span class="o">=</span><span class="s">&#39;w&#39;</span><span class="p">,</span> <span class="n">figsize</span><span class="o">=</span><span class="p">(</span><span class="mi">10</span><span class="p">,</span> <span class="mi">8</span><span class="p">))</span></div><div class='line' id='LC194'><span class="n">pl</span><span class="o">.</span><span class="n">clf</span><span class="p">()</span></div><div class='line' id='LC195'><span class="n">ax</span> <span class="o">=</span> <span class="n">pl</span><span class="o">.</span><span class="n">axes</span><span class="p">([</span><span class="mf">0.</span><span class="p">,</span> <span class="mf">0.</span><span class="p">,</span> <span class="mf">1.</span><span class="p">,</span> <span class="mf">1.</span><span class="p">])</span></div><div class='line' id='LC196'><span class="n">pl</span><span class="o">.</span><span class="n">axis</span><span class="p">(</span><span class="s">&#39;off&#39;</span><span class="p">)</span></div><div class='line' id='LC197'><br/></div><div class='line' id='LC198'><span class="c"># Display a graph of the partial correlations</span></div><div class='line' id='LC199'><span class="n">partial_correlations</span> <span class="o">=</span> <span class="n">edge_model</span><span class="o">.</span><span class="n">precision_</span><span class="o">.</span><span class="n">copy</span><span class="p">()</span></div><div class='line' id='LC200'><span class="n">d</span> <span class="o">=</span> <span class="mi">1</span> <span class="o">/</span> <span class="n">np</span><span class="o">.</span><span class="n">sqrt</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">diag</span><span class="p">(</span><span class="n">partial_correlations</span><span class="p">))</span></div><div class='line' id='LC201'><span class="n">partial_correlations</span> <span class="o">*=</span> <span class="n">d</span></div><div class='line' id='LC202'><span class="n">partial_correlations</span> <span class="o">*=</span> <span class="n">d</span><span class="p">[:,</span> <span class="n">np</span><span class="o">.</span><span class="n">newaxis</span><span class="p">]</span></div><div class='line' id='LC203'><span class="n">non_zero</span> <span class="o">=</span> <span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">abs</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">triu</span><span class="p">(</span><span class="n">partial_correlations</span><span class="p">,</span> <span class="n">k</span><span class="o">=</span><span class="mi">1</span><span class="p">))</span> <span class="o">&gt;</span> <span class="mf">0.02</span><span class="p">)</span></div><div class='line' id='LC204'><br/></div><div class='line' id='LC205'><span class="c"># Plot the nodes using the coordinates of our embedding</span></div><div class='line' id='LC206'><span class="n">pl</span><span class="o">.</span><span class="n">scatter</span><span class="p">(</span><span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">],</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">],</span> <span class="n">s</span><span class="o">=</span><span class="mi">100</span> <span class="o">*</span> <span class="n">d</span> <span class="o">**</span> <span class="mi">2</span><span class="p">,</span> <span class="n">c</span><span class="o">=</span><span class="n">labels</span><span class="p">,</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">cmap</span><span class="o">=</span><span class="n">pl</span><span class="o">.</span><span class="n">cm</span><span class="o">.</span><span class="n">spectral</span><span class="p">)</span></div><div class='line' id='LC208'><br/></div><div class='line' id='LC209'><span class="c"># Plot the edges</span></div><div class='line' id='LC210'><span class="n">start_idx</span><span class="p">,</span> <span class="n">end_idx</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">where</span><span class="p">(</span><span class="n">non_zero</span><span class="p">)</span></div><div class='line' id='LC211'><span class="c">#a sequence of (*line0*, *line1*, *line2*), where::</span></div><div class='line' id='LC212'><span class="c">#            linen = (x0, y0), (x1, y1), ... (xm, ym)</span></div><div class='line' id='LC213'><span class="n">segments</span> <span class="o">=</span> <span class="p">[[</span><span class="n">embedding</span><span class="p">[:,</span> <span class="n">start</span><span class="p">],</span> <span class="n">embedding</span><span class="p">[:,</span> <span class="n">stop</span><span class="p">]]</span></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="n">start</span><span class="p">,</span> <span class="n">stop</span> <span class="ow">in</span> <span class="nb">zip</span><span class="p">(</span><span class="n">start_idx</span><span class="p">,</span> <span class="n">end_idx</span><span class="p">)]</span></div><div class='line' id='LC215'><span class="n">values</span> <span class="o">=</span> <span class="n">np</span><span class="o">.</span><span class="n">abs</span><span class="p">(</span><span class="n">partial_correlations</span><span class="p">[</span><span class="n">non_zero</span><span class="p">])</span></div><div class='line' id='LC216'><span class="n">lc</span> <span class="o">=</span> <span class="n">LineCollection</span><span class="p">(</span><span class="n">segments</span><span class="p">,</span></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">zorder</span><span class="o">=</span><span class="mi">0</span><span class="p">,</span> <span class="n">cmap</span><span class="o">=</span><span class="n">pl</span><span class="o">.</span><span class="n">cm</span><span class="o">.</span><span class="n">hot_r</span><span class="p">,</span></div><div class='line' id='LC218'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">norm</span><span class="o">=</span><span class="n">pl</span><span class="o">.</span><span class="n">Normalize</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="o">.</span><span class="mi">7</span> <span class="o">*</span> <span class="n">values</span><span class="o">.</span><span class="n">max</span><span class="p">()))</span></div><div class='line' id='LC219'><span class="n">lc</span><span class="o">.</span><span class="n">set_array</span><span class="p">(</span><span class="n">values</span><span class="p">)</span></div><div class='line' id='LC220'><span class="n">lc</span><span class="o">.</span><span class="n">set_linewidths</span><span class="p">(</span><span class="mi">15</span> <span class="o">*</span> <span class="n">values</span><span class="p">)</span></div><div class='line' id='LC221'><span class="n">ax</span><span class="o">.</span><span class="n">add_collection</span><span class="p">(</span><span class="n">lc</span><span class="p">)</span></div><div class='line' id='LC222'><br/></div><div class='line' id='LC223'><span class="c"># Add a label to each node. The challenge here is that we want to</span></div><div class='line' id='LC224'><span class="c"># position the labels to avoid overlap with other labels</span></div><div class='line' id='LC225'><span class="k">for</span> <span class="n">index</span><span class="p">,</span> <span class="p">(</span><span class="n">name</span><span class="p">,</span> <span class="n">label</span><span class="p">,</span> <span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">y</span><span class="p">))</span> <span class="ow">in</span> <span class="nb">enumerate</span><span class="p">(</span></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="nb">zip</span><span class="p">(</span><span class="n">names</span><span class="p">,</span> <span class="n">labels</span><span class="p">,</span> <span class="n">embedding</span><span class="o">.</span><span class="n">T</span><span class="p">)):</span></div><div class='line' id='LC227'><br/></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dx</span> <span class="o">=</span> <span class="n">x</span> <span class="o">-</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dx</span><span class="p">[</span><span class="n">index</span><span class="p">]</span> <span class="o">=</span> <span class="mi">1</span></div><div class='line' id='LC230'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dy</span> <span class="o">=</span> <span class="n">y</span> <span class="o">-</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span></div><div class='line' id='LC231'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">dy</span><span class="p">[</span><span class="n">index</span><span class="p">]</span> <span class="o">=</span> <span class="mi">1</span></div><div class='line' id='LC232'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">this_dx</span> <span class="o">=</span> <span class="n">dx</span><span class="p">[</span><span class="n">np</span><span class="o">.</span><span class="n">argmin</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">abs</span><span class="p">(</span><span class="n">dy</span><span class="p">))]</span></div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">this_dy</span> <span class="o">=</span> <span class="n">dy</span><span class="p">[</span><span class="n">np</span><span class="o">.</span><span class="n">argmin</span><span class="p">(</span><span class="n">np</span><span class="o">.</span><span class="n">abs</span><span class="p">(</span><span class="n">dx</span><span class="p">))]</span></div><div class='line' id='LC234'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">this_dx</span> <span class="o">&gt;</span> <span class="mi">0</span><span class="p">:</span></div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">horizontalalignment</span> <span class="o">=</span> <span class="s">&#39;left&#39;</span></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">x</span> <span class="o">=</span> <span class="n">x</span> <span class="o">+</span> <span class="o">.</span><span class="mo">002</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC238'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">horizontalalignment</span> <span class="o">=</span> <span class="s">&#39;right&#39;</span></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">x</span> <span class="o">=</span> <span class="n">x</span> <span class="o">-</span> <span class="o">.</span><span class="mo">002</span></div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="n">this_dy</span> <span class="o">&gt;</span> <span class="mi">0</span><span class="p">:</span></div><div class='line' id='LC241'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">verticalalignment</span> <span class="o">=</span> <span class="s">&#39;bottom&#39;</span></div><div class='line' id='LC242'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">y</span> <span class="o">=</span> <span class="n">y</span> <span class="o">+</span> <span class="o">.</span><span class="mo">002</span></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span><span class="p">:</span></div><div class='line' id='LC244'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">verticalalignment</span> <span class="o">=</span> <span class="s">&#39;top&#39;</span></div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">y</span> <span class="o">=</span> <span class="n">y</span> <span class="o">-</span> <span class="o">.</span><span class="mo">002</span></div><div class='line' id='LC246'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pl</span><span class="o">.</span><span class="n">text</span><span class="p">(</span><span class="n">x</span><span class="p">,</span> <span class="n">y</span><span class="p">,</span> <span class="n">name</span><span class="p">,</span> <span class="n">size</span><span class="o">=</span><span class="mi">10</span><span class="p">,</span></div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">horizontalalignment</span><span class="o">=</span><span class="n">horizontalalignment</span><span class="p">,</span></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">verticalalignment</span><span class="o">=</span><span class="n">verticalalignment</span><span class="p">,</span></div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">bbox</span><span class="o">=</span><span class="nb">dict</span><span class="p">(</span><span class="n">facecolor</span><span class="o">=</span><span class="s">&#39;w&#39;</span><span class="p">,</span></div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">edgecolor</span><span class="o">=</span><span class="n">pl</span><span class="o">.</span><span class="n">cm</span><span class="o">.</span><span class="n">spectral</span><span class="p">(</span><span class="n">label</span> <span class="o">/</span> <span class="nb">float</span><span class="p">(</span><span class="n">n_labels</span><span class="p">)),</span></div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">alpha</span><span class="o">=.</span><span class="mi">6</span><span class="p">))</span></div><div class='line' id='LC252'><br/></div><div class='line' id='LC253'><span class="n">pl</span><span class="o">.</span><span class="n">xlim</span><span class="p">(</span><span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span><span class="o">.</span><span class="n">min</span><span class="p">()</span> <span class="o">-</span> <span class="o">.</span><span class="mi">15</span> <span class="o">*</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span><span class="o">.</span><span class="n">ptp</span><span class="p">(),</span></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span><span class="o">.</span><span class="n">max</span><span class="p">()</span> <span class="o">+</span> <span class="o">.</span><span class="mi">10</span> <span class="o">*</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span><span class="o">.</span><span class="n">ptp</span><span class="p">(),)</span></div><div class='line' id='LC255'><span class="n">pl</span><span class="o">.</span><span class="n">ylim</span><span class="p">(</span><span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span><span class="o">.</span><span class="n">min</span><span class="p">()</span> <span class="o">-</span> <span class="o">.</span><span class="mo">03</span> <span class="o">*</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span><span class="o">.</span><span class="n">ptp</span><span class="p">(),</span></div><div class='line' id='LC256'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span><span class="o">.</span><span class="n">max</span><span class="p">()</span> <span class="o">+</span> <span class="o">.</span><span class="mo">03</span> <span class="o">*</span> <span class="n">embedding</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span><span class="o">.</span><span class="n">ptp</span><span class="p">())</span></div><div class='line' id='LC257'><br/></div><div class='line' id='LC258'><span class="n">pl</span><span class="o">.</span><span class="n">show</span><span class="p">()</span></div></pre></div>
          </td>
        </tr>
      </table>
  </div>

          </div>
        </div>

        <a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
        <div id="jump-to-line" style="display:none">
          <form accept-charset="UTF-8" class="js-jump-to-line-form">
            <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;">
            <button type="submit" class="button">Go</button>
          </form>
        </div>

</div>

<div id="js-frame-loading-template" class="frame frame-loading large-loading-area" style="display:none;">
  <img class="js-frame-loading-spinner" src="https://a248.e.akamai.net/assets.github.com/images/spinners/octocat-spinner-128.gif" height="64" width="64">
</div>


          </div>
        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div>
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">Developer</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>
    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.13710s from fe3.rs.github.com">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/scikit-learn/scikit-learn/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

    
    <span id='server_response_time' data-time='0.13745' data-host='fe3'></span>
    
  </body>
</html>

