



<!DOCTYPE html>
<html lang="en" class="">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    

    <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/frameworks-9ad495a7b5d4473ee5031bc7b3e2d60090320bca47a10a6fa472132f678b6160.css" media="all" rel="stylesheet" />
    <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/github-3ca4d5a0760c7ca10f98748867f95c64b034bd809a90302ab1caf3adf1b7845c.css" media="all" rel="stylesheet" />
    
    
    <link crossorigin="anonymous" href="https://assets-cdn.github.com/assets/site-bfc167aa91c105811fa2848e54f5ae70b6c18ac4a42adf45cfacf9fbb078a71f.css" media="all" rel="stylesheet" />
    

    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta http-equiv="Content-Language" content="en">
    <meta name="viewport" content="width=device-width">
    
    <title>suncasa/main.py at master · sjyu1988/suncasa · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
    <link rel="apple-touch-icon" href="/apple-touch-icon.png">
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-57x57.png">
    <link rel="apple-touch-icon" sizes="60x60" href="/apple-touch-icon-60x60.png">
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-72x72.png">
    <link rel="apple-touch-icon" sizes="76x76" href="/apple-touch-icon-76x76.png">
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114x114.png">
    <link rel="apple-touch-icon" sizes="120x120" href="/apple-touch-icon-120x120.png">
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144x144.png">
    <link rel="apple-touch-icon" sizes="152x152" href="/apple-touch-icon-152x152.png">
    <link rel="apple-touch-icon" sizes="180x180" href="/apple-touch-icon-180x180.png">
    <meta property="fb:app_id" content="1401488693436528">

      <meta content="https://avatars1.githubusercontent.com/u/21953547?v=3&amp;s=400" name="twitter:image:src" /><meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="sjyu1988/suncasa" name="twitter:title" /><meta content="suncasa - CASA scripts and tasks to analysis and visualize solar radio spectroscopic imaging data" name="twitter:description" />
      <meta content="https://avatars1.githubusercontent.com/u/21953547?v=3&amp;s=400" property="og:image" /><meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="sjyu1988/suncasa" property="og:title" /><meta content="https://github.com/sjyu1988/suncasa" property="og:url" /><meta content="suncasa - CASA scripts and tasks to analysis and visualize solar radio spectroscopic imaging data" property="og:description" />
      <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">
    <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">
    <link rel="assets" href="https://assets-cdn.github.com/">
    
    <meta name="pjax-timeout" content="1000">
    
    <meta name="request-id" content="80EBD36B:4889:6A8740C:581A0E14" data-pjax-transient>

    <meta name="msapplication-TileImage" content="/windows-tile.png">
    <meta name="msapplication-TileColor" content="#ffffff">
    <meta name="selected-link" value="repo_source" data-pjax-transient>

    <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
<meta name="google-site-verification" content="ZzhVyEFwb7w3e0-uOTltm8Jsck2F5StVihD0exw2fsA">
    <meta name="google-analytics" content="UA-3769691-2">

<meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" /><meta content="80EBD36B:4889:6A8740C:581A0E14" name="octolytics-dimension-request_id" />
<meta content="/&lt;user-name&gt;/&lt;repo-name&gt;/blob/show" data-pjax-transient="true" name="analytics-location" />



  <meta class="js-ga-set" name="dimension1" content="Logged Out">



        <meta name="hostname" content="github.com">
    <meta name="user-login" content="">

        <meta name="expected-hostname" content="github.com">
      <meta name="js-proxy-site-detection-payload" content="OGQyMmE1M2QwMzU3OGMyZThlYTVkOWFmZDgzMjE5ZGY5ZWQ4NDJmNmMzODkwZDQ0ZjJhYWIzZWRiN2U3ODBkNHx7InJlbW90ZV9hZGRyZXNzIjoiMTI4LjIzNS4yMTEuMTA3IiwicmVxdWVzdF9pZCI6IjgwRUJEMzZCOjQ4ODk6NkE4NzQwQzo1ODFBMEUxNCIsInRpbWVzdGFtcCI6MTQ3ODEwMjU0OCwiaG9zdCI6ImdpdGh1Yi5jb20ifQ==">


      <link rel="mask-icon" href="https://assets-cdn.github.com/pinned-octocat.svg" color="#4078c0">
      <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico">

    <meta name="html-safe-nonce" content="3ddfb8ec5dd413a62ca68d66d9b67fc7454ffb94">
    <meta content="5e1a72eb7281d1c5e014e698969bde666397d86f" name="form-nonce" />

    <meta http-equiv="x-pjax-version" content="f983b7972fc2383d0fabea2ac9a3088f">
    

      
  <meta name="description" content="suncasa - CASA scripts and tasks to analysis and visualize solar radio spectroscopic imaging data">
  <meta name="go-import" content="github.com/sjyu1988/suncasa git https://github.com/sjyu1988/suncasa.git">

  <meta content="21953547" name="octolytics-dimension-user_id" /><meta content="sjyu1988" name="octolytics-dimension-user_login" /><meta content="72456367" name="octolytics-dimension-repository_id" /><meta content="sjyu1988/suncasa" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="72017919" name="octolytics-dimension-repository_parent_id" /><meta content="binchensun/suncasa" name="octolytics-dimension-repository_parent_nwo" /><meta content="72017919" name="octolytics-dimension-repository_network_root_id" /><meta content="binchensun/suncasa" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/sjyu1988/suncasa/commits/master.atom" rel="alternate" title="Recent Commits to suncasa:master" type="application/atom+xml">


      <link rel="canonical" href="https://github.com/sjyu1988/suncasa/blob/master/DataBrowser/QLook/main.py" data-pjax-transient>
  </head>


  <body class="logged-out  env-production  vis-public fork page-blob">
    <div id="js-pjax-loader-bar" class="pjax-loader-bar"><div class="progress"></div></div>
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>

    
    
    



          <header class="site-header js-details-container" role="banner">
  <div class="container-responsive">
    <a class="header-logo-invertocat" href="https://github.com/" aria-label="Homepage" data-ga-click="(Logged out) Header, go to homepage, icon:logo-wordmark">
      <svg aria-hidden="true" class="octicon octicon-mark-github" height="32" version="1.1" viewBox="0 0 16 16" width="32"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path></svg>
    </a>

    <button class="btn-link float-right site-header-toggle js-details-target" type="button" aria-label="Toggle navigation">
      <svg aria-hidden="true" class="octicon octicon-three-bars" height="24" version="1.1" viewBox="0 0 12 16" width="18"><path d="M11.41 9H.59C0 9 0 8.59 0 8c0-.59 0-1 .59-1H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1h.01zm0-4H.59C0 5 0 4.59 0 4c0-.59 0-1 .59-1H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1h.01zM.59 11H11.4c.59 0 .59.41.59 1 0 .59 0 1-.59 1H.59C0 13 0 12.59 0 12c0-.59 0-1 .59-1z"></path></svg>
    </button>

    <div class="site-header-menu">
      <nav class="site-header-nav site-header-nav-main">
        <a href="/personal" class="js-selected-navigation-item nav-item nav-item-personal" data-ga-click="Header, click, Nav menu - item:personal" data-selected-links="/personal /personal">
          Personal
</a>        <a href="/open-source" class="js-selected-navigation-item nav-item nav-item-opensource" data-ga-click="Header, click, Nav menu - item:opensource" data-selected-links="/open-source /open-source">
          Open source
</a>        <a href="/business" class="js-selected-navigation-item nav-item nav-item-business" data-ga-click="Header, click, Nav menu - item:business" data-selected-links="/business /business/partners /business/features /business/customers /business">
          Business
</a>        <a href="/explore" class="js-selected-navigation-item nav-item nav-item-explore" data-ga-click="Header, click, Nav menu - item:explore" data-selected-links="/explore /trending /trending/developers /integrations /integrations/feature/code /integrations/feature/collaborate /integrations/feature/ship /showcases /explore">
          Explore
</a>      </nav>

      <div class="site-header-actions">
            <a class="btn btn-primary site-header-actions-btn" href="/join?source=header-repo" data-ga-click="(Logged out) Header, clicked Sign up, text:sign-up">Sign up</a>
          <a class="btn site-header-actions-btn mr-1" href="/login?return_to=%2Fsjyu1988%2Fsuncasa%2Fblob%2Fmaster%2FDataBrowser%2FQLook%2Fmain.py" data-ga-click="(Logged out) Header, clicked Sign in, text:sign-in">Sign in</a>
      </div>

        <nav class="site-header-nav site-header-nav-secondary mr-md-3">
          <a class="nav-item" href="/pricing">Pricing</a>
          <a class="nav-item" href="/blog">Blog</a>
          <a class="nav-item" href="https://help.github.com">Support</a>
          <a class="nav-item header-search-link" href="https://github.com/search">Search GitHub</a>
              <div class="header-search scoped-search site-scoped-search js-site-search" role="search">
  <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="/sjyu1988/suncasa/search" class="js-site-search-form" data-scoped-search-url="/sjyu1988/suncasa/search" data-unscoped-search-url="/search" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <label class="form-control header-search-wrapper js-chromeless-input-container">
      <div class="header-search-scope">This repository</div>
      <input type="text"
        class="form-control header-search-input js-site-search-focus js-site-search-field is-clearable"
        data-hotkey="s"
        name="q"
        placeholder="Search"
        aria-label="Search this repository"
        data-unscoped-placeholder="Search GitHub"
        data-scoped-placeholder="Search"
        autocapitalize="off">
    </label>
</form></div>

        </nav>
    </div>
  </div>
</header>



    <div id="start-of-content" class="accessibility-aid"></div>

      <div id="js-flash-container">
</div>


    <div role="main">
        <div itemscope itemtype="http://schema.org/SoftwareSourceCode">
    <div id="js-repo-pjax-container" data-pjax-container>
      
<div class="pagehead repohead instapaper_ignore readability-menu experiment-repo-nav">
  <div class="container repohead-details-container">

    

<ul class="pagehead-actions">

  <li>
      <a href="/login?return_to=%2Fsjyu1988%2Fsuncasa"
    class="btn btn-sm btn-with-count tooltipped tooltipped-n"
    aria-label="You must be signed in to watch a repository" rel="nofollow">
    <svg aria-hidden="true" class="octicon octicon-eye" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"></path></svg>
    Watch
  </a>
  <a class="social-count" href="/sjyu1988/suncasa/watchers"
     aria-label="0 users are watching this repository">
    0
  </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fsjyu1988%2Fsuncasa"
    class="btn btn-sm btn-with-count tooltipped tooltipped-n"
    aria-label="You must be signed in to star a repository" rel="nofollow">
    <svg aria-hidden="true" class="octicon octicon-star" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74z"></path></svg>
    Star
  </a>

    <a class="social-count js-social-count" href="/sjyu1988/suncasa/stargazers"
      aria-label="0 users starred this repository">
      0
    </a>

  </li>

  <li>
      <a href="/login?return_to=%2Fsjyu1988%2Fsuncasa"
        class="btn btn-sm btn-with-count tooltipped tooltipped-n"
        aria-label="You must be signed in to fork a repository" rel="nofollow">
        <svg aria-hidden="true" class="octicon octicon-repo-forked" height="16" version="1.1" viewBox="0 0 10 16" width="10"><path d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"></path></svg>
        Fork
      </a>

    <a href="/sjyu1988/suncasa/network" class="social-count"
       aria-label="2 users are forked this repository">
      2
    </a>
  </li>
</ul>

    <h1 class="public ">
  <svg aria-hidden="true" class="octicon octicon-repo-forked" height="16" version="1.1" viewBox="0 0 10 16" width="10"><path d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"></path></svg>
  <span class="author" itemprop="author"><a href="/sjyu1988" class="url fn" rel="author">sjyu1988</a></span><!--
--><span class="path-divider">/</span><!--
--><strong itemprop="name"><a href="/sjyu1988/suncasa" data-pjax="#js-repo-pjax-container">suncasa</a></strong>

    <span class="fork-flag">
      <span class="text">forked from <a href="/binchensun/suncasa">binchensun/suncasa</a></span>
    </span>
</h1>

  </div>
  <div class="container">
    
<nav class="reponav js-repo-nav js-sidenav-container-pjax"
     itemscope
     itemtype="http://schema.org/BreadcrumbList"
     role="navigation"
     data-pjax="#js-repo-pjax-container">

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/sjyu1988/suncasa" aria-selected="true" class="js-selected-navigation-item selected reponav-item" data-hotkey="g c" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /sjyu1988/suncasa" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-code" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path d="M9.5 3L8 4.5 11.5 8 8 11.5 9.5 13 14 8 9.5 3zm-5 0L0 8l4.5 5L6 11.5 2.5 8 6 4.5 4.5 3z"></path></svg>
      <span itemprop="name">Code</span>
      <meta itemprop="position" content="1">
</a>  </span>


  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a href="/sjyu1988/suncasa/pulls" class="js-selected-navigation-item reponav-item" data-hotkey="g p" data-selected-links="repo_pulls /sjyu1988/suncasa/pulls" itemprop="url">
      <svg aria-hidden="true" class="octicon octicon-git-pull-request" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path d="M11 11.28V5c-.03-.78-.34-1.47-.94-2.06C9.46 2.35 8.78 2.03 8 2H7V0L4 3l3 3V4h1c.27.02.48.11.69.31.21.2.3.42.31.69v6.28A1.993 1.993 0 0 0 10 15a1.993 1.993 0 0 0 1-3.72zm-1 2.92c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zM4 3c0-1.11-.89-2-2-2a1.993 1.993 0 0 0-1 3.72v6.56A1.993 1.993 0 0 0 2 15a1.993 1.993 0 0 0 1-3.72V4.72c.59-.34 1-.98 1-1.72zm-.8 10c0 .66-.55 1.2-1.2 1.2-.65 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"></path></svg>
      <span itemprop="name">Pull requests</span>
      <span class="counter">0</span>
      <meta itemprop="position" content="3">
</a>  </span>

  <a href="/sjyu1988/suncasa/projects" class="js-selected-navigation-item reponav-item" data-selected-links="repo_projects new_repo_project repo_project /sjyu1988/suncasa/projects">
    <svg class="octicon" aria-hidden="true" version="1.1" width="15" height="16" viewBox="0 0 15 16">
      <path d="M1 15h13V1H1v14zM15 1v14a1 1 0 0 1-1 1H1a1 1 0 0 1-1-1V1a1 1 0 0 1 1-1h13a1 1 0 0 1 1 1zm-4.41 11h1.82c.59 0 .59-.41.59-1V3c0-.59 0-1-.59-1h-1.82C10 2 10 2.41 10 3v8c0 .59 0 1 .59 1zm-4-2h1.82C9 10 9 9.59 9 9V3c0-.59 0-1-.59-1H6.59C6 2 6 2.41 6 3v6c0 .59 0 1 .59 1zM2 13V3c0-.59 0-1 .59-1h1.82C5 2 5 2.41 5 3v10c0 .59 0 1-.59 1H2.59C2 14 2 13.59 2 13z"></path>
    </svg>
    Projects
    <span class="counter">0</span>
</a>


  <a href="/sjyu1988/suncasa/pulse" class="js-selected-navigation-item reponav-item" data-selected-links="pulse /sjyu1988/suncasa/pulse">
    <svg aria-hidden="true" class="octicon octicon-pulse" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path d="M11.5 8L8.8 5.4 6.6 8.5 5.5 1.6 2.38 8H0v2h3.6l.9-1.8.9 5.4L9 8.5l1.6 1.5H14V8z"></path></svg>
    Pulse
</a>
  <a href="/sjyu1988/suncasa/graphs" class="js-selected-navigation-item reponav-item" data-selected-links="repo_graphs repo_contributors /sjyu1988/suncasa/graphs">
    <svg aria-hidden="true" class="octicon octicon-graph" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path d="M16 14v1H0V0h1v14h15zM5 13H3V8h2v5zm4 0H7V3h2v10zm4 0h-2V6h2v7z"></path></svg>
    Graphs
</a>

</nav>

  </div>
</div>

<div class="container new-discussion-timeline experiment-repo-nav">
  <div class="repository-content">

    

<a href="/sjyu1988/suncasa/blob/ba81fbbf533eb96d2c63828a5def8dd504a31293/DataBrowser/QLook/main.py" class="d-none js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:5f7174fc94ed20de8995c0fbf2b98df1 -->

<div class="file-navigation js-zeroclipboard-container">
  
<div class="select-menu branch-select-menu js-menu-container js-select-menu float-left">
  <button class="btn btn-sm select-menu-button js-menu-target css-truncate" data-hotkey="w"
    
    type="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <i>Branch:</i>
    <span class="js-select-button css-truncate-target">master</span>
  </button>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <svg aria-label="Close" class="octicon octicon-x js-menu-close" height="16" role="img" version="1.1" viewBox="0 0 12 16" width="12"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"></path></svg>
        <span class="select-menu-title">Switch branches/tags</span>
      </div>

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="form-control js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" data-filter-placeholder="Filter branches/tags" class="js-select-menu-tab" role="tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" data-filter-placeholder="Find a tag…" class="js-select-menu-tab" role="tab">Tags</a>
            </li>
          </ul>
        </div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches" role="menu">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <a class="select-menu-item js-navigation-item js-navigation-open selected"
               href="/sjyu1988/suncasa/blob/master/DataBrowser/QLook/main.py"
               data-name="master"
               data-skip-pjax="true"
               rel="nofollow">
              <svg aria-hidden="true" class="octicon octicon-check select-menu-item-icon" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path d="M12 5l-8 8-4-4 1.5-1.5L4 10l6.5-6.5z"></path></svg>
              <span class="select-menu-item-text css-truncate-target js-select-menu-filter-text">
                master
              </span>
            </a>
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div>

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div>

    </div>
  </div>
</div>

  <div class="BtnGroup float-right">
    <a href="/sjyu1988/suncasa/find/master"
          class="js-pjax-capture-input btn btn-sm BtnGroup-item"
          data-pjax
          data-hotkey="t">
      Find file
    </a>
    <button aria-label="Copy file path to clipboard" class="js-zeroclipboard btn btn-sm BtnGroup-item tooltipped tooltipped-s" data-copied-hint="Copied!" type="button">Copy path</button>
  </div>
  <div class="breadcrumb js-zeroclipboard-target">
    <span class="repo-root js-repo-root"><span class="js-path-segment"><a href="/sjyu1988/suncasa"><span>suncasa</span></a></span></span><span class="separator">/</span><span class="js-path-segment"><a href="/sjyu1988/suncasa/tree/master/DataBrowser"><span>DataBrowser</span></a></span><span class="separator">/</span><span class="js-path-segment"><a href="/sjyu1988/suncasa/tree/master/DataBrowser/QLook"><span>QLook</span></a></span><span class="separator">/</span><strong class="final-path">main.py</strong>
  </div>
</div>

<include-fragment class="commit-tease" src="/sjyu1988/suncasa/contributors/master/DataBrowser/QLook/main.py">
  <div>
    Fetching contributors&hellip;
  </div>

  <div class="commit-tease-contributors">
    <img alt="" class="loader-loading float-left" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" />
    <span class="loader-error">Cannot retrieve contributors at this time</span>
  </div>
</include-fragment>
<div class="file">
  <div class="file-header">
  <div class="file-actions">

    <div class="BtnGroup">
      <a href="/sjyu1988/suncasa/raw/master/DataBrowser/QLook/main.py" class="btn btn-sm BtnGroup-item" id="raw-url">Raw</a>
        <a href="/sjyu1988/suncasa/blame/master/DataBrowser/QLook/main.py" class="btn btn-sm js-update-url-with-hash BtnGroup-item">Blame</a>
      <a href="/sjyu1988/suncasa/commits/master/DataBrowser/QLook/main.py" class="btn btn-sm BtnGroup-item" rel="nofollow">History</a>
    </div>


        <button type="button" class="btn-octicon disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg aria-hidden="true" class="octicon octicon-pencil" height="16" version="1.1" viewBox="0 0 14 16" width="14"><path d="M0 12v3h3l8-8-3-3-8 8zm3 2H1v-2h1v1h1v1zm10.3-9.3L12 6 9 3l1.3-1.3a.996.996 0 0 1 1.41 0l1.59 1.59c.39.39.39 1.02 0 1.41z"></path></svg>
        </button>
        <button type="button" class="btn-octicon btn-octicon-danger disabled tooltipped tooltipped-nw"
          aria-label="You must be signed in to make or propose changes">
          <svg aria-hidden="true" class="octicon octicon-trashcan" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path d="M11 2H9c0-.55-.45-1-1-1H5c-.55 0-1 .45-1 1H2c-.55 0-1 .45-1 1v1c0 .55.45 1 1 1v9c0 .55.45 1 1 1h7c.55 0 1-.45 1-1V5c.55 0 1-.45 1-1V3c0-.55-.45-1-1-1zm-1 12H3V5h1v8h1V5h1v8h1V5h1v8h1V5h1v9zm1-10H2V3h9v1z"></path></svg>
        </button>
  </div>

  <div class="file-info">
      <span class="file-mode" title="File mode">executable file</span>
      <span class="file-info-divider"></span>
      426 lines (365 sloc)
      <span class="file-info-divider"></span>
    20.8 KB
  </div>
</div>

  

  <div itemprop="text" class="blob-wrapper data type-python">
      <table class="highlight tab-size js-file-line-container" data-tab-size="8">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code blob-code-inner js-file-line"><span class="pl-c"># The plot server must be running</span></td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code blob-code-inner js-file-line"><span class="pl-c"># Go to http://localhost:5006/bokeh to view this plot</span></td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> json</td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> os</td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> os.path</td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> time</td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code blob-code-inner js-file-line"><span class="pl-k">from</span> sys <span class="pl-k">import</span> platform</td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> matplotlib.cm <span class="pl-k">as</span> cm</td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> matplotlib.colors <span class="pl-k">as</span> colors</td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> numpy <span class="pl-k">as</span> np</td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> pandas <span class="pl-k">as</span> pd</td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> scipy.ndimage <span class="pl-k">as</span> sn</td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code blob-code-inner js-file-line"><span class="pl-k">from</span> bokeh.layouts <span class="pl-k">import</span> row, column, widgetbox</td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code blob-code-inner js-file-line"><span class="pl-k">from</span> bokeh.models <span class="pl-k">import</span> (ColumnDataSource, Button, TextInput, DataTable, TableColumn, BoxSelectTool, TapTool,</td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code blob-code-inner js-file-line">                          HoverTool, Spacer, Div)</td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code blob-code-inner js-file-line"><span class="pl-k">from</span> bokeh.models.widgets <span class="pl-k">import</span> Select</td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code blob-code-inner js-file-line"><span class="pl-k">from</span> bokeh.plotting <span class="pl-k">import</span> figure, curdoc</td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code blob-code-inner js-file-line"><span class="pl-k">import</span> jdutil</td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code blob-code-inner js-file-line">__author__ <span class="pl-k">=</span> [<span class="pl-s"><span class="pl-pds">&quot;</span>Sijie Yu<span class="pl-pds">&quot;</span></span>]</td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code blob-code-inner js-file-line">__email__ <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>sijie.yu@njit.edu<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span>Prepare the ports for FSview<span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code blob-code-inner js-file-line"><span class="pl-k">if</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>linux<span class="pl-pds">&quot;</span></span> <span class="pl-k">or</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>linux2<span class="pl-pds">&quot;</span></span>:</td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code blob-code-inner js-file-line">    <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>Runing QLook in Linux platform<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> ll <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">5100</span>, <span class="pl-c1">5100</span> <span class="pl-k">+</span> <span class="pl-c1">10</span>):</td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code blob-code-inner js-file-line">        os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>fuser -n tcp -k <span class="pl-c1">{}</span><span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code blob-code-inner js-file-line"><span class="pl-k">elif</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>darwin<span class="pl-pds">&quot;</span></span>:</td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code blob-code-inner js-file-line">    <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>Runing QLook in OS X platform<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> ll <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">5100</span>, <span class="pl-c1">5100</span> <span class="pl-k">+</span> <span class="pl-c1">10</span>):</td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code blob-code-inner js-file-line">        os.system(</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code blob-code-inner js-file-line">            <span class="pl-s"><span class="pl-pds">&#39;</span>port=($(lsof -i tcp:<span class="pl-c1">{}</span>|grep python2.7 |cut -f2 -d&quot; &quot;)); [[ -n &quot;$port&quot; ]] &amp;&amp; kill -9 $port<span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code blob-code-inner js-file-line">        os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>port=($(lsof -i tcp:<span class="pl-c1">{}</span>|grep Google |cut -f2 -d&quot; &quot;)); [[ -n &quot;$port&quot; ]] &amp;&amp; kill -9 $port<span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code blob-code-inner js-file-line"><span class="pl-k">elif</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>win32<span class="pl-pds">&quot;</span></span>:</td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code blob-code-inner js-file-line">    <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>Runing QLook in Windows platform<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span>load config file<span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code blob-code-inner js-file-line"><span class="pl-k">with</span> <span class="pl-c1">open</span>(<span class="pl-s"><span class="pl-pds">&#39;</span>../config.json<span class="pl-pds">&#39;</span></span>, <span class="pl-s"><span class="pl-pds">&#39;</span>r<span class="pl-pds">&#39;</span></span>) <span class="pl-k">as</span> fp:</td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code blob-code-inner js-file-line">    config_plot <span class="pl-k">=</span> json.load(fp)</td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code blob-code-inner js-file-line"><span class="pl-k">with</span> <span class="pl-c1">open</span>(<span class="pl-s"><span class="pl-pds">&#39;</span>config_EvtID.json<span class="pl-pds">&#39;</span></span>, <span class="pl-s"><span class="pl-pds">&#39;</span>r<span class="pl-pds">&#39;</span></span>) <span class="pl-k">as</span> fp:</td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code blob-code-inner js-file-line">    config_EvtID <span class="pl-k">=</span> json.load(fp)</td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code blob-code-inner js-file-line">EvtID_list <span class="pl-k">=</span> pd.read_json(<span class="pl-s"><span class="pl-pds">&#39;</span>EvtID_list.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code blob-code-inner js-file-line">EvtID_list <span class="pl-k">=</span> EvtID_list.sort_values([<span class="pl-s"><span class="pl-pds">&#39;</span>date<span class="pl-pds">&#39;</span></span>], <span class="pl-v">ascending</span><span class="pl-k">=</span>[<span class="pl-c1">True</span>])</td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span>define the colormaps<span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code blob-code-inner js-file-line">colormap_jet <span class="pl-k">=</span> cm.get_cmap(<span class="pl-s"><span class="pl-pds">&quot;</span>jet<span class="pl-pds">&quot;</span></span>)  <span class="pl-c"># choose any matplotlib colormap here</span></td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code blob-code-inner js-file-line">bokehpalette_jet <span class="pl-k">=</span> [colors.rgb2hex(m) <span class="pl-k">for</span> m <span class="pl-k">in</span> colormap_jet(np.arange(colormap_jet.<span class="pl-c1">N</span>))]</td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code blob-code-inner js-file-line"><span class="pl-s">-------------------------- panel 1 --------------------------</span></td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code blob-code-inner js-file-line">start_timestamp <span class="pl-k">=</span> time.time()</td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code blob-code-inner js-file-line">database_dir <span class="pl-k">=</span> config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>datadir<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>database<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code blob-code-inner js-file-line">database_dir <span class="pl-k">=</span> os.path.expandvars(database_dir)</td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code blob-code-inner js-file-line">event_id <span class="pl-k">=</span> config_EvtID[<span class="pl-s"><span class="pl-pds">&#39;</span>datadir<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>event_id<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code blob-code-inner js-file-line">specfile <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> config_EvtID[<span class="pl-s"><span class="pl-pds">&#39;</span>datadir<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>event_specfile<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">load_specdata</span>(<span class="pl-smi">specfile</span><span class="pl-k">=</span><span class="pl-c1">None</span>):</td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_specdata, tab1_spec, tab1_tim, tab1_freq, tab1_dtim, tab1_spec_plt, tab1_bl</td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code blob-code-inner js-file-line">    tab1_specdata <span class="pl-k">=</span> np.load(specfile)</td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code blob-code-inner js-file-line">    tab1_bl <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>bl<span class="pl-pds">&#39;</span></span>].item().split(<span class="pl-s"><span class="pl-pds">&#39;</span>;<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code blob-code-inner js-file-line">    tab1_pol <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&#39;</span>I<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code blob-code-inner js-file-line">    bl_index <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code blob-code-inner js-file-line">    tab1_spec <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>spec<span class="pl-pds">&#39;</span></span>][:, :, :, :]</td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code blob-code-inner js-file-line">    tab1_npol <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>npol<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code blob-code-inner js-file-line">    tab1_nbl <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>nbl<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code blob-code-inner js-file-line">    tab1_ntim <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>ntim<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code blob-code-inner js-file-line">    tab1_nfreq <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>nfreq<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code blob-code-inner js-file-line">    tab1_tim <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>tim<span class="pl-pds">&#39;</span></span>][:]</td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code blob-code-inner js-file-line">    tab1_freq <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>] <span class="pl-k">/</span> <span class="pl-c1">1e9</span></td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code blob-code-inner js-file-line">    tab1_spec_sz <span class="pl-k">=</span> tab1_spec.shape</td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code blob-code-inner js-file-line">    spec_sz2, spec_sz1 <span class="pl-k">=</span> <span class="pl-c1">10</span>, <span class="pl-c1">10</span></td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_spec_sz[<span class="pl-c1">3</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">1750</span>:</td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code blob-code-inner js-file-line">        spec_sz2 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">1</span>, <span class="pl-c1">10</span>) <span class="pl-k">if</span> i <span class="pl-k">/</span> <span class="pl-c1">10</span>. <span class="pl-k">*</span> tab1_spec_sz[<span class="pl-c1">3</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">1750</span>)</td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_spec_sz[<span class="pl-c1">2</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">250</span>:</td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code blob-code-inner js-file-line">        spec_sz1 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">1</span>, <span class="pl-c1">10</span>) <span class="pl-k">if</span> i <span class="pl-k">/</span> <span class="pl-c1">10</span>. <span class="pl-k">*</span> tab1_spec_sz[<span class="pl-c1">2</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">250</span>)</td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code blob-code-inner js-file-line">    tab1_spec <span class="pl-k">=</span> sn.interpolation.zoom(tab1_spec, [<span class="pl-c1">1</span>, <span class="pl-c1">1</span>, spec_sz1 <span class="pl-k">/</span> <span class="pl-c1">10.0</span>, spec_sz2 <span class="pl-k">/</span> <span class="pl-c1">10.0</span>], <span class="pl-v">order</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code blob-code-inner js-file-line">    tab1_tim <span class="pl-k">=</span> sn.interpolation.zoom(tab1_tim, spec_sz2 <span class="pl-k">/</span> <span class="pl-c1">10.0</span>, <span class="pl-v">order</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code blob-code-inner js-file-line">    tab1_freq <span class="pl-k">=</span> sn.interpolation.zoom(tab1_freq, spec_sz1 <span class="pl-k">/</span> <span class="pl-c1">10.0</span>, <span class="pl-v">order</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>RR<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code blob-code-inner js-file-line">        tab1_spec_plt <span class="pl-k">=</span> tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :]</td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> tab1_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>LL<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code blob-code-inner js-file-line">        tab1_spec_plt <span class="pl-k">=</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]</td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> tab1_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>I<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code blob-code-inner js-file-line">        tab1_spec_plt <span class="pl-k">=</span> (tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :] <span class="pl-k">+</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]) <span class="pl-k">/</span> <span class="pl-c1">2</span>.</td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> tab1_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>V<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code blob-code-inner js-file-line">        tab1_spec_plt <span class="pl-k">=</span> (tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :] <span class="pl-k">-</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]) <span class="pl-k">/</span> <span class="pl-c1">2</span>.</td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code blob-code-inner js-file-line">    tab1_dtim <span class="pl-k">=</span> tab1_tim <span class="pl-k">-</span> tab1_tim[<span class="pl-c1">0</span>]</td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"># dmax = np.amax(tab1_spec_plt)</span></td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"># dmin = np.amin(tab1_spec_plt)</span></td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"># tab1_spec_plt = (tab1_spec_plt - dmin) / (dmax - dmin) * 255.</span></td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code blob-code-inner js-file-line"><span class="pl-c">#todo lskajdflksadf</span></td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code blob-code-inner js-file-line">load_specdata(specfile)</td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code blob-code-inner js-file-line"><span class="pl-c1">TOOLS</span> <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>pan,wheel_zoom,box_zoom,reset,save<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span>create the dynamic spectrum plot<span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code blob-code-inner js-file-line">tab1_p_dspec <span class="pl-k">=</span> figure(<span class="pl-v">tools</span><span class="pl-k">=</span><span class="pl-c1">TOOLS</span>, <span class="pl-v">webgl</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>WebGL<span class="pl-pds">&#39;</span></span>],</td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">plot_width</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>dspec_wdth<span class="pl-pds">&#39;</span></span>],</td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">plot_height</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>dspec_hght<span class="pl-pds">&#39;</span></span>], <span class="pl-v">x_range</span><span class="pl-k">=</span>(tab1_dtim[<span class="pl-c1">0</span>], tab1_dtim[<span class="pl-k">-</span><span class="pl-c1">1</span>]),</td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">y_range</span><span class="pl-k">=</span>(tab1_freq[<span class="pl-c1">0</span>], tab1_freq[<span class="pl-k">-</span><span class="pl-c1">1</span>]), <span class="pl-v">toolbar_location</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>above<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code blob-code-inner js-file-line">tim0_char <span class="pl-k">=</span> jdutil.jd_to_datetime((tab1_tim[<span class="pl-c1">0</span>] <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>. <span class="pl-k">+</span> <span class="pl-c1">2400000.5</span>) <span class="pl-k">*</span> <span class="pl-c1">86400</span>. <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>.)</td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code blob-code-inner js-file-line">tim0_char <span class="pl-k">=</span> tim0_char.strftime(<span class="pl-s"><span class="pl-pds">&#39;</span>%Y-%b-<span class="pl-c1">%d</span> %H:%M:%S<span class="pl-pds">&#39;</span></span>) <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>.<span class="pl-c1">{}</span><span class="pl-pds">&#39;</span></span>.format(<span class="pl-c1">round</span>(tim0_char.microsecond <span class="pl-k">/</span> <span class="pl-c1">1e3</span>) <span class="pl-k">*</span> <span class="pl-c1">1e3</span>)[<span class="pl-c1">0</span>:<span class="pl-c1">4</span>]</td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.visible <span class="pl-k">=</span> <span class="pl-c1">True</span></td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.title.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>Dynamic spectrum<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L112" class="blob-num js-line-number" data-line-number="112"></td>
        <td id="LC112" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.xaxis.axis_label <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&#39;</span>Seconds since <span class="pl-pds">&#39;</span></span> <span class="pl-k">+</span> tim0_char</td>
      </tr>
      <tr>
        <td id="L113" class="blob-num js-line-number" data-line-number="113"></td>
        <td id="LC113" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.yaxis.axis_label <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&#39;</span>Frequency [GHz]<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L114" class="blob-num js-line-number" data-line-number="114"></td>
        <td id="LC114" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.border_fill_color <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>whitesmoke<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L115" class="blob-num js-line-number" data-line-number="115"></td>
        <td id="LC115" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_out <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L116" class="blob-num js-line-number" data-line-number="116"></td>
        <td id="LC116" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_in <span class="pl-k">=</span> <span class="pl-c1">5</span></td>
      </tr>
      <tr>
        <td id="L117" class="blob-num js-line-number" data-line-number="117"></td>
        <td id="LC117" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_out <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L118" class="blob-num js-line-number" data-line-number="118"></td>
        <td id="LC118" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_in <span class="pl-k">=</span> <span class="pl-c1">3</span></td>
      </tr>
      <tr>
        <td id="L119" class="blob-num js-line-number" data-line-number="119"></td>
        <td id="LC119" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_line_color <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L120" class="blob-num js-line-number" data-line-number="120"></td>
        <td id="LC120" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_line_color <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L121" class="blob-num js-line-number" data-line-number="121"></td>
        <td id="LC121" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L122" class="blob-num js-line-number" data-line-number="122"></td>
        <td id="LC122" class="blob-code blob-code-inner js-file-line">tab1_SRC_dspec <span class="pl-k">=</span> ColumnDataSource(<span class="pl-v">data</span><span class="pl-k">=</span>{<span class="pl-s"><span class="pl-pds">&#39;</span>data<span class="pl-pds">&#39;</span></span>: [tab1_spec_plt], <span class="pl-s"><span class="pl-pds">&#39;</span>xx<span class="pl-pds">&#39;</span></span>: [tab1_dtim], <span class="pl-s"><span class="pl-pds">&#39;</span>yy<span class="pl-pds">&#39;</span></span>: [tab1_freq]})</td>
      </tr>
      <tr>
        <td id="L123" class="blob-num js-line-number" data-line-number="123"></td>
        <td id="LC123" class="blob-code blob-code-inner js-file-line">tab1_r_dspec <span class="pl-k">=</span> tab1_p_dspec.image(<span class="pl-v">image</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>data<span class="pl-pds">&quot;</span></span>, <span class="pl-v">x</span><span class="pl-k">=</span>tab1_dtim[<span class="pl-c1">0</span>], <span class="pl-v">y</span><span class="pl-k">=</span>tab1_freq[<span class="pl-c1">0</span>], <span class="pl-v">dw</span><span class="pl-k">=</span>tab1_dtim[<span class="pl-k">-</span><span class="pl-c1">1</span>] <span class="pl-k">-</span> tab1_dtim[<span class="pl-c1">0</span>],</td>
      </tr>
      <tr>
        <td id="L124" class="blob-num js-line-number" data-line-number="124"></td>
        <td id="LC124" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">dh</span><span class="pl-k">=</span>tab1_freq[<span class="pl-k">-</span><span class="pl-c1">1</span>] <span class="pl-k">-</span> tab1_freq[<span class="pl-c1">0</span>], <span class="pl-v">source</span><span class="pl-k">=</span>tab1_SRC_dspec, <span class="pl-v">palette</span><span class="pl-k">=</span>bokehpalette_jet)</td>
      </tr>
      <tr>
        <td id="L125" class="blob-num js-line-number" data-line-number="125"></td>
        <td id="LC125" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L126" class="blob-num js-line-number" data-line-number="126"></td>
        <td id="LC126" class="blob-code blob-code-inner js-file-line">tab1_spec_sz <span class="pl-k">=</span> tab1_spec.shape</td>
      </tr>
      <tr>
        <td id="L127" class="blob-num js-line-number" data-line-number="127"></td>
        <td id="LC127" class="blob-code blob-code-inner js-file-line">ratio_spec_sz2, ratio_spec_sz1 <span class="pl-k">=</span> <span class="pl-c1">10</span>, <span class="pl-c1">10</span></td>
      </tr>
      <tr>
        <td id="L128" class="blob-num js-line-number" data-line-number="128"></td>
        <td id="LC128" class="blob-code blob-code-inner js-file-line"><span class="pl-k">if</span> tab1_spec_sz[<span class="pl-c1">3</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">520</span>:</td>
      </tr>
      <tr>
        <td id="L129" class="blob-num js-line-number" data-line-number="129"></td>
        <td id="LC129" class="blob-code blob-code-inner js-file-line">    ratio_spec_sz2 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">1</span>, <span class="pl-c1">10</span>) <span class="pl-k">if</span> i <span class="pl-k">/</span> <span class="pl-c1">10</span>. <span class="pl-k">*</span> tab1_spec_sz[<span class="pl-c1">3</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">520</span>)</td>
      </tr>
      <tr>
        <td id="L130" class="blob-num js-line-number" data-line-number="130"></td>
        <td id="LC130" class="blob-code blob-code-inner js-file-line"><span class="pl-k">if</span> tab1_spec_sz[<span class="pl-c1">2</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">75</span>:</td>
      </tr>
      <tr>
        <td id="L131" class="blob-num js-line-number" data-line-number="131"></td>
        <td id="LC131" class="blob-code blob-code-inner js-file-line">    ratio_spec_sz1 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(<span class="pl-c1">1</span>, <span class="pl-c1">10</span>) <span class="pl-k">if</span> i <span class="pl-k">/</span> <span class="pl-c1">10</span>. <span class="pl-k">*</span> tab1_spec_sz[<span class="pl-c1">2</span>] <span class="pl-k">&gt;</span> <span class="pl-c1">75</span>)</td>
      </tr>
      <tr>
        <td id="L132" class="blob-num js-line-number" data-line-number="132"></td>
        <td id="LC132" class="blob-code blob-code-inner js-file-line">tab1_tim_square <span class="pl-k">=</span> sn.interpolation.zoom(tab1_tim, ratio_spec_sz2 <span class="pl-k">/</span> <span class="pl-c1">10</span>., <span class="pl-v">order</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L133" class="blob-num js-line-number" data-line-number="133"></td>
        <td id="LC133" class="blob-code blob-code-inner js-file-line">tab1_freq_square <span class="pl-k">=</span> sn.interpolation.zoom(tab1_freq, ratio_spec_sz1 <span class="pl-k">/</span> <span class="pl-c1">10</span>., <span class="pl-v">order</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L134" class="blob-num js-line-number" data-line-number="134"></td>
        <td id="LC134" class="blob-code blob-code-inner js-file-line">tab1_nfreq_square, <span class="pl-k">=</span> tab1_freq_square.shape</td>
      </tr>
      <tr>
        <td id="L135" class="blob-num js-line-number" data-line-number="135"></td>
        <td id="LC135" class="blob-code blob-code-inner js-file-line">tab1_ntim_square, <span class="pl-k">=</span> tab1_tim_square.shape</td>
      </tr>
      <tr>
        <td id="L136" class="blob-num js-line-number" data-line-number="136"></td>
        <td id="LC136" class="blob-code blob-code-inner js-file-line">tim_map <span class="pl-k">=</span> ((np.tile(tab1_tim_square, tab1_nfreq_square).reshape(tab1_nfreq_square,</td>
      </tr>
      <tr>
        <td id="L137" class="blob-num js-line-number" data-line-number="137"></td>
        <td id="LC137" class="blob-code blob-code-inner js-file-line">    tab1_ntim_square) <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>. <span class="pl-k">+</span> <span class="pl-c1">2400000.5</span>)) <span class="pl-k">*</span> <span class="pl-c1">86400</span>.</td>
      </tr>
      <tr>
        <td id="L138" class="blob-num js-line-number" data-line-number="138"></td>
        <td id="LC138" class="blob-code blob-code-inner js-file-line">freq_map <span class="pl-k">=</span> np.tile(tab1_freq_square, tab1_ntim_square).reshape(tab1_ntim_square, tab1_nfreq_square).swapaxes(<span class="pl-c1">0</span>, <span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L139" class="blob-num js-line-number" data-line-number="139"></td>
        <td id="LC139" class="blob-code blob-code-inner js-file-line">xx <span class="pl-k">=</span> tim_map.flatten()</td>
      </tr>
      <tr>
        <td id="L140" class="blob-num js-line-number" data-line-number="140"></td>
        <td id="LC140" class="blob-code blob-code-inner js-file-line">yy <span class="pl-k">=</span> freq_map.flatten()</td>
      </tr>
      <tr>
        <td id="L141" class="blob-num js-line-number" data-line-number="141"></td>
        <td id="LC141" class="blob-code blob-code-inner js-file-line">tab1_dspecDF_square <span class="pl-k">=</span> pd.DataFrame({<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>: xx <span class="pl-k">-</span> xx[<span class="pl-c1">0</span>], <span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>: yy})</td>
      </tr>
      <tr>
        <td id="L142" class="blob-num js-line-number" data-line-number="142"></td>
        <td id="LC142" class="blob-code blob-code-inner js-file-line">tab1_SRC_dspec_square <span class="pl-k">=</span> ColumnDataSource(tab1_dspecDF_square)</td>
      </tr>
      <tr>
        <td id="L143" class="blob-num js-line-number" data-line-number="143"></td>
        <td id="LC143" class="blob-code blob-code-inner js-file-line"><span class="pl-s"><span class="pl-pds">&#39;&#39;&#39;</span>make the plot lasso selectable<span class="pl-pds">&#39;&#39;&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L144" class="blob-num js-line-number" data-line-number="144"></td>
        <td id="LC144" class="blob-code blob-code-inner js-file-line">tab1_render_square <span class="pl-k">=</span> tab1_p_dspec.square(<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>, <span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>, <span class="pl-v">source</span><span class="pl-k">=</span>tab1_SRC_dspec_square, <span class="pl-v">fill_color</span><span class="pl-k">=</span><span class="pl-c1">None</span>, <span class="pl-v">fill_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>,</td>
      </tr>
      <tr>
        <td id="L145" class="blob-num js-line-number" data-line-number="145"></td>
        <td id="LC145" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">line_color</span><span class="pl-k">=</span><span class="pl-c1">None</span>, <span class="pl-v">line_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>, <span class="pl-v">selection_fill_color</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>red<span class="pl-pds">&#39;</span></span>, <span class="pl-v">selection_fill_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.2</span>, <span class="pl-v">nonselection_fill_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>,</td>
      </tr>
      <tr>
        <td id="L146" class="blob-num js-line-number" data-line-number="146"></td>
        <td id="LC146" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">selection_line_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>, <span class="pl-v">nonselection_line_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>)</td>
      </tr>
      <tr>
        <td id="L147" class="blob-num js-line-number" data-line-number="147"></td>
        <td id="LC147" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L148" class="blob-num js-line-number" data-line-number="148"></td>
        <td id="LC148" class="blob-code blob-code-inner js-file-line"><span class="pl-c">## ----------------baseline &amp; polarization selection------------------------</span></td>
      </tr>
      <tr>
        <td id="L149" class="blob-num js-line-number" data-line-number="149"></td>
        <td id="LC149" class="blob-code blob-code-inner js-file-line">tab1_Select_bl <span class="pl-k">=</span> Select(<span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Baseline:<span class="pl-pds">&quot;</span></span>, <span class="pl-v">value</span><span class="pl-k">=</span>tab1_bl[<span class="pl-c1">0</span>], <span class="pl-v">options</span><span class="pl-k">=</span>tab1_bl, <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>)</td>
      </tr>
      <tr>
        <td id="L150" class="blob-num js-line-number" data-line-number="150"></td>
        <td id="LC150" class="blob-code blob-code-inner js-file-line">tab1_Select_pol <span class="pl-k">=</span> Select(<span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Polarization:<span class="pl-pds">&quot;</span></span>, <span class="pl-v">value</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>I<span class="pl-pds">&quot;</span></span>, <span class="pl-v">options</span><span class="pl-k">=</span>[<span class="pl-s"><span class="pl-pds">&quot;</span>RR<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>LL<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>I<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>V<span class="pl-pds">&quot;</span></span>], <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>)</td>
      </tr>
      <tr>
        <td id="L151" class="blob-num js-line-number" data-line-number="151"></td>
        <td id="LC151" class="blob-code blob-code-inner js-file-line">tab1_Select_colormap <span class="pl-k">=</span> Select(<span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Colormap:<span class="pl-pds">&quot;</span></span>, <span class="pl-v">value</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>linear<span class="pl-pds">&quot;</span></span>, <span class="pl-v">options</span><span class="pl-k">=</span>[<span class="pl-s"><span class="pl-pds">&quot;</span>linear<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>log<span class="pl-pds">&quot;</span></span>], <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>)</td>
      </tr>
      <tr>
        <td id="L152" class="blob-num js-line-number" data-line-number="152"></td>
        <td id="LC152" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L153" class="blob-num js-line-number" data-line-number="153"></td>
        <td id="LC153" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L154" class="blob-num js-line-number" data-line-number="154"></td>
        <td id="LC154" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_dspec</span>(<span class="pl-smi">attrname</span>, <span class="pl-smi">old</span>, <span class="pl-smi">new</span>):</td>
      </tr>
      <tr>
        <td id="L155" class="blob-num js-line-number" data-line-number="155"></td>
        <td id="LC155" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_spec, tab1_dtim, tab1_freq, tab1_bl</td>
      </tr>
      <tr>
        <td id="L156" class="blob-num js-line-number" data-line-number="156"></td>
        <td id="LC156" class="blob-code blob-code-inner js-file-line">    select_pol <span class="pl-k">=</span> tab1_Select_pol.value</td>
      </tr>
      <tr>
        <td id="L157" class="blob-num js-line-number" data-line-number="157"></td>
        <td id="LC157" class="blob-code blob-code-inner js-file-line">    select_bl <span class="pl-k">=</span> tab1_Select_bl.value</td>
      </tr>
      <tr>
        <td id="L158" class="blob-num js-line-number" data-line-number="158"></td>
        <td id="LC158" class="blob-code blob-code-inner js-file-line">    bl_index <span class="pl-k">=</span> tab1_bl.index(select_bl)</td>
      </tr>
      <tr>
        <td id="L159" class="blob-num js-line-number" data-line-number="159"></td>
        <td id="LC159" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> select_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>RR<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L160" class="blob-num js-line-number" data-line-number="160"></td>
        <td id="LC160" class="blob-code blob-code-inner js-file-line">        spec_plt <span class="pl-k">=</span> tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :]</td>
      </tr>
      <tr>
        <td id="L161" class="blob-num js-line-number" data-line-number="161"></td>
        <td id="LC161" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> select_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>LL<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L162" class="blob-num js-line-number" data-line-number="162"></td>
        <td id="LC162" class="blob-code blob-code-inner js-file-line">        spec_plt <span class="pl-k">=</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]</td>
      </tr>
      <tr>
        <td id="L163" class="blob-num js-line-number" data-line-number="163"></td>
        <td id="LC163" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> select_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>I<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L164" class="blob-num js-line-number" data-line-number="164"></td>
        <td id="LC164" class="blob-code blob-code-inner js-file-line">        spec_plt <span class="pl-k">=</span> (tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :] <span class="pl-k">+</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]) <span class="pl-k">/</span> <span class="pl-c1">2</span>.</td>
      </tr>
      <tr>
        <td id="L165" class="blob-num js-line-number" data-line-number="165"></td>
        <td id="LC165" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">elif</span> select_pol <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>V<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L166" class="blob-num js-line-number" data-line-number="166"></td>
        <td id="LC166" class="blob-code blob-code-inner js-file-line">        spec_plt <span class="pl-k">=</span> (tab1_spec[<span class="pl-c1">0</span>, bl_index, :, :] <span class="pl-k">-</span> tab1_spec[<span class="pl-c1">1</span>, bl_index, :, :]) <span class="pl-k">/</span> <span class="pl-c1">2</span>.</td>
      </tr>
      <tr>
        <td id="L167" class="blob-num js-line-number" data-line-number="167"></td>
        <td id="LC167" class="blob-code blob-code-inner js-file-line">        tab1_Select_colormap.value <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&#39;</span>linear<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L168" class="blob-num js-line-number" data-line-number="168"></td>
        <td id="LC168" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_Select_colormap.value <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&#39;</span>log<span class="pl-pds">&#39;</span></span> <span class="pl-k">and</span> select_pol <span class="pl-k">!=</span> <span class="pl-s"><span class="pl-pds">&#39;</span>V<span class="pl-pds">&#39;</span></span>:</td>
      </tr>
      <tr>
        <td id="L169" class="blob-num js-line-number" data-line-number="169"></td>
        <td id="LC169" class="blob-code blob-code-inner js-file-line">        spec_plt <span class="pl-k">=</span> np.log(spec_plt)</td>
      </tr>
      <tr>
        <td id="L170" class="blob-num js-line-number" data-line-number="170"></td>
        <td id="LC170" class="blob-code blob-code-inner js-file-line">    tab1_SRC_dspec.data <span class="pl-k">=</span> {<span class="pl-s"><span class="pl-pds">&#39;</span>data<span class="pl-pds">&#39;</span></span>: [spec_plt], <span class="pl-s"><span class="pl-pds">&#39;</span>xx<span class="pl-pds">&#39;</span></span>: [tab1_dtim], <span class="pl-s"><span class="pl-pds">&#39;</span>yy<span class="pl-pds">&#39;</span></span>: [tab1_freq]}</td>
      </tr>
      <tr>
        <td id="L171" class="blob-num js-line-number" data-line-number="171"></td>
        <td id="LC171" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L172" class="blob-num js-line-number" data-line-number="172"></td>
        <td id="LC172" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L173" class="blob-num js-line-number" data-line-number="173"></td>
        <td id="LC173" class="blob-code blob-code-inner js-file-line">tab1_ctrls <span class="pl-k">=</span> [tab1_Select_bl, tab1_Select_pol, tab1_Select_colormap]</td>
      </tr>
      <tr>
        <td id="L174" class="blob-num js-line-number" data-line-number="174"></td>
        <td id="LC174" class="blob-code blob-code-inner js-file-line"><span class="pl-k">for</span> ctrl <span class="pl-k">in</span> tab1_ctrls:</td>
      </tr>
      <tr>
        <td id="L175" class="blob-num js-line-number" data-line-number="175"></td>
        <td id="LC175" class="blob-code blob-code-inner js-file-line">    ctrl.on_change(<span class="pl-s"><span class="pl-pds">&#39;</span>value<span class="pl-pds">&#39;</span></span>, tab1_update_dspec)</td>
      </tr>
      <tr>
        <td id="L176" class="blob-num js-line-number" data-line-number="176"></td>
        <td id="LC176" class="blob-code blob-code-inner js-file-line"><span class="pl-k">try</span>:</td>
      </tr>
      <tr>
        <td id="L177" class="blob-num js-line-number" data-line-number="177"></td>
        <td id="LC177" class="blob-code blob-code-inner js-file-line">    os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>cp <span class="pl-c1">{}</span>StrID_list.json <span class="pl-c1">{}</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>.format(database_dir <span class="pl-k">+</span> event_id, database_dir <span class="pl-k">+</span> event_id))</td>
      </tr>
      <tr>
        <td id="L178" class="blob-num js-line-number" data-line-number="178"></td>
        <td id="LC178" class="blob-code blob-code-inner js-file-line">    StrIDList <span class="pl-k">=</span> pd.read_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L179" class="blob-num js-line-number" data-line-number="179"></td>
        <td id="LC179" class="blob-code blob-code-inner js-file-line">    StrIDList <span class="pl-k">=</span> StrIDList.sort_values(<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>, <span class="pl-v">ascending</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L180" class="blob-num js-line-number" data-line-number="180"></td>
        <td id="LC180" class="blob-code blob-code-inner js-file-line">    StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>] <span class="pl-k">=</span> [ll <span class="pl-k">-</span> tab1_tim[<span class="pl-c1">0</span>] <span class="pl-k">for</span> ll <span class="pl-k">in</span> StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>]]</td>
      </tr>
      <tr>
        <td id="L181" class="blob-num js-line-number" data-line-number="181"></td>
        <td id="LC181" class="blob-code blob-code-inner js-file-line"><span class="pl-k">except</span>:</td>
      </tr>
      <tr>
        <td id="L182" class="blob-num js-line-number" data-line-number="182"></td>
        <td id="LC182" class="blob-code blob-code-inner js-file-line">    StrIDList <span class="pl-k">=</span> pd.DataFrame({<span class="pl-s"><span class="pl-pds">&#39;</span>date<span class="pl-pds">&#39;</span></span>: [], <span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>: [], <span class="pl-s"><span class="pl-pds">&#39;</span>freqran<span class="pl-pds">&#39;</span></span>: [], <span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>: [], <span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>: [], <span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>: []})</td>
      </tr>
      <tr>
        <td id="L183" class="blob-num js-line-number" data-line-number="183"></td>
        <td id="LC183" class="blob-code blob-code-inner js-file-line">    StrIDList.to_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L184" class="blob-num js-line-number" data-line-number="184"></td>
        <td id="LC184" class="blob-code blob-code-inner js-file-line">    os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>cp <span class="pl-c1">{}</span>StrID_list.json <span class="pl-c1">{}</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>.format(database_dir <span class="pl-k">+</span> event_id, database_dir <span class="pl-k">+</span> event_id))</td>
      </tr>
      <tr>
        <td id="L185" class="blob-num js-line-number" data-line-number="185"></td>
        <td id="LC185" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L186" class="blob-num js-line-number" data-line-number="186"></td>
        <td id="LC186" class="blob-code blob-code-inner js-file-line">tab1_SRC_StrIDPatch <span class="pl-k">=</span> ColumnDataSource(StrIDList)</td>
      </tr>
      <tr>
        <td id="L187" class="blob-num js-line-number" data-line-number="187"></td>
        <td id="LC187" class="blob-code blob-code-inner js-file-line">tab1_render_patch <span class="pl-k">=</span> tab1_p_dspec.patches(<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>, <span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>, <span class="pl-v">source</span><span class="pl-k">=</span>tab1_SRC_StrIDPatch, <span class="pl-v">hover_fill_color</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>firebrick<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L188" class="blob-num js-line-number" data-line-number="188"></td>
        <td id="LC188" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">hover_line_color</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>#cc3333<span class="pl-pds">&#39;</span></span>, <span class="pl-v">hover_fill_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.2</span>, <span class="pl-v">fill_color</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>grey<span class="pl-pds">&quot;</span></span>, <span class="pl-v">fill_alpha</span><span class="pl-k">=</span><span class="pl-c1">0.0</span>, <span class="pl-v">line_color</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L189" class="blob-num js-line-number" data-line-number="189"></td>
        <td id="LC189" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">line_width</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L190" class="blob-num js-line-number" data-line-number="190"></td>
        <td id="LC190" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L191" class="blob-num js-line-number" data-line-number="191"></td>
        <td id="LC191" class="blob-code blob-code-inner js-file-line">tooltips <span class="pl-k">=</span> [(<span class="pl-s"><span class="pl-pds">&quot;</span>StrID<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>@str_id<span class="pl-pds">&quot;</span></span>), (<span class="pl-s"><span class="pl-pds">&quot;</span>Date<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>@date<span class="pl-pds">&quot;</span></span>), (<span class="pl-s"><span class="pl-pds">&quot;</span>TimeRange<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>@timeran<span class="pl-pds">&quot;</span></span>), (<span class="pl-s"><span class="pl-pds">&quot;</span>FreqRange<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>@freqran<span class="pl-pds">&quot;</span></span>)]</td>
      </tr>
      <tr>
        <td id="L192" class="blob-num js-line-number" data-line-number="192"></td>
        <td id="LC192" class="blob-code blob-code-inner js-file-line">tab1_HoverTool_dspec <span class="pl-k">=</span> HoverTool(<span class="pl-v">tooltips</span><span class="pl-k">=</span>tooltips, <span class="pl-v">anchor</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>top_left<span class="pl-pds">&#39;</span></span>, <span class="pl-v">renderers</span><span class="pl-k">=</span>[tab1_render_patch])</td>
      </tr>
      <tr>
        <td id="L193" class="blob-num js-line-number" data-line-number="193"></td>
        <td id="LC193" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.add_tools(tab1_HoverTool_dspec)</td>
      </tr>
      <tr>
        <td id="L194" class="blob-num js-line-number" data-line-number="194"></td>
        <td id="LC194" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.add_tools(BoxSelectTool())</td>
      </tr>
      <tr>
        <td id="L195" class="blob-num js-line-number" data-line-number="195"></td>
        <td id="LC195" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.add_tools(TapTool(<span class="pl-v">renderers</span><span class="pl-k">=</span>[tab1_render_patch]))</td>
      </tr>
      <tr>
        <td id="L196" class="blob-num js-line-number" data-line-number="196"></td>
        <td id="LC196" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.select(BoxSelectTool).select_every_mousemove <span class="pl-k">=</span> <span class="pl-c1">False</span></td>
      </tr>
      <tr>
        <td id="L197" class="blob-num js-line-number" data-line-number="197"></td>
        <td id="LC197" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_out <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L198" class="blob-num js-line-number" data-line-number="198"></td>
        <td id="LC198" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_in <span class="pl-k">=</span> <span class="pl-c1">5</span></td>
      </tr>
      <tr>
        <td id="L199" class="blob-num js-line-number" data-line-number="199"></td>
        <td id="LC199" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_out <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L200" class="blob-num js-line-number" data-line-number="200"></td>
        <td id="LC200" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_in <span class="pl-k">=</span> <span class="pl-c1">3</span></td>
      </tr>
      <tr>
        <td id="L201" class="blob-num js-line-number" data-line-number="201"></td>
        <td id="LC201" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.major_tick_line_color <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L202" class="blob-num js-line-number" data-line-number="202"></td>
        <td id="LC202" class="blob-code blob-code-inner js-file-line">tab1_p_dspec.axis.minor_tick_line_color <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L203" class="blob-num js-line-number" data-line-number="203"></td>
        <td id="LC203" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L204" class="blob-num js-line-number" data-line-number="204"></td>
        <td id="LC204" class="blob-code blob-code-inner js-file-line">tab1_TbCols <span class="pl-k">=</span> [TableColumn(<span class="pl-v">field</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>str_id<span class="pl-pds">&quot;</span></span>, <span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>StrID<span class="pl-pds">&quot;</span></span>), TableColumn(<span class="pl-v">field</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>timeran<span class="pl-pds">&quot;</span></span>, <span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Time Range<span class="pl-pds">&quot;</span></span>),</td>
      </tr>
      <tr>
        <td id="L205" class="blob-num js-line-number" data-line-number="205"></td>
        <td id="LC205" class="blob-code blob-code-inner js-file-line">               TableColumn(<span class="pl-v">field</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>freqran<span class="pl-pds">&quot;</span></span>, <span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Freq Range<span class="pl-pds">&quot;</span></span>), ]</td>
      </tr>
      <tr>
        <td id="L206" class="blob-num js-line-number" data-line-number="206"></td>
        <td id="LC206" class="blob-code blob-code-inner js-file-line">tab1_DataTb_dspec <span class="pl-k">=</span> DataTable(<span class="pl-v">source</span><span class="pl-k">=</span>tab1_render_patch.data_source, <span class="pl-v">columns</span><span class="pl-k">=</span>tab1_TbCols,</td>
      </tr>
      <tr>
        <td id="L207" class="blob-num js-line-number" data-line-number="207"></td>
        <td id="LC207" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">width</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>StrID_DataTb_wdth<span class="pl-pds">&#39;</span></span>],</td>
      </tr>
      <tr>
        <td id="L208" class="blob-num js-line-number" data-line-number="208"></td>
        <td id="LC208" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">height</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>StrID_DataTb_hght<span class="pl-pds">&#39;</span></span>])  <span class="pl-c"># , editable=True)</span></td>
      </tr>
      <tr>
        <td id="L209" class="blob-num js-line-number" data-line-number="209"></td>
        <td id="LC209" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L210" class="blob-num js-line-number" data-line-number="210"></td>
        <td id="LC210" class="blob-code blob-code-inner js-file-line">tab1_Div_Tb <span class="pl-k">=</span> Div(<span class="pl-v">text</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span> <span class="pl-pds">&quot;&quot;&quot;</span></span>, <span class="pl-v">width</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>StrID_DataTb_wdth<span class="pl-pds">&#39;</span></span>])</td>
      </tr>
      <tr>
        <td id="L211" class="blob-num js-line-number" data-line-number="211"></td>
        <td id="LC211" class="blob-code blob-code-inner js-file-line">tab1_Div_exit <span class="pl-k">=</span> Div(<span class="pl-v">text</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L212" class="blob-num js-line-number" data-line-number="212"></td>
        <td id="LC212" class="blob-code blob-code-inner js-file-line"><span class="pl-s">&lt;p&gt;&lt;b&gt;Warning&lt;/b&gt;: 1. Click the &lt;b&gt;Exit QLook&lt;/b&gt; first before closing the tab&lt;/p&gt;</span></td>
      </tr>
      <tr>
        <td id="L213" class="blob-num js-line-number" data-line-number="213"></td>
        <td id="LC213" class="blob-code blob-code-inner js-file-line"><span class="pl-s">&lt;p&gt;&lt;b&gt;Warning&lt;/b&gt;: 2. &lt;b&gt;FSview&lt;/b&gt; or &lt;b&gt;FSview2CASA&lt;/b&gt; tabs will disconnect if &lt;b&gt;Exit QLook is clicked&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L214" class="blob-num js-line-number" data-line-number="214"></td>
        <td id="LC214" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>)</td>
      </tr>
      <tr>
        <td id="L215" class="blob-num js-line-number" data-line-number="215"></td>
        <td id="LC215" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L216" class="blob-num js-line-number" data-line-number="216"></td>
        <td id="LC216" class="blob-code blob-code-inner js-file-line">tab1_selected_SRC_dspec_square <span class="pl-k">=</span> <span class="pl-c1">None</span></td>
      </tr>
      <tr>
        <td id="L217" class="blob-num js-line-number" data-line-number="217"></td>
        <td id="LC217" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L218" class="blob-num js-line-number" data-line-number="218"></td>
        <td id="LC218" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L219" class="blob-num js-line-number" data-line-number="219"></td>
        <td id="LC219" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_SRC_dspec_square_select</span>(<span class="pl-smi">attrname</span>, <span class="pl-smi">old</span>, <span class="pl-smi">new</span>):</td>
      </tr>
      <tr>
        <td id="L220" class="blob-num js-line-number" data-line-number="220"></td>
        <td id="LC220" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_selected_SRC_dspec_square</td>
      </tr>
      <tr>
        <td id="L221" class="blob-num js-line-number" data-line-number="221"></td>
        <td id="LC221" class="blob-code blob-code-inner js-file-line">    tab1_selected_SRC_dspec_square <span class="pl-k">=</span> tab1_SRC_dspec_square.selected[<span class="pl-s"><span class="pl-pds">&#39;</span>1d<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>indices<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L222" class="blob-num js-line-number" data-line-number="222"></td>
        <td id="LC222" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_selected_SRC_dspec_square:</td>
      </tr>
      <tr>
        <td id="L223" class="blob-num js-line-number" data-line-number="223"></td>
        <td id="LC223" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">global</span> dftmp</td>
      </tr>
      <tr>
        <td id="L224" class="blob-num js-line-number" data-line-number="224"></td>
        <td id="LC224" class="blob-code blob-code-inner js-file-line">        dftmp <span class="pl-k">=</span> tab1_dspecDF_square.iloc[tab1_selected_SRC_dspec_square, :]</td>
      </tr>
      <tr>
        <td id="L225" class="blob-num js-line-number" data-line-number="225"></td>
        <td id="LC225" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L226" class="blob-num js-line-number" data-line-number="226"></td>
        <td id="LC226" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L227" class="blob-num js-line-number" data-line-number="227"></td>
        <td id="LC227" class="blob-code blob-code-inner js-file-line">tab1_SRC_dspec_square.on_change(<span class="pl-s"><span class="pl-pds">&#39;</span>selected<span class="pl-pds">&#39;</span></span>, tab1_SRC_dspec_square_select)</td>
      </tr>
      <tr>
        <td id="L228" class="blob-num js-line-number" data-line-number="228"></td>
        <td id="LC228" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L229" class="blob-num js-line-number" data-line-number="229"></td>
        <td id="LC229" class="blob-code blob-code-inner js-file-line">tab1_input_StrID <span class="pl-k">=</span> TextInput(<span class="pl-v">value</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Type in here<span class="pl-pds">&quot;</span></span>, <span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>New StrID:<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L230" class="blob-num js-line-number" data-line-number="230"></td>
        <td id="LC230" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">width</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>StrID_DataTb_BUT_wdth<span class="pl-pds">&#39;</span></span>])</td>
      </tr>
      <tr>
        <td id="L231" class="blob-num js-line-number" data-line-number="231"></td>
        <td id="LC231" class="blob-code blob-code-inner js-file-line">timestart <span class="pl-k">=</span> xx[<span class="pl-c1">0</span>]</td>
      </tr>
      <tr>
        <td id="L232" class="blob-num js-line-number" data-line-number="232"></td>
        <td id="LC232" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L233" class="blob-num js-line-number" data-line-number="233"></td>
        <td id="LC233" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L234" class="blob-num js-line-number" data-line-number="234"></td>
        <td id="LC234" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_addStrID</span>():</td>
      </tr>
      <tr>
        <td id="L235" class="blob-num js-line-number" data-line-number="235"></td>
        <td id="LC235" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> dftmp, tab1_tim_square, timestart, database_dir, tab1_selected_SRC_dspec_squarez</td>
      </tr>
      <tr>
        <td id="L236" class="blob-num js-line-number" data-line-number="236"></td>
        <td id="LC236" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_selected_SRC_dspec_square:</td>
      </tr>
      <tr>
        <td id="L237" class="blob-num js-line-number" data-line-number="237"></td>
        <td id="LC237" class="blob-code blob-code-inner js-file-line">        time0, time1 <span class="pl-k">=</span> dftmp[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>].min() <span class="pl-k">+</span> timestart, dftmp[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>].max() <span class="pl-k">+</span> timestart</td>
      </tr>
      <tr>
        <td id="L238" class="blob-num js-line-number" data-line-number="238"></td>
        <td id="LC238" class="blob-code blob-code-inner js-file-line">        freq0, freq1 <span class="pl-k">=</span> dftmp[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>].min(), dftmp[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>].max()</td>
      </tr>
      <tr>
        <td id="L239" class="blob-num js-line-number" data-line-number="239"></td>
        <td id="LC239" class="blob-code blob-code-inner js-file-line">        date_char <span class="pl-k">=</span> jdutil.jd_to_datetime(timestart <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>.)</td>
      </tr>
      <tr>
        <td id="L240" class="blob-num js-line-number" data-line-number="240"></td>
        <td id="LC240" class="blob-code blob-code-inner js-file-line">        date_char <span class="pl-k">=</span> date_char.strftime(<span class="pl-s"><span class="pl-pds">&#39;</span>%Y-%b-<span class="pl-c1">%d</span><span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L241" class="blob-num js-line-number" data-line-number="241"></td>
        <td id="LC241" class="blob-code blob-code-inner js-file-line">        t0_char <span class="pl-k">=</span> jdutil.jd_to_datetime(time0 <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>.)</td>
      </tr>
      <tr>
        <td id="L242" class="blob-num js-line-number" data-line-number="242"></td>
        <td id="LC242" class="blob-code blob-code-inner js-file-line">        t0_char <span class="pl-k">=</span> t0_char.strftime(<span class="pl-s"><span class="pl-pds">&#39;</span>%H:%M:%S<span class="pl-pds">&#39;</span></span>) <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>.<span class="pl-c1">{<span class="pl-c1">:03d</span>}</span><span class="pl-pds">&#39;</span></span>.format(<span class="pl-c1">int</span>(<span class="pl-c1">round</span>(t0_char.microsecond <span class="pl-k">/</span> <span class="pl-c1">1e3</span>)))</td>
      </tr>
      <tr>
        <td id="L243" class="blob-num js-line-number" data-line-number="243"></td>
        <td id="LC243" class="blob-code blob-code-inner js-file-line">        t1_char <span class="pl-k">=</span> jdutil.jd_to_datetime(time1 <span class="pl-k">/</span> <span class="pl-c1">3600</span>. <span class="pl-k">/</span> <span class="pl-c1">24</span>.)</td>
      </tr>
      <tr>
        <td id="L244" class="blob-num js-line-number" data-line-number="244"></td>
        <td id="LC244" class="blob-code blob-code-inner js-file-line">        t1_char <span class="pl-k">=</span> t1_char.strftime(<span class="pl-s"><span class="pl-pds">&#39;</span>%H:%M:%S<span class="pl-pds">&#39;</span></span>) <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>.<span class="pl-c1">{<span class="pl-c1">:03d</span>}</span><span class="pl-pds">&#39;</span></span>.format(<span class="pl-c1">int</span>(<span class="pl-c1">round</span>(t1_char.microsecond <span class="pl-k">/</span> <span class="pl-c1">1e3</span>)))</td>
      </tr>
      <tr>
        <td id="L245" class="blob-num js-line-number" data-line-number="245"></td>
        <td id="LC245" class="blob-code blob-code-inner js-file-line">        time0, time1 <span class="pl-k">=</span> (time0 <span class="pl-k">/</span> <span class="pl-c1">86400</span>. <span class="pl-k">-</span> <span class="pl-c1">2400000.5</span>) <span class="pl-k">*</span> <span class="pl-c1">24</span>. <span class="pl-k">*</span> <span class="pl-c1">3600</span>., (time1 <span class="pl-k">/</span> <span class="pl-c1">86400</span>. <span class="pl-k">-</span> <span class="pl-c1">2400000.5</span>) <span class="pl-k">*</span> <span class="pl-c1">24</span>. <span class="pl-k">*</span> <span class="pl-c1">3600</span>.</td>
      </tr>
      <tr>
        <td id="L246" class="blob-num js-line-number" data-line-number="246"></td>
        <td id="LC246" class="blob-code blob-code-inner js-file-line">        StrID <span class="pl-k">=</span> pd.DataFrame(<span class="pl-c1">dict</span>(<span class="pl-v">time</span><span class="pl-k">=</span>[[time0, time1, time1, time0]], <span class="pl-v">freq</span><span class="pl-k">=</span>[[freq0, freq0, freq1, freq1]],</td>
      </tr>
      <tr>
        <td id="L247" class="blob-num js-line-number" data-line-number="247"></td>
        <td id="LC247" class="blob-code blob-code-inner js-file-line">            <span class="pl-v">str_id</span><span class="pl-k">=</span>[[tab1_input_StrID.value]], <span class="pl-v">date</span><span class="pl-k">=</span>[[date_char]], <span class="pl-v">timeran</span><span class="pl-k">=</span>[[t0_char <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>~<span class="pl-pds">&#39;</span></span> <span class="pl-k">+</span> t1_char]],</td>
      </tr>
      <tr>
        <td id="L248" class="blob-num js-line-number" data-line-number="248"></td>
        <td id="LC248" class="blob-code blob-code-inner js-file-line">            <span class="pl-v">freqran</span><span class="pl-k">=</span>[[<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-c1">{<span class="pl-c1">:.3f</span>}</span>~<span class="pl-c1">{<span class="pl-c1">:.3f</span>}</span> GHz<span class="pl-pds">&quot;</span></span>.format(freq0, freq1)]]))</td>
      </tr>
      <tr>
        <td id="L249" class="blob-num js-line-number" data-line-number="249"></td>
        <td id="LC249" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> pd.read_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L250" class="blob-num js-line-number" data-line-number="250"></td>
        <td id="LC250" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> pd.concat([StrIDList, StrID])</td>
      </tr>
      <tr>
        <td id="L251" class="blob-num js-line-number" data-line-number="251"></td>
        <td id="LC251" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> StrIDList.sort_values(<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>, <span class="pl-v">ascending</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L252" class="blob-num js-line-number" data-line-number="252"></td>
        <td id="LC252" class="blob-code blob-code-inner js-file-line">        StrIDList.index <span class="pl-k">=</span> <span class="pl-c1">range</span>(StrIDList.index.size)</td>
      </tr>
      <tr>
        <td id="L253" class="blob-num js-line-number" data-line-number="253"></td>
        <td id="LC253" class="blob-code blob-code-inner js-file-line">        StrIDList.to_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L254" class="blob-num js-line-number" data-line-number="254"></td>
        <td id="LC254" class="blob-code blob-code-inner js-file-line">        StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>] <span class="pl-k">=</span> [ll <span class="pl-k">-</span> tab1_tim_square[<span class="pl-c1">0</span>] <span class="pl-k">for</span> ll <span class="pl-k">in</span> StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>]]</td>
      </tr>
      <tr>
        <td id="L255" class="blob-num js-line-number" data-line-number="255"></td>
        <td id="LC255" class="blob-code blob-code-inner js-file-line">        tab1_render_patch.data_source.data <span class="pl-k">=</span> ColumnDataSource(StrIDList).data</td>
      </tr>
      <tr>
        <td id="L256" class="blob-num js-line-number" data-line-number="256"></td>
        <td id="LC256" class="blob-code blob-code-inner js-file-line">        tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;added &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> tab1_input_StrID.value <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;/b&gt;  to the list&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L257" class="blob-num js-line-number" data-line-number="257"></td>
        <td id="LC257" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">else</span>:</td>
      </tr>
      <tr>
        <td id="L258" class="blob-num js-line-number" data-line-number="258"></td>
        <td id="LC258" class="blob-code blob-code-inner js-file-line">        tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;&lt;b&gt;Warning: No time &amp; freq range selected.</span></td>
      </tr>
      <tr>
        <td id="L259" class="blob-num js-line-number" data-line-number="259"></td>
        <td id="LC259" class="blob-code blob-code-inner js-file-line"><span class="pl-s">                                Use the box select tool to select time &amp; freq</span></td>
      </tr>
      <tr>
        <td id="L260" class="blob-num js-line-number" data-line-number="260"></td>
        <td id="LC260" class="blob-code blob-code-inner js-file-line"><span class="pl-s">                                range in the dynamic spectrum first!!!&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L261" class="blob-num js-line-number" data-line-number="261"></td>
        <td id="LC261" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L262" class="blob-num js-line-number" data-line-number="262"></td>
        <td id="LC262" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L263" class="blob-num js-line-number" data-line-number="263"></td>
        <td id="LC263" class="blob-code blob-code-inner js-file-line">tab1_selected_StrID_entry <span class="pl-k">=</span> <span class="pl-c1">None</span></td>
      </tr>
      <tr>
        <td id="L264" class="blob-num js-line-number" data-line-number="264"></td>
        <td id="LC264" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L265" class="blob-num js-line-number" data-line-number="265"></td>
        <td id="LC265" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L266" class="blob-num js-line-number" data-line-number="266"></td>
        <td id="LC266" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_selection_StrID_entry</span>(<span class="pl-smi">attrname</span>, <span class="pl-smi">old</span>, <span class="pl-smi">new</span>):</td>
      </tr>
      <tr>
        <td id="L267" class="blob-num js-line-number" data-line-number="267"></td>
        <td id="LC267" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_selected_StrID_entry</td>
      </tr>
      <tr>
        <td id="L268" class="blob-num js-line-number" data-line-number="268"></td>
        <td id="LC268" class="blob-code blob-code-inner js-file-line">    tab1_selected_StrID_entry <span class="pl-k">=</span> tab1_SRC_StrIDPatch.selected[<span class="pl-s"><span class="pl-pds">&#39;</span>1d<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>indices<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L269" class="blob-num js-line-number" data-line-number="269"></td>
        <td id="LC269" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L270" class="blob-num js-line-number" data-line-number="270"></td>
        <td id="LC270" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L271" class="blob-num js-line-number" data-line-number="271"></td>
        <td id="LC271" class="blob-code blob-code-inner js-file-line">tab1_SRC_StrIDPatch.on_change(<span class="pl-s"><span class="pl-pds">&#39;</span>selected<span class="pl-pds">&#39;</span></span>, tab1_selection_StrID_entry)</td>
      </tr>
      <tr>
        <td id="L272" class="blob-num js-line-number" data-line-number="272"></td>
        <td id="LC272" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L273" class="blob-num js-line-number" data-line-number="273"></td>
        <td id="LC273" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L274" class="blob-num js-line-number" data-line-number="274"></td>
        <td id="LC274" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_deleteStrID</span>():</td>
      </tr>
      <tr>
        <td id="L275" class="blob-num js-line-number" data-line-number="275"></td>
        <td id="LC275" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_tim_square, tab1_selected_StrID_entry, database_dir, event_id</td>
      </tr>
      <tr>
        <td id="L276" class="blob-num js-line-number" data-line-number="276"></td>
        <td id="LC276" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_selected_StrID_entry:</td>
      </tr>
      <tr>
        <td id="L277" class="blob-num js-line-number" data-line-number="277"></td>
        <td id="LC277" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> pd.read_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L278" class="blob-num js-line-number" data-line-number="278"></td>
        <td id="LC278" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> StrIDList.sort_values(<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>, <span class="pl-v">ascending</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L279" class="blob-num js-line-number" data-line-number="279"></td>
        <td id="LC279" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> StrIDList.drop(StrIDList.index[tab1_selected_StrID_entry])</td>
      </tr>
      <tr>
        <td id="L280" class="blob-num js-line-number" data-line-number="280"></td>
        <td id="LC280" class="blob-code blob-code-inner js-file-line">        StrIDList.index <span class="pl-k">=</span> <span class="pl-c1">range</span>(StrIDList.index.size)</td>
      </tr>
      <tr>
        <td id="L281" class="blob-num js-line-number" data-line-number="281"></td>
        <td id="LC281" class="blob-code blob-code-inner js-file-line">        StrIDList.to_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L282" class="blob-num js-line-number" data-line-number="282"></td>
        <td id="LC282" class="blob-code blob-code-inner js-file-line">        StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>] <span class="pl-k">=</span> [ll <span class="pl-k">-</span> tab1_tim_square[<span class="pl-c1">0</span>] <span class="pl-k">for</span> ll <span class="pl-k">in</span> StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>]]</td>
      </tr>
      <tr>
        <td id="L283" class="blob-num js-line-number" data-line-number="283"></td>
        <td id="LC283" class="blob-code blob-code-inner js-file-line">        tab1_render_patch.data_source.data <span class="pl-k">=</span> ColumnDataSource(StrIDList).data</td>
      </tr>
      <tr>
        <td id="L284" class="blob-num js-line-number" data-line-number="284"></td>
        <td id="LC284" class="blob-code blob-code-inner js-file-line">        tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;removed &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> StrIDList.iloc[tab1_selected_StrID_entry[<span class="pl-c1">0</span>]][<span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>][</td>
      </tr>
      <tr>
        <td id="L285" class="blob-num js-line-number" data-line-number="285"></td>
        <td id="LC285" class="blob-code blob-code-inner js-file-line">            <span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;/b&gt; from the list&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L286" class="blob-num js-line-number" data-line-number="286"></td>
        <td id="LC286" class="blob-code blob-code-inner js-file-line">        tab1_selected_StrID_entry <span class="pl-k">=</span> <span class="pl-c1">None</span></td>
      </tr>
      <tr>
        <td id="L287" class="blob-num js-line-number" data-line-number="287"></td>
        <td id="LC287" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">else</span>:</td>
      </tr>
      <tr>
        <td id="L288" class="blob-num js-line-number" data-line-number="288"></td>
        <td id="LC288" class="blob-code blob-code-inner js-file-line">        tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;&lt;b&gt;Warning: No StrID selected. Select one StrID first!!!&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L289" class="blob-num js-line-number" data-line-number="289"></td>
        <td id="LC289" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L290" class="blob-num js-line-number" data-line-number="290"></td>
        <td id="LC290" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L291" class="blob-num js-line-number" data-line-number="291"></td>
        <td id="LC291" class="blob-code blob-code-inner js-file-line">port <span class="pl-k">=</span> <span class="pl-c1">5100</span></td>
      </tr>
      <tr>
        <td id="L292" class="blob-num js-line-number" data-line-number="292"></td>
        <td id="LC292" class="blob-code blob-code-inner js-file-line">ports <span class="pl-k">=</span> []</td>
      </tr>
      <tr>
        <td id="L293" class="blob-num js-line-number" data-line-number="293"></td>
        <td id="LC293" class="blob-code blob-code-inner js-file-line">ports.append(port)</td>
      </tr>
      <tr>
        <td id="L294" class="blob-num js-line-number" data-line-number="294"></td>
        <td id="LC294" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L295" class="blob-num js-line-number" data-line-number="295"></td>
        <td id="LC295" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L296" class="blob-num js-line-number" data-line-number="296"></td>
        <td id="LC296" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_FSviewStrID</span>():</td>
      </tr>
      <tr>
        <td id="L297" class="blob-num js-line-number" data-line-number="297"></td>
        <td id="LC297" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_selected_StrID_entry, tab1_dtim, tab1_specdata, dftmp, timestart, database_dir, event_id</td>
      </tr>
      <tr>
        <td id="L298" class="blob-num js-line-number" data-line-number="298"></td>
        <td id="LC298" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> port, ports, config_plot, tab1_tim_square</td>
      </tr>
      <tr>
        <td id="L299" class="blob-num js-line-number" data-line-number="299"></td>
        <td id="LC299" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> tab1_selected_StrID_entry:</td>
      </tr>
      <tr>
        <td id="L300" class="blob-num js-line-number" data-line-number="300"></td>
        <td id="LC300" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> pd.read_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L301" class="blob-num js-line-number" data-line-number="301"></td>
        <td id="LC301" class="blob-code blob-code-inner js-file-line">        StrIDList <span class="pl-k">=</span> StrIDList.sort_values(<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>, <span class="pl-v">ascending</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L302" class="blob-num js-line-number" data-line-number="302"></td>
        <td id="LC302" class="blob-code blob-code-inner js-file-line">        StrID <span class="pl-k">=</span> StrIDList.iloc[tab1_selected_StrID_entry[<span class="pl-c1">0</span>]]</td>
      </tr>
      <tr>
        <td id="L303" class="blob-num js-line-number" data-line-number="303"></td>
        <td id="LC303" class="blob-code blob-code-inner js-file-line">        out_json <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>.json<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L304" class="blob-num js-line-number" data-line-number="304"></td>
        <td id="LC304" class="blob-code blob-code-inner js-file-line">        StrID.to_json(out_json)</td>
      </tr>
      <tr>
        <td id="L305" class="blob-num js-line-number" data-line-number="305"></td>
        <td id="LC305" class="blob-code blob-code-inner js-file-line">        out_json <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>CurrFS.json<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L306" class="blob-num js-line-number" data-line-number="306"></td>
        <td id="LC306" class="blob-code blob-code-inner js-file-line">        struct_id <span class="pl-k">=</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>/<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L307" class="blob-num js-line-number" data-line-number="307"></td>
        <td id="LC307" class="blob-code blob-code-inner js-file-line">        <span class="pl-c1">FS_config</span> <span class="pl-k">=</span> {<span class="pl-s"><span class="pl-pds">&#39;</span>datadir<span class="pl-pds">&#39;</span></span>: {<span class="pl-s"><span class="pl-pds">&#39;</span>event_id<span class="pl-pds">&#39;</span></span>: event_id, <span class="pl-s"><span class="pl-pds">&#39;</span>struct_id<span class="pl-pds">&#39;</span></span>: struct_id,</td>
      </tr>
      <tr>
        <td id="L308" class="blob-num js-line-number" data-line-number="308"></td>
        <td id="LC308" class="blob-code blob-code-inner js-file-line">                                 <span class="pl-s"><span class="pl-pds">&#39;</span>FS_specfile<span class="pl-pds">&#39;</span></span>: database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> struct_id <span class="pl-k">+</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>_<span class="pl-pds">&#39;</span></span> <span class="pl-k">+</span></td>
      </tr>
      <tr>
        <td id="L309" class="blob-num js-line-number" data-line-number="309"></td>
        <td id="LC309" class="blob-code blob-code-inner js-file-line">                                                StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>date<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>T<span class="pl-pds">&#39;</span></span> <span class="pl-k">+</span> <span class="pl-c1">str</span>(StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>]).translate(<span class="pl-c1">None</span>,</td>
      </tr>
      <tr>
        <td id="L310" class="blob-num js-line-number" data-line-number="310"></td>
        <td id="LC310" class="blob-code blob-code-inner js-file-line">                                     <span class="pl-s"><span class="pl-pds">&#39;</span>:<span class="pl-pds">&#39;</span></span>) <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>.spec.npz<span class="pl-pds">&#39;</span></span>}}</td>
      </tr>
      <tr>
        <td id="L311" class="blob-num js-line-number" data-line-number="311"></td>
        <td id="LC311" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">with</span> <span class="pl-c1">open</span>(out_json, <span class="pl-s"><span class="pl-pds">&#39;</span>w<span class="pl-pds">&#39;</span></span>) <span class="pl-k">as</span> fp:</td>
      </tr>
      <tr>
        <td id="L312" class="blob-num js-line-number" data-line-number="312"></td>
        <td id="LC312" class="blob-code blob-code-inner js-file-line">            json.dump(<span class="pl-c1">FS_config</span>, fp)</td>
      </tr>
      <tr>
        <td id="L313" class="blob-num js-line-number" data-line-number="313"></td>
        <td id="LC313" class="blob-code blob-code-inner js-file-line">        in_json <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>CurrFS.json<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L314" class="blob-num js-line-number" data-line-number="314"></td>
        <td id="LC314" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">with</span> <span class="pl-c1">open</span>(in_json, <span class="pl-s"><span class="pl-pds">&#39;</span>r<span class="pl-pds">&#39;</span></span>) <span class="pl-k">as</span> fp:</td>
      </tr>
      <tr>
        <td id="L315" class="blob-num js-line-number" data-line-number="315"></td>
        <td id="LC315" class="blob-code blob-code-inner js-file-line">            <span class="pl-c1">FS_config</span> <span class="pl-k">=</span> json.load(fp)</td>
      </tr>
      <tr>
        <td id="L316" class="blob-num js-line-number" data-line-number="316"></td>
        <td id="LC316" class="blob-code blob-code-inner js-file-line">        <span class="pl-c1">FS_specfile</span> <span class="pl-k">=</span> <span class="pl-c1">FS_config</span>[<span class="pl-s"><span class="pl-pds">&#39;</span>datadir<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>FS_specfile<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L317" class="blob-num js-line-number" data-line-number="317"></td>
        <td id="LC317" class="blob-code blob-code-inner js-file-line">        <span class="pl-c1">FS_dspecDF</span> <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> struct_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>dspecDF-save<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L318" class="blob-num js-line-number" data-line-number="318"></td>
        <td id="LC318" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">if</span> os.path.exists(<span class="pl-c1">FS_specfile</span>):</td>
      </tr>
      <tr>
        <td id="L319" class="blob-num js-line-number" data-line-number="319"></td>
        <td id="LC319" class="blob-code blob-code-inner js-file-line">            <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>bokeh serve FSview --show --port <span class="pl-c1">{}</span> &amp;<span class="pl-pds">&#39;</span></span>.format(port)</td>
      </tr>
      <tr>
        <td id="L320" class="blob-num js-line-number" data-line-number="320"></td>
        <td id="LC320" class="blob-code blob-code-inner js-file-line">            os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>cd .. &amp; bokeh serve FSview --show --port <span class="pl-c1">{}</span> &amp;<span class="pl-pds">&#39;</span></span>.format(port))</td>
      </tr>
      <tr>
        <td id="L321" class="blob-num js-line-number" data-line-number="321"></td>
        <td id="LC321" class="blob-code blob-code-inner js-file-line">            port <span class="pl-k">+=</span> <span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L322" class="blob-num js-line-number" data-line-number="322"></td>
        <td id="LC322" class="blob-code blob-code-inner js-file-line">            ports.append(port)</td>
      </tr>
      <tr>
        <td id="L323" class="blob-num js-line-number" data-line-number="323"></td>
        <td id="LC323" class="blob-code blob-code-inner js-file-line">            <span class="pl-k">if</span> os.path.exists(<span class="pl-c1">FS_dspecDF</span>):</td>
      </tr>
      <tr>
        <td id="L324" class="blob-num js-line-number" data-line-number="324"></td>
        <td id="LC324" class="blob-code blob-code-inner js-file-line">                tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;sent StrID to &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> database_dir <span class="pl-k">+</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>str_id<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>] <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>.json&lt;/b&gt;&lt;/p&gt;</span></td>
      </tr>
      <tr>
        <td id="L325" class="blob-num js-line-number" data-line-number="325"></td>
        <td id="LC325" class="blob-code blob-code-inner js-file-line"><span class="pl-s">                    &lt;p&gt;sent FS_config to &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>CurrFS.json&lt;/b&gt;&lt;/p&gt;</span></td>
      </tr>
      <tr>
        <td id="L326" class="blob-num js-line-number" data-line-number="326"></td>
        <td id="LC326" class="blob-code blob-code-inner js-file-line"><span class="pl-s">                    &lt;p&gt;Check the &lt;b&gt;FS_view&lt;/b&gt; in the &lt;b&gt;new tab&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L327" class="blob-num js-line-number" data-line-number="327"></td>
        <td id="LC327" class="blob-code blob-code-inner js-file-line">            <span class="pl-k">else</span>:</td>
      </tr>
      <tr>
        <td id="L328" class="blob-num js-line-number" data-line-number="328"></td>
        <td id="LC328" class="blob-code blob-code-inner js-file-line">                tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;Check the &lt;b&gt;FS_clean &lt;/b&gt; in the &lt;b&gt;new tab&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L329" class="blob-num js-line-number" data-line-number="329"></td>
        <td id="LC329" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">else</span>:</td>
      </tr>
      <tr>
        <td id="L330" class="blob-num js-line-number" data-line-number="330"></td>
        <td id="LC330" class="blob-code blob-code-inner js-file-line">            time0, time1 <span class="pl-k">=</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>], StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L331" class="blob-num js-line-number" data-line-number="331"></td>
        <td id="LC331" class="blob-code blob-code-inner js-file-line">            freq0, freq1 <span class="pl-k">=</span> StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>][<span class="pl-c1">0</span>], StrID[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>][<span class="pl-k">-</span><span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L332" class="blob-num js-line-number" data-line-number="332"></td>
        <td id="LC332" class="blob-code blob-code-inner js-file-line">            bl <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>bl<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L333" class="blob-num js-line-number" data-line-number="333"></td>
        <td id="LC333" class="blob-code blob-code-inner js-file-line">            spec <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>spec<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L334" class="blob-num js-line-number" data-line-number="334"></td>
        <td id="LC334" class="blob-code blob-code-inner js-file-line">            npol <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>npol<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L335" class="blob-num js-line-number" data-line-number="335"></td>
        <td id="LC335" class="blob-code blob-code-inner js-file-line">            nbl <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>nbl<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L336" class="blob-num js-line-number" data-line-number="336"></td>
        <td id="LC336" class="blob-code blob-code-inner js-file-line">            ntim <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>ntim<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L337" class="blob-num js-line-number" data-line-number="337"></td>
        <td id="LC337" class="blob-code blob-code-inner js-file-line">            nfreq <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>nfreq<span class="pl-pds">&#39;</span></span>]</td>
      </tr>
      <tr>
        <td id="L338" class="blob-num js-line-number" data-line-number="338"></td>
        <td id="LC338" class="blob-code blob-code-inner js-file-line">            tim <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>tim<span class="pl-pds">&#39;</span></span>][:]</td>
      </tr>
      <tr>
        <td id="L339" class="blob-num js-line-number" data-line-number="339"></td>
        <td id="LC339" class="blob-code blob-code-inner js-file-line">            freq <span class="pl-k">=</span> tab1_specdata[<span class="pl-s"><span class="pl-pds">&#39;</span>freq<span class="pl-pds">&#39;</span></span>] <span class="pl-k">/</span> <span class="pl-c1">1e9</span></td>
      </tr>
      <tr>
        <td id="L340" class="blob-num js-line-number" data-line-number="340"></td>
        <td id="LC340" class="blob-code blob-code-inner js-file-line">            timeidx0 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(ntim) <span class="pl-k">if</span> tim[i] <span class="pl-k">&gt;=</span> time0)</td>
      </tr>
      <tr>
        <td id="L341" class="blob-num js-line-number" data-line-number="341"></td>
        <td id="LC341" class="blob-code blob-code-inner js-file-line">            timeidx1 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(ntim <span class="pl-k">-</span> <span class="pl-c1">1</span>, <span class="pl-k">-</span><span class="pl-c1">1</span>, <span class="pl-k">-</span><span class="pl-c1">1</span>) <span class="pl-k">if</span> tim[i] <span class="pl-k">&lt;=</span> time1) <span class="pl-k">+</span> <span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L342" class="blob-num js-line-number" data-line-number="342"></td>
        <td id="LC342" class="blob-code blob-code-inner js-file-line">            freqidx0 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(nfreq) <span class="pl-k">if</span> freq[i] <span class="pl-k">&gt;=</span> freq0)</td>
      </tr>
      <tr>
        <td id="L343" class="blob-num js-line-number" data-line-number="343"></td>
        <td id="LC343" class="blob-code blob-code-inner js-file-line">            freqidx1 <span class="pl-k">=</span> <span class="pl-c1">next</span>(i <span class="pl-k">for</span> i <span class="pl-k">in</span> <span class="pl-v">xrange</span>(nfreq <span class="pl-k">-</span> <span class="pl-c1">1</span>, <span class="pl-k">-</span><span class="pl-c1">1</span>, <span class="pl-k">-</span><span class="pl-c1">1</span>) <span class="pl-k">if</span> freq[i] <span class="pl-k">&lt;=</span> freq1) <span class="pl-k">+</span> <span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L344" class="blob-num js-line-number" data-line-number="344"></td>
        <td id="LC344" class="blob-code blob-code-inner js-file-line">            spec <span class="pl-k">=</span> spec[:, :, freqidx0:(freqidx1 <span class="pl-k">+</span> <span class="pl-c1">1</span>), timeidx0:(timeidx1 <span class="pl-k">+</span> <span class="pl-c1">1</span>)]</td>
      </tr>
      <tr>
        <td id="L345" class="blob-num js-line-number" data-line-number="345"></td>
        <td id="LC345" class="blob-code blob-code-inner js-file-line">            tim <span class="pl-k">=</span> tim[timeidx0:(timeidx1 <span class="pl-k">+</span> <span class="pl-c1">1</span>)]</td>
      </tr>
      <tr>
        <td id="L346" class="blob-num js-line-number" data-line-number="346"></td>
        <td id="LC346" class="blob-code blob-code-inner js-file-line">            freq <span class="pl-k">=</span> freq[freqidx0:(freqidx1 <span class="pl-k">+</span> <span class="pl-c1">1</span>)] <span class="pl-k">*</span> <span class="pl-c1">1.0e9</span></td>
      </tr>
      <tr>
        <td id="L347" class="blob-num js-line-number" data-line-number="347"></td>
        <td id="LC347" class="blob-code blob-code-inner js-file-line">            ntim <span class="pl-k">=</span> <span class="pl-c1">len</span>(tim)</td>
      </tr>
      <tr>
        <td id="L348" class="blob-num js-line-number" data-line-number="348"></td>
        <td id="LC348" class="blob-code blob-code-inner js-file-line">            nfreq <span class="pl-k">=</span> <span class="pl-c1">len</span>(freq)</td>
      </tr>
      <tr>
        <td id="L349" class="blob-num js-line-number" data-line-number="349"></td>
        <td id="LC349" class="blob-code blob-code-inner js-file-line">            struct_id_dir <span class="pl-k">=</span> database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> struct_id</td>
      </tr>
      <tr>
        <td id="L350" class="blob-num js-line-number" data-line-number="350"></td>
        <td id="LC350" class="blob-code blob-code-inner js-file-line">            <span class="pl-k">if</span> <span class="pl-k">not</span> os.path.exists(struct_id_dir):</td>
      </tr>
      <tr>
        <td id="L351" class="blob-num js-line-number" data-line-number="351"></td>
        <td id="LC351" class="blob-code blob-code-inner js-file-line">                os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>mkdir <span class="pl-c1">{}</span><span class="pl-pds">&#39;</span></span>.format(struct_id_dir))</td>
      </tr>
      <tr>
        <td id="L352" class="blob-num js-line-number" data-line-number="352"></td>
        <td id="LC352" class="blob-code blob-code-inner js-file-line">            np.savez(<span class="pl-c1">FS_specfile</span>, <span class="pl-v">spec</span><span class="pl-k">=</span>spec, <span class="pl-v">tim</span><span class="pl-k">=</span>tim, <span class="pl-v">freq</span><span class="pl-k">=</span>freq, <span class="pl-v">bl</span><span class="pl-k">=</span>bl, <span class="pl-v">npol</span><span class="pl-k">=</span>npol, <span class="pl-v">nbl</span><span class="pl-k">=</span>nbl, <span class="pl-v">nfreq</span><span class="pl-k">=</span>nfreq, <span class="pl-v">ntim</span><span class="pl-k">=</span>ntim)</td>
      </tr>
      <tr>
        <td id="L353" class="blob-num js-line-number" data-line-number="353"></td>
        <td id="LC353" class="blob-code blob-code-inner js-file-line">            <span class="pl-c"># todo save a full resolution dspec within the selected time/freq range</span></td>
      </tr>
      <tr>
        <td id="L354" class="blob-num js-line-number" data-line-number="354"></td>
        <td id="LC354" class="blob-code blob-code-inner js-file-line">            tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;&lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> <span class="pl-c1">FS_specfile</span> <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;/b&gt; saved &gt;&gt;&gt;&gt;&gt;&gt; Click the &lt;b&gt;FS veiw button&lt;/b&gt; again to make aperture synthesis images&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L355" class="blob-num js-line-number" data-line-number="355"></td>
        <td id="LC355" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">else</span>:</td>
      </tr>
      <tr>
        <td id="L356" class="blob-num js-line-number" data-line-number="356"></td>
        <td id="LC356" class="blob-code blob-code-inner js-file-line">        tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;&lt;b&gt;Warning: No StrID selected. Select one StrID first!!!&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L357" class="blob-num js-line-number" data-line-number="357"></td>
        <td id="LC357" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L358" class="blob-num js-line-number" data-line-number="358"></td>
        <td id="LC358" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L359" class="blob-num js-line-number" data-line-number="359"></td>
        <td id="LC359" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_saveStrID</span>():</td>
      </tr>
      <tr>
        <td id="L360" class="blob-num js-line-number" data-line-number="360"></td>
        <td id="LC360" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> database_dir, event_id</td>
      </tr>
      <tr>
        <td id="L361" class="blob-num js-line-number" data-line-number="361"></td>
        <td id="LC361" class="blob-code blob-code-inner js-file-line">    os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>cp <span class="pl-c1">{}</span>StrID_list_tmp.json <span class="pl-c1">{}</span>StrID_list.json<span class="pl-pds">&#39;</span></span>.format(database_dir <span class="pl-k">+</span> event_id, database_dir <span class="pl-k">+</span> event_id))</td>
      </tr>
      <tr>
        <td id="L362" class="blob-num js-line-number" data-line-number="362"></td>
        <td id="LC362" class="blob-code blob-code-inner js-file-line">    tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;StrID data saved to &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span><span class="pl-c1">{}</span>StrID_list.json&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&#39;</span></span>.format(database_dir <span class="pl-k">+</span> event_id)</td>
      </tr>
      <tr>
        <td id="L363" class="blob-num js-line-number" data-line-number="363"></td>
        <td id="LC363" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L364" class="blob-num js-line-number" data-line-number="364"></td>
        <td id="LC364" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L365" class="blob-num js-line-number" data-line-number="365"></td>
        <td id="LC365" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_update_reloadStrID</span>():</td>
      </tr>
      <tr>
        <td id="L366" class="blob-num js-line-number" data-line-number="366"></td>
        <td id="LC366" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">global</span> tab1_tim_square, database_dir, event_id</td>
      </tr>
      <tr>
        <td id="L367" class="blob-num js-line-number" data-line-number="367"></td>
        <td id="LC367" class="blob-code blob-code-inner js-file-line">    os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>cp <span class="pl-c1">{}</span>StrID_list.json <span class="pl-c1">{}</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>.format(database_dir <span class="pl-k">+</span> event_id, database_dir <span class="pl-k">+</span> event_id))</td>
      </tr>
      <tr>
        <td id="L368" class="blob-num js-line-number" data-line-number="368"></td>
        <td id="LC368" class="blob-code blob-code-inner js-file-line">    StrIDList <span class="pl-k">=</span> pd.read_json(database_dir <span class="pl-k">+</span> event_id <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span>StrID_list_tmp.json<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L369" class="blob-num js-line-number" data-line-number="369"></td>
        <td id="LC369" class="blob-code blob-code-inner js-file-line">    StrIDList <span class="pl-k">=</span> StrIDList.sort_values(<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>timeran<span class="pl-pds">&#39;</span></span>, <span class="pl-v">ascending</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L370" class="blob-num js-line-number" data-line-number="370"></td>
        <td id="LC370" class="blob-code blob-code-inner js-file-line">    StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>] <span class="pl-k">=</span> [ll <span class="pl-k">-</span> tab1_tim_square[<span class="pl-c1">0</span>] <span class="pl-k">for</span> ll <span class="pl-k">in</span> StrIDList[<span class="pl-s"><span class="pl-pds">&#39;</span>time<span class="pl-pds">&#39;</span></span>]]</td>
      </tr>
      <tr>
        <td id="L371" class="blob-num js-line-number" data-line-number="371"></td>
        <td id="LC371" class="blob-code blob-code-inner js-file-line">    tab1_render_patch.data_source.data <span class="pl-k">=</span> ColumnDataSource(StrIDList).data</td>
      </tr>
      <tr>
        <td id="L372" class="blob-num js-line-number" data-line-number="372"></td>
        <td id="LC372" class="blob-code blob-code-inner js-file-line">    tab1_Div_Tb.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;StrID data reloaded from &lt;b&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span> <span class="pl-k">+</span> <span class="pl-s"><span class="pl-pds">&#39;</span><span class="pl-c1">{}</span>StrID_list.json&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&#39;</span></span>.format(</td>
      </tr>
      <tr>
        <td id="L373" class="blob-num js-line-number" data-line-number="373"></td>
        <td id="LC373" class="blob-code blob-code-inner js-file-line">        database_dir <span class="pl-k">+</span> event_id)</td>
      </tr>
      <tr>
        <td id="L374" class="blob-num js-line-number" data-line-number="374"></td>
        <td id="LC374" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L375" class="blob-num js-line-number" data-line-number="375"></td>
        <td id="LC375" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L376" class="blob-num js-line-number" data-line-number="376"></td>
        <td id="LC376" class="blob-code blob-code-inner js-file-line">tab1_BUT_OPT <span class="pl-k">=</span> <span class="pl-c1">dict</span>(<span class="pl-v">width</span><span class="pl-k">=</span>config_plot[<span class="pl-s"><span class="pl-pds">&#39;</span>plot_config<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>tab_QLook<span class="pl-pds">&#39;</span></span>][<span class="pl-s"><span class="pl-pds">&#39;</span>StrID_DataTb_BUT_wdth<span class="pl-pds">&#39;</span></span>])</td>
      </tr>
      <tr>
        <td id="L377" class="blob-num js-line-number" data-line-number="377"></td>
        <td id="LC377" class="blob-code blob-code-inner js-file-line">tab1_BUT_addStrID <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>Add to StrID<span class="pl-pds">&#39;</span></span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>success<span class="pl-pds">&#39;</span></span>, <span class="pl-k">**</span>tab1_BUT_OPT)</td>
      </tr>
      <tr>
        <td id="L378" class="blob-num js-line-number" data-line-number="378"></td>
        <td id="LC378" class="blob-code blob-code-inner js-file-line">tab1_BUT_addStrID.on_click(tab1_update_addStrID)</td>
      </tr>
      <tr>
        <td id="L379" class="blob-num js-line-number" data-line-number="379"></td>
        <td id="LC379" class="blob-code blob-code-inner js-file-line">tab1_BUT_deleteStrID <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>Delete StrID<span class="pl-pds">&#39;</span></span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>danger<span class="pl-pds">&#39;</span></span>, <span class="pl-k">**</span>tab1_BUT_OPT)</td>
      </tr>
      <tr>
        <td id="L380" class="blob-num js-line-number" data-line-number="380"></td>
        <td id="LC380" class="blob-code blob-code-inner js-file-line">tab1_BUT_deleteStrID.on_click(tab1_update_deleteStrID)</td>
      </tr>
      <tr>
        <td id="L381" class="blob-num js-line-number" data-line-number="381"></td>
        <td id="LC381" class="blob-code blob-code-inner js-file-line">tab1_BUT_saveStrID <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>Save StrID<span class="pl-pds">&#39;</span></span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>primary<span class="pl-pds">&#39;</span></span>, <span class="pl-k">**</span>tab1_BUT_OPT)</td>
      </tr>
      <tr>
        <td id="L382" class="blob-num js-line-number" data-line-number="382"></td>
        <td id="LC382" class="blob-code blob-code-inner js-file-line">tab1_BUT_saveStrID.on_click(tab1_update_saveStrID)</td>
      </tr>
      <tr>
        <td id="L383" class="blob-num js-line-number" data-line-number="383"></td>
        <td id="LC383" class="blob-code blob-code-inner js-file-line">tab1_BUT_reloadStrID <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>Reload StrID<span class="pl-pds">&#39;</span></span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>warning<span class="pl-pds">&#39;</span></span>, <span class="pl-k">**</span>tab1_BUT_OPT)</td>
      </tr>
      <tr>
        <td id="L384" class="blob-num js-line-number" data-line-number="384"></td>
        <td id="LC384" class="blob-code blob-code-inner js-file-line">tab1_BUT_reloadStrID.on_click(tab1_update_reloadStrID)</td>
      </tr>
      <tr>
        <td id="L385" class="blob-num js-line-number" data-line-number="385"></td>
        <td id="LC385" class="blob-code blob-code-inner js-file-line">tab1_BUT_FSviewStrID <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>FS View<span class="pl-pds">&#39;</span></span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>primary<span class="pl-pds">&#39;</span></span>, <span class="pl-k">**</span>tab1_BUT_OPT)</td>
      </tr>
      <tr>
        <td id="L386" class="blob-num js-line-number" data-line-number="386"></td>
        <td id="LC386" class="blob-code blob-code-inner js-file-line">tab1_BUT_FSviewStrID.on_click(tab1_update_FSviewStrID)</td>
      </tr>
      <tr>
        <td id="L387" class="blob-num js-line-number" data-line-number="387"></td>
        <td id="LC387" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L388" class="blob-num js-line-number" data-line-number="388"></td>
        <td id="LC388" class="blob-code blob-code-inner js-file-line">tab1_SPCR_LFT_DataTb_dspec <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">10</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">100</span>)</td>
      </tr>
      <tr>
        <td id="L389" class="blob-num js-line-number" data-line-number="389"></td>
        <td id="LC389" class="blob-code blob-code-inner js-file-line">tab1_SPCR_RGT_DataTb_dspec <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">20</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">100</span>)</td>
      </tr>
      <tr>
        <td id="L390" class="blob-num js-line-number" data-line-number="390"></td>
        <td id="LC390" class="blob-code blob-code-inner js-file-line">tab1_SPCR_ABV_DataTb_dspec <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">100</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">18</span>)</td>
      </tr>
      <tr>
        <td id="L391" class="blob-num js-line-number" data-line-number="391"></td>
        <td id="LC391" class="blob-code blob-code-inner js-file-line">tab1_SPCR_LFT_But <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">10</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">25</span>)</td>
      </tr>
      <tr>
        <td id="L392" class="blob-num js-line-number" data-line-number="392"></td>
        <td id="LC392" class="blob-code blob-code-inner js-file-line">tab1_SPCR_LFT_DataTb_evt <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">20</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">15</span>)</td>
      </tr>
      <tr>
        <td id="L393" class="blob-num js-line-number" data-line-number="393"></td>
        <td id="LC393" class="blob-code blob-code-inner js-file-line">tab1_SPCR_ABV_DataTb_evt <span class="pl-k">=</span> Spacer(<span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">100</span>, <span class="pl-v">height</span><span class="pl-k">=</span><span class="pl-c1">18</span>)</td>
      </tr>
      <tr>
        <td id="L394" class="blob-num js-line-number" data-line-number="394"></td>
        <td id="LC394" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L395" class="blob-num js-line-number" data-line-number="395"></td>
        <td id="LC395" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L396" class="blob-num js-line-number" data-line-number="396"></td>
        <td id="LC396" class="blob-code blob-code-inner js-file-line"><span class="pl-k">def</span> <span class="pl-en">tab1_exit</span>():</td>
      </tr>
      <tr>
        <td id="L397" class="blob-num js-line-number" data-line-number="397"></td>
        <td id="LC397" class="blob-code blob-code-inner js-file-line">    tab1_Div_exit.text <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;&quot;&quot;</span>&lt;p&gt;&lt;b&gt;You may close the tab anytime you like.&lt;/b&gt;&lt;/p&gt;<span class="pl-pds">&quot;&quot;&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L398" class="blob-num js-line-number" data-line-number="398"></td>
        <td id="LC398" class="blob-code blob-code-inner js-file-line">    <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>You may close the tab anytime you like.<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L399" class="blob-num js-line-number" data-line-number="399"></td>
        <td id="LC399" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> ll <span class="pl-k">in</span> ports:</td>
      </tr>
      <tr>
        <td id="L400" class="blob-num js-line-number" data-line-number="400"></td>
        <td id="LC400" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">if</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>linux<span class="pl-pds">&quot;</span></span> <span class="pl-k">or</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>linux2<span class="pl-pds">&quot;</span></span>:</td>
      </tr>
      <tr>
        <td id="L401" class="blob-num js-line-number" data-line-number="401"></td>
        <td id="LC401" class="blob-code blob-code-inner js-file-line">            os.system(<span class="pl-s"><span class="pl-pds">&#39;</span>fuser -n tcp -k <span class="pl-c1">{}</span><span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L402" class="blob-num js-line-number" data-line-number="402"></td>
        <td id="LC402" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">elif</span> platform <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>darwin<span class="pl-pds">&quot;</span></span>:</td>
      </tr>
      <tr>
        <td id="L403" class="blob-num js-line-number" data-line-number="403"></td>
        <td id="LC403" class="blob-code blob-code-inner js-file-line">            os.system(</td>
      </tr>
      <tr>
        <td id="L404" class="blob-num js-line-number" data-line-number="404"></td>
        <td id="LC404" class="blob-code blob-code-inner js-file-line">                <span class="pl-s"><span class="pl-pds">&#39;</span>port=($(lsof -i tcp:<span class="pl-c1">{}</span>|grep python2.7 |cut -f2 -d&quot; &quot;)); [[ -n &quot;$port&quot; ]] &amp;&amp; kill -9 $port<span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L405" class="blob-num js-line-number" data-line-number="405"></td>
        <td id="LC405" class="blob-code blob-code-inner js-file-line">            os.system(</td>
      </tr>
      <tr>
        <td id="L406" class="blob-num js-line-number" data-line-number="406"></td>
        <td id="LC406" class="blob-code blob-code-inner js-file-line">                <span class="pl-s"><span class="pl-pds">&#39;</span>port=($(lsof -i tcp:<span class="pl-c1">{}</span>|grep Google |cut -f2 -d&quot; &quot;)); [[ -n &quot;$port&quot; ]] &amp;&amp; kill -9 $port<span class="pl-pds">&#39;</span></span>.format(ll))</td>
      </tr>
      <tr>
        <td id="L407" class="blob-num js-line-number" data-line-number="407"></td>
        <td id="LC407" class="blob-code blob-code-inner js-file-line">        <span class="pl-c1">print</span> <span class="pl-s"><span class="pl-pds">&#39;</span>port <span class="pl-c1">{}</span> killed<span class="pl-pds">&#39;</span></span>.format(ll)</td>
      </tr>
      <tr>
        <td id="L408" class="blob-num js-line-number" data-line-number="408"></td>
        <td id="LC408" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">raise</span> <span class="pl-c1">SystemExit</span></td>
      </tr>
      <tr>
        <td id="L409" class="blob-num js-line-number" data-line-number="409"></td>
        <td id="LC409" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L410" class="blob-num js-line-number" data-line-number="410"></td>
        <td id="LC410" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L411" class="blob-num js-line-number" data-line-number="411"></td>
        <td id="LC411" class="blob-code blob-code-inner js-file-line">tab1_BUT_exit <span class="pl-k">=</span> Button(<span class="pl-v">label</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>Exit QLook<span class="pl-pds">&#39;</span></span>, <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>, <span class="pl-v">button_type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>danger<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L412" class="blob-num js-line-number" data-line-number="412"></td>
        <td id="LC412" class="blob-code blob-code-inner js-file-line">tab1_BUT_exit.on_click(tab1_exit)</td>
      </tr>
      <tr>
        <td id="L413" class="blob-num js-line-number" data-line-number="413"></td>
        <td id="LC413" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L414" class="blob-num js-line-number" data-line-number="414"></td>
        <td id="LC414" class="blob-code blob-code-inner js-file-line">panel1 <span class="pl-k">=</span> column(tab1_p_dspec,</td>
      </tr>
      <tr>
        <td id="L415" class="blob-num js-line-number" data-line-number="415"></td>
        <td id="LC415" class="blob-code blob-code-inner js-file-line">    row(widgetbox(tab1_Select_bl, tab1_Select_pol, tab1_Select_colormap, tab1_BUT_exit, tab1_Div_exit, <span class="pl-v">width</span><span class="pl-k">=</span><span class="pl-c1">150</span>),</td>
      </tr>
      <tr>
        <td id="L416" class="blob-num js-line-number" data-line-number="416"></td>
        <td id="LC416" class="blob-code blob-code-inner js-file-line">        tab1_SPCR_LFT_DataTb_evt, tab1_SPCR_LFT_DataTb_dspec, column(tab1_DataTb_dspec, tab1_Div_Tb), tab1_SPCR_LFT_But,</td>
      </tr>
      <tr>
        <td id="L417" class="blob-num js-line-number" data-line-number="417"></td>
        <td id="LC417" class="blob-code blob-code-inner js-file-line">        widgetbox(tab1_BUT_FSviewStrID, tab1_input_StrID, tab1_BUT_addStrID, tab1_BUT_deleteStrID, tab1_BUT_saveStrID,</td>
      </tr>
      <tr>
        <td id="L418" class="blob-num js-line-number" data-line-number="418"></td>
        <td id="LC418" class="blob-code blob-code-inner js-file-line">            tab1_BUT_reloadStrID)))</td>
      </tr>
      <tr>
        <td id="L419" class="blob-num js-line-number" data-line-number="419"></td>
        <td id="LC419" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L420" class="blob-num js-line-number" data-line-number="420"></td>
        <td id="LC420" class="blob-code blob-code-inner js-file-line"><span class="pl-c1">print</span>(<span class="pl-s"><span class="pl-pds">&quot;</span>--- <span class="pl-c1">%s</span> seconds ---<span class="pl-pds">&quot;</span></span> <span class="pl-k">%</span> (time.time() <span class="pl-k">-</span> start_timestamp))</td>
      </tr>
      <tr>
        <td id="L421" class="blob-num js-line-number" data-line-number="421"></td>
        <td id="LC421" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L422" class="blob-num js-line-number" data-line-number="422"></td>
        <td id="LC422" class="blob-code blob-code-inner js-file-line">lout <span class="pl-k">=</span> panel1</td>
      </tr>
      <tr>
        <td id="L423" class="blob-num js-line-number" data-line-number="423"></td>
        <td id="LC423" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L424" class="blob-num js-line-number" data-line-number="424"></td>
        <td id="LC424" class="blob-code blob-code-inner js-file-line">curdoc().add_root(lout)</td>
      </tr>
      <tr>
        <td id="L425" class="blob-num js-line-number" data-line-number="425"></td>
        <td id="LC425" class="blob-code blob-code-inner js-file-line">curdoc().title <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>QLook tool<span class="pl-pds">&quot;</span></span></td>
      </tr>
</table>

  </div>

</div>

<button type="button" data-facebox="#jump-to-line" data-facebox-class="linejump" data-hotkey="l" class="d-none">Jump to Line</button>
<div id="jump-to-line" style="display:none">
  <!-- '"` --><!-- </textarea></xmp> --></option></form><form accept-charset="UTF-8" action="" class="js-jump-to-line-form" method="get"><div style="margin:0;padding:0;display:inline"><input name="utf8" type="hidden" value="&#x2713;" /></div>
    <input class="form-control linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
    <button type="submit" class="btn">Go</button>
</form></div>

  </div>
  <div class="modal-backdrop js-touch-events"></div>
</div>


    </div>
  </div>

    </div>

        <div class="container site-footer-container">
  <div class="site-footer" role="contentinfo">
    <ul class="site-footer-links float-right">
        <li><a href="https://github.com/contact" data-ga-click="Footer, go to contact, text:contact">Contact GitHub</a></li>
      <li><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
      <li><a href="https://shop.github.com" data-ga-click="Footer, go to shop, text:shop">Shop</a></li>
        <li><a href="https://github.com/blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a href="https://github.com/about" data-ga-click="Footer, go to about, text:about">About</a></li>

    </ul>

    <a href="https://github.com" aria-label="Homepage" class="site-footer-mark" title="GitHub">
      <svg aria-hidden="true" class="octicon octicon-mark-github" height="24" version="1.1" viewBox="0 0 16 16" width="24"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
    <ul class="site-footer-links">
      <li>&copy; 2016 <span title="0.05853s from github-fe141-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="https://github.com/site/terms" data-ga-click="Footer, go to terms, text:terms">Terms</a></li>
        <li><a href="https://github.com/site/privacy" data-ga-click="Footer, go to privacy, text:privacy">Privacy</a></li>
        <li><a href="https://github.com/security" data-ga-click="Footer, go to security, text:security">Security</a></li>
        <li><a href="https://status.github.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
        <li><a href="https://help.github.com" data-ga-click="Footer, go to help, text:help">Help</a></li>
    </ul>
  </div>
</div>



    

    <div id="ajax-error-message" class="ajax-error-message flash flash-error">
      <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"></path></svg>
      <button type="button" class="flash-close js-flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
        <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"></path></svg>
      </button>
      You can't perform that action at this time.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/compat-30ce4c86c27f88c3d1b4eb03efda59b45d8d7c871880dee0b8f73d5ef1b25fdf.js"></script>
      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-35b9c541d341309b930b1c790fa1b27a30c7c44ce10c7f8242890e3d83c8adbd.js"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-294739b19304855d2c5077cf4690f1ec89bce7aa73a88859b76dc1f22316a3d8.js"></script>
      
      
      
      
    <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner d-none">
      <svg aria-hidden="true" class="octicon octicon-alert" height="16" version="1.1" viewBox="0 0 16 16" width="16"><path d="M8.865 1.52c-.18-.31-.51-.5-.87-.5s-.69.19-.87.5L.275 13.5c-.18.31-.18.69 0 1 .19.31.52.5.87.5h13.7c.36 0 .69-.19.86-.5.17-.31.18-.69.01-1L8.865 1.52zM8.995 13h-2v-2h2v2zm0-3h-2V6h2v4z"></path></svg>
      <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
      <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
    </div>
    <div class="facebox" id="facebox" style="display:none;">
  <div class="facebox-popup">
    <div class="facebox-content" role="dialog" aria-labelledby="facebox-header" aria-describedby="facebox-description">
    </div>
    <button type="button" class="facebox-close js-facebox-close" aria-label="Close modal">
      <svg aria-hidden="true" class="octicon octicon-x" height="16" version="1.1" viewBox="0 0 12 16" width="12"><path d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48z"></path></svg>
    </button>
  </div>
</div>

  </body>
</html>

