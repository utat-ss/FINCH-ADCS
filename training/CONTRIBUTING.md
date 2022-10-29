# Contributing

`ACS Toolbox` is largely a private effort by experienced Guidance, Navigation and Control (GNC)engineers to make low-precision GNC tools available for junior engineers, academia and the general public.
The development and maintenance is sponsored by Astris Aerospace Inc.
However, we are always grateful for bug reporting and suggestions new features and improvements.

## Table of contents
1. [Ways To Contribute](#ways-to-contribute)
   1. [Report Bugs and Provide Suggestions](#report-bugs-and-provide-suggestions)
   1. [Contribute Bug Fixes and New Features](#contribute-bug-fixes-and-new-features)
1. [Contributing Using Pull Requests](#contributing-using-pull-requests)
   1. [Fork, Branch, Edit and Push](#fork-branch-edit-and-push)
1. [Variable Naming Conventions](#internal-variable-naming-conventions)

## Ways To Contribute

### Report Bugs and Provide Suggestions

Unit tests are implemented ubiquitously. However, new use cases have the possibility of breaking certain utilities.
If you have suggestions for improvement or find something that doesn't work as intended, please open a new ticket in the [issue tracker](https://github.com/ex-astris/acstoolbox/issues).

### Contribute Bug Fixes and New Features

If you have already implemented suggested improvements or made fixes to encountered bugs, please contribute your code by pushing a pull request (PR) to the respository.
One of the developers will review your code, suggest feedback and merge it into the codebase if possible.
See below on how to open a new PR.

<!--TODO: Add Documentation Contributions + TOC-->

## Contributing Using Pull Requests

The best way to improve and expand the `ACS Toolbox` is to push PRs to the repository which will be checked and merged by one of our developers. 
GitHub allows you to make minor source code changes and push a PR directly from your browser and without additional software.
The online [Introduction to Github](https://lab.github.com/githubtraining/introduction-to-github) demonstrates how to do this.
For more complex code changes, please fork the codebase and push a new branch to the `ACS Toolbox` respository.
This can either be done on the command line interface (CLI) by experienced developers or on the [desktop Github application].
Only the latter is outlined below.

<!--Insert CLI instructions-->

### Fork, Branch, Edit and Push

The following steps outline how to push a new PR to the `ACS Toolbox` repository for review and approval. 
Your system will require the [desktop Github application](https://desktop.github.com/) to be installed on your machine and a consistent internet connection.

1. Fork the `ACS Toolbox` repository from your browser while logged into your Github account.
2. Give the forked repository a name and description and click `Create fork` in your browser.
3. In the Github desktop application, select the forked repository under the `Current repository` dropdown menu to perform a local checkout.
4. In the desktop application, select `New branch` under the `Current branch` dropdown menu.
5. Give your branch an appropriate name and click `Create branch`.

<p class="aligncenter">
    <img src="./images/root/1-GitFork.png" alt="centered image" />
</p>

6. Code your feature/bug fix and in the Github desktop application, commit your code to your branch while providing an appropriate description.
7. In the desktop application, select `Publish branch` to push your new code to your online repository

<p class="aligncenter">
    <img src="./images/root/2-GitCommit.png"/>
</p>

8. In your browser, find the branch you just uploaded and select `Compare & pull request` to start a new PR. 
9. Enter a PR title and detailed description. Provide motivation, any references and results. Select `Create pull request` to initiate the continuous integration (CI) pipeline.

<p class="aligncenter">
    <img src="./images/root/3-GitPR.png"/>
</p>

10. Code maintainers will automatically be selected as reviewers and unit tests will be run as part of CI.
If a single test fails, select `Details` to understand why and push an updated commit to the PR.

<p class="aligncenter">
    <img src="./images/root/4-GitCI.png"/>
</p>

<!--TODO: ### Branch Naming Conventions + TOC-->

## Internal Variable Naming Conventions

 This Section requires a [MathJax](https://www.mathjax.org/) browser plugin for the mathematicsto be rendered in between then `$...$`. These are widely available for all the major browsers.

 The `ACS Toolbox` contains algorithms from various branches of estimation, control and astrodynamics theory. The mathematics and notation used in the repository are adapted from vectrix notation developed by Peter Hughes (see [Spacecraft Attitude Dynamics](https://www.amazon.com/Spacecraft-Attitude-Dynamics-Aeronautical-Engineering/dp/0486439259)). Likewise, variable naming conventions must be consistent with the existing code. The following outlines vectrix notation and the corresponding variable naming convention for the `ACS Toolbox`.

 1. Rotation transform from frame $\underrightarrow{\mathcal{F}}_i$ to $\underrightarrow{\mathcal{F}}_b$ ($C_{bi}$)
    - `C_bi`
 1. Satellite position expressed in frame $\underrightarrow{\mathcal{F}}_x$ ($\mathbf{r}_x$)
    - `r_x`
1. Magnitude of satellite position ($r = \|\mathbf{r}_x\|$)
    - `r`
1. Unit vector of Satellite position expressed in frame $\underrightarrow{\mathcal{F}}_x$ ($\hat{\mathbf{r}}_x$)
    - `ur_x`
1. Centuries expressed in frame `TDB` ($T_{TDB}$)
   - `c_tdb`
1. Julian Date in `UT1` ($JD_{UT1}$)
   - `jd_ut1`
1. Julian seconds from J2000 in `TT` ($t_{j2000,TT}$)
   - `jsj2000_tt`
1. Gregorian time (array of yyyy-month-dd-hh-mm-ss) expressed in frame `UTC`
    - `gregorian_utc`

 <!--TODO: ### Move this into a notations file + remove from TOC-->