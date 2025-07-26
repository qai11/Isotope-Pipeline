# Reduction and Isotope pipeline for High resolution Spectra



<!-- ## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)! -->


<!-- ## Integrate with your tools

- [ ] [Set up project integrations](https://eng-git.canterbury.ac.nz/mwega/quin-masters-code/-/settings/integrations) -->

<!-- ## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html) -->

<!-- ## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing (SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html) -->

***

<!-- # Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thanks to [makeareadme.com](https://www.makeareadme.com/) for this template. -->
<!-- 
## Suggestions for a good README

Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information. -->


## Name

Reduction and Isotope pipeline for High resolution Spectra.


## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

The analysis of Magnesium Isotopes in M,F and G stars. 

#### Continuum Normilization and Abundances 
Contains multiple scripts for taking first pass continuum fit data from `Deinira.py` redution pipeline for HERCULES data and performing post redtuction passes along with finding and checking abundances using iSpec.

File: `parameters_pipeline.py` is the main running script for these files. It will first call `rv_combine.py` to perform any barrycentric correction and proper motion adjustment as needed, then `continuum_adjust.py` to correct the normilization of the spectrum, it will then run `find_params.py` to find the stellar parameters for each given star and finally find the abundance for given elements.

All scripts use iSpec and require the ges linelist for finding abundances from the ges_v2 paper.

An extra scipt can interpolatied model atmospheres for use in external scripts. Also found in Useful tools and linelists folder.

#### Useful Tools

Contains scripts for editing iSpec linelist for use with only MOOG and combine iSpec linelist with ones made for MOOG

Contains some code to extract specific flagged lines from the ges line lists. 

#### Isotope pipeline

Isotope pipeline contains four main files which calculate the magneisum isotopic abundance for up to 10 spectral features, finds the uncertainties and combines all files into larger ones for plotting and making of tables. 

`Isotope_pipline_v1.py` is the main file. This requires input of star fits or text files a linelist constructed in the MOOG format and information including approximate Fe, CO2 and C0 abundances in the X/H format along with the vsini. It will also require the use of pymoogi found: https://github.com/madamow/pymoogi/

The code is based off of Madeleine Mckenzie's Ratio.py. https://github.com/madeleine-mckenzie/RAtIO/

Each folder must contain the model atmosphere the linelist and spectrum (as .txt is preferred)

`Isotope_uncertainties.py` must be run before `Isotope_analysis.py` to ensure all files exist. This calculates the uncertainties using the diagonal of the Hessian Matrix for each entry into MOOG.

`Isotope_analysis.py` organises all the files into one location for making tables and plotting. 

`print_params.py` and `print_plots.py` both plot various kinds of plots for checking.

`Isotope_pipline_v2.py` is a work in progress and should not be used currently.


#### Test scripts

Contains various test scripts for finding abundances and running moog. None are complete but might be useful for understanding how things run. 

#### Miscellaneous files 

Examples of linelists csv configuration and atmospheres exist in the repository currently. Some example output files also exist.

<!-- ## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method. -->

## Installation
<!-- Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection. -->

Built to run with MOOG, pymoogi and iSpec. Must have all libraries installed for functionallity. 

Pymoogi: https://github.com/madamow/pymoogi/

iSpec: https://ui.adsabs.harvard.edu/abs/2014ASInC..11...85B

MOOG: https://ui.adsabs.harvard.edu/abs/2012ascl.soft02009S/abstract
## How to clone repositiory

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://eng-git.canterbury.ac.nz/mwega/quin-masters-code.git
git branch -M main
git push -uf origin main
```

## Usage
See miscellaneous files.

<!-- ## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc. -->

## Roadmap

Constantly being updated public releases occasionally. 

<!-- ## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser. -->

<!-- ## Authors and acknowledgment
Show your appreciation to those who have contributed to the project. -->

## License

If using please link the github or credit the paper released on the code. "Placeholder".

## Project status
<!-- If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers. -->

In progress.