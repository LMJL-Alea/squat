## Test environments

* local macOS R installation, R 4.5.0
* continuous integration via GH actions:
  * macOS latest release
  * windows latest release
  * ubuntu 24.04.2 LTS and devel and oldrel-1
* [win-builder](https://win-builder.r-project.org/) (release and devel)
* [macOS-builder](https://mac.r-project.org/macbuilder/submit.html)
* [R-hub](https://r-hub.github.io/rhub/)

## R CMD check results

There was no ERROR and no WARNINGs.

There was 1 NOTE (installed size varies depending on the OS):

    * checking installed package size ... NOTE
      installed size is 13.9Mb
      sub-directories of 1Mb or more:
        data   1.1Mb
        help   1.2Mb
        libs  11.4Mb
