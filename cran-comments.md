## Resubmission

In this resubmission, I have:

- quoted the name of the 'Eigen' C++ library in the description of the package;
- corrected the title of the package to respect title case.

## Test environments

- local macOS R installation, R 4.2.1
- continuous integration via GH actions:
    - macOS latest release
    - windows latest release
    - linux unbuntu release
    - ubuntu 20.04 latest release and devel
- win-builder (release and devel)
- macOS-builder
- R-hub
    - Windows Server 2022, R-devel, 64 bit
    - Ubuntu Linux 20.04.1 LTS, R-release, GCC
    - Fedora Linux, R-devel, clang, gfortran
    - Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

There was no ERROR and no WARNINGs.

There were 2 NOTEs:

    * checking examples ... [28s] NOTE
    Examples with CPU (user + system) or elapsed time > 5s
                    user system elapsed
    plot.prcomp_qts 7.34   0.28    7.77

The size varies according to the system on which the package is installed.

    * checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
      'lastMiKTeXException'

This NOTE appears only on Windows.
