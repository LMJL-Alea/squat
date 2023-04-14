## Test environments

- local macOS R installation, R 4.2.2
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

There were 1 NOTE:

    * checking installed package size ... NOTE
      installed size is 13.0Mb
      sub-directories of 1Mb or more:
        data   1.1Mb
        libs  11.6Mb