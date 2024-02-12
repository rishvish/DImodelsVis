## Test environments

- local Windows install, R 4.3.2
- win-builder (devel and release)
- macOS builder
- rhub platforms
  - debian-clang-devel (Debian Linux, R-devel, clang, ISO-8859-15 locale)
  - debian-gcc-devel (Debian Linux, R-devel, GCC)
  - debian-gcc-patched (Debian Linux, R-patched, GCC)
  - debian-gcc-release (Debian Linux, R-release, GCC)
  - fedora-gcc-devel (Fedora Linux, R-devel, GCC)
  - ubuntu-gcc-devel (Ubuntu Linux 20.04.1 LTS, R-devel, GCC)
  - ubuntu-gcc-release (Ubuntu Linux 20.04.1 LTS, R-release, GCC)
  - ubuntu-rchk (Ubuntu Linux 20.04.1 LTS, R-devel with rchk)
  - windows-x86_64-devel (Windows Server 2022, R-devel, 64 bit)
  - windows-x86_64-oldrel (Windows Server 2022, R-oldrel, 32/64 bit)
  - windows-x86_64-patched (Windows Server 2022, R-patched, 32/64 bit)
  - windows-x86_64-release (Windows Server 2022, R-release, 32/64 bit)
  
## R CMD check results

There were no ERRORS or WARNINGS when building the package on my local machine. I encountered the following three notes.

```
checking CRAN incoming feasibility 
  Maintainer: 'Rishabh Vishwakarma <vishwakr@tcd.ie>'
  New submission
```
- According to this [post](https://mailman.stat.ethz.ch/pipermail/r-devel/2014-March/068497.html) by Uwe Ligges this note should be ignored.
```
Checking HTML version of manual ... NOTE
  Found the following HTML validation problems:
  DImodelsVis-package.html:71:1 (DImodelsVis-package.Rd:52): Warning: inserting implicit <p>
  Skipping checking math rendering: package 'V8' unavailable
```
- I think this is similar to to [R-hub issue #548](https://github.com/r-hub/rhub/issues/548) and is unrelated to my package.

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
- As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

Additionally I also encountered the following notes when building the package on different platforms using `rhub`

```
Possibly misspelled words in DESCRIPTION:
  Compositional (3:9)
  DImodelsVis (11:143)
  compositional (11:40, 11:240)
  ggplot (11:317, 11:383)
```
- Following the suggestion [here](https://github.com/ThinkR-open/prepare-for-cran?tab=readme-ov-file#:~:text=Other%20packages%20should%20be%20between%20%27%20(example%3A%20Lorem%2DIpsum%20Helper%20Function%20for%20%27shiny%27%20Prototyping)), I have added DImodelsVis and ggplot in quotes in the DESCRIPTION to highlight they're R packages and have checked that the spelling for compositional is correct.


```
checking for non-standard things in the check directory
     ''NULL''
```
- As noted in R-hub issue [#560](https://github.com/r-hub/rhub/issues/560),by R-hub maintainer, this should be ignored.

- I also encountered PREPERROR when building the package on Linux systems using rhub but the log files said the build was a success. This is due to the design of the rhub application [#448](https://github.com/r-hub/rhub/issues/448) and is unrelated to my package.

## Reverse dependencies

* None, as this is a new release.

