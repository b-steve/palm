Changes in this update are as follows.
* Fixing missing package anchors in the Rd cross-references, as
  requested in the email from Kurt Hornik.
* Switching the package level documentation to use the special
  sentinel "_PACKAGE".

## Test environments

* Local Xubuntu 24.04 install (release)
* win-builder (both devel and release)

## R CMD check results

With win-builder (both devel and release) there were no ERRORs,
WARNINGs, or NOTEs.

With local Xubuntu 22.04 install (release) there was one NOTE:

  * checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

I have spent time searching for the cause of this NOTE. Based on my
investigation, my best guess is that it is caused by a setting on the
CRAN servers rather than a problem with the package's source code. If
this is incorrect, my apologies, and please return this submission to
me.

## Downstream dependencies

* There are currently no downstream dependencies for this package.
