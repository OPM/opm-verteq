You should have [TeXLive](http://www.tug.org/texlive/) installed; this
is usually already available in Linux distributions, such as
[Ubuntu](https://help.ubuntu.com/community/LaTeX#From_the_Ubuntu_repositories).

Compilation of the technical document is not done by the build system.
Instead you must compile it yourself using a command like this:

    latexmk -pdf verteq.tex

where `verteq.tex` can also be the full path to the document if you want
to build it outside of the source directory.
