The document `verteq.tex` contains a technical description of how the
code implements the vertical equilibrium model detailed in:

"Upscaled models for modeling CO2 injection and migration in geological systems",
by Sarah E. Gasda, Elsa Du Plessis and Helge K. Dahle, 
[DOI:10.1515/9783110282245.1](http://dx.doi.org/10.1515/9783110282245.1)

It does not detail the model, but rather how the model maps to the
particular implementation of the code.

You should have [TeXLive](http://www.tug.org/texlive/) installed; this
is usually already available in Linux distributions, such as
[Ubuntu](https://help.ubuntu.com/community/LaTeX#From_the_Ubuntu_repositories).

Compilation of the technical document is not done by the build system.
Instead you must compile it yourself using a command like this:

    latexmk -pdf verteq.tex

where `verteq.tex` can also be the full path to the document if you want
to build it outside of the source directory.
