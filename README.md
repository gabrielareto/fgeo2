# ForestGEO / CTFS R code review and development

This repository contains working R code for evaluating, updating, and developing analytical workflows for ForestGEO forest plot data, and possibly other compatible forest datasets.

ForestGEO, formerly known as the Center for Tropical Forest Science (CTFS), is an initiative led by the Smithsonian Institution that studies forest diversity via 80+ research sites. Each research site conducts censuses of every tree within the plot at that location. There is substantial convergence in standards and methods within this network. In recent years, ForestGEO has co-led the development of the Alliance for Tropical Science (ATFS) to promote collaboration among tropical forest research networks in the analysis of tropical forest plot data.  

The ForestGEO/CTFS R codebase originated with Richard Condit as the CTFS R Package, a set of scripts with functions designed to analyze the long-term tropical forest plot data collected by the network, totaling around 400 R functions. The package implemented standardized analytical routines for forest dynamics — including growth, mortality, recruitment, biomass, and habitat association analyses — using a consistent database schema. Later, Mauro Lepore and collaborators reorganized and modularized this legacy code into the fgeo suite of R packages, aiming for reproducibility, tidy data workflows, and modular architecture. Unfortunately, and despite its broad interest and use by the ForestGEO network and other tropical forest researchers, the fgeo ecosystem of packages is currently unstable and has been difficult to maintain. Many functions have quickly become outdated, the heavy wrapping and inter-dependency between functions make the code very difficult to debug and improve, and the users frequently experience crashes, namespace conflicts, and opaque error messages, limiting accessibility to researchers and technicians managing these datasets. Because of these problems, many users end up resorting to the old version of the code, which was developed 20 years ago and is at least partially outdated. 

Based on the feedback gathered from the users, ForestGEO has these code-related goals:

- Goal 1: Make the existing body of code more reliable, automated, and easy to use, prioritizing transparency, flexibility, and reproducibility over formal packaging.

- Goal 2: Expand the functionality for QA/QC tasks, adding systematic tools for data checking, error detection, and quality control, which are recurring needs but largely absent from the original code.

## Scope

The repository is used to track code changes, code evaluations, tests, examples, and implementation notes. It is a working development repository, not a formal R package.

The code in this repository may include:

- review of existing CTFS / fgeo functions;
- updated or simplified R functions;
- scripts for common forest plot analyses;
- QA/QC checks for census-style datasets;
- tests and reproducible examples used during development;
- notes documenting code decisions, limitations, and unresolved issues.

The main priorities are:

- transparent code;
- minimal unnecessary dependencies;
- explicit inputs and outputs;
- reproducible examples where available;
- easier debugging and maintenance;
- compatibility with ForestGEO/CTFS-style data structures.

## Repository status

This repository is under active development.

Code may change as functions are reviewed, tested, simplified, replaced, or removed. Users should check script headers, issue notes, and commit history before reusing code outside the current development context.

