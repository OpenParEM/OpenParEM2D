# OpenParEM2D
**Open** **Par**allel **E**lectro**m**agnetic **2D** – A free, open-source electromagnetic simulator for 2D waveguides and transmission lines.

[OpenParEM2D](https://openparem.org) solves waveguides and transmission lines with arbitrary cross sections and conductor counts for the propagation constants, losses (dielectric, conductor, and surface roughness), characteristic impedances, and fields of the dominant and higher-order modes.  It is a full-wave solver simultaneously solving all of Maxwell's equations, so the solutions are good to arbitrarily high frequencys with a problem-dependent limit on lower frequencies.

The tool uses advanced finite elements (FEM) of arbitrarily high order, adaptive mesh refinement, and parallel processing using the Message Passing Interface (MPI).  An extensive set of worked examples demonstrate a variety of setups and solution styles using rectangular waveuide, coax, microstrip, coupled microstrip, and stripline.

A detailed document on the theory, methodology, and accuracy of OpenParEM2D is available along with a comprehensive user's manual.

See the project web site at https://openparem.org.

Installation instructions are provided in the doc directory of the OpenParEMCommon project.  Pre-compiled binaries are available for Ubuntu Linux along with source code for compiling from source.

OpenParEM2D is a free, open-source project licensed under GPLv3 or later.
