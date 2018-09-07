============
 Background
============

:mod:`hop` is a collection of Python modules to analyze solvent
dynamics in molecular dynamics (MD) simulations. It generates a
spatially and temporally coarse grained representation of the dynamics
in terms of a **hopping graph**.

The idea is to first find regions with a density above a given
threshold and catalogue those **sites** (for water, these would be
hydration sites, for other solvent molecules simply high density
locations) . Once this is done, one can analyze water movement in
terms of **hops between those sites**. The complicated solvation
dynamics is thus represented as a **hopping graph** in which hydration
sites are the nodes (or vertices) and movements between sites are the
edges.

However, in principle one is not restricted to using high density
sites. Any geometric partition of space can be used (such as "inside"
and "outside" of a binding site to measure exchange with a binding
site and derive on/off rate constants or "periplasmic" and "cytosolic"
side of a membrane to derive permeation rates through a channel).


