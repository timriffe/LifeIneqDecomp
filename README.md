# This package is now merged into the LifeIneq package, here: https://github.com/timriffe/LifeIneq

You can still separately install and use from here, but it will no longer be maintained. Old readme follows:

# Additive decompositions of lifespan inequality measures: LifeIneqDecomp

## Contact

Tim Riffe, tim.riffe@ehu.eus

Christian Dudel, dudel@demogr.mpg.de

## Description

This package provides a function to decompose a lifespan inequality index into additive components of between-group inequality and within-group inequality. Presently implemented for Theil's index, e-edagger, variance, mean log deviation, the Gini coeficient, mean absolute deviation, lifetable entropy, and the absolute inter-individual difference in lifespan.

For each index, the decomposition is of the form T = W + B, where T is the total lifespan inequality, W is the lifespan inequality within groups, and B is the inequality between groups. For all measures, W is a weighted average of T calculated for each group, T_k, while B is a weighted average of each groups deviation from the population average lifespan or a closely related statistic.


