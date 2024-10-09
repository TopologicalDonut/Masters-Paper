# Master’s Paper

This repository contains the source code and a copy of my paper on
quantifying heterogeneity of racial discrimination in the Canadian labor
market. This was done as the capstone project for my master’s degree at
Simon Fraser University.

# Notes on Structure of Repository

Below are some details of the repository that are not apparent from the
naming.

`code/causal_forest_test.R`: This is a script that compares the
performance of a causal forest against standard linear regression with
interaction terms when the true data generating process is nonlinear.
The results are used in the appendix of my paper.

`code/vimp_causal_forests.R`: This is a slightly modified version of a
function that was created and used by the authors of Bénard and Josse
(2023). It is an implementation of a method to calculate the proportion
of treatment effect variation that is explained by each covariate used
in a causal forest.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-benard2023" class="csl-entry">

Bénard, Clément, and Julie Josse. 2023. “Variable Importance for Causal
Forests: Breaking down the Heterogeneity of Treatment Effects.”
<https://doi.org/10.48550/ARXIV.2308.03369>.

</div>

</div>
