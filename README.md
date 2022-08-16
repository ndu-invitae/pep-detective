# Effect peptide detector

### Introduction

In a hypothetical experiment a scientist dsicovered a small peptide that can potentially enhance enzyme activity, and he/she wishes to know on which enzyme, from a pool of candidates, this enhancer is effectice. To test the idea, a highthroughput experiment was conducted to determine all enzyme activites with and without the addition of this "enhancer". During the experiment, the pH level in the solutions are not perfectly stable, and some of the enzyme activities might be impacted by the pH difference, therefore the variance in repeated avtivity measurements might include contribution from this systematic error. From previous experience we learned that the activity of an enzyme should respond linearly to pH change, and the pH value can be accurately meansured in the testing media. We hypothesize that the enhancer can increase enzyme activity; due to the possible existance of a pH covariate, a mixed effect model plus a one tailed t-test was adpoted for rejecting the hypothesis.  

[![pep-detective tests](https://github.com/ndu-invitae/pep-detective/actions/workflows/test.yaml/badge.svg)](https://github.com/ndu-invitae/pep-detective/actions/workflows/test.yaml)
