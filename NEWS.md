# genomicMateSelectR 0.2.0

Upgrade. Combine **predCrossVar** with other functions for GS pipeline *and* add high level wrapping functions that implement/automatic cross predictions, and more.

# genomicMateSelectR 0.1.0

See [predCrossVar](https://wolfemd.github.io/predCrossVar/), the OG codebase.

# genomicMateSelectR 0.0.0

## Development To Do

-   [ ] Improve input format for SNP-effects into `predictCrosses()`. More flexible SNP-effect list-column naming (instead of opinionated e.g. **allelesubeffectlist**)? Particularly for the **DirDom** model.

-   [ ] Improve phenotype input e.g. for `runGenomicPredictions()` and `runParentWiseCrossVal()`.

    -   Currently demands a number of variables, be named in particular weights, owing to the two-stage genomic prediction process I use for NextGen Cassava. In that Cassava pipeline, I've used de-regressed BLUPs as response variable for genomic prediction and weighted error variances according to a weighting scheme devised by *Garrick et al. 2009*. Users can subvert this, but are currently required to use certain naming: **drgBLUP** for the response data, **BLUP** for data to-be-used as validation values in cross-validation, and a column of weights (**WT**) that could just be set to 1 if users aren't doing weighted anlaysis.

-   [ ] What would make the package BrAPI compliant?

-   [ ] SNP marker ID naming conventions are probably too strict; add flexibility / robustness

-   [ ] More flexibility for the `crosses2predict()` function: reciprocal crosses? don't include selfs?

-   [ ] Change `predCrossMeans()` so that it can accept either `dosages` or a `haploMat`. Currently, `predictCrosses()` requires users to supply both `dosages` *and* `haploMat` . Changing `predCrossMeans()` in this way will allow users of `predictCrosses()` to supply only the `haploMat` if they are predicting both cross means *and* variances, or only `dosages` if predicting only means.
