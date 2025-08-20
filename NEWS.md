# fusedTree 1.1.0

## New features

* Tree fits are now also accepted from  from the `partykit` package
(`constparty` objects), instead of only from the `rpart` package.

* Functions `fusedTree` and `PenOpt` now take a new optional input argument
`symFusion`. Setting this argument to `FALSE` allows for asymmetric fusion,
where omics effects in  nodes having similar predictions of the response
(using only clinical variables) will be shrunken less toward each other compared
to omics effects in nodes having more distinct predictions. 

* For memory efficiency, Matrix algebra is now conduced using sparse matrices
from the `Matrix` package when `X` is high-dimensional.

* Predict function does not require the test response by default anymore. Rather,
it is an optional parameter.

# fusedTree 1.0.1

## New features

* Re-submission with minor changes to make the package more suitable for
CRAN release.


# fusedTree 1.0.0

## New features

* Initial CRAN submission
