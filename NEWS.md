# squat 0.0.1.9000

* Make all functions applicable to a single QTS also applicable to QTS samples,
with appropriate class for the output.
* Enable `as_qts_sample()` to generate a QTS sample of size 1 from a single QTS
as input argument.
* Rename `change_points` argument to the `plot.qts()` function to better reflect
its flexibility.
* Added subset operator for QTS sample objects.
* Added `append` S3 method for QTS sample objects.
* Added `hemispherize()` function to remove any discontinuities in QTS due to
quaternion flips.

# squat 0.0.1

* Added a `NEWS.md` file to track changes to the package.
* Initial version.
