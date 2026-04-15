# vegan can drop dimnames from stratified permat* results; restore them
# before picante matches community and distance matrix taxa.
.restore_permuted_dimnames <- function(permuted, reference) {
  permuted <- as.matrix(permuted)
  dimnames(permuted) <- dimnames(reference)
  permuted
}
