# Stan types used by the semantic analyzer.
#
# A type is list(base, nd) where base is one of "int", "real", "vector",
# "row_vector", "matrix", "void", "string" and nd is the number of array
# dimensions.

ty <- function(base, nd = 0L) list(base = base, nd = as.integer(nd))
T_int <- ty("int")
T_real <- ty("real")
T_vec <- ty("vector")
T_rvec <- ty("row_vector")
T_mat <- ty("matrix")
T_void <- ty("void")

ty_str <- function(t) {
  if (is.null(t)) return("<unknown>")
  if (t$nd > 0L) sprintf("array[%s] %s", strrep(",", t$nd - 1L), t$base) else t$base
}

ty_eq <- function(a, b) identical(a$base, b$base) && identical(a$nd, b$nd)
is_scalar_ty <- function(t) t$nd == 0L && t$base %in% c("int", "real")
is_int_ty <- function(t) t$nd == 0L && t$base == "int"
is_real_scalar_ty <- function(t) t$nd == 0L && t$base == "real"
is_eigen_ty <- function(t) t$nd == 0L && t$base %in% c("vector", "row_vector", "matrix")

# can a value of type `from` be passed where `to` is expected (int to real
# promotion, also inside arrays)
ty_assignable <- function(to, from) {
  if (to$nd != from$nd) return(FALSE)
  if (to$base == from$base) return(TRUE)
  to$base == "real" && from$base == "int"
}

# element type after removing one array dimension, or one vector/matrix index
ty_index_one <- function(t) {
  if (t$nd > 0L) return(ty(t$base, t$nd - 1L))
  switch(t$base,
    vector = T_real, row_vector = T_real, matrix = T_rvec,
    NULL)
}

# number of indexable dimensions (array dims plus vector/matrix dims)
ty_total_dims <- function(t) {
  t$nd + switch(t$base, vector = 1L, row_vector = 1L, matrix = 2L, 0L)
}

# Type of a sized declaration (constrained types map to their base types)
decl_base_type <- function(base) {
  switch(base,
    int = "int", real = "real",
    vector = , simplex = , ordered = , positive_ordered = , unit_vector = ,
    sum_to_zero_vector = "vector",
    row_vector = "row_vector",
    "matrix")
}
