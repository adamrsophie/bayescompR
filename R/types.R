# Stan types used by the semantic analyzer.
#
# A type is list(base, nd) where base is one of "int", "real", "vector",
# "row_vector", "matrix", "complex", "complex_vector", "complex_row_vector",
# "complex_matrix", "void", "string" and nd is the number of array
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

# complex types and their real counterparts
cx_to_real <- c(complex = "real", complex_vector = "vector", complex_row_vector = "row_vector", complex_matrix = "matrix")
real_to_cx <- c(int = "complex", real = "complex", vector = "complex_vector", row_vector = "complex_row_vector", matrix = "complex_matrix")
is_cx_base <- function(b) b %in% names(cx_to_real)
is_cx_ty <- function(t) !is.null(t) && is_cx_base(t$base)
ty_realpart <- function(t) if (is_cx_ty(t)) ty(cx_to_real[[t$base]], t$nd) else t
ty_cx <- function(t) if (is_cx_ty(t)) t else ty(real_to_cx[[t$base]], t$nd)

# can a value of type `from` be passed where `to` is expected (int to real
# promotion, real to complex promotion, also inside arrays)
ty_assignable <- function(to, from) {
  if (to$nd != from$nd) return(FALSE)
  if (to$base == from$base) return(TRUE)
  if (to$base == "real" && from$base == "int") return(TRUE)
  is_cx_base(to$base) && !is_cx_base(from$base) && from$base %in% names(real_to_cx) && real_to_cx[[from$base]] == to$base
}

# element type after removing one array dimension, or one vector/matrix index
ty_index_one <- function(t) {
  if (t$nd > 0L) return(ty(t$base, t$nd - 1L))
  switch(t$base,
    vector = T_real, row_vector = T_real, matrix = T_rvec,
    complex_vector = ty("complex"), complex_row_vector = ty("complex"), complex_matrix = ty("complex_row_vector"),
    NULL)
}

# number of indexable dimensions (array dims plus vector/matrix dims)
ty_total_dims <- function(t) {
  t$nd + switch(t$base, vector = , row_vector = , complex_vector = , complex_row_vector = 1L,
                matrix = , complex_matrix = 2L, 0L)
}

# Type of a sized declaration (constrained types map to their base types)
decl_base_type <- function(base) {
  switch(base,
    int = "int", real = "real",
    vector = , simplex = , ordered = , positive_ordered = , unit_vector = ,
    sum_to_zero_vector = "vector",
    row_vector = "row_vector",
    complex = , complex_vector = , complex_row_vector = , complex_matrix = base,
    "matrix")
}
