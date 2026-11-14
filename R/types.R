# Stan types used by the semantic analyzer.
#
# A type is list(base, nd) where base is one of "int", "real", "vector",
# "row_vector", "matrix", "complex", "complex_vector", "complex_row_vector",
# "complex_matrix", "tuple", "void", "string" and nd is the number of array
# dimensions. A tuple type also has elems, the list of its element types.

ty <- function(base, nd = 0L) list(base = base, nd = as.integer(nd))
ty_tuple <- function(elems, nd = 0L) list(base = "tuple", nd = as.integer(nd), elems = elems)
# the same type with nd array dimensions (keeps tuple elements)
ty_nd <- function(t, nd) { t$nd <- as.integer(nd); t }
is_tuple_ty <- function(t) !is.null(t) && identical(t$base, "tuple")
T_int <- ty("int")
T_real <- ty("real")
T_vec <- ty("vector")
T_rvec <- ty("row_vector")
T_mat <- ty("matrix")
T_void <- ty("void")

ty_str <- function(t) {
  if (is.null(t)) return("<unknown>")
  b <- if (is_tuple_ty(t)) sprintf("tuple(%s)", paste(vapply(t$elems, ty_str, ""), collapse = ", ")) else t$base
  if (t$nd > 0L) sprintf("array[%s] %s", strrep(",", t$nd - 1L), b) else b
}

ty_eq <- function(a, b) {
  if (!identical(a$base, b$base) || !identical(a$nd, b$nd)) return(FALSE)
  if (!is_tuple_ty(a)) return(TRUE)
  length(a$elems) == length(b$elems) && all(mapply(ty_eq, a$elems, b$elems))
}
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
  if (is_tuple_ty(to) || is_tuple_ty(from)) {
    if (!is_tuple_ty(to) || !is_tuple_ty(from) || length(to$elems) != length(from$elems)) return(FALSE)
    return(all(mapply(ty_assignable, to$elems, from$elems)))
  }
  if (to$base == from$base) return(TRUE)
  if (to$base == "real" && from$base == "int") return(TRUE)
  is_cx_base(to$base) && !is_cx_base(from$base) && from$base %in% names(real_to_cx) && real_to_cx[[from$base]] == to$base
}

# number of promotions from type `from` to type `to` (int to real and real
# to complex count one each, int to complex two; tuples add their elements)
promo_cost <- function(to, from) {
  if (is_tuple_ty(to)) return(sum(unlist(mapply(promo_cost, to$elems, from$elems))))
  if (to$base == from$base) return(0L)
  if (to$base == "real") return(1L)
  if (from$base == "int") 2L else 1L
}

# element type after removing one array dimension, or one vector/matrix index
ty_index_one <- function(t) {
  if (t$nd > 0L) return(ty_nd(t, t$nd - 1L))
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
