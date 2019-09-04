##' @title Generate Stan code 
##' @return stan code
##' @author philippe
##' @keywords internal
##' @importFrom glue glue

scode_cov2cor = function(){
    out = "matrix cov2cor(matrix C){
      int dm = rows(C);
      matrix[dm,dm] s;
      matrix[dm,dm] R;
      s = diag_matrix( 1 ./ sqrt(diagonal(C)) );
      R = s*C*s ; //quad_form_diag(C, s);
      return R;
      }"
    out
}
