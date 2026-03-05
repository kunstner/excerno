#' Signature cosine similarity
#'
#' Calculating the cosine similarity between a vector of mutations and a mutational signature
#'
#' @param mutations.vector A vector of mutations. Type string.
#' @param siganture.matrix A signature matrix
#' @return A number between 0 and 1
#'
#' @examples
#'
#' library(MutationalPatterns)
#' library(tidyverse)
#'
#' cosmic.sig4 <- as.matrix(get_known_signatures()[,4])
#' sample.sig4 <- create_signature_sample_vector(cosmic.sig4)
#' signature_cosine_similarity(sample.sig4, cosmic.sig4)
#' @export
signature_cosine_similarity <- function(mutations.vector, signature.matrix){

  # Argument validations
  if (!is.character(mutations.vector)) { stop("argument mutations.vector must be type character") }
  if (!is.matrix(signature.matrix)) { stop("argument signature.matrix must be class matrix") }
  if (length(signature.matrix) != 96) { stop("argument signature.matrix must be length 96") }

  # Create a table of probabilities for the vector of mutations
  tab <- table(mutations.vector)/sum(table(mutations.vector))
  # Transform table to a data frame and assign column names to mutations and frequencies
  mutations.df <- data.frame(tab)
  colnames(mutations.df) <- c("mutations", "frequencies")

  # Defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
  dplyr::left_join(x = x, y = y, by = by) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(col) replace(col, which(is.na(col)), 0)))
}

  # Creates a vector of the 96 mutation types
  mutations <- get_mutation_types()
  # Converts vector of 96 mutation types to a data frame
  mutation.types.df <- data.frame(mutations)
  # Joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  complete.mutations.df <- left_join_NA(mutation.types.df, mutations.df, by = "mutations")

  # Converts data frame of mutation types and frequencies to a data frame and assigns rownames
  mutations.matrix <- as.matrix(complete.mutations.df$frequencies)
  rownames(mutations.matrix) = mutations

  # Calculates and returns cosine similarity to the matrix for the signature of interest
  cosine.similarity.matrix <- cos_sim_matrix(signature.matrix, mutations.matrix)
  return(cosine.similarity.matrix[,1])
}
