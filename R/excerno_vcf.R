#' Exercerno VCF
#'
#' excerno_vcf() produces filtered vcf files. It uses NMF or nonnegative linear combination of mutation signatures to determine contribution of signatures in samples. Then Bayes' Theroem is used to classify each variant.
#'
#' @param files VCF files
#' @param method A string. The method used to determine the signatures (if not given) and calculate the contributions of each signature
#' @param num.signatures Number of signatures. Necessary arugment for when method "linear" is choosen.
#' @param target.sig Matrix of the target signatures.Necessary arugment for when method "linear" is choosen.
#'
#' @return Object containing the vcf objects and classification data frame
#' @examples
#'
#' library(MutationalPatterns)
#' library(tidyverse)
#' library(vcfR)
#' library(Biostrings)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(R.utils)
#'
#' # Load in correct signatures
#' cosmic.sigs <- get_known_signatures()
#' cosmic.sig4 <- as.matrix(cosmic.sigs[, 4])
#' ffpe.sig <- get_ffpe_signature()
#'
#' # Load in files
#' files <- list.files( system.file("extdata", package = "excerno"), pattern = "SIMULATED_SAMPLE_SBS4_\\d.vcf", full.names = TRUE)
#'
#' # Using nmf
#' excerno_vcf(files)
#' excerno_vcf(files, num.signatures = 3)
#'
#' # Using linear method
#'
#' target.sigs <- matrix(nrow = 96, ncol = 2)
#' target.sigs[,1] <- cosmic.sig4
#' target.sigs[,2] <- ffpe.sig
#' rownames(target.sigs) <- get_mutation_types()
#' colnames(target.sigs) <- c("SBS4", "FFPE")
#'
#' excerno_vcf(files, "linear", 2, target.sigs)
#'
#' @export
excerno_vcf <- function(files, method = "nmf", num.signatures = 2, target.sigs = c()) {

  # Arguments validation
  if (!is.character(files)) { stop("argument files must be type character") }
  if (method != "nmf" && method != "linear") { stop("argument method must be \"nmf\" or \"linear\"") }
  if (method == "linear" && is.null(target.sigs)) { stop("argument target.sigs must be non-empty if method equals linear") }
  if (method == "nmf" && length(files) < 2) { stop("argument files must be length greater than 2 if method equals nmf") }

  # Read vcf files
  vcf.data <- list()
  num.samples <- length(files)

  for (i in 1:length(files)) {
    # Hide print
    suppressWarnings(invisible(capture.output(vcf.data[i] <- read.vcfR(files[i]))))
  }

  print_info("Creating mutational vectors")
  samples <- get_mutational_vectors(files)

  # Perform method to get present signatures and contributions of each sample
  if (method == "nmf") {
    print_info("Performing NMF")

    # Argument validation for method = nmf

    samples.df <- data.frame(mutations = get_mutation_types())

    # Add signatures to df
    for (i in 1:num.samples) {
      sample.tbl <- table(samples[[i]])/sum(table(samples[[i]]))
      sample.df <- data.frame(mutations = rownames(sample.tbl))
      sample.df[paste("SAMPLE", toString(i), sep = "")] <- as.vector(sample.tbl)
      samples.df <- left_join(samples.df, sample.df, by = "mutations")
    }

    samples.df[is.na(samples.df)] <- 0

    # Convert data frame to matrix for NMF
    sample.matrix <- as.matrix(samples.df[2:(num.samples + 1)])
    rownames(sample.matrix) <- get_mutation_types()

    # NMF and renaming columns and rows
    nmf.res <- extract_signatures(sample.matrix, rank = num.signatures, nrun = 10)
    sig.names <- find_signature_names(nmf.res$signatures)
    colnames(nmf.res$signatures) <- sig.names
    rownames(nmf.res$contribution) <- sig.names

    signatures <- nmf.res$signatures
    contribution <- nmf.res$contribution
  } else if (method == "linear") {
    print_info("Performing Linear combination method")

    # Defines a function that joins two data frames and replaces NA with 0
    left_join_NA <- function(x, y, by) {
  dplyr::left_join(x = x, y = y, by = by) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), \(col) replace(col, which(is.na(col)), 0)))
}

    # Creates a vector of all 96 snv mutation types
    mutations <- get_mutation_types()
    num.signatures <- length(colnames(target.sigs))
    sig.names <- colnames(target.sigs)

    contribution <- matrix(nrow = num.signatures, ncol = num.samples)
    rownames(contribution) <- sig.names

    for (i in 1:num.samples) {

      # Converts mutations_vector into a table of probabilities
      tab <- table(samples[[i]])/sum(table(samples[[i]]))

      # Converts table of probabilities to a data frame and defines column names
      mutation.probs.df <- data.frame(tab)
      colnames(mutation.probs.df) <- c("mutations", "frequencies")

      # Creates a data.frame of mutation types to join with sample probabilities
      df <- data.frame(mutations)

      # Joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
      mut.profile.df <- left_join_NA(df, mutation.probs.df, by = "mutations")

      # Transform mut.profile.df to mutational profile matrix, add pseudocount to column of probabilities
      mut.profile<- as.matrix(mut.profile.df$frequencies) + 0.0001
      rownames(mut.profile)= mutations

      fits.res <- fit_to_signatures(mut.profile, target.sigs)

      contribution[,i] <- fits.res$contribution
    }
    signatures <- target.sigs
    colnames(contribution) <- paste("SAMPLE", seq(1:num.samples), sep = "")
  }

  print_info("Generating classification data frames")
  # Determine signature probabilities from each sample
  classifications.df <- get_classifications(signatures, contribution)

  print_info("Loading in values to vcf files")
  # Insert probabilities into original vcfR object
  for (i in 1:num.samples) {
    write_classification_to_vcf(files[[i]], classifications.df[[i]])
  }

  # # Create list for outputting results
  # output <- list()
  # output$vcf.data <- vcf.data
  # output$class.df <- classifications.df
  #
  # return (output)
}
