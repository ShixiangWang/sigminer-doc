## Figure: example signatures

library(sigminer)
library(patchwork)

p_sbs = show_cosmic_sig_profile(sig_index = 1:3, sig_db = "SBS", style = "cosmic")
# p_tsb = show_cosmic_sig_profile(sig_index = 1:5, sig_db = "TSB", style = "cosmic")
p_dbs = show_cosmic_sig_profile(sig_index = 1:3, sig_db = "DBS", style = "cosmic")
p_indel = show_cosmic_sig_profile(sig_index = 1:3, sig_db = "ID", style = "cosmic")

load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
                 package = "sigminer", mustWork = TRUE
))
# Show signature profile
p_cn <- show_sig_profile(sig,
                       style = "cosmic",
                       mode = "copynumber",
                       method = "W",
                       normalize = "feature", x_label_angle = 90
)
p_cn

p_all = (p_sbs / p_dbs) | (p_indel / p_cn)

ggplot2::ggsave(filename = "paper/Figure2-example-signature-profile.pdf", plot = p_all, 
                width = 20, height = 8)

## Figure: signature auto-extraction
## Use BRCA data for illustration, 13 signature is detecetd in Nature paper
sbs_sigs <- data.table::fread("paper/PCAWG/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
table(sbs_sigs$`Cancer Types`)
brca_sigs <- sbs_sigs[startsWith(sbs_sigs$`Cancer Types`, "Breast")]

load("paper/PCAWG/lego96.PAN.SNV.091217.RData")
colnames(lego96.SNV) <- gsub(".*__(.+)", "\\1", colnames(lego96.SNV))
brca_lego96 <- lego96.SNV[, colnames(lego96.SNV) %in% brca_sigs$`Sample Names`]

save(brca_sigs, brca_lego96, file = "paper/data/PCAWG_ref.RData")

comps <- colnames(t(brca_lego96))
comps_map <- paste0(
  substr(comps, 3, 3),
  "[",
  substr(comps, 1, 1),
  ">",
  substr(comps, 2, 2),
  "]",
  substr(comps, 4, 4)
)
names(comps_map) <- comps

save(comps_map, file = "paper/data/comps_map.RData")

ref_sigs <- get_sig_db("SBS")$db

## 比较不同方法的重构错误
## SomaticSignatures uses NMF and PCA methods

# https://www.bioconductor.org/packages/release/bioc/vignettes/SomaticSignatures/inst/doc/SomaticSignatures-vignette.html
# can then be chosen such that increasing the number of signatures does not yield
# a significantly better approximation of the data,
# i.e. that the RSS and the explained variance do not change sufficiently for more complex models.
# The first inflection point of the RSS curve has also been proposed as a measure for the number 
# of features in this context
library(SomaticSignatures)
library(MutationalPatterns)
library(sigminer)

## SomaticSignatures
sigs_ss_nmf = identifySignatures(brca_lego96, 13, decomposition = nmfDecomposition)
sigs_ss_pca = identifySignatures(brca_lego96, 13, decomposition = pcaDecomposition)

save(sigs_ss_nmf, sigs_ss_pca, file = "paper/data/Sigs_SomaticSignatures.RData")

## MutationalPatterns
## https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.pdf
mut_mat <- brca_lego96 + 0.0001
library(NMF)
sigs_mp <- extract_signatures(mut_mat, rank = 13)

save(sigs_mp, file = "paper/data/Sigs_MutationalPatterns.RData")

## Sigflow
sigs_sf_nmf <- sig_extract(t(brca_lego96), n_sig = 13, nrun = 200, cores = 4, optimize = TRUE)
sigs_sf_bayes <- sig_auto_extract(t(brca_lego96), 
                                  nrun = 20,
                                  cores = 4,
                                  K0 = 20,
                                  optimize = TRUE,
                                  strategy = "optimal",
                                  result_prefix = "pcawg_brca",
                                  destdir = "paper/bayesianNMF/")

mut_mat2 <- t(brca_lego96)
comps <- colnames(mut_mat2)
colnames(mut_mat2) <- paste0(
  substr(comps, 3, 3),
  "[",
  substr(comps, 1, 1),
  ">",
  substr(comps, 2, 2),
  "]",
  substr(comps, 4, 4)
)
sigprofiler_extract(
  mut_mat2, output = "paper/sigprofiler/PCAWG_BRCA", 
  range = 10:15, nrun = 50, cores = 4, 
  refit = TRUE,
  is_exome = FALSE,
  use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python"
)

sigs_sf_sigprofiler <- sigprofiler_import("paper/sigprofiler/PCAWG_BRCA",
                                          type = "all")
sigs_sf_sigprofiler <- sigs_sf_sigprofiler$solution_list$S13

save(sigs_sf_nmf, sigs_sf_bayes, sigs_sf_sigprofiler, file = "paper/data/Sigs_Sigflow.RData")

load(file = "paper/data/Sigs_Sigflow.RData")
load(file = "paper/data/comps_map.RData")
library(sigminer)

xx = sigs_sf_nmf$Signature.norm
rownames(xx) <- comps_map[rownames(xx)]

yy = get_sig_similarity(xx, sig_db = "SBS")
yy = yy$similarity[, paste0("SBS", c(1, 2, 3, 5, 8, 9, 13, "17a", "17b", 18, 37, 40, 41))]


pheatmap::pheatmap(yy,
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   display_numbers = TRUE,
                   height = 5, width = 7,
                   filename = "paper/cosine_heatmap.png")

pheatmap::pheatmap(yy[paste0("Sig", c(4, 2, 5, 9, 3, 7, 1, 10, 12, 8, 11, 6, 13)), ],
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   display_numbers = TRUE)
pheatmap::pheatmap(yy[paste0("Sig", c(4, 2, 5, 9, 3, 7, 1, 10, 12, 8, 11, 6, 13)), ],
                   cluster_rows = FALSE, cluster_cols = FALSE, 
                   display_numbers = TRUE,
                   height = 5, width = 7,
                   filename = "paper/cosine_heatmap2.png")

sigs_sf_bayes$Raw$summary_run

readr::write_csv(sigs_sf_bayes$Raw$summary_run, "paper/data/bayes_summary_run.csv")


# Comparison
load("paper/data/Sigs_SomaticSignatures.RData")
load("paper/data/Sigs_MutationalPatterns.RData")
load("paper/data/Sigs_Sigflow.RData")

sigs_ss_nmf@fitted %>% sum()
sigs_ss_pca@fitted %>% sum()
sigs_mp$reconstructed %>% sum()
sigs_sf_nmf$Exposure %>% sum()
sigs_sf_bayes$Exposure %>% sum()

rec_ss_nmf <- sigs_ss_nmf@fitted
rec_ss_pca <- sigs_ss_pca@fitted
rec_mp <- sigs_mp$reconstructed
rec_sf_nmf <- sigs_sf_nmf$Signature.norm %*% sigs_sf_nmf$Exposure
rec_sf_bayes <- sigs_sf_bayes$Signature.norm %*% sigs_sf_bayes$Exposure

brca_sig_mat <- brca_sigs %>% 
  dplyr::select(-c("Cancer Types", "Accuracy")) %>% 
  tibble::column_to_rownames("Sample Names") %>% 
  as.matrix()
#rec_sigprofiler <- ref_sigs[, 1:65] %*% t(brca_sig_mat)
rec_sigprofiler <- t(signatures.genome.cosmic.v3.may2019[1:65, ]) %*% t(brca_sig_mat)
xx <- rownames(rec_sigprofiler)
rownames(rec_sigprofiler) <- paste0(
  substr(xx, 3, 3),
  substr(xx, 5, 5),
  substr(xx, 1, 1),
  substr(xx, 7, 7)
)

rec_sigprofiler2 <- sigs_sf_sigprofiler$Signature.norm %*% sigs_sf_sigprofiler$Exposure
xx2 <- rownames(rec_sigprofiler2)
rownames(rec_sigprofiler2) <- paste0(
  substr(xx2, 3, 3),
  substr(xx2, 5, 5),
  substr(xx2, 1, 1),
  substr(xx2, 7, 7)
)

# The ref matrix
brca_lego96

get_measure <- function(mat1, mat2) {
  # Assume component x samples
  # mat1 and mat2 have same dimension and row/col names
  
  samps <- colnames(mat1)
  error <- vector(mode = "numeric", length = length(samps))
  cosine <- vector(mode = "numeric", length = length(samps))
  names(error) <- names(cosine) <- samps
  
  for (i in seq_along(samps)) {
    a <- mat1[, samps[i]]
    a <- a / sum(a)
    b <- mat2[, samps[i]]
    b <- b / sum(b)
    error[[i]] <- sum(abs(a - b))
    cosine[[i]] <- sigminer::cosine(a, b)
  }
  
  message("Mean error: ", mean(error), " Mean cosine: ", mean(cosine))
  return(list(error = error, cosine = cosine))
}

res_ss_nmf <- get_measure(rec_ss_nmf[rownames(brca_lego96), colnames(brca_lego96)],
                          brca_lego96)
res_ss_pca <- get_measure(rec_ss_pca[rownames(brca_lego96), colnames(brca_lego96)],
                          brca_lego96)
res_mp <- get_measure(rec_mp[rownames(brca_lego96), colnames(brca_lego96)],
                      brca_lego96)
res_sf_nmf <- get_measure(rec_sf_nmf[rownames(brca_lego96), colnames(brca_lego96)],
                          brca_lego96)
res_sf_bayes <- get_measure(rec_sf_bayes[rownames(brca_lego96), colnames(brca_lego96)],
                            brca_lego96)
res_sp <- get_measure(rec_sigprofiler[rownames(brca_lego96), colnames(brca_lego96)],
                            brca_lego96)
res_sp2 <- get_measure(rec_sigprofiler2[rownames(brca_lego96), colnames(brca_lego96)],
                      brca_lego96)

save(res_ss_nmf, res_ss_pca, res_mp, res_ss_nmf, res_sf_nmf, res_sf_bayes, res_sp, res_sp2,
     file = "paper/data/extract_comparison.RData")

## Figure: signature fitting
# error calculated as context residue sum 
input_matrix <- t(brca_lego96)
colnames(input_matrix) <- comps_map[colnames(input_matrix)] %>% as.character()

library(deconstructSigs)
library(MutationalPatterns)
library(sigminer)

# deconstructSigs
fit_list_ds <- list() 
fit_list_ds$error <- vector(mode = "numeric", length = nrow(input_matrix))
fit_list_ds$cosine <- vector(mode = "numeric", length = nrow(input_matrix))
names(fit_list_ds$error) <- names(fit_list_ds$cosine) <- rownames(input_matrix)

for (i in seq_len(nrow(input_matrix))) {
  message("Processing ", rownames(input_matrix)[i])
  res_ds <- whichSignatures(tumor.ref = as.data.frame(input_matrix), 
                  signatures.ref = t(ref_sigs) %>% as.data.frame(), 
                  sample.id = rownames(input_matrix)[i], 
                  contexts.needed = TRUE,
                  tri.counts.method = 'genome')
  fit_list_ds$error[[i]] <- res_ds$unknown
  fit_list_ds$cosine[[i]] <- sigminer::cosine(as.numeric(res_ds$tumor), as.numeric(res_ds$product))
  message("Error: ", fit_list_ds$error[[i]], "; Cosine similarity: ", fit_list_ds$cosine[[i]])
}


# MutationalPatterns
fit_mp <- fit_to_signatures(t(input_matrix), ref_sigs)

fit_list_mp <- list() 
fit_list_mp$error <- vector(mode = "numeric", length = nrow(input_matrix))
fit_list_mp$cosine <- vector(mode = "numeric", length = nrow(input_matrix))
names(fit_list_mp$error) <- names(fit_list_mp$cosine) <- rownames(input_matrix)

for (i in seq_len(nrow(input_matrix))) {
  message("Processing ", rownames(input_matrix)[i])
  a <- input_matrix[rownames(input_matrix)[i], ]
  a <- a / sum(a)
  b <- fit_mp$reconstructed[colnames(input_matrix), rownames(input_matrix)[i]]
  b <- b / sum(b)
  fit_list_mp$error[[i]] <- sum(abs(a - b))
  fit_list_mp$cosine[[i]] <- sigminer::cosine(a, b)
  message("Error: ", fit_list_mp$error[[i]], "; Cosine similarity: ", fit_list_mp$cosine[[i]])
}

# sigminer/sigflow

fit_sf <- sig_fit(catalogue_matrix = t(input_matrix),
                  sig = ref_sigs)

sig_fit(catalogue_matrix = t(input_matrix["SP5381", , drop = FALSE]),
        sig = ref_sigs, auto_reduce = TRUE)
fit_sf <- ref_sigs %*% fit_sf

fit_list_sf <- list() 
fit_list_sf$error <- vector(mode = "numeric", length = nrow(input_matrix))
fit_list_sf$cosine <- vector(mode = "numeric", length = nrow(input_matrix))
names(fit_list_sf$error) <- names(fit_list_sf$cosine) <- rownames(input_matrix)

for (i in seq_len(nrow(input_matrix))) {
  message("Processing ", rownames(input_matrix)[i])
  a <- input_matrix[rownames(input_matrix)[i], ]
  a <- a / sum(a)
  b <- fit_sf[colnames(input_matrix), rownames(input_matrix)[i]]
  b <- b / sum(b)
  fit_list_sf$error[[i]] <- sum(abs(a - b))
  fit_list_sf$cosine[[i]] <- sigminer::cosine(a, b)
  message("Error: ", fit_list_sf$error[[i]], "; Cosine similarity: ", fit_list_sf$cosine[[i]])
}

# fit specified signatures
ref_sigs2 <- ref_sigs[, paste0("SBS", c("1", "2", "3", "5", "8", "9", "13",
                                        "17a", "17b", "18", "37", "40", "41"))]
fit_sf_br <- sig_fit(catalogue_matrix = t(input_matrix),
                     sig = ref_sigs2)
fit_sf_br <- ref_sigs2 %*% fit_sf_br

fit_list_sf_br <- list() 
fit_list_sf_br$error <- vector(mode = "numeric", length = nrow(input_matrix))
fit_list_sf_br$cosine <- vector(mode = "numeric", length = nrow(input_matrix))
names(fit_list_sf_br$error) <- names(fit_list_sf_br$cosine) <- rownames(input_matrix)

for (i in seq_len(nrow(input_matrix))) {
  message("Processing ", rownames(input_matrix)[i])
  a <- input_matrix[rownames(input_matrix)[i], ]
  a <- a / sum(a)
  b <- fit_sf_br[colnames(input_matrix), rownames(input_matrix)[i]]
  b <- b / sum(b)
  fit_list_sf_br$error[[i]] <- sum(abs(a - b))
  fit_list_sf_br$cosine[[i]] <- sigminer::cosine(a, b)
  message("Error: ", fit_list_sf_br$error[[i]], "; Cosine similarity: ", fit_list_sf_br$cosine[[i]])
}

save(fit_list_ds, fit_list_mp, fit_list_sf, fit_list_sf_br,
     file = "paper/data/fit_comparison.RData")

# Similar to QP results
# # NNLS
# fit_sf_nnls <- sig_fit(catalogue_matrix = t(input_matrix),
#                   sig = ref_sigs, method = "NNLS")
# fit_sf_nnls <- ref_sigs %*% fit_sf_nnls
# 
# fit_list_sf_nnls <- list() 
# fit_list_sf_nnls$error <- vector(mode = "numeric", length = nrow(input_matrix))
# fit_list_sf_nnls$cosine <- vector(mode = "numeric", length = nrow(input_matrix))
# names(fit_list_sf_nnls$error) <- names(fit_list_sf_nnls$cosine) <- rownames(input_matrix)
# 
# for (i in seq_len(nrow(input_matrix))) {
#   message("Processing ", rownames(input_matrix)[i])
#   a <- input_matrix[rownames(input_matrix)[i], ]
#   a <- a / sum(a)
#   b <- fit_sf_nnls[colnames(input_matrix), rownames(input_matrix)[i]]
#   b <- b / sum(b)
#   fit_list_sf_nnls$error[[i]] <- sum(abs(a - b))
#   fit_list_sf_nnls$cosine[[i]] <- sigminer::cosine(a, b)
#   message("Error: ", fit_list_sf_nnls$error[[i]], "; Cosine similarity: ", fit_list_sf_nnls$cosine[[i]])
# }

# Speed
fit_speed <- microbenchmark::microbenchmark(
  dS = whichSignatures(tumor.ref = as.data.frame(input_matrix), 
                       signatures.ref = t(ref_sigs) %>% as.data.frame(), 
                       sample.id = rownames(input_matrix)[1], 
                       contexts.needed = TRUE,
                       tri.counts.method = 'genome'),
  mP = fit_to_signatures(t(input_matrix[1, , drop = FALSE]), ref_sigs),
  sF = sig_fit(catalogue_matrix = t(input_matrix[1, , drop = FALSE]),
               sig = ref_sigs),
  times = 100L
)

levels(fit_speed$expr) <- c("deconstructSigs", "MutationalPatterns", "Sigflow")
save(fit_speed, file = "paper/data/fit_speed_comparison.RData")

## Plot extract and fit comparison
load("paper/data/extract_comparison.RData")
load("paper/data/fit_comparison.RData")
load("paper/data/fit_speed_comparison.RData")

res_extract <- dplyr::tibble(
  method = rep(c("SomaticSignatures:PCA", "SomaticSignatures:NMF",
                 "MutationalPatterns", "Sigflow:NMF", "Sigflow:BayesianNMF",
                 "SigProfiler"), each = length(res_mp$error)),
  sample = c(names(res_ss_pca$error), names(res_ss_nmf$error),
             names(res_mp$error), names(res_sf_nmf$error),
             names(res_sf_bayes$error), names(res_sp2$error)),
  error = c(res_ss_pca$error, res_ss_nmf$error,
             res_mp$error, res_sf_nmf$error,
             res_sf_bayes$error, res_sp2$error),
  cosine = c(res_ss_pca$cosine, res_ss_nmf$cosine,
             res_mp$cosine, res_sf_nmf$cosine,
             res_sf_bayes$cosine, res_sp2$cosine),
  
)

# fit_list_ds, fit_list_mp, fit_list_sf, fit_list_sf_br
res_fit <- dplyr::tibble(
  method = rep(c("deconstructSigs",
                 "MutationalPatterns",
                 "Sigflow"), each = length(fit_list_mp$error)),
  sample = c(names(fit_list_ds$error), names(fit_list_mp$error),
             names(fit_list_sf$error)),
  error = c(fit_list_ds$error, fit_list_mp$error,
            fit_list_sf$error),
  cosine = c(fit_list_ds$cosine, fit_list_mp$cosine,
             fit_list_sf$cosine),
  
) %>% 
  dplyr::mutate(method = factor(method, 
                                c("MutationalPatterns", "deconstructSigs", "Sigflow")))

library(ggpubr)
library(patchwork)
p1 <- ggboxplot(res_extract, x = "method", y = "error") + rotate_x_text(60) + 
  labs(x = NULL, y = "Sum of absolute residue")
p2 <- ggboxplot(res_extract, x = "method", y = "cosine") + rotate_x_text(60) + 
  labs(x = NULL, y = "Cosine similarity")
p3 <- ggboxplot(res_fit, x = "method", y = "error") + rotate_x_text(60) +
  labs(x = NULL, y = "Sum of absolute residue")
p4 <- ggboxplot(res_fit, x = "method", y = "cosine") + rotate_x_text(60) +
  labs(x = NULL, y = "Cosine similarity")
p5 <- ggplot2::autoplot(fit_speed) + cowplot::theme_cowplot()

pp <- (p1 | p2) /
  (p3 | p4 | p5)

ggplot2::ggsave(filename = "paper/method-comparison.pdf",
                plot = pp, 
                width = 8, height = 7)

## Figure: signature bootstrapping fitting

s <- names(sort(rowSums(input_matrix))[1])
## Just set a sample for examples
bt_res <- suppressMessages(
  sig_fit_bootstrap_batch(input_matrix %>% t(), 
                          sig = ref_sigs,
                          type = "absolute",
                          methods = "QP",
                          p_val_thresholds = c(0.01, 0.05),
                          use_parallel = TRUE,
                          n = 1000,
                          job_id = "pcawg_brca",
                          result_dir = "paper/bootstrap")
)
save(bt_res, file = "paper/data/bt_res.RData")

load(file = "paper/data/bt_res.RData")
# # Pick out the most mutated sample
# which.max(rowSums(maf_tally$SBS_96))
p_stab <- show_sig_bootstrap_stability(bt_res, methods = c("QP"),
                                       add = NULL, ylab = "Signature instability (MRSE)") + 
  ggpubr::rotate_x_text() +
  ggplot2::theme(legend.position = "none")
p_expo <- show_sig_bootstrap_exposure(bt_res, methods = c("QP"), sample = s, highlight = "gold",
                                      highlight_size = 1,
                                      add.params = list(alpha = 0.1)) + ggpubr::rotate_x_text() + 
  ggplot2::theme(legend.position = "none")

input_matrix[s, ] %>% sum()
bt_res$error$errors = bt_res$error$errors / 1203
p_err  <- show_sig_bootstrap_error(bt_res, methods = c("QP"), 
                                   sample = s,
                                   ylab = "Unexplained mutation fraction")

p_bt = p_stab / (p_expo + p_err + patchwork::plot_layout(widths = c(3, 1)))
ggplot2::ggsave(filename = "paper/sig-bootstrap.pdf", plot = p_bt, 
                width = 16, height = 6)

# p value for this sample
