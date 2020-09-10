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

save(bbrca_sigs, brca_lego96, file = "paper/data/PCAWG_ref.RData")

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


## Figure: signature fitting

library(deconstructSigs)
library(MutationalPatterns)

load("paper/maf.RData")
maf = read_maf(maf)
maf

maf_tally <- sig_tally(
  maf,
  mode = "ALL",
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  genome_build = "hg38",
  use_syn = TRUE, 
  add_trans_bias = TRUE
)

mul_sample_all <- sig_fit(catalogue_matrix = maf_tally$SBS_96 %>% t(), 
                          sig_index = 1:30,
                          sig_db = "legacy", type = "relative")
p_fit <- show_sig_fit(mul_sample_all, palette = NULL, add = NULL) +
  ggpubr::rotate_x_text() 
# +
#  ggplot2::ylim(0, 50)

p_fit

bt_res <- sig_fit_bootstrap_batch(maf_tally$SBS_96 %>% t(), 
                                  sig_index = c(1:3, 6:7, 10, 13, 15, 17, 22, 24),
                                  sig_db = "legacy",
                                  type = "relative",
                                  methods = c("LS", "QP"),
                                  p_val_thresholds = c(0.01),
                                  use_parallel = FALSE,
                                  n = 100,
                                  job_id = "tcga_brca",
                                  result_dir = "paper/bootstrap")
save(bt_res, file = "paper/bt_res.RData")

load(file = "paper/bt_res.RData")
# # Pick out the most mutated sample
# which.max(rowSums(maf_tally$SBS_96))
p_stab <- show_sig_bootstrap_stability(bt_res, methods = c("LS", "QP"), add = NULL, ylab = "Signature instability (MRSE)") + ggpubr::rotate_x_text()
p_expo <- show_sig_bootstrap_exposure(bt_res, methods = c("LS", "QP"), add.params = list(alpha = 0.1)) + ggpubr::rotate_x_text()
p_err  <- show_sig_bootstrap_error(bt_res, methods = c("LS", "QP"), ylab = "Reconstruction error (F2 norm)")

p_3 = p_fit / p_stab / (p_expo | p_err)
ggplot2::ggsave(filename = "paper/Figure3-sig-fitting-and-bootstrap.pdf", plot = p_3, 
                width = 8, height = 10)

# p value for this sample
bt_res$p_val[sample == "TCGA-3C-AAAU-01A-11D-A41F-09"]
