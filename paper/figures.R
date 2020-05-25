## Figure 1: workflows

library("DiagrammeR")

x <- grViz("digraph course {
#rankdir = LR
node [shape = box, style=filled]
layout = dot
compound =true 
#color = crimson

subgraph clusterPre{
label = 'Preprocessing'
style = dashed
rank = same
maf [label = 'Mutation annotation format (MAF) data']
mafObj [label = 'Read as MAF object']
tallyMaf [label = 'Generate mutation catalogue matrix,\nincluding SBS, DBS and INDEL etc.']
}

subgraph clusterFind{
label = 'De novo Signatures Discovery'
style = dashed
rank = same
estSig [label = 'Run multiple NMF runs to\na range of signature numbers']
extSig [label = 'Extract specified signatures by NNF']
autoExtSig [label = 'Auto-extract signatures by\nBayesian NMF']
extractRes [label = 'Signature compositions &\nSignature exposures']
sim [label = 'Cosine similarity\nanalysis']
}

subgraph clusterFit{
label = 'Reference Signatures Fitting'
style = dashed
rank = same
refSigs [label = 'Custom signatures or\nCOSMIC reference\nsignatures']

fitting [label = 'Fit specified signatures']
fitOpt [label = 'Optimize signature\nexposures']
fitBootstrap [label = 'Analyze signature instability\nwith bootstrapping']
fitRes [label = 'Optimized exposures']
fitBRes [label = 'Exposure estimation &\nreconstruction errors &\nsignificance']
}

subgraph clusterVis{
label = 'Visualization'
style = dashed
rank = same
numSurvey [label = 'Signature number survey']
catProf [label = 'Total or single sample\ncatalogue profile']
sigProf [label = 'Signature profile']
expoProf [label = 'Exposure profile']
fitProf [label = 'Fitting boxplot']
fitBProf [label = 'Bootstrapping\nresult boxplots']
}


maf -> mafObj
mafObj -> tallyMaf
tallyMaf -> estSig
tallyMaf -> extSig
tallyMaf -> autoExtSig
tallyMaf -> fitting
estSig -> numSurvey
numSurvey -> extSig

extSig -> extractRes
autoExtSig -> extractRes

refSigs -> fitting
fitting -> fitOpt
fitting -> fitBootstrap
fitOpt -> fitRes
fitBootstrap -> fitBRes

tallyMaf -> catProf
extractRes -> sigProf
extractRes -> expoProf
fitRes -> fitProf
fitBRes -> fitBProf

extractRes -> sim
}")

x

yyplot::gv2file(x, file = 'paper/Figure1-sigminer-diagram.pdf')

## Figure 2: signatures

library(sigminer)
library(patchwork)

p_sbs = show_cosmic_sig_profile(sig_index = 1:5, sig_db = "SBS", style = "cosmic")
p_tsb = show_cosmic_sig_profile(sig_index = 1:5, sig_db = "TSB", style = "cosmic")
p_dbs = show_cosmic_sig_profile(sig_index = 1:5, sig_db = "DBS", style = "cosmic")
p_indel = show_cosmic_sig_profile(sig_index = 1:5, sig_db = "ID", style = "cosmic")

p_cosmic = (p_sbs / p_tsb) | (p_dbs / p_indel)

ggplot2::ggsave(filename = "paper/Figure2-cosmic-example-profile.pdf", plot = p_cosmic, 
                width = 20, height = 8)

## Figure 3: signature fitting

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
                          sig_db = "legacy")
p_fit <- show_sig_fit(mul_sample_all, palette = NULL, add = NULL) +
  ggpubr::rotate_x_text() +
  ggplot2::ylim(0, 50)

bt_res <- sig_fit_bootstrap_batch(maf_tally$SBS_96[1:100, ] %>% t(), 
                                  sig_index = c(1:3, 6:7, 10, 13, 15, 17, 22, 24),
                                  sig_db = "legacy",
                                  type = "relative",
                                  methods = c("LS", "QP"),
                                  p_val_thresholds = c(0.05, 0.01), 
                                  n = 100)
