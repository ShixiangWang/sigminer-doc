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