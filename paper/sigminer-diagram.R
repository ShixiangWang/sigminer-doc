library("DiagrammeR")


grViz("digraph course {
rankdir = LR
node [shape = box, style=filled]
layout = dot
compound =true 
#color = crimson
subgraph clusterA{
label = 'Data Import And Transformation'
style = dashed
rank = same
maf [label = 'Mutation Annotation File']
readMaf [label = 'read_maf']
tallyMaf [label = 'sig_tally']
}
subgraph clusterB{
rankdir = LR
label = 'De novo Signatures Discovery'
style = dashed
rank = same
subgraph cluster0 {
label = 'Manual extraction'
sig_estimate
show_sig_number_survey
sig_extract
sig1 [label = 'Signature object']
}
subgraph cluster1 {
label = 'Auto extraction'
sig_auto_extract
sig2 [label = 'Signature object']
}
}
subgraph clusterC{
label = 'Reference Signatures Fitting'
style = dashed
rank = same
subgraph cluster0 {
label = 'Reference signatures'
cosmic [label = 'COSMIC database']
custsig [label = 'Custom signatures']
}
subgraph cluster1 {
label = 'Optimal fitting'
sig_fit
}
subgraph cluster2 {
label = 'Stability analysis\nfor single sample'
sig_fit_bootstrap
report_bootstrap_p_value
}
subgraph cluster3 {
label = 'Stability analysis\nfor multiple samples'
sig_fit_bootstrap_batch
}
}
subgraph clusterD{
label = 'Visualization with ggplot2'
style = dashed
rank = same
subgraph cluster1 {
label = 'Signature fitting'
show_sig_fit
show_sig_bootstrap_exposure
show_sig_bootstrap_error
show_sig_bootstrap_stability
}
subgraph cluster2 {
label = 'Signature/catalogue profile'
show_sig_profile
show_sig_exposure
show_cosmic_sig_profile
show_catalogue
}
}

maf -> readMaf
readMaf -> tallyMaf

tallyMaf -> sig_estimate
tallyMaf -> sig_extract
tallyMaf -> sig_auto_extract
tallyMaf -> sig_fit
tallyMaf -> sig_fit_bootstrap
tallyMaf -> sig_fit_bootstrap_batch
tallyMaf -> show_catalogue

sig_estimate -> show_sig_number_survey
show_sig_number_survey -> sig_extract
sig_extract -> sig1
sig_auto_extract -> sig2

sig1 -> custsig
sig1 -> cosmic
sig2 -> custsig
sig2 -> cosmic

sig_fit_bootstrap -> report_bootstrap_p_value

sig_fit -> show_sig_fit
sig_fit_bootstrap_batch -> show_sig_bootstrap_exposure
sig_fit_bootstrap_batch -> show_sig_bootstrap_error
sig_fit_bootstrap_batch -> show_sig_bootstrap_stability

sig1 -> show_sig_profile
sig1 -> show_sig_exposure
sig2 -> show_sig_profile
sig2 -> show_sig_exposure
cosmic -> show_cosmic_sig_profile

}")  -> x

# Usage
#gene [label = 'Gene of interest\ne.g. from pull-down']
# microarray -> genelist #gseKEGG[lhead=clusterC]
# enrichPathway -> enrichMeSH  [style = invis]
# gseDGN -> cnetplot [ltail=clusterC]

x

#yyplot::gv2file(x,  file = 'sigminer-diagram.pdf')
yyplot::gv2file(x, file = 'sigminer-diagram.png' )