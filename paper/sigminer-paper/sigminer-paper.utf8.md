---
title: Template for Oxford University Press papers
#date: "2020-05-23"
authors:
  - name: Alice Anonymous
    email: alice@example.com
    address: Some Institute of Technology
  - name: Bob Security
    email: bob@example.com
    address: Another University
    corresponding_author: yes
  - name: Cat Memes
    email: cat@example.com
    address: Another University
  - name: Derek Zoolander
    email: derek@example.com
    address: Some Institute of Technology
abstract: |
  This is the abstract.

  It consists of two paragraphs.
acknowledgements: |
  This is an acknowledgement.

  It consists of two paragraphs.
keywords:
  - key
  - dictionary
  - word
#fontsize: 12pt
#spacing: halfline # could also be oneline
#classoptions:
#  - endnotes
bibliography: mybibfile.bib
output: rticles::oup_article
---



# Introduction

This template is based on the generic OUP template available [here](https://academic.oup.com/icesjms/pages/General_Instructions). The original OUP sample tex document, providing more details on prefered formatting for LaTeX documents, is included with the template in the file `ouparticle_sample.tex`.

Here are two sample references: @Feynman1963118 [@Dirac1953888]. Bibliography will appear at the end of the document.

# Materials and methods

An equation with a label for cross-referencing:

\begin{equation}\label{eq:eq1}
\int^{r_2}_0 F(r,\varphi){\rm d}r\,{\rm d}\varphi = [\sigma r_2/(2\mu_0)]
\int^{\infty}_0\exp(-\lambda|z_j-z_i|)\lambda^{-1}J_1 (\lambda r_2)J_0
(\lambda r_i\,\lambda {\rm d}\lambda)
\end{equation}

This equation can be referenced as follows: Eq. \ref{eq:eq1}

## A subsection

A numbered list:

1) First point
2) Second point
    - Subpoint
    
A bullet list:

* First point
* Second point

# Results

Generate a figure.


```r
plot(1:10,main="Some data",xlab="Distance (cm)",ylab="Time (hours)")
```

\begin{figure}[p]
\includegraphics[width=1\linewidth]{sigminer-paper_files/figure-latex/fig1-1} \caption{This is the first figure.}\label{fig:fig1}
\end{figure}

You can reference this figure as follows: Fig. \ref{fig:fig1}.


```r
plot(1:5,pch=19,main="Some data",xlab="Distance (cm)",ylab="Time (hours)")
```

\begin{figure}[p]
\includegraphics[width=1\linewidth]{sigminer-paper_files/figure-latex/fig2-1} \caption{This is the second figure.}\label{fig:fig2}
\end{figure}

Reference to second figure: Fig. \ref{fig:fig2}

Generate a table.


```r
df = data.frame(ID=1:3,code=letters[1:3])
print(xtable(df,caption="This is the table caption",label="tab:tab1"),
      comment=FALSE)
```

\begin{table}[ht]
\centering
\begin{tabular}{rrl}
  \hline
 & ID & code \\ 
  \hline
1 &   1 & a \\ 
  2 &   2 & b \\ 
  3 &   3 & c \\ 
   \hline
\end{tabular}
\caption{This is the table caption} 
\label{tab:tab1}
\end{table}

You can reference this table as follows: Table \ref{tab:tab1}.

# Discussion

You can cross-reference sections and subsections as follows: Section \ref{materials-and-methods} and Section \ref{a-subsection}.

***Note:*** the last section in the document will be used as the section title for the bibliography.

# References
