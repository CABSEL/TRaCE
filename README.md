# TRaCE
<img style = "float: right;" src = "https://github.com/CABSEL/TRaCE/blob/master/traceplus.png" align="right"> 
Transitive Reduction and Closure Ensemble
nferring the structure of gene regulatory networks (GRNs) from expression data is a major topic in systems biology. Despite the large number of methods developed for this purpose, such inference is still an open problem. The network inference is often stated to be underdetermined. The underdetermined nature of the network inference problem implies that there could be multiple network solutions to the inference problem, i.e. the solution consists of an ensemble of networks. TRaCE and its extension TRaCE+ were developed using an ensemble inference approach. TRaCE and TRaCE+ generate the upper-bound and lower-bound directed graphs (digraphs) of the GRN ensemble. The difference between the upper and lower bounds (i.e. edges in the upper bound that are not in the lower bound) defines the set of uncertain edges, representing gene regulations that could not be verified by the data. While TRaCE deals with an ensemble of digraphs, TRaCE+ further considers the signs of the edges. A positive edge represents an activation, while a negative edge represents a repression. These ensemble bounds can be used to optimise gene KO experiments (see our tool: [REDUCE](https://github.com/CABSEL/REDUCE)).

TRaCE and TRaCE+ were implemented and tested on MATLAB 2011a platform.

## Last Update
TRaCE:  27.05.2014

TRaCE+: 11.02.2016  

## License
Redistribution and use in source and binary forms, with or without modification, are permitted provided agreeing to the [Simplified BSD Style License](https://github.com/CABSEL/TRaCE/blob/master/Lincense).

Read about Simplified [BSD Style License](http://www.opensource.org/licenses/bsd-license.php)

## Download & Installation
### TRaCE
Download & unzip the [TRaCE.zip](https://github.com/CABSEL/TRaCE/tree/master/TRaCE) file for codes and data. 

Download & unzip the [TraceWithoutData.zip](https://github.com/CABSEL/TRaCE/blob/master/TraceWithoutData.zip) file for codes without data.

### TRaCE+
Download and unzip the [trace-plus-with-data.zip](https://github.com/CABSEL/TRaCE/blob/master/trace-plus-with-data.zip) for codes and data. 

Download and unzip the [trace-plus.zip](https://github.com/CABSEL/TRaCE/blob/master/trace-plus.zip) for codes without data.

## Please use the following references when using TRaCE and TRaCE+
- Ud-Dean, SMM., Heise, S., Klamt, S. and Gunawan, R. TRaCE+: Ensemble inference of gene regulatory networks from transcriptional expression profiles of gene knock-out experiments. BMC Bioinformatics 17, 252 (2016). [abstract](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1137-z)

- Ud-Dean, SMM., and Gunawan, R. Ensemble Inference and Inferability of Gene Regulatory Networks. PLoS ONE, 9(8), e103812 (2014). [abstract](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0103812)

## Acknowledgement
This work is supported by funding from Swiss National Science Foundation.
