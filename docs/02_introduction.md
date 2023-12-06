<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This workflow is primarily used to assemble genomes from bacterial reads and provide information on features of interest within those assemblies through annotations.

The workflow can provide additional information about the assembly, such as antimicrobial resistance (AMR) analysis and sequence typing through an optional `--isolates` mode. 

In brief, this workflow will perform the following: 

+ De novo (or reference-based) assembly of bacterial genomes 
+ Annotation of regions of interest within the assembly
+ Species identification and sequence typing (`--isolates` mode only)
+ Identify genes and SNVs associated with AMR (`--isolates` mode only)