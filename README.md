# alpseq

Nextflow pipeline for analysing (Illumina) 2x300bp nanobody sequencing data. Performs pre-processing and automatically generates a basic analysis report.

<img src="alpseq_logo.png" alt="alpseq logo" width="1280"/>

## Usage

:star: Before running the pipeline, make sure to edit the nextflow.config file to suit your system and point to the correct input and output directories.

:star: For detailed usage instructions, [check out the documentation website](https://kzeglinski.github.io/alpseq_docs/)

The quickest way to run *alpseq* is by running the command below (after editing the `nextflow.config` file):

```         
nextflow run ./main.nf -profile <docker or singularity> 
```

Or by running it through the [seqera cloud platform's](https://cloud.seqera.io/) GUI.

## Useful links

-   [Our documentation website](https://kzeglinski.github.io/alpseq_docs/)

-   [The *alpseq* preprint](https://www.biorxiv.org/content/10.1101/2025.10.22.661623v1)

-   [The data published alongside *alpseq*](https://www.ebi.ac.uk/ena/browser/view/PRJEB90877)