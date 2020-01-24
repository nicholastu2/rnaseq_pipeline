# RNA-seq analysis pipeline

This RNA-seq analysis pipeline is designed for Brent lab, Washington University. The intention of use is process data obtained from (TF-encoding) gene perturbation, time series, and/or multi-level environmental treatment/stimulation experiments. It includes modules (1) pre-alignment sample QC, (2) RNA-seq reads alignment, (3) transcriptome expression analysis, (4) post-alignment QC, (5) sample quality audit, and (6) differential expression (DE) analysis. Its optional modules provide more thorough quality analysis and help guide future experimental design. The pipeline aim at minimizing human effort in metadata maintenance, sample quality assessment and DE design, while maximizing the versatility in handling complex experiment design.

##
For usage on the cluster:
```
ml rnaseq_pipeline
```
**If this is the first time you have loaded this environment, it will take about 20 minutes to copy the environment into $USER/.conda/envs/rnaseq_pipeline.**
After loading/copying the module once, any subsequent ```ml rnaseq_pipeline``` will take about 5 seconds to load the environment.

See the brent gitlab miniconda repository for instructions on using the lab distribution of conda.

### Installation outside of the brent lab
Clone this repo to your work station. Please see the wiki for dependencies. 
```
git clone https://gitlab.com/brentlab/rnaseq_pipe.git
```

### Wiki
For detailed intallation, usage, and file descriptions, check out our [wiki page](https://gitlab.com/brentlab/rnaseq_pipe/-/wikis/home).

### Contribution
You are more than welcomed to contribute. Feel free to submit any pull request.
