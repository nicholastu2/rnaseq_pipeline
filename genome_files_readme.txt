Genome subdirectories in this directory should be maintained as a common resource.

If you wish to change the gff/gtf, genome, or some other critical aspect of the file,
then the entire codebase should be versioned (i.e. from v1.0 to v1.1). Genome files specific
to a certain version should be maintained in place in order to ensure reproducibility.

Each organism stored in genome_files must contain the following files. The pipeline will pull
files based on their extensions, ie there should be only one gff or gtf in the main
directory structure. Any files with the same extension that are not the files intended for use by
the pipeline should be stored in subdirectories of the organism directory.

Organism (eg KN99)


Annotation file in either .gtf or .gff format (eg crNeo.gtf)
    - this must have metadata denoted with '#' in the top 10 lines describing source, version, 
      and any manipulations performed. If you do manipulate the gtf/gff in any way, include the 
      script in a subdirectory of the organism file along with the original annotation file.

Genome in x format (eg )
    - Same as above, this should have metadata associated with it. It is not required to be in
      in the top 10 lines. However, if metadata (source, version, etc) is not in the genome file
      itself, include it in a README. This should allow anyone to access the file at its original
      source. Similarly, if any manipulations are performed, include the script and the original
      genome in a subdirectory of this particular organism

gene_name_list.txt (eg )
    - Extracted from the annotation file or genome. Include the source of the gene names or 
      the script used to extract

igv_genome (eg )
    - indexed genome for the igv


