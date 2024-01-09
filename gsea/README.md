Full workflow for running gsea with the fgsea package in R. 
1. First, Run Differential expression Analysis for chosen gene sets
2. Create volcano plots from DE data
3. Use that data for the ranked list of gene input to run fgsea
4. Extract common pathways for DE comparisons
5. Extract common leading edge genes from those pathways
6. Count up the common genes' occurences in different pathways for each gene set and for all the common gene sets combined
7. create linked HTML files with gene counts, with the headings linking to the pathway & leading edge gene data for each gene set so
    that you can see which gene set(s) the genes are linked to from the HTML file.
