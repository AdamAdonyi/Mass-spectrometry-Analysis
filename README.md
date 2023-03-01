# Mass-spectrometry-Analysis
Comparison of 3 different groups of mouse-derived heart samples using [Relative label-free quantification (LFQ)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7916846/) or [Redox](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5327050/#:~:text=Redox%20proteomics%20is%20that%20branch,of%20identifying%20the%20target%20proteins.) Mass spectomety methods.
The animal groups are A, B, and C. During comparison all were compared such as AvsB, AvsC, and BvsC. A general background of the animals  [cachexia](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7916846/) then in addational treatment or genetic modification was introduced (coded as A, B and C) but since it is an unpublished data, I am not allowed to reveal the exact background. Upon publication, I am going to share the link to the original paper.

# Main steps

## Test analysis with example dataset
Using pathfindR Enrichment Workflow, I used the RA_input example dataset:
Most important is the dataformat which is recognized by the command (all dataset must be transformed into that format). Gene.symbol, logFC, and adj.P.Val.
The protein-protein interaction network (PIN) pin_name_path can be: “Biogrid”, “STRING”, “GeneMania”, “IntAct”, “KEGG”, “mmu_STRING”. Furthermore, the available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”.

## Adjusting tested code to dataset and use for analysis

## Results
