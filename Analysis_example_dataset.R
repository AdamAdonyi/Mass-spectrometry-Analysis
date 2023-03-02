#EXAMPLE DATA SET of pathfindR 
#DATA PREP
#VISUAL
##worked very well

#original dataset analysis is based on Biogrid
#masspec was tested with biogrid - KEGG is better (more identified genes based on nomenclature)
#additional Java was already installed but update was needed
#additional package update was needed - commented out in line 31 (BiocManager)


#EXAMPLE DATASET:
##load data_RA_input
data("RA_input")
#check_head
head(RA_input)

class(RA_input$Gene.symbol)
class(RA_input$logFC)
class(RA_input$adj.P.Val)


input_processing()

RA_processed <- input_processing(input = RA_input, # the input: in this case, differential expression results
                                 p_val_threshold = 0.05, # p value threshold to filter significant genes
                                 pin_name_path  = "Biogrid", # the name of the PIN to use for active subnetwork search
                                 convert2alias = TRUE)

#already installed because of the error
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")


biocarta_list <- fetch_gene_set(gene_sets = "BioCarta",
                                min_gset_size = 10,
                                max_gset_size = 300)
biocarta_gsets <- biocarta_list[[1]]
biocarta_descriptions <- biocarta_list[[2]]

n_iter <- 10 ## number of iterations
combined_res <- NULL ## to store the result of each iteration

for (i in 1:n_iter) {
  
  ###### Active Subnetwork Search
  snws_file <- paste0("active_snws_", i) # Name of output file
  active_snws <- active_snw_search(input_for_search = RA_processed, 
                                   pin_name_path = "Biogrid", 
                                   snws_file = snws_file,
                                   score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                   sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                   search_method = "GR")
  
  ###### Enrichment Analyses
  current_res <- enrichment_analyses(snws = active_snws,
                                     sig_genes_vec = RA_processed$GENE,
                                     pin_name_path = "Biogrid", 
                                     genes_by_term = biocarta_gsets,
                                     term_descriptions = biocarta_descriptions,
                                     adj_method = "bonferroni",
                                     enrichment_threshold = 0.05,
                                     list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res <- rbind(combined_res, current_res)
}

#View(combined_res)


###### Summarize Combined Enrichment Results
summarized_df <- summarize_enrichment_results(combined_res, 
                                              list_active_snw_genes = TRUE)

View(summarized_df)
###### Annotate Affected Genes Involved in Each Enriched Term
final_res <- annotate_term_genes(result_df = summarized_df, 
                                 input_processed = RA_processed, 
                                 genes_by_term = biocarta_gsets)

##Visualize
visualize_terms(result_df = final_res, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "Biogrid")

head(final_res)
View(final_res)



visualize_active_subnetworks(final_res) #worked - separated folder

term_gene_heatmap(final_res) #worked


UpSet_plot(final_res) #worked
term_gene_graph(final_res) #worked
