
rm(list = ls())
setwd("/Desktop/Analysis_MS/")

#LFQ

#######RUN WITH WHOLE_LFQ DATASETS######
#Comment: 
##AvsB: ok
##AvsC:empty
##BvsC: ok

library(readxl)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(pathfindR)
library(writexl)

#Load_data
#1_Load LFQ data set----
LFQ_AvsB <- read_excel("20220706_Cachexia_LFQ_proteomics(3)_rawdat_mod_AA.xlsx", 
                       sheet = "AvsB")

LFQ_AvsC <- read_excel("20220706_Cachexia_LFQ_proteomics(3)_rawdat_mod_AA.xlsx", 
                       sheet = "AvsC")

LFQ_BvsC <- read_excel("20220706_Cachexia_LFQ_proteomics(3)_rawdat_mod_AA.xlsx", 
                       sheet = "BvsC")


#1.1_Data preparation of AvsB-------
LFQ_AvsB_path_test <- LFQ_AvsB %>% 
  select(Gene_names, log2_FC_A_to_B, negative_log_P_value) %>% 
  rename(Gene.symbol = Gene_names) %>% 
  rename(adj.P.val = negative_log_P_value) %>% 
  rename(logFC = log2_FC_A_to_B) %>% 
  mutate(adj.P.val = 10**(-adj.P.val))
       
View(LFQ_AvsB)    
View(LFQ_AvsB_path_test)

#class(LFQ_AvsB_path_test$Gene.symbol) #char
#class(LFQ_AvsB_path_test$logFC) #num
#class(LFQ_AvsB_path_test$adj.P.val) #num

#1.2_Data preparation of AvsC-------
LFQ_AvsC_path_test <- LFQ_AvsC %>% 
  select(Gene_names, log2_FC_A_to_C, negative_log_P_value) %>% 
  rename(Gene.symbol = Gene_names) %>% 
  rename(adj.P.val = negative_log_P_value) %>% 
  rename(logFC = log2_FC_A_to_C) %>% 
  mutate(adj.P.val = 10**(-adj.P.val))

View(LFQ_AvsC)    
View(LFQ_AvsC_path_test)

#class(LFQ_AvsB_path_test$Gene.symbol) #char
#class(LFQ_AvsB_path_test$logFC) #num
#class(LFQ_AvsB_path_test$adj.P.val) #num

#1.3_Data preparation of BvsC-------
LFQ_BvsC_path_test <- LFQ_BvsC %>% 
  select(Gene_names, log2_FC_B_to_C, negative_log_P_value) %>% 
  rename(Gene.symbol = Gene_names) %>% 
  rename(adj.P.val = negative_log_P_value) %>% 
  rename(logFC = log2_FC_B_to_C) %>% 
  mutate(adj.P.val = 10**(-adj.P.val))

View(LFQ_BvsC)    
View(LFQ_BvsC_path_test)

#class(LFQ_AvsB_path_test$Gene.symbol) #char
#class(LFQ_AvsB_path_test$logFC) #num
#class(LFQ_AvsB_path_test$adj.P.val) #num


############################################################AvsB----
#Comments:1310-->263 because of p values/no significant up or down reg.

RA_processed_AvsB <- input_processing(input = LFQ_AvsB_path_test, # the input: in this case, differential expression results
                                 p_val_threshold = 0.05, # p value threshold to filter significant genes
                                 pin_name_path  = "KEGG", # the name of the PIN to use for active subnetwork search
                                 convert2alias = TRUE)

#already installed because of the error
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

biocarta_list_AvsB <- fetch_gene_set(gene_sets = "KEGG",
                                min_gset_size = 10,
                                max_gset_size = 300)
biocarta_gsets_AvsB <- biocarta_list_AvsB[[1]]
biocarta_descriptions_AvsB <- biocarta_list_AvsB[[2]]

n_iter <- 10 ## number of iterations
combined_res_AvsB <- NULL ## to store the result of each iteration

for (i in 1:n_iter) {
  
  ###### Active Subnetwork Search
  snws_file_AvsB <- paste0("active_snws_", i) # Name of output file
  active_snws_AvsB <- active_snw_search(input_for_search = RA_processed_AvsB, 
                                   pin_name_path = "KEGG", 
                                   snws_file = snws_file_AvsB,
                                   score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                   sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                   search_method = "GR")
  
  ###### Enrichment Analyses
  current_res_AvsB <- enrichment_analyses(snws = active_snws_AvsB,
                                     sig_genes_vec = RA_processed_AvsB$GENE,
                                     pin_name_path = "KEGG", 
                                     genes_by_term = biocarta_gsets_AvsB,
                                     term_descriptions = biocarta_descriptions_AvsB,
                                     adj_method = "bonferroni",
                                     enrichment_threshold = 0.05,
                                     list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res_AvsB <- rbind(combined_res_AvsB, current_res_AvsB)
}

View(combined_res_AvsB)


###### Summarize Combined Enrichment Results
summarized_df_AvsB <- summarize_enrichment_results(combined_res_AvsB, 
                                              list_active_snw_genes = TRUE)

View(summarized_df_AvsB)
###### Annotate Affected Genes Involved in Each Enriched Term
final_res_AvsB <- annotate_term_genes(result_df = summarized_df_AvsB, 
                                 input_processed = RA_processed_AvsB, 
                                 genes_by_term = biocarta_gsets_AvsB)
View(final_res_AvsB)
##Visualize
visualize_terms(result_df = final_res_AvsB, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "KEGG")

head(final_res_AvsB)
View(final_res_AvsB)
enrichment_chart(final_res_AvsB[1:5, ]) #worked


visualize_active_subnetworks(final_res_AvsB) #error

term_gene_heatmap(final_res_AvsB) #worked
visualize_term_interactions(final_res_AvsB) #error
visualize_terms(final_res_AvsB) #error
UpSet_plot(final_res_AvsB) #worked heat and conncetion
term_gene_graph(final_res_AvsB) #worked netz

write_xlsx(summarized_df_AvsB, "/summarized_df_AvsB.xlsx")



############################################################BvsC-----
#2.3_pathfindR in BvsC----1327-->251 because of p, 251-->114 because of PIN


RA_processed_BvsC <- input_processing(input = LFQ_BvsC_path_test, # the input: in this case, differential expression results
                                      p_val_threshold = 0.05, # p value threshold to filter significant genes
                                      pin_name_path  = "KEGG", # the name of the PIN to use for active subnetwork search
                                      convert2alias = TRUE)

#already installed because of the error
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

biocarta_list_BvsC <- fetch_gene_set(gene_sets = "KEGG",
                                     min_gset_size = 10,
                                     max_gset_size = 300)
biocarta_gsets_BvsC <- biocarta_list_BvsC[[1]]
biocarta_descriptions_BvsC <- biocarta_list_BvsC[[2]]

n_iter <- 10 ## number of iterations
combined_res_BvsC <- NULL ## to store the result of each iteration

for (i in 1:n_iter) {
  
  ###### Active Subnetwork Search
  snws_file_BvsC <- paste0("active_snws_", i) # Name of output file
  active_snws_BvsC <- active_snw_search(input_for_search = RA_processed_BvsC, 
                                        pin_name_path = "KEGG", 
                                        snws_file = snws_file_BvsC,
                                        score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                        sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                        search_method = "GR")
  
  ###### Enrichment Analyses
  current_res_BvsC <- enrichment_analyses(snws = active_snws_BvsC,
                                          sig_genes_vec = RA_processed_BvsC$GENE,
                                          pin_name_path = "KEGG", 
                                          genes_by_term = biocarta_gsets_BvsC,
                                          term_descriptions = biocarta_descriptions_BvsC,
                                          adj_method = "bonferroni",
                                          enrichment_threshold = 0.05,
                                          list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res_BvsC <- rbind(combined_res_BvsC, current_res_BvsC)
}

View(combined_res_BvsC)


###### Summarize Combined Enrichment Results
summarized_df_BvsC <- summarize_enrichment_results(combined_res_BvsC, 
                                                   list_active_snw_genes = TRUE)

View(summarized_df_BvsC)
###### Annotate Affected Genes Involved in Each Enriched Term
final_res_BvsC <- annotate_term_genes(result_df = summarized_df_BvsC, 
                                      input_processed = RA_processed_BvsC, 
                                      genes_by_term = biocarta_gsets_BvsC)
View(final_res_BvsC)
##Visualize
visualize_terms(result_df = final_res_BvsC, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "KEGG")

head(final_res_BvsC)
View(final_res_BvsC)
enrichment_chart(final_res_BvsC[1:5, ]) #worked


visualize_active_subnetworks(final_res_BvsC) #error

term_gene_heatmap(final_res_BvsC) #worked
visualize_term_interactions(final_res_BvsC) #error
visualize_terms(final_res_BvsC) #error
UpSet_plot(final_res_BvsC) #worked heat and conncetion
term_gene_graph(final_res_BvsC) #worked netz

write_xlsx(summarized_df_BvsC, "/summarized_df_BvsC.xlsx")

############################################################AvsC----
#empty - no significant gene expression change

RA_processed_AvsC <- input_processing(input = LFQ_AvsC_path_test, # the input: in this case, differential expression results
                                      p_val_threshold = 0.05, # p value threshold to filter significant genes
                                      pin_name_path  = "KEGG", # the name of the PIN to use for active subnetwork search
                                      convert2alias = TRUE)




biocarta_list_AvsC <- fetch_gene_set(gene_sets = "KEGG",
                                     min_gset_size = 10,
                                     max_gset_size = 300)
biocarta_gsets_AvsC <- biocarta_list_AvsC[[1]]
biocarta_descriptions_AvsC <- biocarta_list_AvsC[[2]]

n_iter <- 10 ## number of iterations
combined_res_AvsC <- NULL ## to store the result of each iteration

for (i in 1:n_iter) {
  
  ###### Active Subnetwork Search
  snws_file_AvsC <- paste0("active_snws_", i) # Name of output file
  active_snws_AvsC <- active_snw_search(input_for_search = RA_processed_AvsC, 
                                        pin_name_path = "KEGG", 
                                        snws_file = snws_file_AvsC,
                                        score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                        sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                        search_method = "GR")
  
  ###### Enrichment Analyses
  current_res_AvsC <- enrichment_analyses(snws = active_snws_AvsC,
                                          sig_genes_vec = RA_processed_AvsC$GENE,
                                          pin_name_path = "KEGG", 
                                          genes_by_term = biocarta_gsets_AvsC,
                                          term_descriptions = biocarta_descriptions_AvsC,
                                          adj_method = "bonferroni",
                                          enrichment_threshold = 0.05,
                                          list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
  
  ###### Combine results via `rbind`
  combined_res_AvsC <- rbind(combined_res_AvsC, current_res_AvsC)
}

View(combined_res_AvsC)


###### Summarize Combined Enrichment Results
summarized_df_AvsC <- summarize_enrichment_results(combined_res_AvsC, 
                                                   list_active_snw_genes = TRUE)

View(summarized_df_AvsC)
###### Annotate Affected Genes Involved in Each Enriched Term
final_res_AvsC <- annotate_term_genes(result_df = summarized_df_AvsC, 
                                      input_processed = RA_processed_AvsC, 
                                      genes_by_term = biocarta_gsets_AvsC)
View(final_res_AvsC)
##Visualize
visualize_terms(result_df = final_res_AvsC, 
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "KEGG")

head(final_res_AvsC)
View(final_res_AvsC)
enrichment_chart(final_res_AvsC[1:5, ]) #worked


visualize_active_subnetworks(final_res_AvsC) #error

term_gene_heatmap(final_res_AvsC) #worked
visualize_term_interactions(final_res_AvsC) #error
visualize_terms(final_res_AvsC) #error
UpSet_plot(final_res_AvsC) #worked heat and conncetion
term_gene_graph(final_res_AvsC) #worked netz

write_xlsx(summarized_df_AvsC, "/summarized_df_AvsC.xlsx")
