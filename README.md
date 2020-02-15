# Project Description
This is the repository for Project 1 in BF528. 
This will soon contain files to classify microarray-based tumors.

# Contributors
- Emily Hughes <ehug@bu.edu>
- Simran Makwana <makwana@bu.edu>
- Sumiti Sandhu <sandhu08@bu.edu>
- Michiel Smit <smit2@bu.edu>

# Organization
Recommended organization: project folder with 4 subdirectories
- reference
- samples
- project-1-group-6 repository
- analysis 

# Repository Contents
Code for "Programmer" Role:

Code for "Analyst" Role:
project_1_analyst.Rmd
- assumes results from programmer (expression_data_filtered.csv) and annotations are located in samples directory
- writes output to analysis directory
+ A comma separated file with the filtered results from all three filters from 4.4 (filtered_expression_matrix.csv)
+ A comma separated file with the filtered results from the expression filter from 4.2 (filter2_expression_matrix.csv)
+ A heatmap of the gene-expression of each gene across all samples (heatmap.jpeg)
+ A comma separated file containing the results of the Welch t-test for all genes irrespective of significance for each subtype comparison (welch_results_allgenes.csv)
+ A comma separated file with the t-test results computed on the expression matrix from 4.5. (t_test_filter_2.csv)


