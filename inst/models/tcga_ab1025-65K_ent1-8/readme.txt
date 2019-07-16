This model was built around 95th percentile of tcga sample entropies (0.975-8.1)
over the range of tcga tcr abundances (1024) and below the abundance for each entropy that
would be sufficient to go with uncorrected entropy (44139).

The entropy doesn't drop below 4 for these because 1025 is a high enough abundance to get the exact entropy 
if tht entropy is below 4.

See optomize_diversity_metrics project: 
* exec_elastic_net.R
* plot_fraction_entropy.R
* rscripts_doneas_jobs/subsample_stig_rnaseq_high_abundance_gt.R
* temp_data/subsampled_stig_rnaseq_tcga_distribution_highab6