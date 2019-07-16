Saved a use_mart result from executing:

  mart<-useMart(dataset="hsapiens_gene_ensembl", biomart='ENSEMBL_MART_ENSEMBL')
  BM_results = getBM(
    filters= "ucsc", 
    attributes= c("ucsc", "hgnc_symbol", "entrezgene", "gene_biotype"),
    values= unique_names,
    mart= mart
  )
  saveRDS(BM_results, file = bm_results_path)
  
Where "unique_names" came from a star/salmon run with hg38.

This was run in November of 2018.