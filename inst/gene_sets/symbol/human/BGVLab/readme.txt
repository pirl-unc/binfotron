Repo: human_hgnc_bgvlab v0.1-01

Started with gmt file last modified on 2020-04-09 14:00:01 from source: 
/home/dbortone/R/x86_64-pc-linux-gnu-library/3.5/binfotron/gene_sets/entrez_ids/human/BGVLab/BGVLab_entrez_human.gmt.txt

> BGVLab_IDs contains gene signatures used by the Vincent Lab to assess human
> immune gene signatures. They have been assembled from the following sources:  
> 
>   Iglesia MD et al. Clinical Cancer Research 2014; PMID:24916698
>   Bindea G et al. Immunity 2013; PMID:24138885
>   Fan C et al. BMC Med Genomics 2011; PMID:21214954
>   Palmer C et al. BMC Genomics 2006; PMID:16704732
>   Schmidt M et al Cancer Res 2008; PMID:18593943
>   Rody A et al. Breast Cancer Res 2011; PMID:21978456
>   Prat A et al Breast Cancer Res 2010; PMID:20813035
>   Kardos J and Chai S et al. JCI insights 2016; PMID:27699256
>   Beck et al clin Cancer Research 2009; PMID:19188147
>   Chan et al PNAS 2009; PMID:19666525
>   IPRES derrived Hugo/Lo 2016; PMID26997480
>   Thorsson V et al. Immunity 2018; PMID: 29628290
>   Roufas C et al. Front Oncol 2018
>   And Gene Ontology terms

On Mon Apr 20 07:41:58 PM 2020 converted the gmt file from ENTREZID to SYMBOL 
using binfotron::convert_gmt_file v0.3.1


20210622 - We were able to do 1:1 mapping of transcripts using the gtf file.

We were loosing about 6.7% of our Vincent Lab genes due to missing the genes in 
our transcriptome -> hgnc conversion
After the corrections this dropped to 0.9% (some just couldn't be tracked down)
We were loosing another 5.5% when we restricted ourselves to protein_coding genes
So by going to the new hngc gene signatures we'll go from around 13% of the genes 
missing to 0.9%.

To get the most genes for your signatures here we recommend using Ensembl 
transcript ids and matching to HGNC symbols not to these entrez id signatures. 
Additionally virtually all gene biotypes need to be used to get all of these 
genes. 


Genes biotype found using hgnc ids:

      protein_coding                    3443
      lncRNA                              93
      IG_V_gene                           14
      IG_C_gene                           11
      snoRNA                               9
      transcribed_unprocessed_pseudogene   8
      snRNA                                6
      transcribed_unitary_pseudogene       5
      TR_C_gene                            5
      TR_V_gene                            4
      IG_V_pseudogene                      4
      unprocessed_pseudogene               2
      processed_pseudogene                 2
      polymorphic_pseudogene               2
      transcribed_processed_pseudogene     1
      TR_V_pseudogene                      1
      IG_J_gene                            1
      IG_C_pseudogene                      1