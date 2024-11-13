docker load -i r_with_plink.tar.gz

gunzip bcac_asian_icogs_onco_erneg_meta1.txt.gz
gunzip bcac_asian_icogs_onco_erpos_meta1.txt.gz

dx download PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/w_hm3_updated.snplist

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript HM3_Sumstats_EAS_Latina.R"