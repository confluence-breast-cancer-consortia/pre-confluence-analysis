docker load -i r_with_plink.tar.gz

dx download PreConfluenceAnalysis:/Aim2_Polygenicity/HZ/Results/Clean_summary_data/ -r
dx download PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/snpinfo_hm3_mega.csv

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript HM3_MEGA_Sumstats_EUR_EAS.R"