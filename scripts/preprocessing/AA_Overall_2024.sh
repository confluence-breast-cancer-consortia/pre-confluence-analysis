docker load -i r_with_plink.tar.gz

# mount PWD and run bigsnpr using my_r_script.r
docker run -v $PWD:/data -w /data --entrypoint /bin/bash r_with_plink -l -c "Rscript AA_Overall_2024.R"