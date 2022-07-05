echo "PREPROCESS" && TE_detective preprocess -i test_sim.bam -r ref_fofn \
    && echo 'DISCOVER' && TE_detective discover -i test_sim.bam -r ref_fofn \
    && echo 'NADISCOVER' && TE_detective nadiscover -i test_sim.bam -r ref_fofn --polyA --discord_cluster_dens 5  --polyA \
    && echo 'ANALYZE' && TE_detective analyze -i test_sim.bam -r ref_fofn --inp initial_predictions.txt \
    && echo 'CLUSTER2D' && TE_detective cluster2d -i test_sim.bam -r ref_fofn --discord_cluster_dens 5 --all \
    && echo 'FILTER' && TE_detective filter -i final_results.tsv --bed rmsk_ucsc_mm10.bed
