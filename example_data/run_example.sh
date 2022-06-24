echo "PREPROCESS" && TE_detective preprocess --bam test_sim.bam --ref ref_fofn \
    && echo 'DISCOVER' && TE_detective discover --bam test_sim.bam --ref ref_fofn \
    && echo 'NADISCOVER' && TE_detective nadiscover --bam test_sim.bam --ref ref_fofn --pat --drd 5  \
    && echo 'ANALYZE' && TE_detective analyze --bam test_sim.bam --ref ref_fofn --inp initial_predictions.txt \
    && echo 'CLUSTER2D' && TE_detective cluster2d --bam test_sim.bam --ref ref_fofn \
    && echo 'FILTER' && TE_detective filter --ofa final_results.tsv --bed rmsk_ucsc_mm10.bed
