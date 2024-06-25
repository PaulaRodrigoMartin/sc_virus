#Snakefile, place in 

CELLTYPES = ["Secretory","Macro","T","DC","Ciliated","Neu", "NK","Plasma","Squamous", "Mono"]
HYPOTHESES = ["hypo1","hypo2", "hypo3"]
EPITHELIAL = ["Secretory","Ciliated","Squamous"]


rule all:
    input:
        # expand("log/longermer_upreg_filtered_{celltype}.rds", celltype = CELLTYPES),
        # expand("log/longermer_downreg_filtered_{celltype}.rds", celltype = CELLTYPES),
        # expand("res/filtered_7mers_hypo1_{celltype}.rds", celltype = CELLTYPES),
        # expand("res/filtered_7mers_hypo2_{celltype}.rds", celltype = CELLTYPES),
        # expand("res/filtered_7mers_hypo3_{celltype}.rds", celltype = CELLTYPES),
        # expand("res/filtered_table_hypo1_{celltype}.png", celltype = CELLTYPES),
        # expand("res/filtered_table_hypo2_{celltype}.png", celltype = CELLTYPES),
        # expand("res/filtered_table_hypo3_{celltype}.png", celltype = CELLTYPES),
        # expand("res/common_7mers_{hypo}.rds", hypo=HYPOTHESES),
        # expand("summary_plots/plots_{hypo}_{celltype}.pdf", hypo=HYPOTHESES, celltype = CELLTYPES),
        # expand("res/venn_7mers_{hypo}.rds", hypo = HYPOTHESES),
        # expand("log/stats_{celltype}.pdf", celltype = CELLTYPES)
        # expand("summary_plots/plots_{hypo}_{celltype}.pdf", hypo = HYPOTHESES, celltype = CELLTYPES),
        # expand("summary_plots/plots_{hypo}_{celltype}.png", hypo = HYPOTHESES, celltype = CELLTYPES)
        expand("res/filtered_7mers_hypo3_{celltype}.rds", celltype = CELLTYPES),
        expand("res/filtered_table_hypo3_{celltype}.png", celltype = CELLTYPES)

##########################
## ANNOTATION OF 7MERS
##########################

# Annotate variables that are not cell type specific (fex miRNA)
rule general_annot:
    input:
        activity_mireact = "res/mirnaActivity_th.rds",
        annotations =  "data/annotation_patients_more.csv",
        mir_annotations =  "data/humir.rds",
        rbp_annot = "data/rbps_geneyeo.RDS",
        occurrences = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
        sars = "data/EPI_ISL_402124.fasta",
        script_init_annot = "code/flow/init_annot.R"
    output:
        celltypes = expand("log/df_{celltype}.h5ad", celltype = CELLTYPES)
        #folder = "res"
    threads:
        1
    resources:
        mem_mb = 8000,
        runtime = 540
    shell:"""
	Rscript {input.script_init_annot} {input.activity_mireact} {input.annotations} {input.mir_annotations} {input.rbp_annot} {input.occurrences} {input.sars} {output.celltypes}
	"""

# Annotate variables that are not cell type specific (fex differential activity between infecetd and healthy cells)
rule cell_specific_annot:
    input:
        celltype = "log/df_{celltype}.h5ad",
        occurrences = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
        script_cell_annot = "code/flow/cell_annot.R"
    output:
        cellannot1 = "log/annot1_{celltype}.h5ad"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 540
    shell:"""
	Rscript {input.script_cell_annot} {input.celltype} {input.occurrences} {output.cellannot1}
	"""

# Annotate correlation with viral load
# run if there's data for viral load in cells, else, comment
rule viralload:
    input:
        celltype = "log/annot1_{celltype}.h5ad",
        script_fur_cell_annot = "code/flow/fur_cell_annot.R"
    output:
        cellannot2 = "log/annot2_{celltype}.h5ad"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 540
    shell:"""
	Rscript {input.script_fur_cell_annot} {input.celltype} {output.cellannot2}
	"""

# Annotate gene ontology for all cell types
rule goterms:
    input:
        celltypesss = "log/annot2_{celltype}.h5ad",
        script_filtering = "code/flow/filtering.r"
    output:
        filtered = "log/annot3_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {output.filtered}
	"""

# Annotate hits with SARS-CoV-2 genome (hypothesis )
rule sarshit:
    input:
        celltypesss = "log/annot3_{celltype}.rds",
        sars = "data/EPI_ISL_402124.fasta",
        script_filtering = "code/flow/sarshit.r"
    output:
        filtered = "log/annot4_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {input.sars} {output.filtered}
	"""

# Hits with the reverse complementary sequence of the virus
rule sarshit_rev:
    input:
        celltypesss = "log/annot4_{celltype}.rds",
        sars = "data/EPI_ISL_402124.fasta",
        script_filtering = "code/flow/sarshit_rev.r"
    output:
        filtered = "log/annot5_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {input.sars} {output.filtered}
	"""

# Ranks kmers based on mean rank of above calculated variables
rule rank:
    input:
        celltypesss = "log/annot5_{celltype}.rds",
        script_filtering = "code/flow/rank.r"
    output:
        rank1 = "res/ranked_7mers_hypo1_{celltype}.rds",
        rank2 = "res/ranked_7mers_hypo2_{celltype}.rds",
        rank3 = "res/ranked_7mers_hypo3_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {output.rank1} {output.rank2} {output.rank3}
	"""


## Summary statistics of the ranked and summary lists
rule rank_indiv_1:
    input:
        rank1 = "res/ranked_7mers_hypo1_{celltype}.rds",
        script_filtering = "code/flow/rank_indiv_1.r"
    output:
        rank_corr = "res/filtered_7mers_hypo1_{celltype}.rds",
        rank_pval = "res/filtered_table_hypo1_{celltype}.png"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.rank1} {output.rank_corr} {output.rank_pval} 
	"""

rule rank_indiv_2:
    input:
        rank1 = "res/ranked_7mers_hypo2_{celltype}.rds", 
        script_filtering = "code/flow/rank_indiv_2.r"
    output:
        rank_corr = "res/filtered_7mers_hypo2_{celltype}.rds", 
        rank_pval = "res/filtered_table_hypo2_{celltype}.png"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
    Rscript {input.script_filtering} {input.rank1} {output.rank_corr} {output.rank_pval} 
    """

rule rank_indiv_3:
    input:
        rank1 = "res/ranked_7mers_hypo3_{celltype}.rds",
        script_filtering = "code/flow/rank_indiv_3.r"
    output:
        rank_corr = "res/filtered_7mers_hypo3_{celltype}.rds", 
        rank_pval = "res/filtered_table_hypo3_{celltype}.png"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.rank1} {output.rank_corr} {output.rank_pval} 
	"""

rule process_hypothesis:
    input:
        lambda hypo: expand("res/filtered_7mers_{hypo}_{celltype}.rds", hypo=hypo, celltype=CELLTYPES),
        script_filtering = "code/flow/common.R"
    output:
        "res/common_7mers_{hypo}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:
        """
        Rscript {input.script_filtering} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {output}
        """

rule general_annot_fig:
    input:
        rank = "res/filtered_7mers_{hypo}_{celltype}.rds",
        occurrencess = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
        mir_annotationss =  "data/humir.rds",
        seqlist = "data/hs.utr3.seqlist.rds",
        activity_mireact = "res/mirnaActivity_th.rds",
        annot = "data/annotation_patients_more.csv",
        script_init_fig = "code/flow/plot_gen.R"
    output:
        figures_general = "summary_plots/plots_{hypo}_{celltype}.pdf",
        fig_png = "summary_plots/plots_{hypo}_{celltype}.png"
    threads:
        1
    resources:
        mem_mb = 15000,
        runtime = 180
    shell:"""
	Rscript {input.script_init_fig} {input.rank} {input.occurrencess} {input.mir_annotationss} {input.seqlist} {input.activity_mireact} {input.annot} {output.figures_general} {output.fig_png}
	"""

# Venn diagram for epithelial cells
# rule venn:
#     input:
#         lambda hypo: expand("res/filtered_7mers_{hypo}_{epithelial}.rds", hypo=hypo, epithelial=EPITHELIAL),
#         annot = "data/annotation_patients_more.csv",
#         mir_annotationss =  "data/humir.rds",
#         script_filtering = "code/flow/venn.r"
#     output:
#         "res/venn_7mers_{hypo}.rds"
#     threads:
#         8
#     resources:
#         mem_mb = 100000,
#         runtime = 600
#     shell:
#         """
#         Rscript {input.script_filtering} {input[0]} {input[1]} {input[2]} {input.annot} {input.mir} {output}
#         """

##########################
## CONSTRUCT LONGERMERS
##########################

rule longermer_up:
    input:
        celltypesss = "log/annot4_{celltype}.rds",
        script_filtering = "code/flow/longermer_up.R"
    output:
        filtered_up = "log/longermer_upreg_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 10000,
        runtime = 5000
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {output.filtered_up}
	"""

rule longermer_down:
    input:
        celltypesss = "log/annot4_{celltype}.rds",
        script_filtering = "code/flow/longermer_down.R"
    output:
        filtered_down = "log/longermer_downreg_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 10000,
        runtime = 5000
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {output.filtered_down}
	"""

# Check biological relevance of the longermers
rule check_longermer_up:
    input:
        upreg = "log/longermer_upreg_{celltype}.rds",
        celltypesss = "log/annot4_{celltype}.rds",
        activity_mireact = "res/mirnaActivity_th.rds",
        annotations =  "data/annotation_patients_more.csv",
        script_filtering = "code/flow/check_longermers_up.R"
    output:
        filtered_up = "log/longermer_upreg_filtered_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 10000,
        runtime = 2100
    shell:"""
	Rscript {input.script_filtering} {input.upreg} {input.celltypesss} {input.activity_mireact} {input.annotations} {output.filtered_up} 
	"""

rule check_longermer_down:
    input:
        downreg = "log/longermer_downreg_{celltype}.rds",
        celltypesss = "log/annot4_{celltype}.rds",
        activity_mireact = "res/mirnaActivity_th.rds",
        annotations =  "data/annotation_patients_more.csv",
        script_filtering = "code/flow/check_longermers_down.R"
    output:
        filtered_down = "log/longermer_downreg_filtered_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 10000,
        runtime = 2100
    shell:"""
	Rscript {input.script_filtering} {input.downreg} {input.celltypesss} {input.activity_mireact} {input.annotations} {output.filtered_down}
	"""

# Check hits with SARS-CoV-2 of the longmers (both strands)
rule annot_long:
    input:
        celltypesss = "log/longermer_downreg_filtered_{celltype}.rds",
        sars = "data/EPI_ISL_402124.fasta",
        script_filtering = "code/flow/sarshit_long.r"
        ##### need more input around here
    output:
        filtered = "log/annot_longmer_{celltype}.rds"
    threads:
        8
    resources:
        mem_mb = 100000,
        runtime = 600
    shell:"""
	Rscript {input.script_filtering} {input.celltypesss} {input.sars} {output.filtered}
	"""
#########################################################

# ## Merge all ranked lists and rank them based on abundancy among celltypes and among mean rank
# rule analysis_hypo1:
#     input:
#         celltypes_f = expand("res/ranked_7mers_hypo1_{celltype}.rds", celltype = CELLTYPES),
#         mir_annotations =  "data/humir.rds",
#         rbp_annot = "data/rbps_geneyeo.RDS",
#         occurrences = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
#         sars = "data/EPI_ISL_402124.fasta",
#         script_analysis = "code/flow/analysis1.R"
#     output:
#         analysiss = "res/candidates_hypo1.rds",
#         analysiss1 = "res/candidates_hypo1_abundant.rds"
#     threads:
#         8
#     resources:
#         mem_mb = 100000,
#         runtime = 540
#     shell:"""
# 	Rscript {input.script_analysis} {input.celltypes_f} {input.mir_annotations} {input.rbp_annot} {input.occurrences} {input.sars} {output.analysiss} {output.analysiss1}
# 	"""

# rule analysis_hypo2:
#     input:
#         celltypes_f = expand("res/ranked_7mers_hypo2_{celltype}.rds", celltype = CELLTYPES),
#         mir_annotations =  "data/humir.rds",
#         rbp_annot = "data/rbps_geneyeo.RDS",
#         occurrences = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
#         sars = "data/EPI_ISL_402124.fasta",
#         script_analysis = "code/flow/analysis2.R"
#     output:
#         analysiss = "res/candidates_hypo2.rds",
#         analysiss1 = "res/candidates_hypo2_abundant.rds"
#     threads:
#         8
#     resources:
#         mem_mb = 100000,
#         runtime = 540
#     shell:"""
# 	Rscript {input.script_analysis} {input.celltypes_f} {input.mir_annotations} {input.rbp_annot} {input.occurrences} {input.sars} {output.analysiss} {output.analysiss1}
# 	"""

# rule analysis_hypo3:
#     input:
#         celltypes_f = expand("res/ranked_7mers_hypo3_{celltype}.rds", celltype = CELLTYPES),
#         mir_annotations =  "data/humir.rds",
#         rbp_annot = "data/rbps_geneyeo.RDS",
#         occurrences = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
#         sars = "data/EPI_ISL_402124.fasta",
#         script_analysis = "code/flow/analysis3.R"
#     output:
#         analysiss = "res/candidates_hypo3.rds",
#         analysiss1 = "res/candidates_hypo3_abundant.rds"
#     threads:
#         8
#     resources:
#         mem_mb = 100000,
#         runtime = 540
#     shell:"""
# 	Rscript {input.script_analysis} {input.celltypes_f} {input.mir_annotations} {input.rbp_annot} {input.occurrences} {input.sars} {output.analysiss} {output.analysiss1}
# 	"""




# rule to make celltype specific summary stats --> distrib pvals, etc
rule summ_stats:
    input:
        celltype = "log/annot2_{celltype}.h5ad",
        script_sumstats_annot = "code/flow/plot_temp.R"
    output:
        cellannot3 = "log/stats_{celltype}.pdf"
    threads:
        1
    resources:
        mem_mb = 30000,
        runtime = 120
    shell:"""
	Rscript {input.script_sumstats_annot} {input.celltype} {output.cellannot3}
	"""

# rule analyze_sep:
#     input:
#         rank1 = "log/ranked_7mers_hypo1_{celltype}.rds",
#         rank2 = "log/ranked_7mers_hypo2_{celltype}.rds",
#         rank3 = "log/ranked_7mers_hypo3_{celltype}.rds",
#         script_filtering = "code/flow/ana.r"
#     output:
#         fig1 = "log/fig_7mers_hypo1_{celltype}.pdf",
#         fig2 = "log/fig_7mers_hypo2_{celltype}.pdf",
#         fig3 = "log/fig_7mers_hypo3_{celltype}.pdf"
#     threads:
#         8
#     resources:
#         mem_mb = 100000,
#         runtime = 600
#     shell:"""
# 	Rscript {input.script_filtering} {input.celltypesss} {output.rank1} {output.rank2} {output.rank3}
# 	"""