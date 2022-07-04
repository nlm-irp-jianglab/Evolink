suppressPackageStartupMessages({
    suppressWarnings(library(tidyverse))
    suppressWarnings(library(treeio))
    suppressWarnings(library(ggtree))
    suppressWarnings(library(optparse))
})
 
# Usage: Rscript --vanilla Evolink_plot.R -g gene.tsv -t trait.tsv -n tree.nwk -r output_dir/result.tsv -c 3.0 -a 9 -b 8 -m 1 -o output_dir

option_list = list(
    make_option(c("-g", "--gene"), type="character", default="gene.tsv", 
              help="gene file", metavar="character"),
    make_option(c("-t", "--trait"), type="character", default="trait.tsv", 
              help="trait file", metavar="character"),
    make_option(c("-n", "--tree"), type="character", default="tree.nwk", 
              help="tree file", metavar="character"),
    make_option(c("-r", "--result"), type="character", default="results.tsv", 
              help="Evolink results file", metavar="character"),
    make_option(c("-c", "--cutoff"), type="double", default=NULL, 
              help="threshold for Evolink index", metavar="float"),
    make_option(c("-a", "--top_pos"), type="integer", default=5, 
              help="Number of top up genes, default is 5", metavar="int"),
    make_option(c("-b", "--top_neg"), type="integer", default=5, 
              help="Number of top down genes, default is 5", metavar="int"),
    make_option(c("-m", "--display_mode"), type="integer", default=1, 
              help="Tree display mode/layout. [1: circular, 2: rectangular, default is 1]", metavar="int"),              
    make_option(c("-o", "--outdir"), type="character", default="Evolink_plot", 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree <- read.tree(opt$tree)
gene <- read.table(opt$gene, header=T, row.names=1, comment.char="", sep="\t")
trait <- read.table(opt$trait, header=T)
res <- read.table(opt$result, header=T, comment.char="", sep="\t")
cutoff = opt$cutoff
top_pos = opt$top_pos
top_neg = opt$top_neg
outdir = opt$outdir
display_mode = opt$display_mode

dir.create(outdir, showWarnings = FALSE)

# 1. ggtree plot
plot_tree <- function(p, gene, genes_list, top_gene_ct=5, color="#FE6D73", name="positive_genes", layout="circular"){
    dt = gene %>% filter(rownames(gene) %in% genes_list[1:min(length(genes_list), top_gene_ct)])
    dt_n = rownames(dt)
    dt <- sapply(dt, as.character)
    rownames(dt) <- dt_n
    dt <- t(dt) %>% as.data.frame()
    p = gheatmap(p, dt, offset = 0, color="#D2D9D4", width=0.8,
            colnames_position="top", 
            colnames_angle=90, colnames_offset_y=0,
            hjust=0, font.size=1.5) +
    scale_fill_manual(values=c("white", color), breaks=0:1, name=name)
    if(layout=="circular"){
        return(p + layout_circular())
    }else{
        return(p + coord_cartesian(clip = 'off') + 
        theme_tree(plot.margin=margin(30, 6, 10, 6)) +
        layout_rectangular())
    }
}

trait$Status = factor(trait$Status, levels=c(0, 1))
pos_genes = res %>% filter(significance=="sig", Evolink_index>0) %>% arrange(desc(Evolink_index)) %>% .[["orthoID"]]
neg_genes = res %>% filter(significance=="sig", Evolink_index<0) %>% arrange(Evolink_index) %>% .[["orthoID"]]

if(display_mode==1){
    layout="circular"
}else if(display_mode==2){
    layout="rectangular"
}else{
    print("Warnings: Display_mode is neither 1 or 2. Force display_mode=1")
    layout="circular"
}

ggtree_p <- ggtree(tree, size=0.15, open.angle=30) %<+% trait +
     geom_tippoint(mapping=aes(color=Status), 
                   size=2, show.legend=TRUE) +
     scale_color_manual(values=c("#ADA9B7", "#FFB844"), name="Phenotype")

if(length(pos_genes)>=1){
    suppressMessages({p1 = plot_tree(ggtree_p, gene, pos_genes, top_gene_ct=top_pos, color="#FE6D73", name="positive_genes", layout=layout)})
    tree_pos_path = file.path(outdir, "Tree_mapping_pos.pdf")
    ggsave(tree_pos_path, device="pdf", width=20, height=20, units="cm")
    print(paste("Generate", tree_pos_path, "done."))
}

if(length(neg_genes)>=1){
    suppressMessages({p2 = plot_tree(ggtree_p, gene, neg_genes, top_gene_ct=top_neg, color="#227C9D", name="negative_genes", layout=layout)})
    tree_neg_path = file.path(outdir, "Tree_mapping_neg.pdf")
    ggsave(tree_neg_path, device="pdf", width=20, height=20, units="cm")
    print(paste("Generate", tree_neg_path, "done."))
}

# 2. Evolink plot
bound = max(abs(res$Evolink_index))
plot_data = res %>% mutate(genes = case_when((Evolink_index > 0 & significance == "sig") ~ "up",
                        (Evolink_index < 0 & significance == "sig") ~ "down",
                        TRUE ~ "background"))
p = ggplot(plot_data, aes(x=Prevalence_index, y=Evolink_index, color=genes)) + 
    geom_hline(yintercept=0, color="black", size=0.2, linetype="dotted") +
    geom_vline(xintercept=0, color="black", size=0.2, linetype="dotted") +
    geom_point(alpha=0.9) +
    scale_color_manual(breaks = c("up","down","background"), values=c("#FE6D73", "#227C9D", "grey")) +
    ylim(c(-1*bound, bound)) +
    theme(text = element_text(size = 10)) +
    theme_classic()

evolink_plot_path = file.path(outdir, "Evolink.pdf")
ggsave(evolink_plot_path, device="pdf", width=15, height=10, units="cm")
print(paste("Generate", evolink_plot_path, "done."))

# 3. Manhattan plot
plot_data = res %>% mutate(genes=case_when((Evolink_index > 0 & significance == "sig") ~ "up",
                            (Evolink_index < 0 & significance == "sig") ~ "down",
                            TRUE ~ "background"))
plot_data$x = 1:length(plot_data$Evolink_index)

p = ggplot(plot_data, aes(x=x, y=abs(Evolink_index), color=genes)) + 
    geom_point(alpha=0.9) + 
    scale_color_manual(breaks=c("up", "down", "background"), values=c("#FE6D73", "#227C9D", "grey")) +
    geom_hline(yintercept=cutoff, color="black") + 
    xlab("genes") +
    ylab("absolute(Evolink_index)") + 
    theme_classic()

manhattan_plot_path = file.path(outdir, "Manhattan.pdf")
ggsave(manhattan_plot_path, device="pdf", width=20, height=10, units="cm")
print(paste("Generate", manhattan_plot_path, "done."))
