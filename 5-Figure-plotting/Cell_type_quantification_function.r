## Cell type quantification function:
# factorize cell type before running the functions, so that the cell types are in order:

## Example
## plotlist <- quantify_cell_types_per_sample(metadata=metadata, 
##                                             ident='samples', 
##                                             group='group', 
##                                             group_levels=NA,  
##                                             cell_type_var='Broad_cell_identity', 
##                                             stats=TRUE, 
##                                             group_comparisons=NULL, 
##                                             p_output_format="p.signif", 
##                                             group_colors=NULL)


library(dplyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


quantify_cell_types_per_sample <- function(metadata=metadata, ident='samples', min_cell_ident=0, group='group', group_levels=NA,  cell_type_var='Broad_cell_identity', stats=TRUE, group_comparisons=NA, p_output_format="p.signif", group_colors=NA, save_percentages=FALSE, file_name=NA){

    grouping_meta <- metadata %>%
                        dplyr::select(!!sym(ident), !!sym(group)) %>%
                        remove_rownames() %>%
                        unique()

    count <- metadata %>% 
                group_by(!!sym(ident)) %>% 
                dplyr::count(!!sym(cell_type_var)) %>%
                mutate(percent=(n/sum(n)*100), number_cells=sum(n)) %>%
                filter(number_cells>=min_cell_ident) %>%
                select(-n)%>%
                select(-number_cells)

    #sum_plot <- reshape::cast(count, !!sym(ident)~!!sym(cell_type_var), value="percent")
    sum_plot <- count %>% tidyr::spread(!!sym(cell_type_var), percent, fill=0 )
    sum_plot <- merge(sum_plot, grouping_meta, by=eval(ident))

    if(save_percentages){
        write.csv(sum_plot, file_name, row.names=FALSE)
    }


    if(!is.na(group_levels)){
        sum_plot[,eval(group)] <- factor(sum_plot[,eval(group)], levels=group_levels)
    }

    ## # remove NAs from df and replace with 0:
    ## sum_plot[is.na(sum_plot)] <- 0
    # remove odd sign from cell type names:
    colnames(sum_plot) <- gsub("\\+", "pos", colnames(sum_plot))
    colnames(sum_plot) <- gsub("-", "_", colnames(sum_plot))
    colnames(sum_plot) <- gsub(" ", "_", colnames(sum_plot))

    cell_type_names <- colnames(sum_plot)[2:(length(colnames(sum_plot))-1)]
    plotlist <- list()

    if(is.na(group_colors)){
        group_colors = colorRampPalette(RColorBrewer::brewer.pal(11, "Paired"))(length(unique(sum_plot[,eval(group)])))
    }

    if(stats){
        for (c in cell_type_names){
            
            if(max(sum_plot[,eval(c)]) < 1 ){
                    label_pos <- max(sum_plot[,eval(c)])+0.5
                }else if(max(sum_plot[,eval(c)]) < 5 ){
                    label_pos <- max(sum_plot[,eval(c)])+2
                } else if(max(sum_plot[,eval(c)]) < 12  ){
                    label_pos <- max(sum_plot[,eval(c)])+5 
                } else if(max(sum_plot[,eval(c)]) < 30 ){
                    label_pos <- max(sum_plot[,eval(c)])+10
                } else if (max(sum_plot[,eval(c)]) < 60){
                    label_pos <- max(sum_plot[,eval(c)])+25
                } else {
                    label_pos <- max(sum_plot[,eval(c)])+50
            }

            if(length(unique(sum_plot[,group]))>2){
                my_comparisons <- group_comparisons
                plotlist[[eval(c)]] <- ggplot(sum_plot, aes_string(x=group, y=c, fill=group)) +
                                                    geom_boxplot(outlier.size=0.1, outlier.color="white") + 
                                                    scale_fill_manual(values=group_colors) +
                                                    ylab(paste0("% ",c," in sample")) +
                                                    expand_limits(y=0) +
                                                    theme_classic() +
                                                    ggtitle(eval(c))+
                                                    geom_jitter(width = 0.1, height = 0.01) +
                                                    theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))+
                                                    stat_compare_means(comparisons = my_comparisons, label = p_output_format) + 
                                                    stat_compare_means(label.y = label_pos,label ="p.format" )
            } else{
                plotlist[[eval(c)]] <- ggplot(sum_plot, aes_string(x=group, y=c, fill=group)) +
                                        geom_boxplot(outlier.size=0.1, outlier.color="white") + 
                                        scale_fill_manual(values=group_colors) +
                                        ylab(paste0("% ",c," in sample")) +
                                        expand_limits(y=0) +
                                        theme_classic() +
                                        ggtitle(eval(c))+
                                        geom_jitter(width = 0.1, height = 0.01) +
                                        theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))+
                                        stat_compare_means(label = p_output_format)
                
            }

        }

    } else {
        for (c in cell_type_names){
            plotlist[[eval(c)]] <- ggplot(sum_plot, aes_string(x=group, y=c, fill=group)) +
                                                geom_boxplot(outlier.size=0.1, outlier.color="white") + 
                                                scale_fill_manual(values=group_colors) +
                                                ylab(paste0("% ",c," in sample")) +
                                                expand_limits(y=0) +
                                                theme_classic() +
                                                ggtitle(eval(c))+
                                                geom_jitter(width = 0.1, height = 0.01) +
                                                theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))
                               
        }
    }

        return(plotlist)

}


