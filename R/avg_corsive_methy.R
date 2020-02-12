average_corsiv_methylation <- function(DF){

    probe_loc <- read.csv("./CoRSIV_ESS_SIV_CG_sites_clusters_hg38.csv",header = T,stringsAsFactors = F)
    probe_loc <- probe_loc[order(probe_loc$chr, probe_loc$pos),]
    probe_loc <- probe_loc[!duplicated(probe_loc$CG),]

    temp_clinical <- DF[ , !grepl( "cg[0-9]" , names( DF ) ) ]
    probe_data <-DF[ , grepl( "cg[0-9]" , names( DF ) )]

    avg_prob_df <- data.frame(matrix(ncol = 0, nrow = dim(probe_data)[1]))

    for(probe_cluster_id in unique(probe_loc$cluster_id)){
        #print(probe_cluster_id)
        #probe_cluster_id <- "single14"
        probs_temp <- probe_loc[probe_loc$cluster_id==probe_cluster_id,]$CG
        if(length(intersect(colnames(probe_data),probs_temp)) > 0){
            data_prob_cluster <- probe_data[intersect(colnames(probe_data),probs_temp)]
            avg_prob_df[[probs_temp[1]]] <- apply(data_prob_cluster, 1,mean)
        }
    }
    return(cbind(temp_clinical,avg_prob_df))
}
#asdfdas
