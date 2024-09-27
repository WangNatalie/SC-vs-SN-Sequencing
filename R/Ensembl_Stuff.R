create_maps <- function() {
	#require("biomaRt")		
	ensg_name_map <- read.table("Hsap_ensg2symbol.txt", header=T);
	musg_name_map <- read.table("Mmus_ensg2symbol.txt");
	ensg2musg <- read.table("human-mouse_orthologs.txt");
	
	ensg2musg_all <- ensg2musg[,1:2]
	ensg2musg_one2one <- ensg2musg[ensg2musg[,4]=="ortholog_one2one",1:2]
	ensg_map_obj_list <- list(ensg_name_map=ensg_name_map, ensg2musg=ensg2musg_all, musg_name_map=musg_name_map, ensg2musg_one2one=ensg2musg_one2one);
	saveRDS(ensg_map_obj_list, "/cluster/projects/macparland/TA/ExternalData/Ensembl/ID_Maps.rds");
}

#ensg_map_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/ExternalData/ID_Maps_fixed.rds")
ensg_map_obj <- readRDS("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Genes_111_my_ensembl_name_maps.rds") # Changed 26-03-2024
ensg_name_map <- ensg_map_obj[["ensg2symbol"]]
ensg2musg <- ensg_map_obj[["ensg2musg"]]
musg_name_map <- ensg_map_obj[["musg2symbol"]]

#### Rectify duplicate ####
if (!file.exists("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Genes_111_deduped_name_maps.rds")) {
	ensg_name_map$has.ortho <- ensg_name_map[,1] %in% ensg2musg[,1]
	has.dups <- unique(ensg_name_map[duplicated(ensg_name_map[,2]),2])
	for (g in has.dups) {
		ensg_name_map_trimmed <- ensg_name_map[-1*which(ensg_name_map[,2] == g),]
		if (sum(ensg_name_map[ensg_name_map[,2]==g,"has.ortho"]) == 0) {
			keep <- ensg_name_map[ensg_name_map[,2]==g,][1,]
			ensg_name_map_trimmed <- rbind(ensg_name_map_trimmed, keep)
		} else {
			keep <- ensg_name_map[ensg_name_map[,2]==g & ensg_name_map$has.ortho==TRUE,][1,]
			ensg_name_map_trimmed <- rbind(ensg_name_map_trimmed, keep)
		}
		ensg_name_map <- ensg_name_map_trimmed
	}

	musg_name_map$has.ortho <- musg_name_map[,1] %in% ensg2musg[,1]
	has.dups <- unique(musg_name_map[duplicated(musg_name_map[,2]),2])
	for (g in has.dups) {
		musg_name_map_trimmed <- musg_name_map[-1*which(musg_name_map[,2] == g),]
		if (sum(musg_name_map[musg_name_map[,2]==g,"has.ortho"]) == 0) {
			keep <- musg_name_map[musg_name_map[,2]==g,][1,]
			musg_name_map_trimmed <- rbind(musg_name_map_trimmed, keep)
		} else {
			keep <- musg_name_map[musg_name_map[,2]==g & musg_name_map$has.ortho==TRUE,][1,]
			musg_name_map_trimmed <- rbind(musg_name_map_trimmed, keep)
		}
		musg_name_map <- musg_name_map_trimmed
	}
	saveRDS(list(ensg_name_map=ensg_name_map, musg_name_map=musg_name_map), file="/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Genes_111_deduped_name_maps.rds")
} else {
	tmp <- readRDS("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Genes_111_deduped_name_maps.rds")
	ensg_name_map <- tmp$ensg_name_map
	musg_name_map <- tmp$musg_name_map
	rm(tmp)
}

h_dups = unique(ensg2musg[(duplicated(ensg2musg[,1])),1])
m_dups = unique(ensg2musg[(duplicated(ensg2musg[,2])),2])
ensg2musg_one2one <- ensg2musg[!(ensg2musg[,1] %in% h_dups) & !(ensg2musg[,2] %in% m_dups),]

map_symbol_ensg <- function(genes, is.org=c("Hsap","Mmus"), is.name=c("symbol","ensg")) {
	if (is.org[1] == "Hsap") {
		if(is.name[1]=="symbol") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Mmus") {
		if(is.name[1]=="symbol") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

map_Hsap_Mmus <- function(genes, is.org=c("Hsap","Mmus"), one2one=FALSE) {
	if (one2one) {
		local_map <- ensg2musg_one2one
	} else {
		local_map <- ensg2musg
	}
	if (is.org[1] == "Hsap") {
		new = as.character(local_map[match(genes, local_map[,1]),2])
	} else if (is.org[1] == "Mmus") {
		new = as.character(local_map[match(genes, local_map[,2]),1])
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

General_Map <- function(genes, in.org=c("Hsap","Mmus"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus"), out.name=c("symbol","ensg")) {
	if (in.org == out.org & in.name == out.name) {
		# No change
		return(genes)
	}
	if (in.org == out.org) {
		# change names not spp
		return(map_symbol_ensg(genes, is.org=in.org, is.name=in.name))

	} else {
		if (in.name == "symbol") {
			tmp <- map_symbol_ensg(genes, is.org=in.org, is.name=in.name)
		} else {
			tmp <- genes
		}
		tmp <- map_Hsap_Mmus(tmp, is.org=in.org)
		if (out.name =="symbol") {
			out <- map_symbol_ensg(tmp, is.org=out.org, is.name="ensg")
		} else {
			out <- tmp
		}
		return(out)	
	}
}

M_symbol_2_H_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Mmus", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Mmus")
        tmp <- map_symbol_ensg(tmp, is.org="Hsap", is.name="ensg")
        return(tmp)
}

H_symbol_2_M_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Hsap", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Hsap")
        tmp <- map_symbol_ensg(tmp, is.org="Mmus", is.name="ensg")
        return(tmp)
}

