library(RJSONIO)
library(RCurl)

source("~/git/Ronald-and-Mert/preprocessing_tool.r")
source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r")
source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")
source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
source("~/git/Ronald-and-Mert/calc_stats.r")
source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")
library(DESeq)

CompareNormalCancerPipeline <- function(projectsVec = F, projectNamesFile = F) {
  # Entire pipeline for analyzing the normal vs cancer HTSeq counts data in a single or 
  # multiple projects from the GDC
  #
  # Args:
  #   projectsVec: A project name you wish to analyze alone or a vector of 
  #                many project names you wish analyze together
  #   projectNamesFile: A file or a list of files that contain project
  #                     names separated by line
  #
  # Returns:
  # (None)
  # Downloads all HTSeq Counts files
  # Generates PCoA colored by whether the sample is normal or tumor
  # Generates Heatmap from top 1% of significant genes

  # builds projectsVec from a newline separated projects id file
  if (projectNamesFile) {
    projectsVec <- vector()
    for(i in projects) {
      projectsVec <- c(projectsVec , scan(i, what = "character", sep = "\n"))
    }
  }
  # downloads normal samples for every project from GDC
  for(i in projectsVec) {
    download_normal_from_GDC(i)
  }
  # unzips files
  system("gunzip *.gz")

  # download matching cancer samples from GDC
  system("ls | grep .counts$ > normal_filenames")
  vec.normal.filenames <- scan("normal_filenames", what = "character", sep = "\n")
  for (f in vec.normal.filenames) {
    file.name <- strsplit(f, "\\.")[[1]][1]
    my_call <- paste0("https://gdc-api.nci.nih.gov/cases?fields=case_id&size=99999&pretty=true&filters=%7b%0d%0a+++%22op%22+%3a+%22%3d%22+%2c%0d%0a+++%22content%22+%3a+%7b%0d%0a+++++++%22field%22+%3a+%22files.file_id%22+%2c%0d%0a+++++++%22value%22+%3a+%5b+%22", f, "%22+%5d%0d%0a+++%7d%0d%0a%7d")
    my_call.json <- fromJSON(getURL(my_call))
    case.id <- unlist(my_call.json$data$hits)
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%0D%0A%09%22op%22%3A%22and%22%2C%0D%0A%09%22content%22%3A%5B%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22cases.samples.sample_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22Primary+Tumor%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22analysis.workflow_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22files.data_format%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22TXT%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22%3D%22%2C%0D%0A%09++++%22content%22%3A%7B%0D%0A%09++++%09%22field%22%3A%22cases.case_id%22%2C%0D%0A%09++++%09%22value%22%3A%5B%22", case.id, "%22%5D%0D%0A%09++++%7D%0D%0A%09%7D%5D%0D%0A%7D")
    my_call.json <- fromJSON(getURL(my_call))
    matching.cancer.id <- unlist(my_call.json$data$hits)
    system(paste0("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
          matching.cancer.id, "'"))
  }

  # unzip cancer files
  system("gunzip *.gz")
  system("ls | grep .counts$ > all_filenames")

  # merges abundance data together
  GDC_raw_count_merge(id_list="all_filenames")

  # gets UUIDs of data merged together
  export_listof_UUIDs(tsv = "all_file_UUIDs")

  # gets metadata file
  get_GDC_metadata("all_file_UUIDs.merged_data_file_UUIDs", my_rot = "yes")

  # normalizing the data
  preprocessing_tool(data_in = "all_filenames.merged_data.txt", produce_boxplots = TRUE)

  # add normOrCancer col
  metadata.table <- read.table(file="all_filenames.merged_data_file_UUIDs.GDC_METADATA.txt",row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE)
  normOrCancer <- c()
  for (r in row.names(metadata.table)) {
    if (r %in% vec.normal.filenames) {
      normOrCancer <- c(normOrCancer, "Normal")
    } else {
      normOrCancer <- c(normOrCancer, "Cancer")
    }
  }
  added.col.metadata.table <- cbind(metadata.table, normOrCancer)
  write.table(added.col.metadata.table, file = "all_filenames.merged_data_file_UUIDs.GDC_METADATA.txt.addedcol", 
      sep = "\t", quote = F, col.names = NA, row.names = T)

  # generate PCoA
  render_calculated_pcoa(PCoA_in= "all_filenames.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
      metadata_table = "all_filenames.merged_data_file_UUIDs.GDC_METADATA.txt.addedcol", 
      metadata_column_index = "normOrCancer")

  # calc stats
  calc_stats(data_table = "all_filenames.merged_data.txt.DESeq_blind.PREPROCESSED.txt", 
      metadata_table = "all_filenames.merged_data_file_UUIDs.GDC_METADATA.txt.addedcol", 
      metadata_column = "normOrCancer")

  # take top 1% of genes
  stats.table <- read.table(file = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Kruskal-Wallis.normal.or.tumor.col.STATS_RESULTS.txt",
      row.names = 1,header = TRUE,sep = "\t",
      colClasses = "character", check.names = FALSE,comment.char = "", 
      fill = TRUE, blank.lines.skip = FALSE)
  filtered.stats.table <- stats.table[1:floor(nrow(stats.table)/100), 1:(ncol(stats.table) - 7)]
  write.table(filtered.stats.table, file = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Kruskal-Wallis.normal.or.tumor.col.STATS_REMOVED.filtered.txt", 
      sep = "\t", quote = F, col.names = NA, row.names = T)

  # generate heatmap
  heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.Kruskal-Wallis.normal.or.tumor.col.STATS_REMOVED.filtered.txt", 
      metadata_table = "all_filenames.merged_data_file_UUIDs.GDC_METADATA.txt.addedcol", 
      metadata_column = "normOrCancer")
}

download_normal_from_GDC <- function(projects) {
  # download normal sample HTSeq Counts data for a given project
  # 
  # Args:
  #   project: vector of project names
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to be unzippped afterward
  for (p in projects) {
    print(p)
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22Blood+Derived+Normal%22%2C%22Solid+Tissue+Normal%22%2C%22Bone+Marrow+Normal%22%2C%22Fibroblasts+from+Bone+Marrow+Normal%22%2C%22Buccal+Cell+Normal%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22TXT%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22", p, "%22%5D%7D%7D%5D%7D")
    my_call.json <- fromJSON(getURL(my_call))
    UUID.list <- unlist(my_call.json$data$hits)
    for(j in UUID.list) {
      print(paste0(j, ": ", j))
      system(paste0("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                  j, "'"))
    }
  }
}

download_cancer_from_GDC <- function(projects) {
  # download cancer sample HTSeq Counts data for a given project
  # 
  # Args:
  #   project: vector of project names
  #
  # Returns:
  #   (None)
  #   Downloads all HTSeq Counts files as *.gz, need to be unzippped afterward
  for (p in projects) {
    my_call <- paste0("https://gdc-api.nci.nih.gov/files?fields=file_id&size=99999&pretty=true&filters=%7B%0D%0A%09%22op%22%3A%22and%22%2C%0D%0A%09%22content%22%3A%5B%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22cases.samples.sample_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22Primary+Tumor%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22analysis.workflow_type%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22HTSeq+-+Counts%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22in%22%2C%0D%0A%09%09%22content%22%3A%7B%0D%0A%09%09%09%22field%22%3A%22files.data_format%22%2C%0D%0A%09%09%09%22value%22%3A%5B%22TXT%22%5D%0D%0A%09%09%09%7D%0D%0A%09%09%7D%2C%7B%0D%0A%09%09%22op%22%3A%22%3D%22%2C%0D%0A%09++++%22content%22%3A%7B%0D%0A%09++++%09%22field%22%3A%22cases.project.project_id%22%2C%0D%0A%09++++%09%22value%22%3A%5B%22", p, "%22%5D%0D%0A%09++++%7D%0D%0A%09%7D%5D%0D%0A%7D")
    my_call.json <- fromJSON(getURL(my_call))
    UUID.list <- unlist(my_call.json$data$hits)
    for(j in UUID.list) {
      print(paste0(j, ": ", j))
      system(paste("curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/", 
                   j,
                   "'",
                   sep=""))
    }
  }
}
