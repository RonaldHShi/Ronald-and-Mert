library(RJSONIO)
library(RCurl)

whole_process <- function(projects, # use files with project_ids or a vector
                          # of project_ids
                          # parameters to skip certain parts of script
                          skip_to_merging = F,
                          skip_to_metadata = F, 
                          skip_rendering = F, 
                          skip_to_corr_matrix = F, 
                          # size_lim for every project: only downloads XX 
                          # number of files per project if XX files are
                          # available
                          size_lim = F, 
                          # option to render heatmap dendrogram
                          heatmap_option = F) {	

  if (skip_to_metadata == T) {
    skip_to_merging = T
  }

  # for logging purposes
  time <- Sys.time()
  time <- gsub(" ", "", time, fixed = TRUE)

  for(i in seq_len(sink.number())){
    sink(NULL)
  }
  log <- file(paste(paste(projects,collapse="."),".full.log",sep=""), open="wt")
  sink(log)

  to_log("START")

  # from Kevin_R_scripts
  source("~/git/Ronald-and-Mert/install_r_prereqs.r")       
  source("~/git/Ronald-and-Mert/preprocessing_tool.r")
  source("~/git/Ronald-and-Mert/heatmap_dendrogram.r")
  source("~/git/Ronald-and-Mert/calculate_pco.r")
  source("~/git/Ronald-and-Mert/render_calculated_pcoa.r")
  # from Ronald-and-Mert
  source("~/git/Ronald-and-Mert/GDC_raw_count_merge.vRonald.r")
  source("~/git/Ronald-and-Mert/get_listof_UUIDs.r")
  source("~/git/Ronald-and-Mert/GDC_metadata_download.RandM.r")
  source("~/git/Ronald-and-Mert/correlate_pcoa.r")
  source("~/git/Ronald-and-Mert/improve_cor_matrix.r")
  library(DESeq)

  if (!skip_to_merging) {
    to_log("Downloading data")
    # if user inputed a file name with project_ids
    if (any(file.exists(projects))) {
      vec_of_projects <- vector()
      # loop through every filename and read in projects to vec_of_projects
      for(i in projects) {
        vec_of_projects <- c(vec_of_projects, scan(i, what = "character"))
      }
      # loop through vec_of_projects and download all HTSeq counts files
      for(i in vec_of_projects) {
        if (size_lim) {
          download_project_data(i, size_lim = size_lim)
        } else {
          download_project_data(i)
        }
      }
    } else {
      # else if user inputed a vector of project_ids
      for(i in projects) {
        if (size_lim) {
          download_project_data(i, size_lim = size_lim)
        } else {
          download_project_data(i)
        }
      }
    }

    to_log("Download Completed")

    # unzips data into .counts files
    to_log("Unzipping")
    system("gunzip *.gz")
    to_log("Unzipping completed")
  }
  if (!skip_to_metadata) {
    # spits all of .counts file names in counts_files
    system("ls | grep .counts$ > counts_files")

    # merges abundance data together
    to_log("Merging data")
    GDC_raw_count_merge(id_list="counts_files")
    to_log("Merge completed")

    # gets UUIDs of data merged together
    to_log("Exporting UUID's")
    export_listof_UUIDs(tsv = "counts_files.merged_data.txt")
    to_log("Export Completed")
  }
  # gets metadata file
  to_log("Getting Metadata")
  get_GDC_metadata("counts_files.merged_data_file_UUIDs", my_rot = "yes")

  # gets metadata_filename 
  # (obsolete now but useful for interacting with old metadata filenames)
  files <- list.files()[grep("GDC_METADATA.txt", list.files())]
  details <- file.info(list.files()[grep("GDC_METADATA.txt", list.files())])
  details <- details[with(details, order(as.POSIXct(mtime))), ]
  metadata_filename <- rownames(details)[length(files)] 
  print(metadata_filename)
  to_log("Metadata Download Completed")

  # normalizing the data
  to_log("Preprocessing")
  preprocessing_tool(data_in = "counts_files.merged_data.txt", produce_boxplots = TRUE)
  to_log("Preprocessing Completed")
  # visualize the data
  # generate heatmap dendrogram
  if (heatmap_option == T) {
    heatmap_dendrogram(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  }

  to_log("Calculating PCoA")
  calculate_pco(file_in = "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt")
  print(metadata_filename)
  to_log("Calculated PCoA")
  if (skip_rendering == F) {
    to_log("Rendering PCoA")
    # renders a PCoA for every piece of metadata
    render_calcualted_pcoa("counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA", 
                           metadata_table = metadata_filename, 
                           use_all_metadata_columns = T)
    to_log("Rendering Completed")
  }	

  to_log("Creating Correlation Matrix")
  correlate_pcoa("counts_files.merged_data_file_UUIDs.GDC_METADATA.txt",
                 "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")
  improve_corr_matrix(corr_matrix="cor_summary_spearman.txt", metadata_matrix="counts_files.merged_data_file_UUIDs.GDC_METADATA.txt")
  improve_corr_matrix(corr_matrix="cor_summary_pearson.txt", metadata_matrix="counts_files.merged_data_file_UUIDs.GDC_METADATA.txt")
  improve_corr_matrix(corr_matrix="cor_summary_pearson_test.txt", metadata_matrix="counts_files.merged_data_file_UUIDs.GDC_METADATA.txt")
  improve_corr_matrix(corr_matrix="cor_summary_spearman_test.txt", metadata_matrix="counts_files.merged_data_file_UUIDs.GDC_METADATA.txt")
  to_log("Created Correlation Matrix")
  to_log("END")
  sink()

}

# logging function: shows time and memory usage
to_log <- function(message){
  print(message)
  print(Sys.time())
  print(gc(verbose = T,reset = F))
}

individual_analysis <- function() {
  # perform the whole_process on every project
  #
  # Args:
  #   None
  # Returns:
  #   Produces a folder named the project name in the current directory

  # API call to get list of projects
  project_list <- fromJSON(getURL("https://gdc-api.nci.nih.gov/projects?fields=project_id&from=1&size=50&sort=project.project_id:asc&pretty=true"))$data$hits
  for(p in project_list) {
    print(p)
    # check if folder exists
    if (!dir.exists(file.path("/mnt/single_projects", p))) {
      # create folder if it doesn't exist
      dir.create(file.path("/mnt/single_projects", p))
      # go into folder and run whole_process()
      setwd(file.path("/mnt/single_projects", p))
      whole_process(p, skip_rendering = T)
      # return back to /mnt
      setwd("/mnt/single_projects")
    } else if (file.exists(file.path("/mnt/single_projects", p, 
                                     "counts_files.merged_data_file_UUIDs")) & 
               !file.exists(file.path("/mnt/single_projects", p, 
                                      "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA"))) {
      # check if it has already downloaded all the data
      setwd(file.path("/mnt/single_projects", p))
      whole_process(p, skip_rendering = T, skip_to_metadata = T)
    }

    # obsolete now (was used to generate correlation matrices for the projects
    # that didn't have them previously when running earlier version of code
    temp <- list.files()
    metadataFile <- temp[grepl("GDC_METADATA.txt$", temp)]
    if (!is.null(metadataFile) & 
        file.exists(file.path("/mnt/single_projects", p, 
                              "counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")) &
        !file.exists(file.path("/mnt/single_projects", p, "cor_summary.txt"))) {
      for(i in metadataFile) {
        source("~/git/Ronald-and-Mert/correlate_pcoa.r")
        setwd(file.path("/mnt/single_projects", p))
        correlate_pcoa(i,"counts_files.merged_data.txt.DESeq_blind.PREPROCESSED.txt.euclidean.PCoA")
      }
    }
  }
}
