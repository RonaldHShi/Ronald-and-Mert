# Metadata retrieval test

# source libraries
library(RJSONIO)
library(RCurl)
library(matlab)


crop_metadata <- function(metadata,exclude,unique_rows=F){
  #input: metadata = metadata file, exclude = vector containing substrings that are to be 
  #                                           excluded from metadata
  #output: cropped metadata file named cropped_metadata
  metadata <- read.delim(metadata,header = F)
  metadata_orginal <<-metadata
  thrash <- c()
  changeddata <- metadata_orginal
  
  if(unique_rows){
  for (i in 1:(length(metadata[1,])-1)){
    if(length(unique(metadata[2:nrow(metadata),i]))==1){
      thrash <- c(thrash,i)
    }
  }
  changeddata <- metadata[-thrash]
  }
  for (x in exclude){
  new <- unlist(changeddata[1,])
  new <- as.character(new)
  thrash2 <- (grep(x,new))
  changeddata <<- changeddata[-thrash2]
  }
  #print(length(changeddata[1,]))
  write.table(changeddata,"cropped_metadata",sep="\t",row.names = F,col.names = F,quote = F)
  
}



get_mapping_options <- function(){
  #returns a vector of endpoints under mapping/expand
  system("curl 'https://gdc-api.nci.nih.gov/cases/_mapping' > mapping_JSON.txt")
  endpoints <- fromJSON("mapping_JSON.txt")["expand"]
  #print(endpoints)
  nodes <- c()
  for (i in 1:length(endpoints$expand)){
    a <- unlist(as.vector(strsplit(endpoints$expand[i],"\\.")))
    nodes <- c(nodes,a[1])
  }
  nodes <- unique(nodes)
  return(nodes)
}

endpoints <- get_mapping_options()

############################
### ### ### MAIN ### ### ###
############################
get_GDC_metadata <- function(id_list, my_rot="no", output="file",  order_rows=TRUE,  order_columns=TRUE, verbose=FALSE, debug=FALSE){

    # import list of ids
    my_ids <- flatten_list(as.list(scan(file=id_list, what="character")))

    metadata_matrix <- matrix()
    
    for (i in 1:length(my_ids)){
        
        #if(verbose==TRUE){print(paste("Processing sample (", i, ")"))}
        print(paste("Processing sample (", i, ":", my_ids[i], ")"))

        raw_metadata_vector <- vector()
        unique_metadata_vector <- vector
        
        if(debug==TRUE){print("Made it here (1)")}
        
        ###########CHANGED

        endpoints <- get_mapping_options() 
		# this line is under the library commands inorder to perform it once
        for (endpoint in endpoints){
			temp <- vector()
			for(j in c(1:50)){
				temp <- metadata_cases(my_id=my_ids[i], dict=endpoint)
				if(is.atomic(temp))
					Sys.sleep(10)
				else
					break
			}
            raw_metadata_vector <- c(raw_metadata_vector,
								 	flatten_list(temp$data))
        }

		###########END CHANGED

        if(debug==TRUE){print("Made it here (2)")}
        
        for (j in 1:length(unique(names(raw_metadata_vector)))){
            my_name <- unique(names(raw_metadata_vector))[j]
            #print(my_name)
            unique_name_value_pairs <- unique(  raw_metadata_vector[  which( names(raw_metadata_vector) == my_name )  ]  )
            name_list <- vector()
            #names(unique_name_value_pairs) <- rep(my_name, length(unique_name_value_pairs))
            for (k in 1:length(unique_name_value_pairs)){
                name_list <- c(name_list, (paste( my_name, ".", k ,sep="")))
            }
            names(unique_name_value_pairs) <- name_list
            unique_metadata_vector <- c(unique_metadata_vector, unique_name_value_pairs)
        }

        if(debug==TRUE){print("Made it here (3)")}
    
        if( i==1 ){ # on first sample, create the data matrix
            metadata_matrix <- matrix(unique_metadata_vector, ncol=1)
            rownames(metadata_matrix) <- names(unique_metadata_vector)
            colnames(metadata_matrix) <- convert_fileUUID_to_name(my_ids[i])
            if(debug==TRUE){print("Made it here (3.1)")}
        }else{ # for all additional samples add on to the existing matrix
            if(debug==TRUE){print("Made it here (3.2)")}
            sample_metadata_matrix <- matrix(unique_metadata_vector, ncol=1)
            if(debug==TRUE){print("Made it here (3.3)")}
            rownames(sample_metadata_matrix) <- names(unique_metadata_vector)
            if(debug==TRUE){print("Made it here (3.4)")}
            colnames(sample_metadata_matrix) <- convert_fileUUID_to_name(my_ids[i])
            if(debug==TRUE){
                print("Made it here (3.5)")
                matrix1 <<- metadata_matrix
                matrix2 <<- sample_metadata_matrix
            }
            if(debug==TRUE){print("Made it here (3.6)")}
            # place merge code here
			

            # Note - merging changes class of metadata_matrix from "matrix" to "data frame"; it's converted back below
            metadata_matrix <- combine_matrices_by_column(matrix1=metadata_matrix, matrix2=sample_metadata_matrix, func_order_rows=order_rows, func_order_columns=order_columns, func_debug=debug)
            if(debug==TRUE){
                print("Made it here (3.7)")
                print(paste("DATA_CLASS:", class(metadata_matrix)))
            }
        }

        
        
    }
    #ADDED
    #print(head(metadata_matrix))
    metadata_matrix <- metadata_matrix[-1,]
      # temporary edit to get rid of junk metadata that was screwing up the table
    #END ADDED

    # covert data product from data.frame back to matrix
        metadata_matrix <- as.matrix(metadata_matrix)
        
        #ADDED
        #sorts matrix by name
        #metadata_matrix <- metadata_matrix[, order(colnames(metadata_matrix))]
        
        if(debug==TRUE){print("Made it here (4)")}
        
        # rotate the matrix if that option is selected
        if( identical(my_rot, "yes")==TRUE ){
            metadata_matrix <- rot90(rot90(rot90(metadata_matrix)))
        }

        if(debug==TRUE){print("Made it here (5)")}
    
        # output the matrix as a flat file (default) or return as matrix
        # create name for the output file        
        if( identical(output, "file")==TRUE ){
            date_tag <- gsub(" ", "_", date())
            date_tag <- gsub("__", "_", date_tag)
            date_tag <- gsub(":", "-", date_tag)
            #output_name =paste(id_list, ".", date_tag, ".GDC_METADATA.txt", sep="")
            output_name = gsub(" ", "", paste(id_list,".GDC_METADATA.txt"))
            if(debug==TRUE){ print(paste("FILE_OUT: ", output_name)) }
            export_data(metadata_matrix, output_name)
        }else{
            return(metadata_matrix)
        }
        if(debug==TRUE){print("Made it here (6)")}
    
}
############################
############################
############################                                                                                           


############################
### ### ### SUBS ### ### ###
############################
metadata_cases <- function(before_id="https://gdc-api.nci.nih.gov/cases?filters=%7b%0d%0a+++%22op%22+%3a+%22%3d%22+%2c%0d%0a+++%22content%22+%3a+%7b%0d%0a+++++++%22field%22+%3a+%22files.file_id%22+%2c%0d%0a+++++++%22value%22+%3a+%5b+%22",
                           after_id="%22+%5d%0d%0a+++%7d%0d%0a%7d&pretty=true&fields=case_id&expand=",
                           my_id="07218202-2cd3-4db1-93e7-071879e36f27", 
                           dict="")
{ 
  my_call <- gsub(" ", "", paste(before_id, my_id, after_id, dict))
  my_call.json <- fromJSON(getURL(my_call))
  #print(my_call.json)
  return(my_call.json)
  #my_call.list <- flatten_list(my_call.json)
  #print(my_call.list)
}



export_data <- function(data_object, file_name){
    write.table(data_object, file=file_name, sep="\t", col.names = NA, row.names = TRUE, quote = FALSE, eol="\n")
}



flatten_list <- function(some_list){
    flat_list <- unlist(some_list)
    flat_list <- gsub("\r","",flat_list)
    flat_list <- gsub("\n","",flat_list)
    flat_list <- gsub("\t","",flat_list)
}



combine_matrices_by_column <- function(matrix1, matrix2, func_order_rows=TRUE, func_order_columns=TRUE, func_debug=debug){

    # perform the merge
    comb_matrix<- merge(data.frame(matrix1), data.frame(matrix2), by="row.names", all=TRUE)
    if(func_debug==TRUE){
        print("Made it here (3.6.1)")
        print(paste("MATRIX_1", dim(matrix1)))
        print(paste("MATRIX_2",dim(matrix2)))
        print(paste("MATRIX_C",dim(comb_matrix)))
        #print(colnames(matrix1))
        #print(colnames(matrix2))
        matrix3 <<- comb_matrix
    }
    
    # undo garbage formatting that merge introduces
    rownames(comb_matrix) <- comb_matrix$Row.names
    comb_matrix$Row.names <- NULL

    matrix4 <<- comb_matrix

    #if(func_debug==TRUE){print(paste("MATRIX DIM:", dim(comb_matrix)))}
    colnames(comb_matrix) <- c(colnames(matrix1), colnames(matrix2))
    #if(func_debug==TRUE){print("Made it here (3.6.2)")}
    
    # order columns
    if( func_order_rows==TRUE){
        ordered_rownames <- order(rownames(comb_matrix))
        comb_matrix <- comb_matrix[ordered_rownames,]
    }
    #if(func_debug==TRUE){print("Made it here (3.6.3)")}

    # order rows
    if( func_order_columns==TRUE){
        ordered_colnames <- order(colnames(comb_matrix))
        comb_matrix <- comb_matrix[,ordered_colnames]
    }
    #if(func_debug==TRUE){print("Made it here (3.6.4)")}

    #if(func_debug==TRUE){ matrix5 <<- comb_matrix }
    
    return(comb_matrix)
}



combine_matrices_by_row <- function(matrix1, matrix2, pseudo_fudge=10000, func_order_rows=TRUE, func_order_columns=TRUE, func_debug=debug){

    # perform the merge
    comb_matrix<- merge(matrix1, matrix2, by="col.names", all=TRUE)

    # undo garbage formatting that merge introduces
    rownames(comb_matrix) <- comb_matrix$Col.names
    comb_matrix$Col.names <- NULL
    colnames(comb_matrix) <- c(colnames(matrix1), colnames(matrix2))
    
    # order columns
    if( func_order_rows==TRUE){
        ordered_rownames <- order(rownames(comb_matrix))
        comb_matrix <- comb_matrix[ordered_rownames,]
    }

    # order rows
    if( func_order_columns==TRUE){
        ordered_colnames <- order(colnames(comb_matrix))
        comb_matrix <- comb_matrix[,ordered_colnames]
    }
    
    return(comb_matrix)
}

convert_fileUUID_to_name <- function(UUID)
{
  before_id <- "https://gdc-api.nci.nih.gov/files/"
  after_id <- "?&fields=file_name&pretty=true"
  my_call <- gsub(" ", "", paste(before_id, UUID, after_id))
  my_call.json <- fromJSON(getURL(my_call))
  return(substr(unname(my_call.json$data),1,36))
}

############################
############################
############################



############################
### ### ## NOTES ### ### ###
############################
## length(my_metadata)
## length(unique(my_metadata))
## length(names(my_metadata))
## length(unique(names(my_metadata)))
     
## length(unique_metadata)
## length(unique(unique_metadata))
## length(names(unique_metadata))
## length(unique(names(unique_metadata)))
############################
############################
############################
