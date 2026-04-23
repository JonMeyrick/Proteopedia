########### Import Libraries ####################################################################################################################################
#' @import BiocGenerics
#' @import clusterProfiler
#' @import data.table
#' @import DescTools
#' @import DOSE
#' @import dplyr
#' @import enrichplot
#' @import forcats
#' @import ggplot2
#' @import ggrepel
#' @import ggside
#' @import ggupset
#' @import grid
#' @import highcharter
#' @import htmlwidgets
#' @import imputeLCMD
#' @import Proteopedia
#' @import limma
#' @import mzR
#' @import org.Hs.eg.db
#' @import patchwork
#' @import Peptides
#' @import RColorBrewer
#' @import ReactomePA
#' @import stats
#' @import stringr
#' @import tidyr

########### Package/Basic Function Loading ####################################################################################################################################
.onLoad <- function(libname = "~/Desktop/Positron_Scripts/R_Packages", pkgname = "Proteopedia") {
  ggplot2::theme_set(ggplot2::theme_classic(base_size = 20) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 26), 
                                                                             strip.text.x = ggplot2::element_text(size = 26),
                                                                             strip.text.y = ggplot2::element_text(size = 26),
                                                                             strip.background.x = ggplot2::element_blank(),
                                                                             strip.background.y = ggplot2::element_blank()))
  set.seed(123)
}
.onAttach <- function(libname = "~/Desktop/Positron_Scripts/R_Packages", pkgname = "Proteopedia") {
  packageStartupMessage("Proteopedia Loaded: Welcome Back...")
}
#' @export
"%!in%" <- function(x,y)!("%in%"(x,y))
#' @export
"%like%" <- function(x,y){like(x,y)}
#' @export
End_Timer <- function(Start){
  time.taken <- format(difftime(Sys.time(), Start, unit= "secs"), scientific = F)
  Minutes <- floor(as.numeric(gsub("(.*)\\..*", "\\1", time.taken))/60)
  Seconds <- as.numeric(gsub("(.*)\\..*", "\\1", time.taken)) - floor(as.numeric(gsub("(.*)\\..*", "\\1", time.taken))/60)*60
  Milliseconds <- gsub("(..).*", "\\1", as.numeric(gsub(" secs", "", gsub(".*\\.(.*)", "\\1", time.taken))))
  message(paste0("Time Elapsed: ", Minutes, " min ", Seconds, " s ", Milliseconds, " ms"))
}
#' @export
Update_Proteopedia <- function(){
  library(devtools)
  library(roxygen2)
  setwd("/Users/JM/Desktop/Positron_Scripts/R_Packages/Proteopedia")
  devtools::document()
  devtools::build()
  devtools::install()
}
#' @export
Check_AggregateRows <- function(InputData, IDColumns){
  InputData[, .N, IDColumns][N > 1]
}
#' @export
Reset_Dev <- function(){
  done <- FALSE
  while(!done){
    result <- try(dev.off(),silent = TRUE)
    done <- class(result) == "try-error"
  }
}

########### ZenoTOF Setup Functions ####################################################################################################################################
#' @export
Generate_Metadata <- function(Cell = 0, Drug = 0, Conc = 0, Time = 0, Reps = 3, OrderedBy = c("Cell", "Drug", "Conc", "Time", "Rep"), Output_Directory = "~/Desktop"){
  
  Sample_Metadata <- tidyr::crossing(Cell, Drug, Conc, Time, 1:Reps) |> data.table::data.table() |>  data.table::setnames("1:Reps", "Rep") |> 
                        data.table::setorderv(OrderedBy) |> tibble::rowid_to_column("ID")
  suppressWarnings(
    Sample_Metadata <- Sample_Metadata[, Sample := paste0(Cell,"_",Drug,"_",Conc,"_",Time,"h_R",Rep)]
  )
  data.table::fwrite(Sample_Metadata, paste0(Output_Directory,"/Sample_Metadata.csv"))
}
#' @export
Generate_ZenoTOFBatch <- function(Input_Metadata, Run_Date = paste0("E",substr(gsub("-","",lubridate::today()), 3, 8)), 
                                  Sample1Position, SampleNPosition = 0, LC_Run = "30method_Berlin_Nov_2024_CCP_depth_0.4mm", 
                                  Rack_Position = 1, InjectionVolume = 2, Blank_Position = "H12", Blank_Every = 6, InjectorID = "JM"){
  Sample_Metadata <- data.table::fread(Input_Metadata)
  N_Samples <- nrow(Sample_Metadata)
  
  # Define Sample Positions in Well
  PlateMap <- tidyr::crossing(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) |> data.table::data.table() |> 
    data.table::setnames(c("RowID", "ColID")) |> tibble::rowid_to_column("SampleN")
  suppressWarnings(
  PlateMap[, Well := paste(RowID, ColID)]
  )
  if(SampleNPosition == 0){
    Sample_Positions <- PlateMap[SampleN <= N_Samples, Well]
  } else {
    Sample_Positions <- PlateMap[ColID <= as.numeric(gsub(".(\\d+)", "\\1", SampleNPosition)), Well][1:N_Samples]
  }
  
  N_Injections <- 1 + N_Samples + (N_Samples/Blank_Every)
  Blank_Positions <- seq(from = 1, to = N_Injections, by = Blank_Every+1)
  
  Batch_File <- data.table::data.table("Sample_Name" = rep(0, times = N_Injections), "Sample_ID" = seq(1, to = N_Injections), "Barcode_ID" = "",
                           "MS_Method" = 0, "Processing_Method" = "", "LC_Method" = 0, "Rack_Type" = "Sample Manager", "Rack_Position" = Rack_Position,
                           "Plate_Type" = "Custom-96-Position", "Vial_Position" = 0, "Sample_Type" = "Unknown", "Dilution_Factor" = 1, "Weight_Volume" = 0,
                           "Data_File" = 0, "Results_File" = "", "Comment" = "", "Injection_Volume" = 0,
                           "Marker_Well" = FALSE)
  
  Batch_File[, `:=`(MS_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_12min", "uFlow_ZenoSWATH_85VW_11ms_30T_ZENO_ON_Berlin"),
                    LC_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_5grad_12method_depth_0.5mm", LC_Run),
                    Injection_Volume = data.table::fifelse(Sample_ID %in% Blank_Positions, 5, InjectionVolume))]
  # Add Blank ID
  count <- 1
  for(i in Blank_Positions){
    Batch_File$Sample_Name[i] <- paste0(Run_Date,"_",Batch_File$Sample_ID[i],"_Zeno_",InjectorID,"_QC_12_Blank",count)
    Batch_File$Vial_Position[i] <- Blank_Position
    count <- count + 1
  }

  # Add Sample ID  & Vial Positions
  Sample_IDs <- 1:N_Injections
  Sample_IDs <- Sample_IDs[Sample_IDs %!in% Blank_Positions]
  Sample_Metadata$Run <- 0
  count <- 1
  for(i in Sample_IDs){
    Batch_File$Sample_Name[i] <- paste0(Run_Date,"_",Batch_File$Sample_ID[i],"_Zeno_",InjectorID,"_IN_30_",Sample_Metadata$Sample[count])
    Sample_Metadata$Run[count] <- paste0(Run_Date,"_",Batch_File$Sample_ID[i],"_Zeno_",InjectorID,"_IN_30_",Sample_Metadata$Sample[count])
    Batch_File$Vial_Position[i] <- Sample_Positions[count]
    count <- count + 1
  }

  setwd(gsub("/Sample_Metadata.csv","", paste0(Input_Metadata)))
  data.table::fwrite(Sample_Metadata, paste0(Input_Metadata))
  
  Batch_File <- Batch_File[, Data_File := paste0("\\", substr(Run_Date, 1,3), "\\", substr(Run_Date, 1,5), "\\", Run_Date, "\\", Sample_Name)]
  suppressWarnings(
    Batch_File[, Sample_ID := ""]
  )
  colnames(Batch_File) <- gsub("Weight_Volume", "Weight/Volume", colnames(Batch_File))
  colnames(Batch_File) <- gsub("_", " ", colnames(Batch_File))
  data.table::fwrite(Batch_File, paste0(Run_Date,"_JMBatch.csv"))
}
#' @export
Generate_QCScaleBatch <- function(InjectionMasses = c(50, 75, 100, 150, 200), 
                                  Run_Date = paste0("E",substr(gsub("-","",lubridate::today()), 3, 8)), 
                                  Conc25nguLPosition, Conc50nguLPosition, LC_Run = "30method_Berlin_Nov_2024_CCP_depth_0.4mm", 
                                  Rack_Position = 1, Blank_Position = "H12", InjectorID = "JM"){
  
  N_Samples = length(InjectionMasses)
  N_Injections <- 1 + N_Samples + ceiling(N_Samples/Blank_Every)
  Blank_Positions <- seq(from = 1, to = N_Injections, by = 6)
  
  Batch_File <- data.table::data.table("Sample_Name" = rep(0, times = N_Injections), "Sample_ID" = seq(1, to = N_Injections), "Barcode_ID" = "",
                                       "MS_Method" = 0, "Processing_Method" = "", "LC_Method" = 0, "Rack_Type" = "Sample Manager", "Rack_Position" = Rack_Position,
                                       "Plate_Type" = "Custom-96-Position", "Vial_Position" = 0, "Sample_Type" = "Unknown", "Dilution_Factor" = 1, "Weight_Volume" = 0,
                                       "Data_File" = 0, "Results_File" = "", "Comment" = "", "Injection_Volume" = 0, "Marker_Well" = FALSE)
  
  Batch_File[, `:=`(MS_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_12min", "uFlow_ZenoSWATH_85VW_11ms_30T_ZENO_ON_Berlin"),
                    LC_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_5grad_12method_depth_0.5mm", LC_Run),
                    Injection_Volume = data.table::fifelse(Sample_ID %in% Blank_Positions, 5, 0))]
  # Add Blank ID
  count <- 1
  for(i in Blank_Positions){
    Batch_File$Sample_Name[i] <- paste0(Run_Date,"_",Batch_File$Sample_ID[i],"_Zeno_",InjectorID,"_QC_12_Blank",count)
    Batch_File$Vial_Position[i] <- Blank_Position
    count <- count + 1
  }
  
  # Add Sample ID  & Vial Positions
  Sample_IDs <- 1:N_Injections
  Sample_IDs <- Sample_IDs[Sample_IDs %!in% Blank_Positions]
  count <- 1
  for(i in Sample_IDs){
    Batch_File$Sample_Name[i] <- paste0(Run_Date,"_",Batch_File$Sample_ID[i],"_Zeno_",InjectorID,"_IN_30_K562_", InjectionMasses[count], "ng")
    count <- count + 1
  }
  Batch_File[!grepl("Blank", Sample_Name), Vial_Position := data.table::fifelse(as.numeric(gsub(".*K562_(\\d+)ng", "\\1", Sample_Name)) < 100, 
                                                    Conc25nguLPosition, Conc50nguLPosition)]
  Batch_File[!grepl("Blank", Sample_Name), Injection_Volume := data.table::fifelse(as.numeric(gsub(".*K562_(\\d+)ng", "\\1", Sample_Name)) < 100, 
                                                       as.numeric(gsub(".*K562_(\\d+)ng", "\\1", Sample_Name))/25, 
                                                       as.numeric(gsub(".*K562_(\\d+)ng", "\\1", Sample_Name))/50)]
  
  Batch_File <- Batch_File[, Data_File := paste0("\\", substr(Run_Date, 1,3), "\\", substr(Run_Date, 1,5), "\\", Run_Date, "\\", Sample_Name)]
  suppressWarnings(
    Batch_File[, Sample_ID := ""]
  )
  colnames(Batch_File) <- gsub("Weight_Volume", "Weight/Volume", colnames(Batch_File))
  colnames(Batch_File) <- gsub("_", " ", colnames(Batch_File))
  data.table::fwrite(Batch_File, paste0(Run_Date,"_QCScaleBatch.csv"))
}

########### MS Analysis: Data QC Functions ####################################################################################################################################
#' @export
Calculate_AUC <- function(InputDirectory){
  setwd(InputDirectory)
  mzML_Files = list.files(pattern = "mzML")
  AUCs = data.table::data.table()
  for(i in mzML_Files){tic_data <- mzR::tic(openMSfile(paste0(InputDirectory, "/",i)))
    AUCs = rbind(AUCs, data.table(mzML = i, AUC = DescTools::AUC(tic_data$rtime, tic_data$intensity)))}
  fwrite(AUCs, "AUC_Data.csv")
}
#' @export
Simplify_Data <- function(x){data.table::data.table(simplify(x, cutoff = 0.7))}
#' @export
Calculate_LFQ <- function(InputData, LFQ_colname, SILAC = F){    
  if(SILAC == T){
    tmp <- iq::fast_MaxLFQ(InputData[, .(protein_list = ProteinGroup, sample_list = Sample, id = Precursor.Id, quant = log2(Precursor.Quantity))]) 
  } else {
    tmp <- iq::fast_MaxLFQ(InputData[, .(protein_list = ProteinGroup, sample_list = Sample, id = Precursor.Id, quant = log2(Precursor.Normalised))]) 
  }
  tmp <- data.table::data.table(tmp$estimate, Precursor_group = tmp$annotation, keep.rownames = "ProteinGroup") 
  tmp <- data.table::melt.data.table(tmp, id.vars = c("ProteinGroup", "Precursor_group"), variable.name = "Sample", value.name = LFQ_colname) 
  tmp <- tmp[ Precursor_group == "" & !is.na( get(LFQ_colname) ) ][, -"Precursor_group" ] 
  tmp <- tmp[, (LFQ_colname) := 2^get(LFQ_colname)] 
}
#' @export
Merge_PrecursorData <- function(x, y){merge(x, y, by = c("ProteinGroup", "Sample"), all = T)}
#' @export
Count_Proteins <- function(dt, var_name){ 
  tmp <- dt[, .(N_samples = .N), ProteinGroup][, .(N_proteins = .N), N_samples]
  tmp[ order(-N_samples), cumulative_protein_N := cumsum(N_proteins)]
  data.table::data.table(tmp, Precursors = var_name)
}
#' @export
Calculate_CV <- function(x){sd(x)/mean(x)*100} 
#' @export
Calculate_PearsonsSkew <- function(x){3*((mean(x, na.rm = T) - median(x, na.rm = T))/sd(x, na.rm = T))}
#' @export
Calculate_RelativeTemp <- function(CurrentTemp, Temp = c(37, 42, 47, 52, 57), Shift = 1) {
  NextPos = which(Temp == unique(CurrentTemp)) + Shift
  if(data.table::between(NextPos, 1, length(Temp))){return(Temp[NextPos])}else{return(NA_character_)}
}
#' @export
Calculate_DensityPeak <- function(x){max(density(x, na.rm = T)$y)}
#' @export
Calculate_VolcanoLog2FC <- function(x){log2(abs(x))*sign(x)}

########### MS Analysis: Plot Aesthetics ####################################################################################################################################
#' @export
Test_ColourPalette <- function(ColourPalette){
  data.table::data.table(Colour = ColourPalette)[, ID := .I] |> ggplot2::ggplot(ggplot2::aes(x = ID, y = 1, fill = Colour)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_identity() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                                    axis.line = ggplot2::element_blank())
}

#' @export
Add_AbundanceAxes <- function(){ggplot2::labs(x = expression("Log"[2]~"FC in Protein Abundance"), y = expression("-Log"[10]~"P-Value"))}
#' @export
Add_KlossAxes <- function(scale = "Difference"){if(scale == "Log2FC"){ggplot2::labs(x = expression("Log"[2]~"FC in Protein Turnover Rate (k"[loss]~")"), y = expression("-Log"[10]~"P-Value"))} else 
{ggplot2::labs(x = expression("Difference in Protein Turnover Rate (k"[loss]~")"), y = expression("-Log"[10]~"P-Value"))}}
#' @export
Add_KsynAxes <- function(scale = "Difference"){if(scale == "Log2FC"){ggplot2::labs(x = expression("Log"[2]~"FC in Protein Synthesis Rate (k"[syn]~")"), y = expression("-Log"[10]~"P-Value"))} else 
{ggplot2::labs(x = expression("Difference in Protein Synthesis Rate (k"[syn]~")"), y = expression("-Log"[10]~"P-Value"))}}
#' @export
Add_IsotopeRatioAxes <- function(){ggplot2::labs(x = expression("Log"[2]~"FC Heavy:Light Protein LFQ Ratio"), y = expression("-Log"[10]~"P-Value"))}
#' @export
Add_NotSigBox <- function(){ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1)}
#' @export
Add_Rsq <- function(Subgroups = T){if(Subgroups == F){ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(rr.label),  group = 1), geom = "text")
                                    } else {ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(rr.label)), geom = "text")}}
#' @export
Add_NotSigBox <- function(){ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1)}
#' @export
Add_XYLine <- function(colour = "black"){ggplot2::geom_abline(linetype = "dashed", colour = colour)}
#' @export
Clean_SideDensities <- function(x = T, y = T){ggplot2::theme(ggside.axis.line = ggplot2::element_blank(), ggside.axis.ticks = ggplot2::element_blank(), 
                                                             ggside.axis.text = ggplot2::element_blank())}
#' @export
Add_GSEAAxes <- function(){ggplot2::labs(x = "Normalised Enrichment Score (NES)", y = expression("-Log"[10]~"Adj. P-Value"))}
#' @export
Add_Isotope_Colour <- function(){ggplot2::scale_colour_manual(values = c("Heavy" = "#90F", "Light" = "#F09"))}
#' @export
Add_Isotope_Fill <- function(){ggplot2::scale_fill_manual(values = c("Heavy" = "#90F", "Light" = "#F09"))}

########### MS Analysis: LF DIA-NN To Protein Data Functions ####################################################################################################################################
#' @export
Export_SampleRenameFile <- function(InputDirectory){
  setwd(InputDirectory)
  if(length(list.files(pattern = "report.tsv")) > 0){
    Sample_Rename <- data.table::fread(list.files(pattern = "report.tsv")[1])[, Run] |> unique() |> data.table(keep.rownames = T) |> data.table::setnames("Run")
    Sample_Rename[, Position := as.numeric(gsub("E\\d+_(\\d+)_.*", "\\1", Run))]
    Sample_Rename <- Sample_Rename |> data.table::setorder(Position)
    Sample_Rename[, `:=`(Renamed = NA, Position = NULL)]
    data.table::fwrite(Sample_Rename, "Sample_Rename.csv")
  } else {message("ERROR: No DIA-NN Report File Found")}
}
#' @export
Process_LF_DIANN <- function(InputDirectory, CtlGroup, ProteotypicFiltering = F){
  start.time <- Sys.time()
  set.seed(123)
  message("Importing DIA-NN Report File")
  {
    setwd(InputDirectory)
    if(length(list.files(pattern = "report.tsv")) > 0){
      InputFile <- list.files(pattern = "report.tsv")[1]
      PrecursorData <- data.table::fread(InputFile)[, .(Run, Protein.Group, Protein.Ids, First.Protein.Description, Genes, 
                                             Stripped.Sequence, Precursor.Id, Proteotypic, Precursor.Normalised,    
                                             Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value, 
                                             Lib.PG.Q.Value)] |> data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), 
                                                                                      c("ProteinGroup", "Gene", "ProteinDescription"))
    } else {return(message("ERROR: No Input File Found"))}
    
    if(file.exists("Sample_Rename.csv") == TRUE){
      Sample_Metadata <- data.table::fread("Sample_Rename.csv") 
      PrecursorData$Sample <- Sample_Metadata$Renamed[match(unlist(PrecursorData$Run), Sample_Metadata$Run)]
      PrecursorData <- PrecursorData[!is.na(Sample)]
    }
    
    if(any(grepl("SILAC-", PrecursorData$Precursor.Id))){
      return(message("ERROR: Label-Free Processing on SILAC Data"))
    }
  }
  message("Defining Metadata")
  {
    PrecursorData[, Cell := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\1", Sample)]
    PrecursorData[, Conc := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\2", Sample)]
    PrecursorData[, Drug := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\3", Sample)]
    PrecursorData[, Time := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\4", Sample)]
    PrecursorData[, Replicate := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\5", Sample)]
    PrecursorData[, Sample := gsub("_0", "", Sample)]
    PrecursorData[, Condition := gsub("(.*)_R\\d", "\\1", Sample)]

    # Order Samples
    Sample_Order <- unique(PrecursorData[, .(Sample, Condition, Replicate)]) |> dplyr::arrange(!grepl(CtlGroup, Condition), Condition, Replicate)
    PrecursorData$Sample <- factor(PrecursorData$Sample, levels = Sample_Order$Sample)
  }
  message("Filtering Precursors")
  {
    if(ProteotypicFiltering == T){PrecursorData <- PrecursorData[Proteotypic >= 1]} else {PrecursorData <- PrecursorData}
    PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01] 
    PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
  }
  message("Compiling & Exporting Data")
  {
    data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
    LFQ <- Calculate_LFQ(PrecursorData, "LFQ")
    tot_intensities <- PrecursorData[, .(Intensity = sum(Precursor.Normalised)), .(ProteinGroup, Sample)]
    PrecursorCounts <- PrecursorData[, .(N_precursors = data.table::uniqueN(Precursor.Id)), .(ProteinGroup, Sample)]
    proteotypic_counts <- PrecursorData[, .(N_precursors_proteotypic = sum(Proteotypic)), .(ProteinGroup, Sample)]
    annotations <- unique(PrecursorData[, .(ProteinGroup, Sample, Condition, Replicate, ProteinDescription, Gene)])
    ProteinData <- Reduce(Merge_PrecursorData, list(LFQ, tot_intensities, PrecursorCounts, proteotypic_counts, annotations))
    data.table::fwrite(ProteinData, "LF_DIANN_Output.csv.gz")
  }
  message("Plotting Intensities")
  {
    Intensities_Data <- data.table::rbindlist(list(
      PrecursorData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Precursor.Normalised), Type = "Precursor Quantity")],
      ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ), Type = "Protein MaxLFQ")],
      ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity), Type = "Protein Intensity")]
    ))
    Intensities_Data[, Type := factor(Type, levels = c("Precursor Quantity", "Protein MaxLFQ", "Protein Intensity"))]   # Set plotting order
    
    Intensity_Plot <- Intensities_Data |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = `log2 quantity`, colour = Condition)) + 
      ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
      ggplot2::facet_wrap("Type", scales = "free_x") + ggplot2::ylab("Log2 Value") + ggplot2::coord_flip() + 
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }
  message("Plotting Precursor, Peptide & Protein Counts")
  {
    all_counts <- PrecursorData[, lapply(.SD, data.table::uniqueN), .(Sample, Condition), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")]
    facet_labels <- ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups")) # Create plot labellers for facets
    
    Count_Plot <- data.table::melt.data.table(all_counts, id.vars = c("Sample","Condition"), value.name = "IDs") |> 
      ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = IDs/1000, fill = Condition, label = format(IDs, big.mark = ",", scientific = FALSE))) +
      ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.2) +
      ggplot2::facet_wrap("variable", scales = "free_x", labeller = facet_labels) + ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + 
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }
  message("Plotting Data Completeness")
  {
    data_completeness <- rbind(Count_Proteins(ProteinData, "All"),
                               Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                               Count_Proteins(ProteinData[N_precursors_proteotypic >= 2], "≥ 2 Proteotypic"))
    
    NAs_Plot <- data_completeness |> ggplot2::ggplot(ggplot2::aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors))+
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("All" = "black", "≥ 2" = "darkgrey", "≥ 2 Proteotypic" = "orange3")) +
      ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") +
      ggplot2::scale_x_continuous( breaks = seq(1, 1000, 1)) +
      ggplot2::scale_y_continuous( limits = c(0, max(data_completeness$cumulative_protein_N )/1000))+
      ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     legend.position = "inside", legend.position.inside = c(0.25, 0.25))
  }
  message("Plotting Data Skewing (No Output)")
  {
    SkewData <- data.table::data.table(PrecursorData[, .(Sample, Condition, Replicate)] |> unique())
    SkewData[, `:=`(PearsonsSkew = 0, Median = 0, Mean = 0)]
    for(RowIndex in 1:nrow(SkewData)){
      SampleID <- SkewData$Sample[RowIndex]
      SkewData$PearsonsSkew[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ])
      SkewData$Median[RowIndex] <- median(ProteinData[Sample == SampleID, LFQ], na.rm = T)
      SkewData$Mean[RowIndex] <- mean(ProteinData[Sample == SampleID, LFQ], na.rm = T)
    }
    
    SkewPlot <- ProteinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = LFQ, colour = Condition)) + 
      ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
      ggplot2::ylab("Protein LFQ") + ggplot2::coord_flip() + ggplot2::scale_y_log10() +
      ggplot2::geom_text(data = SkewData, aes(y = Median*2.5, label = paste0("Pearson's Skew\n", round(PearsonsSkew, 2))), size = 6) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Calculating Missed Trypsinisation Sites")
  {
    PrecursorData <- PrecursorData[, MissedCleavage := grepl("[RK][^P]", Stripped.Sequence)]
    PrecursorCount <- PrecursorData[, N := .N, by = Sample]
    TrypsinData <- PrecursorData[MissedCleavage == TRUE, .(N = .N), by = .(Sample, Condition, MissedCleavage)]
    
    TrypsinData$PercentPrecursors <- 0
    for(i in 1:nrow(TrypsinData)){
      TrypsinData$PercentPrecursors[i] <- TrypsinData$N[i]/PrecursorCount$N[which(PrecursorCount$Sample == TrypsinData$Sample[i])]
    }
    
    MissedCleavage_Plot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentPrecursors*100, fill = Condition)) + 
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::labs(x = "", y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::coord_flip() + 
      ggplot2::lims(y = c(0, max(TrypsinData$PercentPrecursors*200)))
  }
  message("Plotting Precursor & Protein Variation")
  {
    precursor_CVs <- PrecursorData[, .(CV = Calculate_CV(Precursor.Normalised), N = .N), .(Precursor.Id, Condition)]
    precursor_CVs <- precursor_CVs[, Rank := data.table::frank(CV), Condition] 
    precursor_CVs$ID <- "Precursors"
    
    protein_CVs <- ProteinData[, .(CV = Calculate_CV(LFQ), N = .N), .(ProteinGroup, Condition)]  
    protein_CVs <- protein_CVs[, Rank := data.table::frank(CV), Condition] 
    protein_CVs$ID <- "Protein Groups"
    
    all_CVs <- precursor_CVs[, Precursor.Id := NULL] |> rbind(protein_CVs[, ProteinGroup := NULL])
    
    Variation_Plot <- all_CVs |> ggplot2::ggplot(ggplot2::aes(x = Rank/1000, y = CV, colour = Condition))+
      ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Coeff. of Variation [%]") +
      ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + ggplot2::coord_cartesian(ylim = c(0,50)) + 
      ggplot2::facet_wrap(~ID, scales = "free") +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     legend.position = "inside", legend.position.inside = c(0.6, 0.6))
  }
  message("Exporting QC Plots")
  {
    plot_design <- "AAAA
                    BBBB
                    CCDD
                    EEEE"
    
    DIANN_QC_Plot <- Intensity_Plot + Count_Plot + patchwork::free(NAs_Plot, type = "label") + 
      patchwork::free(MissedCleavage_Plot, type = "label") + patchwork::free(Variation_Plot, type = "label") +
      patchwork::plot_layout(design = plot_design) + patchwork::plot_annotation(tag_levels = "A")
    
    pdf("PrecursorQC_DIANN_Plot.pdf", width = 18, height = 20)
    print(DIANN_QC_Plot)
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_LF_PrecursorBiochemistry <- function(InputDirectory){
  start.time <- Sys.time()
  set.seed(123)
  message("Loading Precursor Data")
  {
    setwd(InputDirectory)
    PrecursorData <- data.table::fread(list.files(pattern = ".*PrecursorData.csv.*"))[, .(Run, ProteinGroup, ProteinDescription, 
                                                                                      Gene, Stripped.Sequence, Precursor.Normalised,
                                                                                      Sample, Condition, Cell, Drug, Time, Replicate,
                                                                                      Precursor.Length)] |> 
      setnames(c("Stripped.Sequence", "Precursor.Length"), c("Sequence", "Length"))
  }
  message("Annotating with Biochemical Measures")
  {
    PrecursorData[, `:=`(Aliphatic_Score = Peptides::aIndex(Sequence), Boman_Interaction_Score = Peptides::boman(Sequence),
                         Hydrophobicity_Score = Peptides::hydrophobicity(Sequence, scale = "KyteDoolittle"),
                         Instability_Score = Peptides::instaIndex(Sequence), MW = Peptides::mw(Sequence), 
                         pI = Peptides::pI(Sequence, pKscale = "Dawson"))]
    
    for(ColIndex in which(colnames(PrecursorData) == "Length"):ncol(PrecursorData)){
      if(is.numeric(PrecursorData[, get(colnames(PrecursorData)[ColIndex])])){
        message(paste0("Analysing ", colnames(PrecursorData)[ColIndex]), " Trend")
        SubsetData = PrecursorData[, .(Sequence, Length, Condition, Replicate, Precursor.Normalised, get(colnames(PrecursorData)[ColIndex]))]
        SubsetData |> data.table::setnames("V6", "Subset")
        
        pval <- summary(stats::lm(Precursor.Normalised ~ Subset, data = SubsetData))$coefficients[2,4]
        
        TrendPlot <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Subset, y = Precursor.Normalised, colour = Condition, shape = factor(Replicate))) + 
          ggplot2::geom_smooth(method = "lm", alpha = 0.1) + 
          ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Precursor.Normalised], na.rm = T)*0.93, 
                            label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) +
          ggplot2::labs(x = paste0("Precursor ", gsub("_", " ", colnames(PrecursorData)[ColIndex])), y = "Normalised Precursor Intensity") + 
          ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
        
        pdf(paste0(colnames(PrecursorData)[ColIndex],"_Trend.pdf"), width = 12, height = 10)
          print(TrendPlot)
          print(TrendPlot + ggplot2::labs(x = "", y = ""))
        Proteopedia::Reset_Dev()
      }
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Analyse_LF_Proteins <- function(InputDirectory, ExpGroups, CtlGroups, MinPrecursors = 2, ImputationQ = 0.01, ImputationSigma = 1){
  # Auto-Read Directory & Produce Metadata
  start.time <- Sys.time()
  set.seed(123)
  for(ExpGroup in ExpGroups){
    for(CtlGroup in CtlGroups){
      message("Loading Input File")
      {
        setwd(InputDirectory)
        spectra_read <- data.table::fread(list.files(pattern = "DIANN_Output.csv"))[, Log2LFQ := log2(LFQ)]
        
        spectra_read |> data.table::setnames(
          c(colnames(spectra_read)[grepl("protein.*group", ignore.case = T, colnames(spectra_read))], 
            colnames(spectra_read)[grepl("protein.*desc", ignore.case = T, colnames(spectra_read))], 
            colnames(spectra_read)[grepl("Gene", ignore.case = T, colnames(spectra_read)) & !grepl("group", ignore.case = T, colnames(spectra_read))]),
          c("ProteinGroup", "ProteinDescription", "Gene"))
        meta <- spectra_read[grepl(ExpGroup, Condition)|grepl(CtlGroup, Condition), .(Sample, Condition, Replicate)] |> dplyr::distinct() |> 
          dplyr::arrange(c(which(grepl(CtlGroup, Condition)), which(grepl(ExpGroup, Condition))))
        protein_name_conversion <- spectra_read[, .(ProteinGroup, ProteinDescription, Gene)] |> dplyr::distinct()
        
        spectra_read <- spectra_read[Sample %in% meta$Sample]
        spectra_read[, Condition := factor(Condition, levels = c(unique(spectra_read[grepl(CtlGroup, Condition), Condition]), 
                                                                 unique(spectra_read[grepl(ExpGroup, Condition), Condition])))]
        
        if(dir.exists(paste0(getwd(),"/",ExpGroup,"_vs_",CtlGroup,"_Output")) == TRUE){
          unlink(paste0(getwd(),"/",ExpGroup,"_vs_",CtlGroup,"_Output"), recursive = TRUE)
        }
        dir.create(paste0(getwd(),"/",ExpGroup,"_vs_",CtlGroup,"_Output"), showWarnings = TRUE)
        setwd(paste0(getwd(),"/",ExpGroup,"_vs_",CtlGroup,"_Output"))
        data.table::fwrite(meta, file = "Sample_Metadata.csv")
      }
      message("Performing PCA")
      {
        prePCA_data <- spectra_read[, .(ProteinGroup, Gene, Sample, Log2LFQ)] |> 
          tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Log2LFQ, names_from = Sample, values_fill = NA) |> 
          tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
        
        prePCA_summary <- summary(prePCA_data)$importance
        prePCA_data <- data.table::data.table(prePCA_data$x, keep.rownames = "Sample") |> data.table::merge.data.table(meta)
        
        prePCA12_plot <- ggplot2::ggplot(prePCA_data, ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = as.factor(Replicate))) +
          ggplot2::geom_point(size = 4) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
          ggplot2::labs(x = paste("PC1 [", round(prePCA_summary[rownames(prePCA_summary) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
               y = paste("PC2 [", round(prePCA_summary[rownames(prePCA_summary) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
          ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank(), legend.position = "none")
        
        prePCA34_plot <- ggplot2::ggplot(prePCA_data, ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = as.factor(Replicate))) +
          ggplot2::geom_point(size = 4) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
          ggplot2::labs(x = paste("PC3 [", round(prePCA_summary[rownames(prePCA_summary) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""), 
                        y = paste("PC4 [", round(prePCA_summary[rownames(prePCA_summary) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
          ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
      }
      message("Filtering Proteins")
      {
        filtering_data <- spectra_read[,.(N_Samples = .N, Min_Precursors = min(N_precursors), N_conditions = data.table::uniqueN(Condition)), ProteinGroup]
        
        filtering_data[, UniProtID := gsub("[;-].*", "", ProteinGroup)]
        
        Retained_1 <- filtering_data[N_Samples == data.table::uniqueN(spectra_read$Sample) & Min_Precursors >= MinPrecursors]
        Retained_2 <- filtering_data[N_Samples == data.table::uniqueN(spectra_read$Sample)-1 & Min_Precursors >= (MinPrecursors+1)]
        Retained_3 <- filtering_data[N_Samples == floor(data.table::uniqueN(spectra_read$Sample)/2) & Min_Precursors >= (MinPrecursors+1) & N_conditions >= 1]
        
        Retained_Proteins <- c(Retained_1$ProteinGroup, Retained_2$ProteinGroup, Retained_3$ProteinGroup)
        Filtered_Proteins <- spectra_read[ProteinGroup %!in% Retained_Proteins]
        spectra_read <- spectra_read[ProteinGroup %in% Retained_Proteins]
        
        data.table::fwrite(Filtered_Proteins, file = "Filtered_Proteins.csv")
        data.table::fwrite(Retained_3, file = "Imputed_Proteins.csv")
        
        #Create Upset of Retained Proteins
        Filtering_Plot_Data_A <- spectra_read[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")]
        Filtering_Plot_Data_B <- Filtered_Proteins[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")]
        Filtering_Plot_Data <- Filtering_Plot_Data_A |> rbind(Filtering_Plot_Data_B) |> 
          dplyr::group_by(Sample, Condition, Replicate, Inclusion) |> dplyr::summarise(n = dplyr::n()) |> 
          dplyr::mutate(Colour = ifelse(Condition == CtlGroup, RColorBrewer::brewer.pal(3, "Set1")[2], RColorBrewer::brewer.pal(3, "Set1")[1]))
        
        Protein_Counts_Bar <- Filtering_Plot_Data |> ggplot2::ggplot(ggplot2::aes(x = Sample, y = n, fill = Colour, alpha = Inclusion)) + 
          ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_identity(guide = "none") + 
          ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") + 
          ggplot2::geom_text(ggplot2::aes(label = n), size = 6, colour = ifelse(Filtering_Plot_Data$Inclusion == "Retained", "white","black"), position = ggplot2::position_stack(), vjust = ifelse(Filtering_Plot_Data$Inclusion == "Retained", 1.5, -0.5)) +                                  
          ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x") +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) + ggplot2::labs(x = NULL, y = "Count", fill = NULL) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), strip.text.x = ggplot2::element_blank(), 
                                              strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0,"line"))
        
        Upset_Plot <- spectra_read[, .(Sample = list(Sample)), by = ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
          ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) + 
          ggupset::scale_x_upset(order_by = "degree", reverse = TRUE, sets = spectra_read[order(Condition), unique(Sample)]) + 
          ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) 
      }
      message("Performing Median Normalisation")
      {
        spectra_mean <- median(spectra_read$Log2LFQ, na.rm = TRUE)
        spectra_read[, Log2LFQ_Norm := Log2LFQ - median(Log2LFQ, na.rm = TRUE) + spectra_mean, by = Sample]  
        spectra_read[, Colour := data.table::fifelse(Condition == CtlGroup, RColorBrewer::brewer.pal(3, "Set1")[2], RColorBrewer::brewer.pal(3, "Set1")[1])]
        
        Normalisation_Plot <- ggplot2::ggplot(data.table::melt.data.table(spectra_read, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                     ggplot2::aes(x = Condition, fill = Colour, y = value, group = Sample))+
          ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
          ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::ggtitle("Normalisation") +
          ggplot2::scale_fill_identity(guide = "none") + ggplot2::labs(x = NULL, y = "Log2LFQ")
      }
      message("Imputing NA Values")
      {
        spectra_read <- spectra_read[,.(ProteinGroup, Gene, Sample, Log2LFQ_Norm)] |> 
          tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Log2LFQ_Norm, names_from = Sample, values_fill = NA) |> 
          data.table::data.table()
        
        spectra_imp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(spectra_read, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma), 
                                  keep.rownames = "ProteinGroup")
        spectra_all <- data.table::merge.data.table(data.table::melt.data.table(spectra_read, id.vars = "ProteinGroup", value.name = "Measured_LFQ", variable.name = "Sample"),
                                         data.table::melt.data.table(spectra_imp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ", variable.name = "Sample"))
        spectra_all[, Data := data.table::fifelse(is.na(Measured_LFQ), "Imputed","Measured")]
    
        ImputedDensity <- spectra_all |> ggplot2::ggplot(ggplot2::aes(x = Imputed_LFQ, fill = Data)) +
          ggplot2::geom_density(adjust = 2, alpha = 0.5) + ggplot2::scale_fill_manual(values = c("Measured" = "black", "Imputed" = "magenta3")) +
          ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))
        
        imputedNAs <- data.table::melt.data.table(spectra_read[ProteinGroup %in% Retained_3$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
        imputedNAs <- imputedNAs[is.na(value)]
        imputedNAs[, value := NULL] 
        imputedNAs <- imputedNAs |> dplyr::group_by(ProteinGroup) |> dplyr::summarise(vector = paste(Sample, collapse = ", "))
        data.table::fwrite(imputedNAs, "Imputed_LFQs.csv")
    
        spectra_all[, Log2LFQ := data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ)]
        spectra_all <- spectra_all[, .(ProteinGroup, Sample, Log2LFQ)]
      }
      message("Performing Paired T-Testing")
      {
        spectra_all <- spectra_all |> data.table::merge.data.table(meta[, .(Condition, Sample)], by = "Sample")
        spectra_all[, LFQ := 2^(Log2LFQ)]
        
        spectra_Ttest <- spectra_all |> dplyr::group_by(Condition, ProteinGroup) |> dplyr::summarise(
          N = dplyr::n(), Log2MeanLFQ = log2(mean(LFQ)), CV = (sd(LFQ)/mean(LFQ))*100) |> 
          tidyr::pivot_wider(id_cols = ProteinGroup, names_from = Condition, values_from = c(N, Log2MeanLFQ, CV)) 
        
        ctl_col_index <- which(grepl("Log2MeanLFQ_", colnames(spectra_Ttest)) & grepl(CtlGroup, colnames(spectra_Ttest)))
        exp_col_index <-  which(grepl("Log2MeanLFQ_", colnames(spectra_Ttest)) & grepl(ExpGroup, colnames(spectra_Ttest)))
        spectra_Ttest$Log2FC <- spectra_Ttest[exp_col_index] - spectra_Ttest[ctl_col_index]
        
        Ttest_Output <- spectra_all[, P.Value := stats::t.test(Log2LFQ ~ Condition)$p.value, by = ProteinGroup]
        spectra_Ttest <- data.table::merge.data.table(spectra_Ttest |> data.table(), Ttest_Output, by = "ProteinGroup")
        spectra_Ttest[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
        spectra_Ttest <- spectra_Ttest |> data.table::merge.data.table(protein_name_conversion) |> data.table()
        spectra_Ttest[, Log2FC := as.numeric(Log2FC)]
        data.table::fwrite(spectra_Ttest |> dplyr::distinct(), "Paired_T-Test_Output.csv")
      }
      message("Fitting to Linear Model & Exporting Data")
      {
        design <- stats::model.matrix(~0 + Condition, data = meta)
        colnames(design) <- gsub("Condition", "", colnames(design))
        rownames(design) <- meta$Sample
        design <- design[, (c(which(grepl(CtlGroup, colnames(design))), which(grepl(ExpGroup, colnames(design)))))]
        
        contr.matrix <- matrix(nrow = 2, ncol = 1)
        contr.matrix[,1] <- c(-1,1)
        dimnames(contr.matrix) <- list("Levels" = colnames(design), "Contrasts" = "comp")
        
        limma_input <- spectra_all |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |> 
          data.table::setcolorder(neworder = c("ProteinGroup", meta$Sample))
    
        efit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(limma_input, design), contr.matrix))
        
        if(!is.finite(efit$df.prior)){
          message("Warning: Limma df.prior is Infinite")
        }
        
        Mean_Var_Data <- data.table(efit$genes, "Mean" = efit$Amean,
                                    "Variance" = sqrt(efit$sigma))
        Mean_Var_Data[, Data := data.table::fifelse(ProteinGroup %in% imputedNAs$ProteinGroup, "Imputed", "Measured")]
        
        Mean_Var_TrendPlot <- Mean_Var_Data |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance, colour = Data)) +
          ggplot2::geom_point() + ggplot2::scale_colour_manual(values = c("Measured" = "black", "Imputed" = "magenta3")) + ggplot2::labs(x = "Mean Log2LFQ", y = "Variance", colour = NULL)
        
        LimmaOutput <- limma::topTable(efit, coef=1, adjust.method = "BH", n=Inf) |> data.table() |> data.table::setnames("logFC", "Log2FC")
        LimmaOutput <- LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = TRUE),] |> data.table::merge.data.table(protein_name_conversion, all.x = T)
        LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
        LimmaOutput <- data.table::merge.data.table(LimmaOutput, limma_input, by = "ProteinGroup", all.x = TRUE) 
        LimmaOutput$Isoforms <- 1
        for(i in 1:nrow(LimmaOutput)){
          if(length(stringr::str_extract_all(LimmaOutput$ProteinGroup[i], "-\\d", simplify = T)) > 0){
            LimmaOutput$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaOutput$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          } else {LimmaOutput$Isoforms[i] <- 1}
        }    
        LimmaOutput[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
        LimmaOutput[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
        LimmaOutput[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
        data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
      }
      message("Generating Volcano Plots")
      {
        MeanLog2FC <- round(mean(LimmaOutput$Log2FC, na.rm = T), digits = 3)
        
        LimmaOutput <- LimmaOutput |> dplyr::arrange(desc(abs(t)))
        
        all_limma_volcano <- LimmaOutput |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
          ggplot2::scale_colour_manual("black") + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) +
          ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") +
          ggplot2::annotate("text", x = min(LimmaOutput$Log2FC)*0.9, y = 0, label = paste0("Mean Log2FC\n",MeanLog2FC), size = 5) +
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        imputed_limma_data <- LimmaOutput |> dplyr::mutate(Imputed = ifelse(ProteinGroup %in% imputedNAs$ProteinGroup, "Yes","No"))
        
        imputed_limma_volcano <- imputed_limma_data |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual(values = c("No" = "black", "Yes" = "magenta3"), 
                                                                               labels = c("No" = "Measured", "Yes" = "Imputed"), guide = "none") +
          ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) +
          ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") +
          ggplot2::annotate("text", x = min(LimmaOutput$Log2FC), y = 0, label = paste0("Mean Log2FC\n",MeanLog2FC), size = 5) +
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        pdf("Limma_Volcanoes.pdf", height = 10, width = 14)  
        print(all_limma_volcano)
        print(imputed_limma_volcano)
        Proteopedia::Reset_Dev()
        
        MeanLog2FC <- round(mean(spectra_Ttest$Log2FC), digits = 3)
        
        all_Ttest_volcano <- spectra_Ttest |> 
          ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
          ggplot2::scale_colour_manual("black") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) +
          ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") +
          ggplot2::annotate("text", x = min(spectra_Ttest$Log2FC)*0.9, y = 0, label = paste0("Mean Log2FC\n",MeanLog2FC), size = 5) +
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        imputed_Test_data <- spectra_Ttest |> dplyr::mutate(Imputed = ifelse(ProteinGroup %in% imputedNAs$ProteinGroup, "Yes","No"))
        
        ImputedTtestVolcano <- imputed_Test_data |> 
          ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
          ggplot2::scale_colour_manual(values = c("No" = "black", "Yes" = "magenta3"), labels = c("No" = "Measured", "Yes" = "Imputed"), guide = "none") +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) +
          ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") +
          ggplot2::annotate("text", x = min(imputed_Test_data$Log2FC), y = 0, label = paste0("Mean Log2FC\n",MeanLog2FC), size = 5) +
          Proteopedia::Add_AbundanceAxes() + ggtitle("Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        pdf("PairedTtest_Volcanoes.pdf", height = 10, width = 14)
        print(all_Ttest_volcano)
        print(ImputedTtestVolcano)
        Proteopedia::Reset_Dev()
      }
      message("Exporting QC Plot")
      {
        plot_design <- "AABBCCCC\nDDEEFFFF\nGGGGGGGG"
        Limma_QC_Plot <-  patchwork::free(prePCA12_plot, type = "label") + patchwork::free(prePCA34_plot, type= "label") + 
          Protein_Counts_Bar + patchwork::free(Upset_Plot, type = "label") + 
          patchwork::free(Normalisation_Plot) + patchwork::free(ImputedDensity) + patchwork::free(Mean_Var_TrendPlot, type = "label") + #free(all_volcano, type = "label") +
          patchwork::plot_layout(design = plot_design) + patchwork::plot_annotation(tag_levels = "A")
        
        pdf("ProteinQC_Limma_Plot.pdf", width = 18, height = 20)
        print(Limma_QC_Plot)
        Proteopedia::Reset_Dev()
      }
      message("Exporting HTML Volcano Plot")
      {
        InteractiveData <- data.table("Protein" = LimmaOutput$Gene, 
                                      "Log2FC" = round(LimmaOutput$Log2FC, digits = 2),
                                      "PValue" = LimmaOutput$P.Value,
                                      "Significance" = factor(LimmaOutput$Significance),
                                      "Scien_PValue" = formatC(LimmaOutput$P.Value, format = "e", digits = 2),
                                      "GeneGroup" = LimmaOutput$GeneGroup)
        
        InteractiveData[, URL := paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup))]
        
        pHC <- highcharter::hchart(InteractiveData, "scatter", 
                                   highcharter::hcaes(x = Log2FC, y = -log10(PValue), group = Significance)) |> 
          highcharter::hc_chart(zoomType = "xy") |>  # Enables zooming on both x and y axes
          highcharter::hc_xAxis(title = list(text = paste0(unique(meta$Condition)[1]," vs ", unique(meta$Condition)[2]," Log2 Fold-Difference")), 
                   lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black", gridLineWidth = 0 ) |>
          highcharter::hc_yAxis(title = list(text = "-Log10 P-Value"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                   tickColor = "black", gridLineWidth = 0 ) |>
          highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}") |>
          highcharter::hc_plotOptions(scatter = list(marker = list(radius = 3),
                                        states = list(hover = list(enabled = TRUE), inactive = list(enabled = FALSE)),   # Don"t dim inactive series
                                        point = list(events = list( click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))) |>
          highcharter::hc_colors(c("#999999", "#880000","#0033FF"))
        htmlwidgets::saveWidget(pHC, "InteractiveVolcanoPlot.html")
      }
      message("Exporting Analysis Parameters")
      Analysis_Parameters <- data.table("Experimental_Condition" = paste0(ExpGroup), "Control_Condition" = paste0(CtlGroup),
                                        "Min_Precursors" = paste0(MinPrecursors), "Imputation_Q-value" = ImputationQ,
                                        "Imputation_Sigma" = ImputationSigma)
      data.table::fwrite(Analysis_Parameters, "Analysis_Parameters.csv")
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
########### MS Analysis: Static SILAC DIA-NN To Protein Data Functions ####################################################################################################################################
#' @export
Process_StaticSILAC_DIANN <- function(InputDirectory, CtlGroup, ProteotypicFiltering = F){
    set.seed(123)
    start.time <- Sys.time()
    message("Importing DIA-NN Report File")
    {
      setwd(InputDirectory)
      if(length(list.files(pattern = "report.tsv")[1]) > 0){
        InputFile <- list.files(pattern = "report.tsv")[1]
        PrecursorData <- data.table::fread(InputFile)[, .(Run, Protein.Group, Protein.Ids, First.Protein.Description, Genes, 
                                              Stripped.Sequence, Precursor.Id, Proteotypic, Precursor.Quantity,    
                                              Precursor.Translated, Channel.Q.Value,
                                              Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value, 
                                              Lib.PG.Q.Value)] |> data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), 
                                                                            c("ProteinGroup", "Gene", "ProteinDescription"))
      } else {return(message("ERROR: No InputFile Found"))}
      
      if(file.exists("Sample_Rename.csv") == TRUE){
        PrecursorData <- PrecursorData |> data.table::merge.data.table(fread("Sample_Rename.csv") |> data.table::setnames("Renamed", "Sample"))
      }
    }
    message("Defining Metadata")
    {
      PrecursorData[, Cell := gsub("(.*)_(.*)_(.*)_(.*)_R(.*)","\\1", Sample)]
      PrecursorData[, Drug := gsub("(.*)_(.*)_(.*)_(.*)_R(.*)","\\2", Sample)]
      PrecursorData[, Conc := gsub("(.*)_(.*)_(.*)_(.*)_R(.*)","\\3", Sample)]
      PrecursorData[, Time := gsub("h", "", gsub("(.*)_(.*)_(.*)_(.*)_R(.*)","\\4", Sample))]
      PrecursorData[, Replicate := gsub("(.*)_R(.*)", "\\2", Sample)]
      PrecursorData[, Sample := gsub("0_", "", Sample)]
      PrecursorData[, Condition := gsub("(.*)_R(.*)","\\1", Sample)]
      
      Sample_Order <- dplyr::arrange(unique(PrecursorData[, .(Sample, Condition, Replicate)]), !grepl(CtlGroup, Condition), Condition, Replicate)
      PrecursorData[, Sample := factor(Sample, levels = Sample_Order$Sample)]
    }
    message("Filtering Precursors")
    {
      if(ProteotypicFiltering == T){PrecursorData <- PrecursorData[Proteotypic >= 1]} else {PrecursorData <- PrecursorData}
      PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.01]
      PrecursorData[Precursor.Quantity == 0, Precursor.Quantity   := NA]
      PrecursorData[Precursor.Translated == 0, Precursor.Translated := NA]
      PrecursorData <- PrecursorData[!is.na(Precursor.Quantity)]
    }
    message("Normalising Precursor Quantities")
    {
      PrecursorData[, Precursor.Quantity   := Precursor.Quantity/sum(Precursor.Quantity)*PrecursorData[, sum(Precursor.Quantity), Run][, median(V1)], by = Run]
      PrecursorData[, Precursor.Translated := Precursor.Translated/sum(Precursor.Translated, na.rm = TRUE)*PrecursorData[, sum(Precursor.Quantity), Run][, median(V1)], Run]
    }
    message("Processing SILAC Labels")
    {
      if(PrecursorData[Precursor.Id %like% "SILAC-.-L" & Precursor.Id %like% "SILAC-.-H", .N] != 0){
        message("ERROR: Multi-Label Precursors Detected")
        MultiLabelDetection = T} else {MultiLabelDetection = F}
      
      PrecursorData[data.table::like(Precursor.Id, "SILAC-.-L"), Label := "L"]
      PrecursorData[data.table::like(Precursor.Id, "SILAC-.-H"), Label := "H"]
      PrecursorData[,Precursor.Id.nolabels := gsub("SILAC-.-.", "SILAC", Precursor.Id)]
      
      FullIsotopeRatio <- PrecursorData[!is.na(Precursor.Translated) & !is.na(Label), .(LabelIntensity = sum(Precursor.Translated)), by = list(Sample, Condition, Replicate, Label)] |> 
        data.table::merge.data.table(PrecursorData[!is.na(Precursor.Translated) & !is.na(Label), .(TotalIntensity = sum(Precursor.Translated)), by = list(Sample, Condition, Replicate)])
      FullIsotopeRatio[, Prop := LabelIntensity/TotalIntensity]
      
      IsotopeIncorporation <- FullIsotopeRatio |> ggplot2::ggplot(ggplot2::aes(x = gsub("_", " ", Sample), y = Prop, fill = Condition, alpha = Label)) +
        ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::scale_alpha_manual(values = c("H" = 0.5, "L" = 1), guide = "none") + ggplot2::ylab("Label Incorporation Ratio") +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      IntensitiesRatioData <- PrecursorData |> data.table::dcast(Run+ProteinGroup+Protein.Ids+ProteinDescription+Gene+Stripped.Sequence+
                                                                 Proteotypic+Sample+Cell+Drug+Conc+Time+Replicate+Condition+Precursor.Id.nolabels ~Label, 
                                                                 value.var = "Precursor.Translated")
      
      IntensitiesRatioData <- IntensitiesRatioData[!is.na(H) & H > 0 & !is.na(L) & L > 0]
      IntensitiesRatioData <- IntensitiesRatioData[,.(H = sum(H, na.rm = T), L = sum(L, na.rm = T), .N), .(ProteinGroup, ProteinDescription, Gene, Sample)]
      IntensitiesRatioData[, Ratio := H/L]
      data.table::fwrite(IntensitiesRatioData, "IntensitiesRatioData.csv")
    }
    message("Compiling & Exporting Data")
    {
      PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
      data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
      
      LFQ_T <- Proteopedia::Calculate_LFQ(PrecursorData, "LFQ", SILAC = T)
      Intensity_T <- PrecursorData[,.(Intensity = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_T <- PrecursorData[, .(N_precursors = uniqueN(Precursor.Id), N_precursors_proteotypic = sum(Proteotypic)), .(ProteinGroup, Sample)]
      
      LFQ_L <- Proteopedia::Calculate_LFQ(PrecursorData[Label == "L"], "LFQ_L", SILAC = T)
      Intensity_L <- PrecursorData[Label == "L",.(Intensity_L = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_L <- PrecursorData[Label == "L" , .(N_precursors_L = uniqueN(Precursor.Id), N_precursors_proteotypic_L = sum(Proteotypic)), .(ProteinGroup, Sample)]
      
      LFQ_H <- Proteopedia::Calculate_LFQ(PrecursorData[Label == "H"], "LFQ_H", SILAC = T)
      Intensity_H <- PrecursorData[Label == "H",.(Intensity_H = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_H <- PrecursorData[Label == "H" , .(N_precursors_H = uniqueN(Precursor.Id), N_precursors_proteotypic_H = sum(Proteotypic)), .(ProteinGroup, Sample)]

      
      #LFQ_Ratio <- PrecursorData |> data.table::dcast(Run+ProteinGroup+Protein.Ids+ProteinDescription+Gene+Stripped.Sequence+
      #                                     Proteotypic+Sample+Cell+Drug+Conc+Time+Replicate+Condition+Precursor.Id.nolabels ~Label, 
      #                                   value.var = "Precursor.Translated") 
      #LFQ_Ratio[, Precursor.Quantity := H/L]
      #LFQ_Ratio |> data.table::setnames("Precursor.Id.nolabels", "Precursor.Id")
      #LFQ_Ratio <- Proteopedia::Calculate_LFQ(LFQ_Ratio[!is.na(Precursor.Quantity)], "LFQ_Ratio", SILAC = T)
      
      Annotations <- unique(PrecursorData[, .(ProteinGroup, Run, Sample, Condition, Cell, Drug, Conc, Time, Replicate, ProteinDescription, Gene)])
      ProteinData <- Reduce(Proteopedia::Merge_PrecursorData, list(LFQ_T, LFQ_H, LFQ_L, Intensity_T, Intensity_H, Intensity_L, Counts_T, Counts_H, Counts_L, Annotations))                               
      ProteinData[, LFQ_Ratio := LFQ_H/LFQ_L]
      data.table::fwrite(ProteinData, "SILAC_DIANN_Output.csv.gz")
    }
    message("Plotting Intensities")
    {  
      PrecursorData[, Label := gsub("L", "Light", Label)]
      PrecursorData[, Label := gsub("H", "Heavy", Label)]
      
      IntensitiesData <- data.table::rbindlist(list(
        PrecursorData[!is.na(Label), .(Sample, Condition, Replicate, Label, `log2 quantity` = log2(Precursor.Quantity), Type = "Precursor Quantity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ_L), Label = "Light", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity_L), Label = "Light", Type = "Protein Intensity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ_H), Label = "Heavy", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity_H), Label = "Heavy", Type = "Protein Intensity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ), Label = "Total", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity), Label = "Total", Type = "Protein Intensity")]
      ), use.names = TRUE)
      IntensitiesData[, Type := factor(Type, levels = c("Precursor Quantity", "Max. Protein LFQ", "Protein Intensity"))]
      
      IntensityPlot <- IntensitiesData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = `log2 quantity`, colour = Condition, alpha = Label)) + 
        ggplot2::geom_boxplot(outliers = FALSE)+
        ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1, "Total" = 1), guide = "none") + 
        ggplot2::facet_grid(cols = ggplot2::vars(Type), rows = ggplot2::vars(Label), scales = "free_x") + 
        ggplot2::ylab(expression("Log"[2]~"Value")) + ggplot2::coord_flip() + ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Plotting Precursor, Peptide & Protein Counts")
    {
      CountPlot <- data.table::melt(PrecursorData[!is.na(Label), lapply(.SD, uniqueN), .(Sample, Condition, Label), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")], 
                        id.vars = c("Sample","Condition", "Label"), value.name = "IDs") |> 
        ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = IDs/1000, fill = Condition, alpha = Label, 
                                     label = format(IDs, big.mark = ",", scientific = FALSE))) +
        ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.1) +
        ggplot2::facet_grid(ggplot2::vars(Label), ggplot2::vars(variable), scales = "free_x", 
                   labeller = ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups", 
                                            Light = "Light", Heavy = "Heavy"))) + 
        ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                                                        strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Calculating Data Completeness")
    {
      CompletenessData <- rbind(Proteopedia::Count_Proteins(ProteinData, "All"), Proteopedia::Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                                            Proteopedia::Count_Proteins(ProteinData[N_precursors_proteotypic >= 2], "≥ 2 Proteotypic"))
      
      NAsPlot <- CompletenessData |> ggplot2::ggplot(ggplot2::aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors))+
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("All" = "black", "≥ 2" = "darkgrey", "≥ 2 Proteotypic" = "orange3")) +
        ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") + ggplot2::scale_x_continuous( breaks = seq(1, 1000, 1)) +
        ggplot2::scale_y_continuous( limits = c(0, max(CompletenessData$cumulative_protein_N)/1000))+
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
              legend.position = "inside", legend.position.inside = c(0.25, 0.25))
    }
    message("Calculating Channel Skewing")
    {
      SkewData <- data.table::data.table(meta)
      SkewData[, `:=`(PearsonsSkewRatio = 0, PearsonsSkew_L = 0, PearsonsSkew_H = 0, Median = 0, Mean = 0)]
      for(RowIndex in 1:nrow(SkewData)){
        SampleID <- SkewData$Sample[RowIndex]
        SkewData$PearsonsSkewRatio[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ_Ratio])
        SkewData$PearsonsSkew_L[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ_L])
        SkewData$PearsonsSkew_H[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ_H])
        SkewData$Median[RowIndex] <- median(ProteinData[Sample == SampleID, LFQ_Ratio], na.rm = T)
        SkewData$Mean[RowIndex] <- mean(ProteinData[Sample == SampleID, LFQ_Ratio], na.rm = T)
      }
  
      SkewPlot <- ProteinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = LFQ_Ratio, colour = Condition)) + 
        ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
        ggplot2::ylab("Heavy:Light Protein LFQ Ratio") + ggplot2::coord_flip() + 
        ggplot2::geom_text(data = SkewData, aes(y = Median*2, label = paste0("Pearson's Skew\n", round(PearsonsSkewRatio, 2))), size = 6) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Plotting Missed Trypsinisation Sites")
    {
      TrypsinData <- PrecursorData |> data.table::copy()
      TrypsinData[, MissedTrypsin := grepl("[RK][^P]", Stripped.Sequence)]
      TrypsinData[, N_Trypsin := .N, by = .(Sample, MissedTrypsin, Label)]
      TrypsinData[, N_Sample := .N, by = .(Sample, Label)]
      TrypsinData <- TrypsinData[MissedTrypsin == T, .(Sample, Condition, Replicate, Label, N_Trypsin, N_Sample)] |> dplyr::distinct()
      TrypsinData[, PercentTrypsin := (N_Trypsin/N_Sample)*100]
      
      MissedCleavagePlot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentTrypsin, fill = Condition, alpha = Label)) + 
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
        ggplot2::labs(y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                       strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Calculating Precursor & Protein Variation")
    {
      PrecursorCVs <- PrecursorData[, .(CV = Proteopedia::Calculate_CV(Precursor.Translated), N = .N), .(Precursor.Id, Condition, Label)]
      PrecursorCVs <- PrecursorCVs[N >= 3]  
      PrecursorCVs[, rank := data.table::frank(CV), .(Condition, Label)] 
      PrecursorCVs[, ID := "Precursors"]
      
      ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_L), N = .N), .(ProteinGroup, Condition)][, Label := "Light"] |> 
                      rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_H), N = .N), .(ProteinGroup, Condition)][, Label := "Heavy"])
      ProteinCVs <- ProteinCVs[N >= 3]
      ProteinCVs <- ProteinCVs[, rank := data.table::frank(CV), .(Condition, Label)] 
      ProteinCVs[, ID := "Protein Groups"]
      
      AllCVs <- PrecursorCVs[, Precursor.Id := NULL] |> rbind(ProteinCVs[, ProteinGroup := NULL])
      
      VariationPlot <- AllCVs |> ggplot2::ggplot(ggplot2::aes(x = rank/1000, y = CV, colour = gsub("_", " ", Condition), alpha = Label)) +
        ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Variation [%]") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, name = "Condition") +
        ggplot2::coord_cartesian(ylim = c(0,50)) + ggplot2::facet_wrap(~ID, scales = "free") + 
        ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
              legend.position = "inside", legend.position.inside = c(0.7, 0.7), strip.background = ggplot2::element_blank(), 
              strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Exporting QC Plots")
    {
      pdf("PrecursorQC_DIANN_Plot.pdf", width = 18, height = 20)
        print(IntensityPlot + CountPlot + free(NAsPlot) + free(IsotopeIncorporation) + SkewPlot + MissedCleavagePlot + free(VariationPlot) + 
                plot_layout(design = "AAAAAA\nBBBBCC\nDEEEFF\nGGGGGG") + plot_annotation(tag_levels = "A"))
      Proteopedia::Reset_Dev()
    }
    Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_StaticSILAC_PrecursorBiochemistry <- function(InputDirectory){
  start.time <- Sys.time()
  set.seed(123)
  message("Loading Precursor Data")
  {
    setwd(InputDirectory)
    PrecursorData <- data.table::fread(list.files(pattern = ".*PrecursorData.csv.*"))[, .(Run, ProteinGroup, ProteinDescription, 
                                                                                          Gene, Stripped.Sequence, Precursor.Translated,
                                                                                          Sample, Condition, Cell, Drug, Time, Replicate,
                                                                                          Precursor.Length, Label)] |> 
      setnames(c("Stripped.Sequence", "Precursor.Length"), c("Sequence", "Length"))
  }
  message("Annotating with Biochemical Measures")
  {
    PrecursorData[, `:=`(Aliphatic_Score = Peptides::aIndex(Sequence), Boman_Interaction_Score = Peptides::boman(Sequence),
                         Hydrophobicity_Score = Peptides::hydrophobicity(Sequence, scale = "KyteDoolittle"),
                         Instability_Score = Peptides::instaIndex(Sequence), MW = Peptides::mw(Sequence), 
                         pI = Peptides::pI(Sequence, pKscale = "Dawson"))]
    
    for(ColIndex in which(colnames(PrecursorData) == "Length"):ncol(PrecursorData)){
      if(is.numeric(PrecursorData[, get(colnames(PrecursorData)[ColIndex])])){
        message(paste0("Analysing ", colnames(PrecursorData)[ColIndex]), " Trend")
        SubsetData = PrecursorData[, .(Sequence, Length, Condition, Replicate, Label, Precursor.Translated, get(colnames(PrecursorData)[ColIndex]))]
        SubsetData |> data.table::setnames("V7", "Subset")
        
        pval_L <- summary(stats::lm(Precursor.Translated ~ Subset*Condition, data = SubsetData[Label == "Light"]))$coefficients[4,4]
        pval_H <- summary(stats::lm(Precursor.Translated ~ Subset*Condition, data = SubsetData[Label == "Heavy"]))$coefficients[4,4]
        
        TrendPlot <- SubsetData[Label %in% c("Light", "Heavy")] |> 
          ggplot2::ggplot(ggplot2::aes(x = Subset, y = Precursor.Translated, colour = Condition, alpha = Label)) + 
          ggplot2::geom_smooth(method = "lm") + scale_alpha_manual(values = c("Heavy" = 0.1, "Light" = 0.3)) +
          ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Precursor.Translated], na.rm = T)*0.93, 
                            label = paste0("Light P-Value: ", ifelse(pval_L < 0.01, formatC(pval_L, format = "e", digits = 2), round(pval_L, digits = 2)), "\n",
                                           "Heavy P-Value: ", ifelse(pval_H < 0.01, formatC(pval_H, format = "e", digits = 2), round(pval_H, digits = 2))), size = 6) +
          ggplot2::labs(x = paste0("Precursor ", gsub("_", " ", colnames(PrecursorData)[ColIndex])), y = "Normalised Precursor Intensity") + 
          ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
        
        pdf(paste0(colnames(PrecursorData)[ColIndex],"_Trend.pdf"), width = 12, height = 10)
        print(TrendPlot)
        print(TrendPlot + ggplot2::labs(x = "", y = ""))
        Proteopedia::Reset_Dev()
      }
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Analyse_StaticSILAC_Proteins <- function(InputDirectory, ExpGroup, CtlGroup, MinPrecursors = 2, ImputationQ = 0.01, ImputationSigma = 1){
    set.seed(123)
    start.time <- Sys.time()
    message("Loading & Formatting Data")
    {
      setwd(InputDirectory)
      if(length(list.files(pattern = "SILAC_DIANN_Output.csv")) > 0){
        InputFile <- list.files(pattern = "SILAC_DIANN_Output.csv")[1]
        meta <- data.table::fread(InputFile)[grepl(ExpGroup, Condition)|grepl(CtlGroup, Condition), .(Sample, Condition, Replicate)] |> unique() |> dplyr::arrange(!grepl(CtlGroup, Condition), Condition, Replicate)
        
        SpectraRead <- data.table::fread(InputFile)[, Log2LFQ_H := log2(LFQ_H)][, Log2LFQ_L := log2(LFQ_L)][, Log2LFQ := log2(LFQ)][Sample %in% meta$Sample]
        ProteinInfo <- data.table::fread(InputFile)[,.(ProteinGroup, ProteinDescription, Gene)] |> dplyr::distinct()
        RatioData <- data.table::fread(InputFile)[Sample %in% meta$Sample]
      } else {return(message("ERROR: No InputFile Found"))}
      
      design <- stats::model.matrix(~0 + Condition, data = meta)
      colnames(design) <- gsub("Condition", "", colnames(design))
      rownames(design) <- meta$Sample
      design <- design[, (c(which(grepl(CtlGroup, colnames(design))), which(grepl(ExpGroup, colnames(design)))))]
      contr.matrix <- matrix(nrow = 2, ncol = 1)
      contr.matrix[, 1] <- c(-1, 1)
      dimnames(contr.matrix) <- list(Levels = colnames(design), Contrasts = "comp")
      
      data.table::fwrite(meta, file = "Sample_Metadata.csv")
    }
    message("Analysing Both Channels As Label-Free Proxy")
    {
      setwd(InputDirectory)
      if(dir.exists("Total_Analysis") == TRUE){
        unlink("Total_Analysis", recursive = TRUE)
      }
      dir.create("Total_Analysis", showWarnings = TRUE)
      setwd("Total_Analysis")
      SpectraReadT <- SpectraRead |> data.table::copy()
      message("Label-Free Proxy: Performing PCA")
      {
        PCAData <- stats::prcomp(t(data.frame(tidyr::drop_na(data.table::dcast(SpectraReadT, ProteinGroup ~ Sample, value.var = "Log2LFQ", values_fill = NA)), row.names = "ProteinGroup")), scale. = TRUE)
        PCASummary <- summary(PCAData)$importance
        PCAData <- data.table::merge.data.table(data.table::data.table(PCAData$x, keep.rownames = "Sample"), meta)
        PCAData[, Replicate := paste0("R", Replicate)]
        PCAPlot <- (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate)) + 
                      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
                      ggplot2::labs(x = paste("PC1 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC1"] * 100, 0), "%]", sep = ""), 
                                    y = paste("PC2 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC2"] * 100, 0), "%]", sep = "")) + 
                      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank(), 
                                     legend.position = "none")) +
          (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate)) + 
             ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
             ggplot2::labs(x = paste("PC3 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC3"] * 100, 0), "%]", sep = ""), 
                           y = paste("PC4 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC4"] * 100, 0), "%]", sep = "")) + 
             ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
      }
      message("Label-Free Proxy: Filtering Proteins")
      {
        FilteringData <- SpectraReadT[, .(N_Samples = .N, Min_Precursors = min(N_precursors), N_conditions = data.table::uniqueN(Condition)), ProteinGroup]
        Retained1 <- FilteringData[N_Samples == data.table::uniqueN(SpectraReadT$Sample) & Min_Precursors >= MinPrecursors]
        Retained2 <- FilteringData[N_Samples == data.table::uniqueN(SpectraReadT$Sample) - 1 & Min_Precursors >= (MinPrecursors + 1)]
        Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraReadT$Sample)/2) & Min_Precursors >= (MinPrecursors + 1) & N_conditions >= 1]
        RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
        FilteredProteins <- SpectraReadT[ProteinGroup %!in% RetainedProteins]
        SpectraReadT <- SpectraReadT[ProteinGroup %in% RetainedProteins]
        data.table::fwrite(FilteredProteins, file = "Filtered_Proteins.csv")
        data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")
        FilteringPlotDataA <- SpectraReadT[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")]
        FilteringPlotDataB <- FilteredProteins[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")]
        FilteringPlotData <- dplyr::mutate(dplyr::summarise(dplyr::group_by(rbind(FilteringPlotDataA, FilteringPlotDataB), Sample, Condition, 
                                                                            Replicate, Inclusion), n = dplyr::n()), 
                                           Colour = ifelse(grepl(CtlGroup, Condition), RColorBrewer::brewer.pal(3, "Set1")[2], RColorBrewer::brewer.pal(3, "Set1")[1]))
        
        CountsBar <- ggplot2::ggplot(FilteringPlotData, ggplot2::aes(x = gsub("_", " ", Sample), y = n, fill = Colour, alpha = Inclusion)) + 
          ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_identity(guide = "none") + 
          ggplot2::scale_alpha_manual(values = c(Excluded = 0.4, Retained = 1), guide = "none") + 
          ggplot2::geom_text(ggplot2::aes(label = n), size = 6, colour = ifelse(FilteringPlotData$Inclusion == "Retained", "white", "black"), 
                                 position = ggplot2::position_stack(), vjust = ifelse(FilteringPlotData$Inclusion == "Retained", 1.5, -0.5)) + 
          ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x") + 
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) + 
          ggplot2::labs(x = NULL, y = "No. Proteins", fill = NULL) + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = ggplot2::element_blank(), 
                         strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0, "line"))
        
        UpsetPlot <- ggplot2::ggplot(SpectraReadT[,.(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup], 
                                     ggplot2::aes(x = Sample)) + ggplot2::geom_bar() + 
          ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) + 
          ggupset::scale_x_upset(order_by = "degree", reverse = TRUE, sets = gsub("_", " ", SpectraReadT[order(Condition), unique(Sample)])) + 
          ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) + 
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))
      }
      message("Label-Free Proxy: Performing Median Normalisation")
      {
        SpectraReadT[, `:=`(Log2LFQ_Norm, Log2LFQ - median(Log2LFQ, na.rm = TRUE) + median(SpectraReadT$Log2LFQ, na.rm = TRUE)), by = Sample]
        NormalisationPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraReadT, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")), 
                                             ggplot2::aes(x = Condition, colour = Condition,  y = value, group = Sample)) + 
          ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c(Log2LFQ = "Pre", Log2LFQ_Norm = "Post"))) + 
          ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
          ggplot2::labs(x = NULL, y = expression("Log"[2]~ "Protein LFQ")) + ggplot2::theme(strip.background = ggplot2::element_blank())
      }
      message("Label-Free Proxy: Imputing NA Values")
      {
        SpectraReadT <- data.table::data.table(tidyr::pivot_wider(SpectraReadT[, .(ProteinGroup, Gene, Sample, Log2LFQ_Norm)], 
                                                                  id_cols = ProteinGroup, values_from = Log2LFQ_Norm, 
                                                                  names_from = Sample, values_fill = NA))
        SpectraImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraReadT, rownames = "ProteinGroup"), q = ImputationQ, 
                                                                        tune.sigma = ImputationSigma), keep.rownames = "ProteinGroup")
        SpectraAll <- data.table::merge.data.table(data.table::melt.data.table(SpectraReadT, id.vars = "ProteinGroup", value.name = "Measured_LFQ", 
                                                                                          variable.name = "Sample"), 
                                                   data.table::melt.data.table(SpectraImp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ", 
                                                                                          variable.name = "Sample"))
        SpectraAll[, `:=`(Data, data.table::fifelse(is.na(Measured_LFQ), "Imputed", "Measured"))]
        
        ImputedDensity <- ggplot2::ggplot(SpectraAll, ggplot2::aes(x = Imputed_LFQ, fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + 
          ggplot2::scale_fill_manual(values = c(Measured = "black", Imputed = "magenta3")) + 
          ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) + 
          ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))
        
        ImputedNAs <- data.table::melt.data.table(SpectraReadT[ProteinGroup %in% Retained3$ProteinGroup], 
                                                             id.vars = "ProteinGroup", variable.name = "Sample")[is.na(value)]
        ImputedNAs[, `:=`(value, NULL)]
        ImputedNAs <- dplyr::summarise(dplyr::group_by(ImputedNAs, ProteinGroup), vector = paste(Sample, collapse = ", "))
        data.table::fwrite(ImputedNAs, "Imputed_LFQs.csv")
        SpectraAll[, `:=`(Log2LFQ, data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ))]
        SpectraAll <- SpectraAll[, .(ProteinGroup, Sample, Log2LFQ)]
      }
      message("Label-Free Proxy: Performing Paired T-Testing")
      {
        SpectraAll <- data.table::merge.data.table(SpectraAll, meta[, .(Condition, Sample)], by = "Sample")[, `:=`(LFQ, 2^(Log2LFQ))]
        SpectraTtest <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(SpectraAll, Condition, ProteinGroup), N = dplyr::n(), 
                                                            Log2MeanLFQ = log2(mean(LFQ)), CV = (sd(LFQ)/mean(LFQ)) * 100), 
                                           id_cols = ProteinGroup, names_from = Condition, values_from = c(N, Log2MeanLFQ, CV))
        ctl_col_index <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(CtlGroup, colnames(SpectraTtest)))
        exp_col_index <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(ExpGroup, colnames(SpectraTtest)))
        SpectraTtest$Log2FC <- SpectraTtest[exp_col_index] - SpectraTtest[ctl_col_index]
        Ttest_Output <- SpectraAll[, `:=`(P.Value, stats::t.test(Log2LFQ ~ Condition)$p.value), by = ProteinGroup]
        SpectraTtest <- data.table::merge.data.table(data.table::data.table(SpectraTtest), Ttest_Output, by = "ProteinGroup")
        SpectraTtest[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        SpectraTtest <- data.table::data.table(data.table::merge.data.table(SpectraTtest, ProteinInfo))
        SpectraTtest[, `:=`(Log2FC, as.numeric(Log2FC))]
        data.table::fwrite(dplyr::distinct(SpectraTtest), "Paired_T-Test_Output.csv")
      }
      message("Label-Free Proxy: Fitting to Linear Model & Exporting Data")
      {
        LimmaInput <- data.table::dcast(SpectraAll, formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
          data.table::setcolorder(neworder = c("ProteinGroup", meta$Sample))
        
        efit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, design), contr.matrix))
        if (!is.finite(efit$df.prior)) {
          message("Warning: Limma df.prior is Infinite")
        }
        MeanVarData <- data.table::data.table(efit$genes, Mean = efit$Amean, Variance = sqrt(efit$sigma))
        MeanVarData[, `:=`(Data, data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured"))]
        
        MeanVarPlot <- ggplot2::ggplot(MeanVarData, ggplot2::aes(x = Mean, y = Variance, colour = Data)) + 
          ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_manual(values = c(Measured = "black", Imputed = "magenta3"), guide = "none") + 
          ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "magenta3", stroke = NA) +
          ggplot2::labs(x = expression("Mean Log"[2]~"Protein LFQ"), y = "Variance", colour = NULL)
        
        LimmaOutput <- data.table::setnames(data.table::data.table(limma::topTable(efit, coef = 1, adjust.method = "BH", n = Inf)), "logFC", "Log2FC")
        LimmaOutput <- data.table::merge.data.table(LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = TRUE), ], ProteinInfo, all.x = T)
        LimmaOutput[, `:=`(Significance, data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", 
                                                             data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None")))]
        LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = TRUE)
        LimmaOutput$Isoforms <- 1
        for(i in 1:nrow(LimmaOutput)){
          if(length(stringr::str_extract_all(LimmaOutput$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
            LimmaOutput$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaOutput$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          } else {LimmaOutput$Isoforms[i] <- 1}
        }
        LimmaOutput[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", GeneGroup))]
        data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
      }
      message("Label-Free Proxy: Generating Volcano Plots")
      {
        MeanLog2FC <- round(mean(LimmaOutput$Log2FC, na.rm = T), digits = 3)
        LimmaOutput <- dplyr::arrange(LimmaOutput, desc(abs(t)))
        LimmaVolcano <- ggplot2::ggplot(LimmaOutput, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("black") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput$Log2FC) * 0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        ImputedLimmaData <- dplyr::mutate(LimmaOutput, Imputed = ifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes", "No"))
        ImputedLimmaVolcano <- ggplot2::ggplot(ImputedLimmaData, ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          ggplot2::scale_colour_manual(values = c(No = "black", Yes = "magenta3"), labels = c(No = "Measured", Yes = "Imputed"), guide = "none") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput$Log2FC), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Imputed Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        pdf("Limma_Volcanoes.pdf", height = 10, width = 14)
        print(LimmaVolcano)
        print(ImputedLimmaVolcano)
        Proteopedia::Reset_Dev()
        
        SpectraTtestVolcanoData <- SpectraTtest[, .(Gene, ProteinGroup, Log2FC, P.Value)] |> dplyr::arrange(desc(abs(Log2FC*-log10(P.Value)))) |> dplyr::distinct()
        MeanLog2FC <- round(mean(SpectraTtestVolcanoData$Log2FC), digits = 3)
        
        TtestVolcano <- ggplot2::ggplot(SpectraTtestVolcanoData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("black") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(SpectraTtest$Log2FC) * 0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        ImputedTtestVolcano <- ggplot2::ggplot(SpectraTtestVolcanoData[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes", "No")], 
                                               ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          ggplot2::scale_colour_manual(values = c(No = "black", Yes = "magenta3"), labels = c(No = "Measured", Yes = "Imputed"), guide = "none") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(SpectraTtestVolcanoData$Log2FC), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + Proteopedia::Add_AbundanceAxes() + 
          ggplot2::ggtitle("Imputed Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        pdf("PairedTtest_Volcanoes.pdf", height = 10, width = 14)
        print(TtestVolcano)
        print(ImputedTtestVolcano)
        Proteopedia::Reset_Dev()
      }
      message("Label-Free Proxy: Exporting QC Plot")
      {
        pdf("ProteinQC_Limma_Plot.pdf", width = 18, height = 20)
        print(PCAPlot + patchwork::free(CountsBar) + patchwork::free(UpsetPlot) + patchwork::free(NormalisationPlot) + patchwork::free(ImputedDensity) + MeanVarPlot + 
                patchwork::plot_layout(design = "ABCC\nDDEE\nFFGG") + patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F"))))
        Proteopedia::Reset_Dev()
      }
    }
    message("Analysing Each Channel For Stability/Synthesis Measures")
    {
      setwd(InputDirectory)
      if(dir.exists("Channel_Analysis") == TRUE){
        unlink("Channel_Analysis", recursive = TRUE)
      }
      dir.create("Channel_Analysis", showWarnings = TRUE)
      setwd("Channel_Analysis")
      SpectraRead_L <- SpectraRead |> data.table::copy()
      SpectraRead_H <- SpectraRead |> data.table::copy()
      message("SILAC Channels: Performing PCA")
      {
        PCAData_L <- SpectraRead_L[, .(ProteinGroup, Gene, Sample, Log2LFQ_L)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_L") |> 
          tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
        PCASummary_L <- summary(PCAData_L)$importance
        PCAData_L <- data.table::data.table(PCAData_L$x, keep.rownames = "Sample") |> data.table::merge.data.table(meta)
        PCAData_L[, Dataset := "Light"]
        
        PCAData_H <- SpectraRead_H[, .(ProteinGroup, Gene, Sample, Log2LFQ_H)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_H") |> 
          tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
        PCASummary_H <- summary(PCAData_H)$importance
        PCAData_H <- data.table::data.table(PCAData_H$x, keep.rownames = "Sample") |> data.table::merge.data.table(meta)
        PCAData_H[, Dataset := "Heavy"]
        
        PCAData <- PCAData_H |> rbind(PCAData_L)
        PCAData[, Replicate := paste0("R", Replicate)]
        PCASummary <- PCASummary_H |> rbind(PCASummary_L)
        
        PCAPlot <- (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, alpha = Dataset)) + 
                      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_alpha_manual(values = c("Heavy" = 1, "Light" = 0.5)) +
                      ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
                      ggplot2::labs(x = paste("PC1 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC1"] * 100, 0), "%]", sep = ""), 
                                    y = paste("PC2 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC2"] * 100, 0), "%]", sep = "")) + 
                      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
                                     legend.title = ggplot2::element_blank(), legend.position = "none")) + 
          (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, alpha = Dataset)) + 
             ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_alpha_manual(values = c("Heavy" = 1, "Light" = 0.5)) +
             ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
             ggplot2::labs(x = paste("PC3 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC3"] * 100, 0), "%]", sep = ""), 
                           y = paste("PC4 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC4"] * 100, 0), "%]", sep = "")) + 
             ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
                            legend.title = ggplot2::element_blank()))
      }
      message("SILAC Channels: Filtering Proteins")
      {
        FilteringData <- SpectraRead |> dplyr::group_by(ProteinGroup) |> dplyr::summarise(N_Samples = dplyr::n(), 
                                                                                            Min_Precursors_L = min(N_precursors_L), 
                                                                                            Min_Precursors_H = min(N_precursors_H), 
                                                                                            N_conditions = uniqueN(Condition)) |> data.table::data.table()
        
        Retained1_L <- FilteringData[N_Samples == uniqueN(SpectraRead_L$Sample) & Min_Precursors_L >= MinPrecursors]
        Retained2_L <- FilteringData[N_Samples == uniqueN(SpectraRead_L$Sample)-1 & Min_Precursors_L >= (MinPrecursors+1)]
        Retained3_L <- FilteringData[N_Samples == floor(uniqueN(SpectraRead_L$Sample)/2) & Min_Precursors_L >= (MinPrecursors+1) & N_conditions == 1]
        RetainedProteins_L <- c(Retained1_L$ProteinGroup, Retained2_L$ProteinGroup, Retained3_L$ProteinGroup)
        FilteredProteins_L <- SpectraRead_L[ProteinGroup %!in% RetainedProteins_L]
        SpectraRead_L <- SpectraRead_L[ProteinGroup %in% RetainedProteins_L] |> dplyr::select(!ends_with("_H"))
        data.table::fwrite(FilteredProteins_L, file = "Light_Filtered_Proteins.csv")
        data.table::fwrite(Retained3_L, file = "Light_Imputed_Proteins.csv")
        
        Retained1_H <- FilteringData[N_Samples == uniqueN(SpectraRead_H$Sample) & Min_Precursors_H >= MinPrecursors]
        Retained2_H <- FilteringData[N_Samples == uniqueN(SpectraRead_H$Sample)-1 & Min_Precursors_H >= (MinPrecursors+1)]
        Retained3_H <- FilteringData[N_Samples == floor(uniqueN(SpectraRead_H$Sample)/2) & Min_Precursors_H >= (MinPrecursors+1) & N_conditions == 1]
        RetainedProteins_H <- c(Retained1_H$ProteinGroup, Retained2_H$ProteinGroup, Retained3_H$ProteinGroup)
        FilteredProteins_H <- SpectraRead_H[ProteinGroup %!in% RetainedProteins_H]
        SpectraRead_H <- SpectraRead_H[ProteinGroup %in% RetainedProteins_H] |> dplyr::select(!ends_with("_L"))
        data.table::fwrite(FilteredProteins_H, file = "Heavy_Filtered_Proteins.csv")
        data.table::fwrite(Retained3_H, file = "Heavy_Imputed_Proteins.csv")
        
        #Create Upset of Retained Proteins
        FilteringPlotData_L <- SpectraRead_L[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")] |> 
          rbind(FilteredProteins_L[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")]) |> 
          dplyr::group_by(Sample, Condition, Replicate, Inclusion) |> dplyr::summarise(N = dplyr::n()) |> data.table::data.table()
        FilteringPlotData_L[, Dataset := "Light"]
        
        FilteringPlotData_H <- SpectraRead_L[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")] |> 
          rbind(FilteredProteins_H[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")]) |> 
          dplyr::group_by(Sample, Condition, Replicate, Inclusion) |> dplyr::summarise(N = dplyr::n()) |> data.table::data.table()
        FilteringPlotData_H[, Dataset := "Heavy"]
        
        FilteringPlotData <- FilteringPlotData_L |> rbind(FilteringPlotData_H)
        
        CountsBar <- FilteringPlotData |> ggplot2::ggplot(ggplot2::aes(x = gsub("_", " ", Sample), y = N, fill = Condition, alpha = Inclusion)) + 
          ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
          ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") + 
          ggplot2::geom_text(ggplot2::aes(label = N), size = 6, colour = ifelse(FilteringPlotData$Inclusion == "Retained", "white","black"), 
                             position = ggplot2::position_stack(), vjust = ifelse(FilteringPlotData$Inclusion == "Retained", 1.5, -0.5)) +                                  
          ggplot2::facet_wrap(~Dataset, scales = "free_x") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) + 
          ggplot2::labs(x = NULL, y = "No. Proteins") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = ggplot2::element_blank(), 
                                                                       strip.text.y = ggplot2::element_text(size = 26))
        
        UpsetPlot <- (SpectraRead_H[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
                        ggplot2::geom_text(stat="count", ggplot2::aes(label =ggplot2::after_stat(count)), vjust = -0.5, size = 3) + ggplot2::ggtitle("Light") +
                        ggupset::scale_x_upset(order_by = "degree", reverse = TRUE, sets = gsub("_", " ", SpectraRead_H[order(Condition), unique(Sample)])) + 
                        ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) +
                        ggplot2::theme(plot.title = ggplot2::element_text(size = 24))) + 
          (SpectraRead_L[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
             ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) + ggplot2::ggtitle("Heavy") +
             ggupset::scale_x_upset(order_by = "degree", reverse = TRUE, sets = gsub("_", " ", SpectraRead_L[order(Condition), unique(Sample)])) + 
             ggplot2::labs(x = NULL, y = NULL) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) +
             ggplot2::theme(plot.title = ggplot2::element_text(size = 24), axis.text.x = ggplot2::element_blank())) + patchwork::plot_layout(guides = "collect")
      }
      message("SILAC Channels: Performing Median Normalisation")
      {
        SpectraRead_L[, Log2LFQ_Norm_L := Log2LFQ_L - median(Log2LFQ_L, na.rm = TRUE) + median(SpectraRead_L$Log2LFQ_L, na.rm = TRUE), by = Sample]  
        colnames(SpectraRead_L) <- gsub("_L", "", colnames(SpectraRead_L))
        SpectraRead_L[ , Dataset := "Light"]
        
        SpectraRead_H[, Log2LFQ_Norm_H := Log2LFQ_H - median(Log2LFQ_H, na.rm = TRUE) + median(SpectraRead_H$Log2LFQ_H, na.rm = TRUE), by = Sample]  
        colnames(SpectraRead_H) <- gsub("_H", "", colnames(SpectraRead_H))
        SpectraRead_H[ , Dataset := "Heavy"]
        
        NormalisationPlot <- ggplot2::ggplot(melt(SpectraRead_L |> rbind(SpectraRead_H), measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                             ggplot2::aes(x = Condition, colour = Condition, y = value, group = Sample))+
          ggplot2::facet_grid(c("Dataset","variable"), labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
          ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
          ggplot2::labs(x = NULL, y = expression("Log"[2]~"Protein LFQ")) + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
      }
      message("SILAC Channels: Imputing NA Values")
      {
        SpectraRead_L <- SpectraRead_L[, .(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
        SpectraImpL <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead_L, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma), 
                                              keep.rownames = "ProteinGroup")
        SpectraAllL <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead_L, id.vars = "ProteinGroup", value.name = "MeasuredLFQ", variable.name = "Sample"),
                                         data.table::melt.data.table(SpectraImpL, id.vars = "ProteinGroup", value.name = "ImputedLFQ", variable.name = "Sample"))
        SpectraAllL[, Data := data.table::fifelse(is.na(MeasuredLFQ), "Imputed","Measured")]
        SpectraAllL[, Dataset := "Light"]
        
        SpectraRead_H <- SpectraRead_H[, .(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
        SpectraImpH <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead_L, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma), 
                                              keep.rownames = "ProteinGroup")
        SpectraAllH <- data.table::merge.data.table(melt(SpectraRead_H, id.vars = "ProteinGroup", value.name = "MeasuredLFQ", variable.name = "Sample"),
                                         data.table::melt.data.table(SpectraImpH, id.vars = "ProteinGroup", value.name = "ImputedLFQ", variable.name = "Sample"))
        SpectraAllH[, Data := data.table::fifelse(is.na(MeasuredLFQ), "Imputed","Measured")]
        SpectraAllH[, Dataset := "Heavy"]
        
        SpectraAll <- SpectraAllL |> rbind(SpectraAllH)
        
        ImputedDensity <- SpectraAll |> ggplot2::ggplot(ggplot2::aes(x = ImputedLFQ, fill = Data)) +
          ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_fill_manual(values = c("Measured" = "black", "Imputed" = "magenta3")) +
          ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) + ggplot2::facet_wrap(~Dataset, nrow = 2, strip.position =  "right") + 
          ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8), strip.background = ggplot2::element_blank(), 
                         strip.text.y = ggplot2::element_text(size = 26))
        
        ImputedNAs_L <- data.table::melt.data.table(SpectraRead_L[ProteinGroup %in% Retained3_L$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
        ImputedNAs_L <- ImputedNAs_L |> dplyr::filter(is.na(value)) |> dplyr::select(!value) |> dplyr::group_by(ProteinGroup) |> dplyr::summarise(vector = paste(Sample, collapse=", "))
        ImputedNAs_L$Dataset <- "Light"
        ImputedNAs_H <- data.table::melt.data.table(SpectraRead_H[ProteinGroup %in% Retained3_H$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
        ImputedNAs_H <- ImputedNAs_H |> dplyr::filter(is.na(value)) |> dplyr::select(!value) |> dplyr::group_by(ProteinGroup) |> dplyr::summarise(vector = paste(Sample, collapse=", "))
        ImputedNAs_H$Dataset <- "Heavy"
        data.table::fwrite(ImputedNAs_L |> rbind(ImputedNAs_H), "Imputed_LFQs.csv")
        
        SpectraAll[, Log2LFQ := data.table::fifelse(is.na(MeasuredLFQ), ImputedLFQ, MeasuredLFQ)]
      }
      message("SILAC Channels: Performing Paired T-Testing")
      {
        SpectraAll[, Sample := gsub("Log2LFQ_._(.*)", "\\1", Sample)]
        SpectraAll[, Condition := gsub("(.*)_R\\d", "\\1", Sample)]
        SpectraAll[, LFQ := 2^(Log2LFQ)]
        # Light T-Test
        {
          SpectraTtest_L <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(SpectraAll[Dataset == "Light"], Condition, ProteinGroup), 
                                                                N = dplyr::n(), Log2MeanLFQ = log2(mean(LFQ)), 
                                                                CV = (sd(LFQ)/mean(LFQ)) * 100), 
                                               id_cols = ProteinGroup, names_from = Condition, values_from = c(N, Log2MeanLFQ, CV))
          CtlIndex <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_L)) & grepl(CtlGroup, colnames(SpectraTtest_L)))
          ExpIndex <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_L)) & grepl(ExpGroup, colnames(SpectraTtest_L)))
          SpectraTtest_L$Log2FC <- SpectraTtest_L[ExpIndex] - SpectraTtest_L[CtlIndex]
          Ttest_Output <- SpectraAll[Dataset == "Light", `:=`(P.Value, stats::t.test(Log2LFQ ~ Condition)$p.value), by = ProteinGroup]
          SpectraTtest_L <- data.table::merge.data.table(data.table(SpectraTtest_L), Ttest_Output, by = "ProteinGroup")
          SpectraTtest_L[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
          SpectraTtest_L <- data.table::data.table(data.table::merge.data.table(SpectraTtest_L, ProteinInfo))
          SpectraTtest_L[, `:=`(Log2FC, as.numeric(Log2FC))]
        }
        # Heavy T-Test
        {
          SpectraTtest_H <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(SpectraAll[Dataset == "Heavy"], Condition, ProteinGroup), 
                                                                N = dplyr::n(), Log2MeanLFQ = log2(mean(LFQ)), 
                                                                CV = (sd(LFQ)/mean(LFQ)) * 100), 
                                               id_cols = ProteinGroup, names_from = Condition, values_from = c(N, Log2MeanLFQ, CV))
          CtlIndex <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_H)) & grepl(CtlGroup, colnames(SpectraTtest_H)))
          ExpIndex <- which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_H)) & grepl(ExpGroup, colnames(SpectraTtest_H)))
          SpectraTtest_H$Log2FC <- SpectraTtest_H[ExpIndex] - SpectraTtest_H[CtlIndex]
          Ttest_Output <- SpectraAll[Dataset == "Heavy", `:=`(P.Value, stats::t.test(Log2LFQ ~ Condition)$p.value), by = ProteinGroup]
          SpectraTtest_H <- data.table::merge.data.table(data.table(SpectraTtest_H), Ttest_Output, by = "ProteinGroup")
          SpectraTtest_H[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
          SpectraTtest_H <- data.table::data.table(data.table::merge.data.table(SpectraTtest_H, ProteinInfo))
          SpectraTtest_H[, `:=`(Log2FC, as.numeric(Log2FC))]
        }
        SpectraTtest <- SpectraTtest_L |> rbind(SpectraTtest_H)
        data.table::fwrite(SpectraTtest, "Paired_T-Test_Output.csv")
      }
      message("SILAC Channels: Fitting to Linear Model & Exporting Data")
      {
        LimmaInput_L <- SpectraAll[Dataset == "Light"] |> copy() |> data.table::dcast(ProteinGroup+Dataset ~ Sample, value.var = "Log2LFQ") |> data.table::setcolorder(c("ProteinGroup", rownames(design)))
        LimmaInput_H <- SpectraAll[Dataset == "Heavy"] |> copy() |> data.table::dcast(ProteinGroup+Dataset ~ Sample, value.var = "Log2LFQ") |> data.table::setcolorder(c("ProteinGroup", rownames(design)))

        efit_L <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput_L[, Dataset := NULL], design), contrasts = contr.matrix))
        efit_H <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput_H[, Dataset := NULL], design), contrasts = contr.matrix))
        
        MeanVarData <- data.table::data.table(efit_L$genes, "Mean" = efit_L$Amean, "Variance" = sqrt(efit_L$sigma), "Dataset" = "Light") |> 
          rbind(data.table::data.table(efit_H$genes, "Mean" = efit_H$Amean, "Variance" = sqrt(efit_H$sigma), "Dataset" = "Heavy"))
        
        MeanVarPlot <- MeanVarData[, Data := data.table::fifelse(Dataset == "Light" & ProteinGroup %in% ImputedNAs_L$ProteinGroup, "Imputed", 
                                                                 data.table::fifelse(Dataset == "Heavy" & ProteinGroup %in% ImputedNAs_H$ProteinGroup, "Imputed", "Measured"))] |> 
          ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) + ggplot2::geom_point(stroke = NA) + 
          ggplot2::scale_colour_manual(values = c(Measured = "black", Imputed = "magenta3"), guide = "none") + 
          ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "magenta3", stroke = NA) +
          ggplot2::labs(x = expression("Mean Log"[2]~"Protein LFQ"), y = "Variance") + ggplot2::facet_wrap(~Dataset, scales = "free_y", nrow = 2, strip.position = "right") +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size = 20), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26), 
                         panel.spacing.x = ggplot2::unit(0,"line"))
        
        LimmaOutput_L <- limma::topTable(efit_L, coef=1, adjust.method = "BH", n=Inf) |> data.table::data.table() |> data.table::setnames("logFC", "Log2FC")
        LimmaOutput_L <- LimmaOutput_L[order(abs(LimmaOutput_L$t), decreasing = T),] |> data.table::merge.data.table(ProteinInfo, sort = F)
        LimmaOutput_L[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
        LimmaOutput_L <- LimmaOutput_L |> data.table::merge.data.table(LimmaInput_L, by = "ProteinGroup", sort = F) |> unique() 
        LimmaOutput_L$Isoforms <- 1
        for(i in 1:nrow(LimmaOutput_L)){
          if (length(stringr::str_extract_all(LimmaOutput_L$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
            LimmaOutput_L$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaOutput_L$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          }else {LimmaOutput_L$Isoforms[i] <- 1}
        }
        LimmaOutput_L[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput_L[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
        LimmaOutput_L[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs_L$ProteinGroup, "Imputed", "Measured")]
        data.table::fwrite(LimmaOutput_L, file = "Light_Limma_Output.csv")
        
        LimmaOutput_H <- limma::topTable(efit_H, coef=1, adjust.method = "BH", n=Inf) |> data.table::data.table() |> data.table::setnames("logFC", "Log2FC")
        LimmaOutput_H <- LimmaOutput_H[order(abs(LimmaOutput_H$t), decreasing = TRUE),] |> data.table::merge.data.table(ProteinInfo, sort = F)
        LimmaOutput_H[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
        LimmaOutput_H <- LimmaOutput_H |> data.table::merge.data.table(LimmaInput_H, by = "ProteinGroup", sort = F) |> unique() 
        LimmaOutput_H$Isoforms <- 1
        for(i in 1:nrow(LimmaOutput_H)){
          if (length(stringr::str_extract_all(LimmaOutput_H$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
            LimmaOutput_H$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaOutput_H$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          }else {LimmaOutput_H$Isoforms[i] <- 1}
        }
        LimmaOutput_H[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput_H[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
        LimmaOutput_H[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs_H$ProteinGroup, "Imputed", "Measured")]
        data.table::fwrite(LimmaOutput_H, file = "Heavy_Limma_Output.csv")
        
        data.table::fwrite(LimmaOutput_L[, .(ProteinGroup, Isoforms, Gene, ProteinDescription, Log2FC, P.Value, AveExpr, t, adj.P.Val, B, Significance, Imputed, GeneGroup)][, Dataset := "Light"] |> 
                             rbind(LimmaOutput_H[, .(ProteinGroup, Isoforms, Gene, ProteinDescription, Log2FC, P.Value, AveExpr, t, adj.P.Val, B, Significance, Imputed, GeneGroup)][, Dataset := "Heavy"]), file = "Limma_Output.csv")
      }
      message("SILAC Channels: Generating Volcano Plots")
      {
        MeanLog2FC_L <- round(mean(LimmaOutput_L[, Log2FC]), digits = 3)
        MeanLog2FC_H <- round(mean(LimmaOutput_H[, Log2FC]), digits = 3)
        
        LimmaVolcano_L <- ggplot2::ggplot(LimmaOutput_L, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC_L, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput_L[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_L), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Light Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        ImputedLimmaVolcano_L <- ggplot2::ggplot(LimmaOutput_L, ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          ggplot2::scale_colour_manual(values = c(Measured = "black", Imputed = "magenta3"), guide = "none") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC_L, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput_L[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_L), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Imputed Light Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        LimmaVolcano_H <- ggplot2::ggplot(LimmaOutput_H, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC_H, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput_H[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_H), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Heavy Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        ImputedLimmaVolcano_H <- ggplot2::ggplot(LimmaOutput_H, ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + 
          ggplot2::scale_colour_manual(values = c(Measured = "black", Imputed = "magenta3"), guide = "none") + 
          Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
          ggplot2::geom_vline(xintercept = MeanLog2FC_H, linetype = "dashed", colour = "black") + 
          ggplot2::annotate("text", x = min(LimmaOutput_H[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_H), size = 5) + 
          Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Imputed Heavy Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())
        
        pdf("Limma_Volcanoes.pdf", height = 10, width = 14)
          print(LimmaVolcano_L)
          print(ImputedLimmaVolcano_L)
          print(LimmaVolcano_H)
          print(ImputedLimmaVolcano_H)
        Proteopedia::Reset_Dev()
        
        # Light Interactive Plot
        InteractiveData_L <- data.table(Protein = LimmaOutput_L[, Gene], 
                                        Log2FC = round(LimmaOutput_L[, Log2FC], digits = 2), 
                                        PValue = LimmaOutput_L[, P.Value], 
                                        Significance = factor(LimmaOutput_L[, Significance]), 
                                        Scien_PValue = formatC(LimmaOutput_L[, P.Value], format = "e", digits = 2), 
                                        GeneGroup = LimmaOutput_L[, GeneGroup])
        InteractiveData_L[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup)))]
        pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData_L, 
                                                                                                                                                                                "scatter", highcharter::hcaes(x = Log2FC, 
                                                                                                                                                                                                              y = -log10(PValue), group = Significance)), 
                                                                                                                                                            zoomType = "xy"), title = list(text = paste0(ExpGroup, " vs ", CtlGroup, " Log2 Fold-Difference")), 
                                                                                                                                      lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                                      tickColor = "black", gridLineWidth = 0), title = list(text = "-Log10 P-Value"), 
                                                                                                                lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                tickColor = "black", gridLineWidth = 0), headerFormat = "", 
                                                                                          pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}"), 
                                                                  scatter = list(marker = list(radius = 3), 
                                                                                 states = list(hover = list(enabled = TRUE), 
                                                                                               inactive = list(enabled = FALSE)), point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))), 
                                      c("#999999", "#880000", "#0033FF"))
        htmlwidgets::saveWidget(pHC, "LightInteractiveVolcanoPlot.html")
        
        # Heavy Interactive Plot
        InteractiveData_H <- data.table(Protein = LimmaOutput_H[, Gene], 
                                        Log2FC = round(LimmaOutput_H[, Log2FC], digits = 2), 
                                        PValue = LimmaOutput_H[, P.Value], 
                                        Significance = factor(LimmaOutput_H[, Significance]), 
                                        Scien_PValue = formatC(LimmaOutput_H[, P.Value], format = "e", digits = 2), 
                                        GeneGroup = LimmaOutput_H[, GeneGroup])
        InteractiveData_H[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup)))]
        pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData_H, 
                                                                                                                                                                                "scatter", highcharter::hcaes(x = Log2FC, 
                                                                                                                                                                                                              y = -log10(PValue), group = Significance)), 
                                                                                                                                                            zoomType = "xy"), title = list(text = paste0(ExpGroup, " vs ", CtlGroup, " Log2 Fold-Difference")), 
                                                                                                                                      lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                                      tickColor = "black", gridLineWidth = 0), title = list(text = "-Log10 P-Value"), 
                                                                                                                lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                tickColor = "black", gridLineWidth = 0), headerFormat = "", 
                                                                                          pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}"), 
                                                                  scatter = list(marker = list(radius = 3), 
                                                                                 states = list(hover = list(enabled = TRUE), 
                                                                                               inactive = list(enabled = FALSE)), point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))), 
                                      c("#999999", "#880000", "#0033FF"))
        htmlwidgets::saveWidget(pHC, "HeavyInteractiveVolcanoPlot.html")
      }
      message("SILAC Channels: Exporting QC Plot")
      {
        pdf("ProteinQC_Limma_Plot.pdf", width = 18, height = 20)
        print(PCAPlot + patchwork::free(CountsBar) + patchwork::free(UpsetPlot) + patchwork::free(NormalisationPlot) + 
                patchwork::free(ImputedDensity) + MeanVarPlot + patchwork::plot_layout(design = "AABBCCCCCCCC\nDDDDDDEEEEEE\nFFFFFFGGGGGG") + 
                patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "", "D", "E", "F"))))
        Proteopedia::Reset_Dev()
      }
    }
    message("Analysing Channel Ratio For Turnover Measures")
    {
      setwd(InputDirectory)
      if(dir.exists("Ratio_Analysis") == TRUE){unlink("Ratio_Analysis", recursive = TRUE)}
      dir.create("Ratio_Analysis", showWarnings = TRUE)
      setwd("Ratio_Analysis")
      message("Channel Ratio: Performing PCA")
      {
        PCAData <- stats::prcomp(t(data.frame(tidyr::drop_na(data.table::dcast(RatioData, ProteinGroup ~ Sample, value.var = "LFQ_Ratio", 
                                                                               values_fill = NA)), row.names = "ProteinGroup")), scale. = TRUE)
        PCASummary <- summary(PCAData)$importance
        PCAData <- data.table::merge.data.table(data.table::data.table(PCAData$x, keep.rownames = "Sample"), meta)
        PCAData[, Replicate := paste0("R", Replicate)]
        
        PCAPlot <- (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate)) + 
                      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
                      ggplot2::labs(x = paste("PC1 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC1"] * 100, 0), "%]", sep = ""), 
                                    y = paste("PC2 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC2"] * 100, 0), "%]", sep = "")) + 
                      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank(), 
                                     legend.position = "none")) + 
          (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate)) + 
             ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
             ggplot2::labs(x = paste("PC3 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC3"] * 100, 0), "%]", sep = ""), 
                           y = paste("PC4 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC4"] * 100, 0), "%]", sep = "")) + 
             ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
      }
      message("Channel Ratio: Filtering Proteins")
      {
        FilteringData <- RatioData[, .(N_Samples = .N, Min_Precursors = min(N_precursors), N_conditions = data.table::uniqueN(Condition)), ProteinGroup]
        Retained1 <- FilteringData[N_Samples == data.table::uniqueN(RatioData$Sample) & Min_Precursors >= MinPrecursors]
        Retained2 <- FilteringData[N_Samples == data.table::uniqueN(RatioData$Sample) - 1 & Min_Precursors >= (MinPrecursors + 1)]
        Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(RatioData$Sample)/2) & Min_Precursors >= (MinPrecursors + 1) & 
                                     N_conditions >= 1]
        RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
        FilteredProteins <- RatioData[ProteinGroup %!in% RetainedProteins]
        RatioData <- RatioData[ProteinGroup %in% RetainedProteins]
        data.table::fwrite(FilteredProteins, file = "Filtered_Proteins.csv")
        data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")
        
        FilteringPlotDataA <- RatioData[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")]
        FilteringPlotDataB <- FilteredProteins[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")]
        FilteringPlotData <- dplyr::mutate(dplyr::summarise(dplyr::group_by(rbind(FilteringPlotDataA, FilteringPlotDataB), 
                                                                            Sample, Condition, Replicate, Inclusion), N = dplyr::n()))
        
        CountsBar <- ggplot2::ggplot(FilteringPlotData, ggplot2::aes(x = gsub("_", " ", Sample), y = N, fill = Condition, alpha = Inclusion)) + 
          ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
          ggplot2::scale_alpha_manual(values = c(Excluded = 0.4, Retained = 1), guide = "none") + 
          ggplot2::geom_text(ggplot2::aes(label = N), size = 6, colour = ifelse(FilteringPlotData$Inclusion == "Retained", "white", "black"), 
                                 position = ggplot2::position_stack(), vjust = ifelse(FilteringPlotData$Inclusion == "Retained", 1.5, -0.5)) + 
          ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x") + 
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) + ggplot2::labs(x = NULL, y = "No. Proteins", fill = NULL) + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = ggplot2::element_blank(), 
                         strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0, "line"))
        
        UpsetPlot <- ggplot2::ggplot(RatioData[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup], ggplot2::aes(x = Sample)) + 
          ggplot2::geom_bar() + ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) + 
          ggupset::scale_x_upset(order_by = "degree", reverse = TRUE, sets = gsub("_", " ", RatioData[order(Condition), unique(Sample)])) + 
          ggplot2::labs(x = NULL, y = "Post-Filtering Count") + 
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))
        
        data.table::fwrite(RatioData, "FilteredRatioData.csv")
      }
      message("Channel Ratio: Plotting Experimental vs Control Ratio Data")
      {
        RatioDataWide <- RatioData[, .(Log2MeanRatio = log2(mean(LFQ_Ratio))), by = list(ProteinGroup, ProteinDescription, Condition)] |>
          data.table::dcast(ProteinGroup+ProteinDescription ~ Condition, value.var = "Log2MeanRatio") |> 
          data.table::merge.data.table(ProteinInfo) |> setnames(c(CtlGroup, ExpGroup), c("CtlGroup", "ExpGroup"))
        RatioDataWide[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        RatioDataWide[, Histone := fifelse(grepl("Histone H", ProteinDescription), T, F)]
        
        RatioDataLong <- RatioDataWide |> data.table::melt.data.table(id.vars = c("ProteinGroup", "ProteinDescription", "Gene"), 
                                                                      measure.vars = c("CtlGroup", "ExpGroup"), value.name = "SILACRatio", 
                                                                      variable.name = "Condition") 
        
        RatioDataLong[, Condition := data.table::fifelse(Condition == "CtlGroup", CtlGroup, ExpGroup)]
        
        DistributionPval <- wilcox.test(RatioDataLong[Condition == CtlGroup, SILACRatio], RatioDataLong[Condition == ExpGroup, SILACRatio])$p.value
        
        DistributionPval <- data.table::fifelse(DistributionPval < 0.001, "p < 0.001 ***", data.table::fifelse(DistributionPval < 0.01, "p < 0.01 **",
                                                data.table::fifelse(DistributionPval < 0.05, "p < 0.05 *", "p > 0.05")))
        
        pdf(paste0(ExpGroup, "_", CtlGroup, "_RatioComp.pdf"), height = 10, width = 14)
        print(RatioDataWide |> ggplot2::ggplot(ggplot2::aes(x = CtlGroup, y = ExpGroup)) + ggplot2::geom_point(stroke = NA) +
                Proteopedia::Add_XYLine("grey") + Proteopedia::Add_Rsq(Subgroups = F) + 
                ggplot2::annotate("text", x = max(RatioDataWide$CtlGroup, na.rm = T)*0.9, y = min(RatioDataWide$ExpGroup, na.rm = T)*0.9, label = "Histones", colour = "#F0F") +
                ggrepel::geom_text_repel(ggplot2::aes(label = Gene)) + ggplot2::geom_point(data = RatioDataWide[Histone ==T], stroke = NA, colour = "#F0F") +
                labs(x = expression("Log"[2]~"Heavy:Light Protein LFQ Ratio"), y = expression("Log"[2]~"Heavy:Light Protein LFQ Ratio"),
                     title = paste0(ExpGroup, " vs. ", CtlGroup)))
        print(RatioDataLong |> ggplot2::ggplot(ggplot2::aes(x = SILACRatio, fill = factor(Condition, levels = c(ExpGroup, CtlGroup)))) + ggplot2::geom_density(alpha = 0.7) + 
                ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, name = "") + ggplot2::labs(x = expression("Log"[2]~"Heavy:Light Protein LFQ Ratio"), y = "Protein Frequency") +
                ggplot2::annotate("segment", x = stats::quantile(RatioDataLong$SILACRatio, 0.01, na.rm = T), xend = stats::quantile(RatioDataLong$SILACRatio, 0.01, na.rm = T),
                                  y = -Inf, yend= Inf, linetype = "dashed") +
                ggplot2::annotate("segment", x = stats::quantile(RatioDataLong$SILACRatio, 0.99, na.rm = T), xend = stats::quantile(RatioDataLong$SILACRatio, 0.99, na.rm = T),
                                  y = -Inf, yend= Inf, linetype = "dashed") + 
                ggplot2::annotate("text", x = median(RatioDataLong$SILACRatio, na.rm = T), y = Proteopedia::Calculate_DensityPeak(RatioDataLong$SILACRatio)*1.1, label = DistributionPval) +
                ggplot2::annotate("label", label = "1%", y = 0, hjust = 0.25, x = stats::quantile(RatioDataLong$SILACRatio, 0.01, na.rm = T), angle = 90) +
                ggplot2::annotate("label", label = "99%", y = 0, hjust = 0.25, x = stats::quantile(RatioDataLong$SILACRatio, 0.99, na.rm = T), angle = 90) +
                ggplot2::coord_cartesian(xlim = c(stats::quantile(RatioDataLong$SILACRatio, 0.01, na.rm = T), 
                                                  stats::quantile(RatioDataLong$SILACRatio, 0.99, na.rm = T))) +
                theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)))
        Proteopedia::Reset_Dev()
        
        InteractiveData <- data.table(Protein = RatioDataWide$Gene, 
                                      Log2FC_CtlGroup = round(RatioDataWide$CtlGroup, digits = 3), 
                                      Log2FC_ExpGroup = round(RatioDataWide$ExpGroup, digits = 3), 
                                      GeneGroup = RatioDataWide$GeneGroup)
        InteractiveData[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup)))]
        pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(
          highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC_CtlGroup, y = Log2FC_ExpGroup)), zoomType = "xy"), 
                                title = list(text = paste0("Log2 Heavy:Light Protein LFQ Ratio (", CtlGroup, ")")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black", gridLineWidth = 0), 
          title = list(text = paste0("Log2 Heavy:Light Protein LFQ Ratio (", ExpGroup, ")")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black", gridLineWidth = 0), headerFormat = "", 
          pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC (Ctl): {point.Log2FC_CtlGroup:.3f}<br>Log2FC (Exp): {point.Log2FC_ExpGroup:.3f}"), 
          scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = TRUE), inactive = list(enabled = FALSE)), 
                         point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))), 
          c("#999999", "#880000", "#0033FF")) |> hc_title(text = paste0(ExpGroup, " vs. ", CtlGroup))
        htmlwidgets::saveWidget(pHC, "InteractiveCorrPlot.html")
      }
      message("Channel Ratio: Performing Median Normalisation")
      {
        RatioData[, `:=`(LFQ_Ratio_Norm, LFQ_Ratio - median(LFQ_Ratio, na.rm = TRUE) + median(RatioData$LFQ_Ratio, na.rm = TRUE)), by = Sample]
        NormalisationPlot <- ggplot2::ggplot(data.table::melt.data.table(RatioData, measure.vars = c("LFQ_Ratio", "LFQ_Ratio_Norm")), 
                                             ggplot2::aes(x = Condition, colour = Condition, y = value, group = Sample)) + 
          ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c(LFQ_Ratio = "Pre", SILACRatio_Norm = "Post"))) + 
          ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
          ggplot2::labs(x = NULL, y = expression("Log"[2]~"Heavy:Light Protein LFQ Ratio")) + ggplot2::theme(strip.background = ggplot2::element_blank())
      }
      message("Channel Ratio: Imputing NA Values")
      {
        RatioData <- data.table::data.table(tidyr::pivot_wider(RatioData[, .(ProteinGroup, Gene, Sample, LFQ_Ratio_Norm)], 
                                                               id_cols = ProteinGroup, values_from = LFQ_Ratio_Norm, 
                                                               names_from = Sample, values_fill = NA))
        RatioImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(RatioData, rownames = "ProteinGroup"), q = ImputationQ, 
                                                                      tune.sigma = ImputationSigma), keep.rownames = "ProteinGroup")
        RatioAll <- data.table::merge.data.table(data.table::melt.data.table(RatioData, id.vars = "ProteinGroup", value.name = "MeasuredRatio", 
                                                                                        variable.name = "Sample"), 
                                                 data.table::melt.data.table(RatioImp, id.vars = "ProteinGroup", value.name = "ImputedRatio", 
                                                                                        variable.name = "Sample"))
        RatioAll[, `:=`(Data, data.table::fifelse(is.na(MeasuredRatio), "Imputed", "Measured"))]
        
        ImputedDensity <- ggplot2::ggplot(RatioAll, ggplot2::aes(x = log10(ImputedRatio), fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + 
          ggplot2::scale_fill_manual(values = c(Measured = "black", Imputed = "magenta3")) + 
          ggplot2::labs(x = expression("Normalised Log"[10]~ "Heavy:Light Protein LFQ Ratio"), y = "Density", fill = NULL) + 
          ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))
        
        ImputedNAs <- data.table::melt.data.table(RatioData[ProteinGroup %in% Retained3$ProteinGroup], id.vars = "ProteinGroup", 
                                                             variable.name = "Sample")[is.na(value)][, `:=`(value, NULL)]
        ImputedNAs <- dplyr::summarise(dplyr::group_by(ImputedNAs, ProteinGroup), vector = paste(Sample, collapse = ", "))
        data.table::fwrite(ImputedNAs, "Imputed_Ratios.csv")
        RatioAll[, `:=`(Ratio, data.table::fifelse(is.na(MeasuredRatio), ImputedRatio, MeasuredRatio))]
        RatioAll <- RatioAll[, .(ProteinGroup, Sample, Ratio)]
      }
      message("Channel Ratio: Performing Paired T-Testing")
      {
        RatioAll <- data.table::merge.data.table(RatioAll, meta[, .(Condition, Sample)], by = "Sample")
        Ratio_Ttest <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(RatioAll, Condition, ProteinGroup), N = dplyr::n(), 
                                                           Log2MeanRatio = log2(mean(Ratio)), CV = (sd(Ratio)/mean(Ratio)) * 100), 
                                          id_cols = ProteinGroup, names_from = Condition, values_from = c(N, Log2MeanRatio, CV))
        CtlIndex <- which(grepl("Log2MeanRatio_", colnames(Ratio_Ttest)) & grepl(CtlGroup, colnames(Ratio_Ttest)))
        ExpIndex <- which(grepl("Log2MeanRatio_", colnames(Ratio_Ttest)) & grepl(ExpGroup, colnames(Ratio_Ttest)))
        Ratio_Ttest$Log2FC <- Ratio_Ttest[ExpIndex] - Ratio_Ttest[CtlIndex]
        Ttest_Output <- RatioAll[, `:=`(P.Value, stats::t.test(Ratio ~ Condition)$p.value), by = ProteinGroup]
        Ratio_Ttest <- data.table::merge.data.table(data.table(Ratio_Ttest), Ttest_Output, by = "ProteinGroup")
        Ratio_Ttest[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        Ratio_Ttest <- data.table::data.table(data.table::merge.data.table(Ratio_Ttest, ProteinInfo))
        Ratio_Ttest[, `:=`(Log2FC, as.numeric(Log2FC))]
        data.table::fwrite(dplyr::distinct(Ratio_Ttest), "Ratio_Paired_T-Test_Output.csv")
      }
      message("Channel Ratio: Fitting to Linear Model & Exporting Data")
      {
        LimmaInput <- data.table::dcast(RatioAll[Ratio > 0, Log2Ratio := log2(Ratio)], formula = ProteinGroup ~ Sample, value.var = "Log2Ratio") |>
                          data.table::setcolorder(neworder = c("ProteinGroup", meta$Sample))
        efit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, design), contr.matrix))
        if (!is.finite(efit$df.prior)) {
          message("Warning: Limma df.prior is Infinite")
        }
        MeanVarData <- data.table::data.table(efit$genes, Mean = efit$Amean, Variance = sqrt(efit$sigma))
        MeanVarData[, `:=`(Data, data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured"))]
        
        MeanVarPlot <- ggplot2::ggplot(MeanVarData, ggplot2::aes(x = Mean, y = Variance, colour = Data)) + 
          ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_manual(values = c(Measured = "black", Imputed = "magenta3"), guide = "none") + 
          ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "magenta3", stroke = NA) +
          ggplot2::labs(x = expression("Mean Log"[2]~"Heavy:Light Protein LFQ Ratio"), y = "Variance", colour = NULL) + 
          ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.9, 0.9))
        
        LimmaOutput <- data.table::setnames(data.table(limma::topTable(efit, coef = 1, adjust.method = "BH", n = Inf)), "logFC", "Log2FC")
        LimmaOutput <- data.table::merge.data.table(LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = TRUE), ], ProteinInfo, all.x = T)
        LimmaOutput[, `:=`(Significance, data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", 
                                                             data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None")))]
        LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = TRUE)
        LimmaOutput$Isoforms <- 1
        for (i in 1:nrow(LimmaOutput)) {
          if (length(stringr::str_extract_all(LimmaOutput$ProteinGroup[i],  "-\\d", simplify = T)) > 0) {
            LimmaOutput$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaOutput$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          } else { LimmaOutput$Isoforms[i] <- 1}
        }
        LimmaOutput[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
        LimmaOutput[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", GeneGroup))]
        data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
      }
      message("Channel Ratio: Generating Volcano Plots")
      {
        MeanLog2FC <- round(mean(LimmaOutput$Log2FC,  na.rm = T), digits = 3)
        LimmaOutput <- dplyr::arrange(LimmaOutput, desc(abs(t)))

        Ratio_TtestVolcanoData <- dplyr::arrange(Ratio_Ttest[, .(Gene, ProteinGroup, P.Value, Log2FC)], desc(abs(-log10(P.Value)*Log2FC))) |> 
          dplyr::distinct()
        
        pdf("Ratio_Limma_Volcanoes.pdf", height = 10, width = 14)
        print(ggplot2::ggplot(LimmaOutput, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
                ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("black") + 
                ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) + 
                ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
                ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
                ggplot2::annotate("text", x = min(LimmaOutput$Log2FC, na.rm = T)*0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
                ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Protein LFQ Ratio"), y = expression("-Log"[10]~"P-Value")) + 
                ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank()))
        print(ggplot2::ggplot(LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes", "No")], 
                              ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
                ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual(values = c(No = "black", Yes = "magenta3"), labels = c(No = "Measured", Yes = "Imputed"), guide = "none") + 
                ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) + 
                ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
                ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
                ggplot2::annotate("text", x = min(LimmaOutput$Log2FC, na.rm = T)*0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
                ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Protein LFQ Ratio"), y = expression("-Log"[10]~"P-Value")) + 
                ggplot2::ggtitle("Imputed Limma") + ggplot2::theme(legend.title = ggplot2::element_blank()))
        print(ggplot2::ggplot(Ratio_TtestVolcanoData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
                ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("black") + 
                ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) + 
                ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
                ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
                ggplot2::annotate("text", x = min(Ratio_Ttest$Log2FC, na.rm = T)*0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
                ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Protein LFQ Ratio"), y = expression("-Log"[10]~"P-Value")) + 
                ggplot2::ggtitle("Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank()))
        print(ggplot2::ggplot(Ratio_TtestVolcanoData[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes", "No")], 
                              ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) + 
                ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual(values = c(No = "black", Yes = "magenta3"), labels = c(No = "Measured", Yes = "Imputed"), guide = "none") + 
                ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) + 
                ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) + 
                ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black") + 
                ggplot2::annotate("text", x = min(Ratio_Ttest$Log2FC, na.rm = T)*0.9, y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC), size = 5) + 
                ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Protein LFQ Ratio"), y = expression("-Log"[10]~"P-Value")) + 
                ggplot2::ggtitle("Imputed Paired T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank()))
        Proteopedia::Reset_Dev()
        
        InteractiveData <- data.table(Protein = LimmaOutput[, Gene], 
                                        Log2FC = round(LimmaOutput[, Log2FC], digits = 2), 
                                        PValue = LimmaOutput[, P.Value], 
                                        Significance = factor(LimmaOutput[, Significance]), 
                                        Scien_PValue = formatC(LimmaOutput[, P.Value], format = "e", digits = 2), 
                                        GeneGroup = LimmaOutput[, GeneGroup])
        InteractiveData[, `:=`(URL, paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup)))]
        pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData, 
                                                                                                                                                                                "scatter", highcharter::hcaes(x = Log2FC, 
                                                                                                                                                                                                              y = -log10(PValue), group = Significance)), 
                                                                                                                                                            zoomType = "xy"), title = list(text = paste0(ExpGroup, " vs ", CtlGroup, " Log2 Fold-Difference")), 
                                                                                                                                      lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                                      tickColor = "black", gridLineWidth = 0), title = list(text = "-Log10 P-Value"), 
                                                                                                                lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                                                                                                                tickColor = "black", gridLineWidth = 0), headerFormat = "", 
                                                                                          pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}"), 
                                                                  scatter = list(marker = list(radius = 3), 
                                                                                 states = list(hover = list(enabled = TRUE), 
                                                                                               inactive = list(enabled = FALSE)), point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))), 
                                      c("#999999", "#880000", "#0033FF"))
        htmlwidgets::saveWidget(pHC, "InteractiveVolcanoPlot.html")
      }
      message("Channel Ratio: Exporting QC Plot")
      {
        pdf("ProteinQC_Limma_Plot.pdf", width = 18, height = 20)
        print(PCAPlot + patchwork::free(CountsBar) + patchwork::free(UpsetPlot) + patchwork::free(NormalisationPlot) + patchwork::free(ImputedDensity) + MeanVarPlot + 
                patchwork::plot_layout(design = "ABCC\nDDEE\nFFGG") + 
                patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E", "F"))))
        Proteopedia::Reset_Dev()
      }
    }
    # Output Analysis Parameters
    setwd(InputDirectory)
    data.table::fwrite(data.table::data.table("Input Directory" = paste0(InputDirectory), "Experimental Condition" = paste0(ExpGroup),
                                              "Control Condition" = paste0(CtlGroup), "Min. Precursors" = paste0(MinPrecursors, " Precursors"),
                                              "Imputation q-value" = ImputationQ, "Imputation Sigma-value" = ImputationSigma),
                                              "Analysis_Parameters.csv")
    Proteopedia::End_Timer(Start = start.time)
} 
########### MS Analysis: Timecourse SILAC DIA-NN To Protein Data Functions ####################################################################################################################################
#' @export
Process_TimecourseSILAC_DIANN <- function(InputDirectory, ProteotypicFiltering = F){
  set.seed(123)
  start.time <- Sys.time()
  setwd(InputDirectory)
  if(length(list.files(pattern = "AUC_Data.csv")) != 1){
    message("Check AUC_Data.csv File Availability\nUsing Summed Quantity Normalisation")
  }
  message("Loading Precursor Data & Generating Metadata")
  {
    if(length(list.files(pattern = "report.parquet")) > 0){
      PrecursorData <- data.table::data.table(arrow::read_parquet(list.files(pattern = "report.parquet")[1]))[, .(Run, Protein.Group, Protein.Ids, 
                                                                                                                  #First.Protein.Description, 
                                                                                                                  Genes, Stripped.Sequence, Channel,
                                                                                                                  Precursor.Id, Proteotypic, 
                                                                                                                  Precursor.Quantity, Quantity.Quality, 
                                                                                                                  #Precursor.Translated, 
                                                                                                                  Empirical.Quality, Channel.Q.Value, Q.Value, 
                                                                                                                  Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, 
                                                                                                                  Lib.Q.Value, Lib.PG.Q.Value)] |> 
        data.table::setnames(c("Protein.Group", "Genes", "Channel"), c("ProteinGroup", "Gene", "Label"))
      DIANN_Version <- 2.2
    } else if(length(list.files(pattern = "report.tsv")) > 0){
      PrecursorData <- data.table::fread(list.files(pattern = "report.tsv")[1])[, .(Run, Protein.Group, Protein.Ids, First.Protein.Description, Genes, 
                                                                                    Stripped.Sequence, Precursor.Id, Proteotypic, Precursor.Quantity,    
                                                                                    Quantity.Quality, Precursor.Translated, Channel.Q.Value,
                                                                                    Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value, 
                                                                                    Lib.PG.Q.Value)] |> data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), 
                                                                                                                             c("ProteinGroup", "Gene", "ProteinDescription"))
      DIANN_Version <- 1.8
    } else {return(message("ERROR: No InputFile Found"))}
    
    if(file.exists("Sample_Rename.csv") == TRUE){
      Sample_Rename <- data.table::fread("Sample_Rename.csv") 
      PrecursorData$Sample <- Sample_Rename$Renamed[match(unlist(PrecursorData$Run), Sample_Rename$Run)]
      PrecursorData <- PrecursorData[!is.na(Sample)]
    } else {
      PrecursorData[, Sample := gsub(".*IN_30_(.*)", "\\1", Run)]
    }
    
    PrecursorData[, Cell := gsub("(.*)_(.*)_(.*)_(.*)h_R(\\d)","\\1", Sample)]
    PrecursorData[, Drug := gsub("(.*)_(.*)_(.*)_(.*)h_R(\\d)","\\2", Sample)]
    PrecursorData[, Conc := gsub("(.*)_(.*)_(.*)_(.*)h_R(\\d)","\\3", Sample)]
    PrecursorData[, Time := gsub("(.*)_(.*)_(.*)_(.*)h_R(\\d)","\\4", Sample)]
    PrecursorData[, Replicate := gsub("(.*)_(.*)_(.*)_(.*)h_R(\\d)","\\5", Sample)]
    PrecursorData[, Sample := paste(Cell, data.table::fifelse(Drug == 0, "-", Drug), data.table::fifelse(Drug == 0, "-", Conc), Time, Replicate, sep = "_")]
    PrecursorData[, Sample := gsub("_-", "", Sample)]
    PrecursorData[, Condition := gsub("(.*)_\\d", "\\1", Sample)]
    PrecursorData[, Cluster := gsub("(.*)_.*_\\d","\\1", Sample)]
    if(DIANN_Version == 2.2){
      PrecursorData[ , Label := gsub("L$", "Light", Label)]
      PrecursorData[ , Label := gsub("H$", "Heavy", Label)]
    }
    PrecursorData[, Sample := gsub("^0_", "",  gsub("_0$", "", Sample))]
    PrecursorData[, Condition := gsub("^0_", "",  gsub("_0$", "", Condition))]
    PrecursorData[, Cluster := gsub("^0_", "",  gsub("_0$", "", Cluster))]
    
    Metadata <- PrecursorData[, .(Run, Sample, Cell, Drug, Conc, Time, Replicate, Condition, Cluster)] |> dplyr::distinct() |> dplyr::arrange(Cell, Drug, Conc, as.numeric(gsub("h", "", Time)), Replicate)
    SampleLevels <- unique(Metadata$Sample)
    CellLevels <- unique(Metadata$Cell)
    DrugLevels <- unique(Metadata$Drug)
    ConcLevels <- unique(Metadata$Conc)
    TimeLevels <- unique(Metadata$Time)
    ConditionLevels <- unique(Metadata$Condition)
    ClusterLevels <- unique(Metadata$Cluster)
    
    Metadata[, Sample := factor(Sample, levels = SampleLevels)]
    Metadata[, Cell := factor(Cell, levels = CellLevels)]
    Metadata[, Drug := factor(Drug, levels = DrugLevels)]
    Metadata[, Conc := factor(Conc, levels = ConcLevels)]
    Metadata[, Time := factor(Time, levels = TimeLevels)]
    Metadata[, Condition := factor(Condition, levels = ConditionLevels)]
    Metadata[, Cluster := factor(Cluster, levels = ClusterLevels)]
    
    # Plot PCA on Original Data
    if(DIANN_Version != 2.2){
      All_PCA <- PrecursorData[Precursor.Quantity == 0, Precursor.Quantity := NA] |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
      
      OriginalData_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate)) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel(ggplot2::aes(label = Time)) +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate)) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel(ggplot2::aes(label = Time)) +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Original Data")
    } else {
      All_PCA <- PrecursorData |> data.table::copy()
      All_PCA <- All_PCA[Precursor.Quantity == 0, Precursor.Quantity := NA][, Precursor.Id := paste0(Precursor.Id, "_", Label)] |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
      
      OriginalData_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate)) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel(ggplot2::aes(label = paste0(Time, "h"))) +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate)) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel(ggplot2::aes(label = paste0(Time, "h"))) +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Original Data")
    }
  }
  message("Filtering Precursors")
  {
    PrecursorData[Precursor.Quantity == 0, Precursor.Quantity := NA]
    if(DIANN_Version != 2.2){
      PrecursorData[Precursor.Translated == 0, Precursor.Translated := NA]
      N_Samples <- PrecursorData[, .(Original = sum(Precursor.Quantity, na.rm = T), Translated = sum(Precursor.Translated, na.rm = T)), .(Run, Sample)] |> 
        data.table::melt.data.table(id.vars = c("Run", "Sample"), measure.vars = c("Original", "Translated"), variable.name = "Measure", value.name = "Total_Intensity") |>
        data.table::merge.data.table(Metadata)
    } else {
      N_Samples <- PrecursorData[, .(Original = sum(Precursor.Quantity, na.rm = T)), .(Run, Sample)] |> 
        data.table::melt.data.table(id.vars = c("Run", "Sample"), measure.vars = "Original", variable.name = "Measure", value.name = "Total_Intensity") |>
        data.table::merge.data.table(Metadata)
    }
    
    if(length(list.files(pattern = "AUC_Data.csv")) == 1){
      AUCs <- data.table::fread("AUC_Data.csv")[, Run := mzML]
      Loading_Plot <- data.table::merge.data.table(N_Samples, AUCs[, Run := stringr::str_remove(mzML,".mzML.*")], by = "Run") |> 
        ggplot2::ggplot(ggplot2::aes(x = AUC, y = Total_Intensity, colour = Cluster, label = paste0(Time, "h"), shape = Replicate)) + 
        ggplot2::geom_point() + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + Proteopedia::Add_Rsq(Subgroups = F) + 
        ggrepel::geom_text_repel() + ggplot2::facet_wrap(~Measure) + ggplot2::labs(x = "Area Under Curve (AUC)", y = "Total Intensity") +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = ggplot2::aes(label = NA))) + 
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
    } else {
      SummedPrecursorQuants <- PrecursorData[, .(SummedQuants = sum(Precursor.Quantity, na.rm = T)), by = Run]
      Loading_Plot <- data.table::merge.data.table(N_Samples, SummedPrecursorQuants, by = "Run") |> 
        ggplot2::ggplot(ggplot2::aes(x = SummedQuants, y = Total_Intensity, colour = Cluster, label = paste0(Time, "h"), shape = Replicate)) + 
        ggplot2::geom_point() + ggplot2::scale_colour_brewer(palette = "Set1") + Proteopedia::Add_Rsq(Subgroups = F) + 
        ggrepel::geom_text_repel() + ggplot2::facet_wrap(~Measure) + ggplot2::labs(x = "Summed Precursor Quantities", y = "Total Intensity") +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = ggplot2::aes(label = NA))) + 
        ggplot2::theme(strip.background = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
    }
    
    if(ProteotypicFiltering == T){PrecursorData <- PrecursorData[Proteotypic >= 1]}
    if(DIANN_Version != 2.2){
      PreIonCV_Data <- PrecursorData[, .(CV = sd(Precursor.Quantity, na.rm = T)/mean(Precursor.Quantity, na.rm =T), 
                                         N_Points = .N, Avg_Quant = mean(Precursor.Quantity, na.rm = T), 
                                         Quality = mean(Quantity.Quality, na.rm = T)), by = .(Precursor.Id, Condition)]
      PrecursorData <- PrecursorData[!is.na(Precursor.Quantity) & Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & 
                                       Channel.Q.Value < max(c(0.01, min(PrecursorData[grepl("SILAC-[RK]-H", Precursor.Id) & grepl(".*_0_\\d", Sample)]$Channel.Q.Value, na.rm = T)))]    
      PostIonCV_Data <- PrecursorData[, .(CV = sd(Precursor.Quantity, na.rm = T)/mean(Precursor.Quantity, na.rm =T), 
                                          N_Points = .N, Avg_Quant = mean(Precursor.Quantity, na.rm = T),
                                          Quality = mean(Quantity.Quality, na.rm = T)), by = .(Precursor.Id, Condition)]
    } else {
      PreIonCV_Data <- PrecursorData[, .(CV = sd(Precursor.Quantity, na.rm = T)/mean(Precursor.Quantity, na.rm =T), 
                                         Empirical.Quality = mean(Empirical.Quality, na.rm = T), 
                                         N_Points = .N, Avg_Quant = mean(Precursor.Quantity, na.rm = T), 
                                         Quality = mean(Quantity.Quality, na.rm = T)), by = .(Precursor.Id, Label, Condition)]
      PrecursorData <- PrecursorData[!is.na(Precursor.Quantity) & Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & 
                                       Channel.Q.Value < max(c(0.01, min(PrecursorData[grepl("(SILAC)", Precursor.Id) & Label == "H" & grepl(".*_0_\\d", Sample)]$Channel.Q.Value, na.rm = T)))]
      PostIonCV_Data <- PrecursorData[, .(CV = sd(Precursor.Quantity, na.rm = T)/mean(Precursor.Quantity, na.rm =T), 
                                          Empirical.Quality = mean(Empirical.Quality, na.rm = T),
                                          N_Points = .N, Avg_Quant = mean(Precursor.Quantity, na.rm = T),
                                          Quality = mean(Quantity.Quality, na.rm = T)), by = .(Precursor.Id, Label, Condition)]
    }
    
    IonCV_Data <- PreIonCV_Data[, Dataset := "Pre-Filtering"] |> rbind(PostIonCV_Data[, Dataset := "Post-Filtering"])
    IonCV_Data <- IonCV_Data[N_Points == length(unique(PrecursorData$Replicate))]
    IonCV_Data[, CutQuality := gsub("]", "", gsub(".*,", "≤ ", cut(Quality, 10)))]
    
    IonCV_Plot <- IonCV_Data |> ggplot2::ggplot(ggplot2::aes(x = CutQuality, y = CV)) + 
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "#F00", alpha = 0.3) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.5, fill = "#F90", alpha = 0.3) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.2, fill = "#0F0", alpha = 0.3) +
      ggplot2::geom_boxplot(outliers = F, fill = NA) + ggplot2::facet_wrap(~forcats::fct_rev(Dataset)) + ggplot2::labs(x = "Quality", y = "Ion-Level CV") + 
      ggplot2::theme(strip.background = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    # Plot PCA on Filtered Data
    All_PCA <- PrecursorData |> data.table::copy()
    if(DIANN_Version == 2.2){
      All_PCA <- All_PCA[, Precursor.Id := paste0(Precursor.Id, "_", Label)] |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    } else {
      All_PCA <- All_PCA |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    }
    
    FilteredData_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "Filtered Data")
  }
  message("Normalising Data")
  {
    if(length(list.files(pattern = "AUC_Data.csv")) == 1){
      if(DIANN_Version != 2.2){
        PrecursorData <- PrecursorData |> data.table::merge.data.table(AUCs[, .(Run, AUC)], by = "Run")
        PrecursorData[, Precursor.Quantity   := Precursor.Quantity / AUC * median(AUCs$AUC), by = Sample]
        PrecursorData[, Precursor.Translated := Precursor.Translated / AUC *  median(AUCs$AUC), Sample]
      } else {
        PrecursorData <- PrecursorData |> data.table::merge.data.table(AUCs[, .(Run, AUC)], by = "Run")
        PrecursorData[, Precursor.Quantity   := Precursor.Quantity / AUC * median(AUCs$AUC), by = Sample]
      }
    } else {
      if(DIANN_Version != 2.2){
        PrecursorData[, `:=`(Precursor.Quantity, Precursor.Quantity/sum(Precursor.Quantity, na.rm = T) * PrecursorData[, sum(Precursor.Quantity, na.rm = T), Run][, median(V1, na.rm = T)]), by = Run]
        PrecursorData[, `:=`(Precursor.Translated, Precursor.Translated/sum(Precursor.Translated, na.rm = T) * PrecursorData[, sum(Precursor.Translated, na.rm = T), Run][, median(V1, na.rm = T)]), by = Run]
      }
      PrecursorData[, `:=`(Precursor.Quantity, Precursor.Quantity/sum(Precursor.Quantity, na.rm = T) * PrecursorData[, sum(Precursor.Quantity, na.rm = T), Run][, median(V1, na.rm = T)]), by = Run]
    }
    PrecursorData <- PrecursorData[!is.na(Precursor.Quantity)]
    
    QuantityQuality_Density <- PrecursorData |> ggplot2::ggplot(ggplot2::aes(x = Quantity.Quality, fill = Cluster)) + 
      ggplot2::geom_histogram() + ggplot2::labs(x = "Quantity-Quality", y = "Density") + 
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
    
    Log2Intensity_Density <- PrecursorData |> ggplot2::ggplot(ggplot2::aes(x = log2(Precursor.Quantity), fill = Cluster)) + 
      ggplot2::geom_histogram() + ggplot2::labs(x = expression("Log"[2]~"Ion Precursor Quantity"), y = "Density") + 
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
      ggplot2::scale_y_continuous(expand = c(0, 0))
    
    # Plot PCA on AUC-Normalised Data
    All_PCA <- PrecursorData |> data.table::copy() 
    if(DIANN_Version == 2.2){
      All_PCA <- All_PCA[, Precursor.Id := paste0(Precursor.Id, "_", Label)] |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    } else {
      All_PCA <- All_PCA |> tidyr::pivot_wider(id_cols = Precursor.Id, values_from = Precursor.Quantity, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = "Precursor.Id") |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    }
    NormalisedData_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "Normalised Data")    
    
    CorrData <- PrecursorData |> data.table::copy() 
    if(DIANN_Version == 2.2){
      CorrSamples <- stats::cor(CorrData[, Precursor.Id := paste0(Precursor.Id, "_", Label)] |> data.table::dcast(Precursor.Id~Sample , value.var =  "Precursor.Quantity") |> tibble::column_to_rownames("Precursor.Id"),
                                use = "pairwise.complete.obs", method = "spearman")
    } else {
      CorrSamples <- stats::cor(CorrData |> data.table::dcast(Precursor.Id~Sample, value.var =  "Precursor.Quantity") |> tibble::column_to_rownames("Precursor.Id"),
                                use = "pairwise.complete.obs", method = "spearman")
    }
    
    AnnotationCol <- data.table::data.table(ID = rownames(CorrSamples), Cluster = gsub("(.*)_.*_\\d", "\\1", rownames(CorrSamples)),
                                            Time = factor(gsub(".*_(.*)_\\d", "\\1", rownames(CorrSamples)), 
                                                          levels = sort(as.numeric(unique(gsub("h","", gsub(".*_(.*)_\\d", "\\1", 
                                                                                                            rownames(CorrSamples))))))),
                                            Replicate = gsub(".*_(\\d)", "\\1", rownames(CorrSamples))) |> tibble::column_to_rownames("ID") 
    
    AnnotationColours <- list(Cluster = Proteopedia::NiceColourPalette[1:length(ClusterLevels)],
                              Time = rev(RColorBrewer::brewer.pal(length(TimeLevels), "Spectral")),
                              Replicate = RColorBrewer::brewer.pal(length(unique(Metadata$Replicate)), "Pastel1"))
    names(AnnotationColours$Cluster) <- unique(AnnotationCol$Cluster)
    names(AnnotationColours$Time) <- levels(AnnotationCol$Time)
    names(AnnotationColours$Replicate) <- unique(AnnotationCol$Replicate)
    
    Corr_Heatmap <- pheatmap::pheatmap(CorrSamples, main = "Spearman Complete Observation Post-Filtering", show_rownames = T, show_colnames = F, 
                                       annotation_col = AnnotationCol, annotation_colors = AnnotationColours)
  }
  message("Processing SILAC Labels & Calculating Ratios")
  {
    if(DIANN_Version != 2.2){
      if(PrecursorData[stringr::str_detect(Precursor.Id , "SILAC-.-L") & stringr::str_detect(Precursor.Id ,"SILAC-.-H"), .N ] != 0){
        message("ERROR: Multi-Label Precursors Detected")
        MultiLabelDetection <- T
      } else {MultiLabelDetection <- F} 
      
      PrecursorData[stringr::str_detect(Precursor.Id ,"SILAC-.-L"), Label := "Light"]
      PrecursorData[stringr::str_detect(Precursor.Id ,"SILAC-.-H"), Label := "Heavy"]
      PrecursorData[, Precursor.Id.nolabels := gsub("SILAC-.-.", "SILAC", Precursor.Id)]
      PrecursorData[, Log2Precursor.Value := log2(Precursor.Translated)]
    } else {
      PrecursorData[, Precursor.Id.nolabels := Precursor.Id]
      PrecursorData[, Log2Precursor.Value := log2(Precursor.Quantity)]
    }
    PrecursorData <- PrecursorData[, Label := data.table::fifelse(nchar(Label) == 0, NA, Label)]
    PrecursorData <- PrecursorData[!is.na(Label)]
    
    ChannelRatio_Data <- PrecursorData|> data.table::dcast(Precursor.Id.nolabels+ProteinGroup+Sample ~ Label, value.var = "Log2Precursor.Value")
    ChannelRatio_Data[, Log2Ratio := Heavy/Light]
    ChannelRatio_Data <- ChannelRatio_Data[!is.na(Log2Ratio), .(ProteinGroup, Log2Ratio, Sample, ProteinGroup)]
    ChannelRatio_Data <- ChannelRatio_Data[, .(Log2Ratio = median(Log2Ratio, na.rm = TRUE), N_Precursors = .N), .(Sample, ProteinGroup)]
    
    # Channel Ratio PCA
    All_PCA <- ChannelRatio_Data |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Log2Ratio, names_from = Sample, values_fill = NA) |> 
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(All_PCA)$importance
    All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    
    ChannelRatio_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "Heavy:Light Channel Ratio Data")    
    
    ChannelCounts <- PrecursorData[, .(Count = sum(Precursor.Quantity)), by = .(Sample, Label)] |> 
      data.table::merge.data.table(PrecursorData[, .(TotalCount = sum(Precursor.Quantity)), by = Sample]) |> 
      data.table::merge.data.table(Metadata)
    ChannelCounts[, LabelProp := Count/TotalCount] 
    
    ChannelPropPlot <- ChannelCounts |> ggplot2::ggplot(ggplot2::aes(x = Time, y = LabelProp, alpha = Label, fill = Cluster)) + 
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::facet_grid(ggplot2::vars(Cluster), ggplot2::vars(paste0("Rep. ", Replicate))) +
      ggplot2::scale_fill_brewer(palette = "Set1", guide = "none") + 
      ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), labels = c("Heavy", "Light"), name = "Channel") + 
      ggplot2::labs(x = "Time (hours)", y = "Proportion of Channel Intensities") + 
      ggplot2::theme(strip.background = ggplot2::element_blank(),strip.text.y = ggplot2::element_text(size = 26), strip.text.x = ggplot2::element_text(size = 22))
  }
  message("Calculating Protein Weights")
  {
    ProteinWeights <- PrecursorData[,.(Sum_Intensities  = sum(log2(Precursor.Quantity), na.rm = T)), by = .(Sample, Label, ProteinGroup)] |> 
      data.table::merge.data.table(Metadata)
    ProteinWeights[, MeanSum := mean(Sum_Intensities, na.rm =T ), by = .(ProteinGroup, Label)]
    ProteinWeights[, RowID := paste0(ProteinGroup, "_", Label)]
    ProteinWeights[, Weight := Sum_Intensities/MeanSum]
    
    # Protein Weights PCA
    All_PCA <- ProteinWeights |> tidyr::pivot_wider(id_cols = RowID, values_from = Weight, names_from = Sample, values_fill = NA) |> 
      tidyr::drop_na() |> data.frame(row.names = "RowID") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(All_PCA)$importance
    All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    
    ProteinWeights_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "Protein Weights")    
    
    ProteinWeight_Density <- ProteinWeights |> ggplot2::ggplot(ggplot2::aes(x = Weight, fill = Time)) + 
      ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap("Label", nrow = 2, strip.position = "right", scales = "free_x") + 
      ggplot2::scale_x_log10() + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) +
      ggplot2::labs(title = "Protein Weights for Limma", x = "Protein Weight", y = "No. Proteins", fill = "Time") +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Compiling & Exporting Data")
  {
    LFQ_T <- Calculate_LFQ(PrecursorData, "LFQ", SILAC = T)  # For total sample (all precursors)
    LFQ_H <- Calculate_LFQ(PrecursorData[Label == "Heavy"], "LFQ_H", SILAC = T)  # For heavy fraction
    LFQ_L <- Calculate_LFQ(PrecursorData[Label == "Light"], "LFQ_L", SILAC = T)  # For light fraction
    Intensities_T <- PrecursorData[, .(Intensity = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]  # For total sample (all precursors)
    Intensities_H <- PrecursorData[Label == "Heavy", .(Intensity_H = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]  # For heavy fraction 
    Intensities_L <- PrecursorData[Label == "Light", .(Intensity_L = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]  # For light fraction 
    Counts_H <- PrecursorData[Label == "Heavy", .(N_precursors_H = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_H = sum(Proteotypic)), .(ProteinGroup, Sample) ]
    Counts_L <- PrecursorData[Label == "Light", .(N_precursors_L = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_L = sum(Proteotypic)), .(ProteinGroup, Sample) ]
    if(DIANN_Version != 2.2){
      Annotations <- unique(PrecursorData[, .(ProteinGroup, Run, Sample, Condition, Cluster, Cell, Drug, Conc, Time, Replicate, ProteinDescription, Gene)])
    } else {
      Annotations <- unique(PrecursorData[, .(ProteinGroup, Run, Sample, Condition, Cluster, Cell, Drug, Conc, Time, Replicate, Gene)])
    }
    ProteinData <- Reduce(Proteopedia::Merge_PrecursorData, list(LFQ_T, LFQ_H, LFQ_L, Intensities_T, Intensities_H, Intensities_L, 
                                                                 Counts_H, Counts_L, Annotations))                               
    ProteinData <- ProteinData[!is.na(Intensity)] |> data.table::setcolorder(neworder = c(2, 1, 16, 17))
    ProteinData[, SILACRatio := LFQ_H/LFQ_L]
    
    data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
    data.table::fwrite(ProteinWeights, "ProteinWeights.csv.gz")
    data.table::fwrite(ProteinData, "ProteinQuantities.csv.gz")
    
    LFQIntensity_Plot <- ProteinData[order(ProteinGroup)][1:1000,] |> ggplot2::ggplot(ggplot2::aes(x = LFQ_L, y = Intensity_L, colour = ProteinGroup)) + 
      ggplot2::geom_point(stroke = NA, alpha = 0.5) + ggplot2::labs(x = "Light Protein LFQ", y = "Light Protein Intensity") +
      ggplot2::scale_x_log10() +ggplot2::scale_y_log10() + Proteopedia::Add_XYLine("grey") + Proteopedia::Add_Rsq(F) +
      ggplot2::theme(legend.position = "none")
    
    # Protein Light LFQ PCA
    All_PCA <- ProteinData |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = LFQ_L, names_from = Sample, values_fill = NA) |> 
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(All_PCA)$importance
    All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    
    LightLFQ_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "Light LFQ Data")    
    
    # Protein Heavy LFQ PCA
    #All_PCA <- ProteinData[Time != "0h"] |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = LFQ_H, names_from = Sample, values_fill = NA) |> 
    #  tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    #SummaryPCA <- summary(All_PCA)$importance
    #All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    
    #HeavyLFQ_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cell, shape = Replicate, label = paste0(Time, "h"))) +
    #  ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
    #  ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
    #                y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
    #  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
    #  All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cell, shape = Replicate, label = paste0(Time, "h"))) +
    #  ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
    #  ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
    #                y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
    #  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
    #  patchwork::plot_annotation(title = "Heavy LFQ Data")  
    
    # Protein All LFQ PCA
    All_PCA <- ProteinData |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = LFQ, names_from = Sample, values_fill = NA) |> 
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(All_PCA)$importance
    All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    
    AllLFQ_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                    y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
      All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
      ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
      ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                    y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
      patchwork::plot_annotation(title = "All LFQ Data")     
  }
  message("Plotting Intensities") 
  {  
    IntensitiesData <- data.table::rbindlist(list(
      PrecursorData[!is.na(Label), .(Sample, Label, `log2 quantity` = log2(Precursor.Quantity), Type = "Precursor Quantity")],
      ProteinData[, .(Sample, `log2 quantity` = log2(LFQ_L), Label = "Light", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, `log2 quantity` = log2(Intensity_L), Label = "Light", Type = "Protein Intensity")],
      ProteinData[, .(Sample, `log2 quantity` = log2(LFQ_H), Label = "Heavy", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, `log2 quantity` = log2(Intensity_H), Label = "Heavy", Type = "Protein Intensity")],
      ProteinData[, .(Sample, `log2 quantity` = log2(LFQ), Label = "Total", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, `log2 quantity` = log2(Intensity), Label = "Total", Type = "Protein Intensity")]), use.names = TRUE) |>
      data.table::merge.data.table(Metadata)
    IntensitiesData[, Type := factor(Type, levels = c("Precursor Quantity", "Max. Protein LFQ", "Protein Intensity"))]
    
    IntensityPlot <- IntensitiesData |> ggplot2::ggplot(ggplot2::aes(x = factor(gsub(".*_(\\d+)_(\\d)", "\\1h R\\2", Sample),
                                                                                levels = unique(gsub(".*_(\\d+)_(\\d)", "\\1h R\\2", SampleLevels))), 
                                                                     y = `log2 quantity`, colour = Cluster)) + 
      ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + 
      ggplot2::facet_grid(cols = ggplot2::vars(Type), rows = ggplot2::vars(Label), scales = "free_x") + 
      ggplot2::ylab(expression("Log"[2]~"Value")) + ggplot2::coord_flip() + 
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                     strip.text.y = ggplot2::element_text(size = 26), legend.position = "inside", 
                     legend.position.inside = c(0.2, 0.2), legend.title = ggplot2::element_blank())
  }
  message("Plotting Precursor, Peptide & Protein Counts") 
  {
    CountData <- data.table::melt.data.table(PrecursorData[!is.na(Label), lapply(.SD, data.table::uniqueN), .(Sample, Condition, Label), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")], 
                                             id.vars = c("Sample","Condition", "Label"), value.name = "IDs") |> data.table::merge.data.table(Metadata)
    
    CountPlot <- CountData |> ggplot2::ggplot(ggplot2::aes(x = factor(gsub("_", " ", Sample), levels = gsub("_", " ", SampleLevels)), y = IDs/1000, fill = Cluster, alpha = Label, 
                                                           label = format(IDs, big.mark = ",", scientific = FALSE))) +
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.1) +
      ggplot2::facet_grid(ggplot2::vars(Label), ggplot2::vars(variable), scales = "free_x", 
                          labeller = ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups", 
                                                            Light = "Light", Heavy = "Heavy"))) + 
      ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                     strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Calculating Data Completeness") 
  {
    CompletenessData <- rbind(Proteopedia::Count_Proteins(ProteinData, "All"), 
                              Proteopedia::Count_Proteins(ProteinData[N_precursors_H >= 2 | N_precursors_L >= 2], "≥ 2"), 
                              Proteopedia::Count_Proteins(ProteinData[N_precursors_proteotypic_H >= 2 | N_precursors_proteotypic_L >= 2], "≥ 2 Proteotypic"))
    NAsPlot <- ggplot2::ggplot(CompletenessData, ggplot2::aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors)) + 
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c(All = "black", `≥ 2` = "darkgrey", `≥ 2 Proteotypic` = "orange3")) + 
      ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") + ggplot2::scale_x_continuous(breaks = seq(1, 1000, 1)) + 
      ggplot2::scale_y_continuous(limits = c(0, max(CompletenessData$cumulative_protein_N)/1000)) + 
      ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"), 
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"), 
                     legend.position = "inside", legend.position.inside = c(0.25, 0.25))
  }
  message("Calculating Channel Skewing")
  {
    SkewData <- data.table::data.table(Metadata)
    SkewData[, `:=`(PearsonsSkewRatio = 0, PearsonsSkew_L = 0, PearsonsSkew_H = 0, Median = 0, Mean = 0)]
    for(RowIndex in 1:nrow(SkewData)){
      SampleID <- SkewData$Sample[RowIndex]
      SkewData$PearsonsSkewRatio[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, SILACRatio])
      SkewData$PearsonsSkew_L[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ_L])
      SkewData$PearsonsSkew_H[RowIndex] <- Proteopedia::Calculate_PearsonsSkew(ProteinData[Sample == SampleID, LFQ_H])
      SkewData$Median[RowIndex] <- median(ProteinData[Sample == SampleID, SILACRatio], na.rm = T)
      SkewData$Mean[RowIndex] <- mean(ProteinData[Sample == SampleID, SILACRatio], na.rm = T)
    }
    
    SkewPlot <- ProteinData |> ggplot2::ggplot(ggplot2::aes(x = factor(gsub("_", " ", Sample), levels = gsub("_", " ", SampleLevels)), 
                                                            y = SILACRatio, colour = Cluster)) + 
      ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
      ggplot2::ylab("Heavy:Light Protein LFQ Ratio") + ggplot2::coord_flip() + 
      ggplot2::geom_text(data = SkewData[!is.na(PearsonsSkewRatio)], ggplot2::aes(y = 0, label = paste0(round(PearsonsSkewRatio, 2))), size = 4, hjust = 1) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Plotting Missed Trypsinisation Sites")
  {
    TrypsinData <- data.table::copy(PrecursorData)
    TrypsinData[, `:=`(MissedTrypsin, grepl("[RK][^P]", Stripped.Sequence))]
    TrypsinData[, `:=`(N_Trypsin, .N), by = .(Sample, MissedTrypsin, Label)]
    TrypsinData[, `:=`(N_Sample, .N), by = .(Sample, Label)]
    TrypsinData <- dplyr::distinct(TrypsinData[MissedTrypsin == T, .(Sample, Condition, Replicate, Label, N_Trypsin, N_Sample)]) |>
      data.table::merge.data.table(Metadata)
    TrypsinData[, `:=`(PercentTrypsin, (N_Trypsin/N_Sample) * 100)]
    MissedCleavagePlot <- ggplot2::ggplot(TrypsinData, ggplot2::aes(x = Replicate, y = PercentTrypsin, fill = Cluster, alpha = Label)) + 
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
      ggplot2::facet_grid(ggplot2::vars(factor(Cluster, levels = ClusterLevels)), ggplot2::vars(factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h"))), scales = "free_y") +
      ggplot2::labs(y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::scale_alpha_manual(values = c(Heavy = 0.5, Light = 1), guide = "none") + 
      ggplot2::coord_flip() + ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                                             strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Calculating Precursor & Protein Variation")
  {
    PrecursorCVs <- PrecursorData[, .(CV = Proteopedia::Calculate_CV(Precursor.Quantity), N = .N), .(Precursor.Id, Condition, Cluster, Cell, Drug, Conc, Time, Label)]
    PrecursorCVs <- PrecursorCVs[N >= 3]
    PrecursorCVs[, `:=`(rank, data.table::frank(CV)), .(Condition, Cluster, Cell, Drug, Conc, Time, Label)]
    PrecursorCVs[, `:=`(ID, "Precursors")]
    ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_L), N = .N), .(ProteinGroup, Cell, Drug, Conc, Time, Condition, Cluster)][, `:=`(Label, "Light")] |>
      rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_H), N = .N), .(ProteinGroup, Cell, Drug, Conc, Time, Condition, Cluster)][, `:=`(Label, "Heavy")])
    ProteinCVs <- ProteinCVs[N >= 3]
    ProteinCVs <- ProteinCVs[, `:=`(rank, data.table::frank(CV)), .(Condition, Cluster, Cell, Drug, Conc, Time, Label)]
    ProteinCVs[, `:=`(ID, "Protein Groups")]
    AllCVs <- rbind(PrecursorCVs[, `:=`(Precursor.Id, NULL)], ProteinCVs[, `:=`(ProteinGroup, NULL)])
    VariationPlot <- AllCVs |> ggplot2::ggplot(ggplot2::aes(x = rank/1000, y = CV, colour = Cluster, alpha = Label)) + 
      ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Variation [%]") + 
      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
      ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") + ggplot2::coord_cartesian(ylim = c(0, 50)) + 
      ggplot2::facet_grid(ggplot2::vars(ID), ggplot2::vars(factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h"))), scales = "free") +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"), 
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"), 
                     strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Exporting QC Plots")
  {
    pdf("QC_Processing_Plots.pdf", width = 18, height = 20)
    print(patchwork::free(Loading_Plot) + patchwork::free(IonCV_Plot) + patchwork::free(QuantityQuality_Density) + patchwork::free(Log2Intensity_Density) + 
            ChannelPropPlot + ProteinWeight_Density + patchwork::free(LFQIntensity_Plot) + 
            patchwork::plot_layout(design = "AABB\nCDEE\nFFGG") + patchwork::plot_annotation(tag_levels = "A"))
    print(IntensityPlot)
    print(CountPlot)
    print(NAsPlot + SkewPlot + MissedCleavagePlot + patchwork::free(VariationPlot) +
            patchwork::plot_layout(design = "AB\nCD") + patchwork::plot_annotation(tag_levels = "A"))
    #print(Corr_Heatmap)
    Proteopedia::Reset_Dev()
    
    pdf("PCA_Plots.pdf", width = 18, height = 10)
    print(OriginalData_PCA)
    print(FilteredData_PCA)
    print(NormalisedData_PCA)
    print(ChannelRatio_PCA)
    print(ProteinWeights_PCA)
    print(LightLFQ_PCA)
    #print(HeavyLFQ_PCA)
    print(AllLFQ_PCA)
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(start.time)
}
#' @export
Analyse_TimecourseSILAC_Proteins <- function(InputDirectory, CtlGroup, ExpGroups, GenerateDataPlots = F, SameInitialAbundance = T, LightMinSamples = 0.5, 
                                              HeavyModel = "NLS", HeavyMinMonotonicity = 0.5, HeavyMinSamples = 0.5, HeavyMaxCV = 0.3, 
                                              UseLightKloss = F, ReplicatesUsed = c(1, 2, 3)){
  set.seed(123)
  start.time <- Sys.time()
  for(ExpGroup in ExpGroups){
    setwd(InputDirectory)
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Loading Data"))
    {
      ProtLFQsInput <- data.table::fread(list.files(pattern = "ProteinQuant"))[grepl(CtlGroup, Condition) | grepl(ExpGroup, Condition)] |> data.table::setnames("Protein_group", "ProteinGroup", skip_absent = T)
      ProtLFQsInput |> data.table::setnames(c(colnames(ProtLFQsInput)[grepl("Gene.*", colnames(ProtLFQsInput))], colnames(ProtLFQsInput)[grepl("Protein.*group", colnames(ProtLFQsInput), ignore.case = T)]), 
                                            c("Gene", "ProteinGroup"))
      ProtLFQsInput <- ProtLFQsInput[Replicate %in% ReplicatesUsed]
      
      for(ColumnID in c("Cell", "Drug", "Conc")){
        if(length(unique(grepl(ColumnID, colnames(ProtLFQsInput)))) == 1){
          ProtLFQsInput[, PlaceholderName := 0]
          ProtLFQsInput |> data.table::setnames("PlaceholderName", ColumnID)
        }
      }
      if(length(unique(grepl("Cluster", colnames(ProtLFQsInput)))) == 1){
        ProtLFQsInput[, Cluster := gsub("(.*)_\\d+h", "\\1", Condition)]
      }
      if(length(unique(grepl("Time", colnames(ProtLFQsInput)))) == 1){
        ProtLFQsInput[, Time := as.numeric(gsub(".*_(\\d+)h", "\\1", Condition))]
      }
      
      Metadata <- ProtLFQsInput[, .(Run, Sample, Cell, Drug, Conc, Time, Replicate, Condition, Cluster)] |> unique() |> 
        data.table::setorderv(c("Cell", "Drug", "Conc", "Time", "Replicate"))
      SampleLevels <- unique(Metadata$Sample)
      CellLevels <- unique(Metadata$Cell)
      DrugLevels <- unique(Metadata$Drug)
      ConcLevels <- unique(Metadata$Conc)
      TimeLevels <- unique(Metadata$Time)
      ConditionLevels <- unique(Metadata$Condition)
      ClusterLevels <- unique(Metadata$Cluster)
      Metadata[, Sample := factor(Sample, levels = SampleLevels)]
      Metadata[, Cell := factor(Cell, levels = CellLevels)]
      Metadata[, Drug := factor(Drug, levels = DrugLevels)]
      Metadata[, Conc := factor(Conc, levels = ConcLevels)]
      Metadata[, Time := factor(Time, levels = TimeLevels)]
      Metadata[, Condition := factor(Condition, levels = ConditionLevels)]
      Metadata[, Replicate := factor(Replicate)]
      Metadata[, Cluster := factor(Cluster, levels = ClusterLevels)]
      
      ProteinInfo <- ProtLFQsInput[,.(ProteinGroup, #ProteinDescription, 
                                      Gene)] |> unique()
      
      ProtWeights <- data.table::fread(list.files(pattern = "ProteinWeight"))[grepl(CtlGroup, Condition) | grepl(ExpGroup, Condition)] |> data.table::setnames(c("Protein.Group", "weights"), c("ProteinGroup", "Weight"), skip_absent = T)
      ProtWeights |> data.table::setnames(c("label", colnames(ProtWeights)[grepl("Protein.*group", colnames(ProtWeights), ignore.case = T)], colnames(ProtWeights)[grepl("weight", colnames(ProtWeights), ignore.case = T)],
                                            colnames(ProtWeights)[grepl("sum.*int", colnames(ProtWeights), ignore.case = T)], colnames(ProtWeights)[grepl("mean.*sum", colnames(ProtWeights), ignore.case = T)]), 
                                          c("Label", "ProteinGroup", "Weight", "SumIntensities", "MeanSum"), skip_absent = T)
      ProtWeights <- ProtWeights[Replicate %in% ReplicatesUsed]
      
      
      if(length(unique(grepl("Cluster", colnames(ProtWeights)))) == 1){
        ProtWeights[, Cluster := gsub("(.*)_\\d+h", "\\1", Condition)]
      }
      if(length(unique(grepl("Time", colnames(ProtWeights)))) == 1){
        ProtWeights[, Time := as.numeric(gsub(".*_(\\d+)h", "\\1", Condition))]
      }
      
      ProtWeights[, Label := gsub("L$", "Light", Label)]
      ProtWeights[, Label := gsub("H$", "Heavy", Label)]
      ProtWeights <- ProtWeights[, .(Run, Sample, Label, ProteinGroup, SumIntensities, MeanSum, Weight)] |> data.table::merge.data.table(Metadata)
      
      if (dir.exists(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output")) == TRUE) {
        unlink(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"), recursive = TRUE)
      }
      dir.create(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"), showWarnings = TRUE)
      setwd(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      data.table::fwrite(Metadata, file = "Sample_Metadata.csv")
    }
    # Analyse Light Channel
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Isolating Light Channel Proteins"))
    {
      ProtLFQsInput_L <- ProtLFQsInput |> data.table::dcast(ProteinGroup ~ Sample, value.var = "LFQ_L") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
      ProtWeights_L <- ProtWeights[Label == "Light"] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Weight") |> 
        tibble::column_to_rownames("ProteinGroup") |> as.matrix()
      
      ProteinLFQs_L <- ProtLFQsInput_L |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
      ProteinLFQs_L <- data.table::melt.data.table(ProteinLFQs_L, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |> 
        data.table::merge.data.table(Metadata)
      ProteinLFQs_L[, N_Values := sum(is.finite(Abundance)), by = .(ProteinGroup, Cluster, Time)]
      ProteinLFQs_L[, MeanAbundance := mean(Abundance, na.rm = T), by = .(ProteinGroup)]
      ProteinLFQs_L[, NormAbundance := Abundance - MeanAbundance]
      ProteinLFQs_L[, Sum_Values := sum(is.finite(Abundance)), by = .(ProteinGroup, Time)]
      ProteinLFQs_L[, Diff_Detected := Sum_Values - N_Values]
      ProteinLFQs_L[, Time := factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels))]
      
      LightAbunBoxplot <- ProteinLFQs_L |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")), 
                                                                        y = Abundance, colour = Cluster)) +
        ggplot2::geom_boxplot(outliers = F) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::labs(x = "Time", y = "Light Protein Abundance") + 
        ggplot2::theme(legend.title = ggplot2::element_blank())
      
      Relative_Time <- function(CurrentTime, Times = gsub("h", "", TimeLevels), Shift = -1){
        NextPos = which(Times == unique(CurrentTime)) + Shift
        if(data.table::between(NextPos,1, length(TimeLevels))){return(Times[NextPos])}else{return(NA_character_)}
      }
      ProteinLFQs_L[, NextTime := Relative_Time(Time, Times = gsub("h", "", TimeLevels), Shift = 1), by = .(Time)]
      NextTime <- ProteinLFQs_L[,.(NextTimeSamples = mean(N_Values)),by = .(Time, ProteinGroup, Cluster)] |> data.table::setnames("Time","NextTime")
      ProteinLFQs_L <- ProteinLFQs_L |> data.table::merge.data.table(NextTime, by = c("NextTime","ProteinGroup","Cluster"), all.x = T)
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Imputing Light Abundances & Weights"))
    {
      ProteinLFQs_L[NextTimeSamples == 0 & N_Values == 0 & Diff_Detected == 3, Impute :=  T]
      ProteinLFQs_L[is.na(Impute), Impute := F]
      
      ModelProteins_L <-  imputeLCMD::impute.MinProb(log2(ProtLFQsInput_L), q = 0.01) |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> 
        data.table::data.table() |> data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Sample", value.name = "ImpAbundance")
      ModelProteins_L[, ImpAbundance := ImpAbundance-1]
      ModelProteins_L <- ProteinLFQs_L |> data.table::merge.data.table(ModelProteins_L, by = c("ProteinGroup", "Sample"))
      ModelProteins_L[Impute == T, Abundance := 2^(ImpAbundance)]
      ModelProteins_L[, ImpAbundance := NULL]
      
      AbundanceMatrix_L <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Abundance') |> 
        tibble::column_to_rownames('ProteinGroup') |> as.matrix()
      AbundanceMatrix_L <- AbundanceMatrix_L[matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) <= 1-LightMinSamples,]
      
      ProtWeights_L_Imp <- ProtWeights_L[rownames(AbundanceMatrix_L), colnames(AbundanceMatrix_L)]
      
      ImputationMatrix <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Impute') |> 
        tibble::column_to_rownames('ProteinGroup') |> as.matrix()
      ImputationMatrix <- ImputationMatrix[rownames(ProtWeights_L_Imp), colnames(ProtWeights_L_Imp)]
      
      ProtWeights_L_Imp[ImputationMatrix==T] <- (ProtWeights_L_Imp |> min(na.rm = T))*0.1
      
      N_NAs <- matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) |> tibble::enframe(name = "ProteinGroup", value = "Prop_NA") |> data.table::data.table()
      
      # Plot PCA on Processed Light Data
      All_PCA <- log(AbundanceMatrix_L) |> data.frame() |> tidyr::drop_na() |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
      
      LightProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Processed Light Data")
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Modelling Light Data"))
    {
      Targets <- Metadata[, .(Sample, Cluster, Time)][ , Time := as.numeric(gsub("h", "", Time))]
      Targets <- Targets[match(ProtWeights_L_Imp |> colnames(), Sample)]
      Time <- Targets$Time
      Cluster <- factor(Targets$Cluster)
      if(SameInitialAbundance == T){design <- model.matrix(~Time + Time:Cluster)} else {design <- model.matrix(~Time + Cluster + Time:Cluster)}
      colnames(design) <- gsub(paste0("Cluster", ExpGroup), "Cluster_Exp", colnames(design))
      colnames(design) <- gsub(paste0("Cluster", CtlGroup), "Cluster_Ctl", colnames(design))
      
      LimmaOutput <- limma::eBayes(limma::lmFit(Biobase::ExpressionSet(assayData = log(AbundanceMatrix_L)), design, method = "robust", 
                                                weights = ProtWeights_L_Imp))
      Limma_Slopes <- LimmaOutput$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup")
      Limma_Slopes_SD <- LimmaOutput$stdev.unscaled |> data.table::data.table(keep.rownames = T) |> data.table::copy() 
      Limma_Slopes_SD <- Limma_Slopes_SD |> data.table::setnames(c("rn", "(Intercept)", "Time", "Cluster_Exp", "Time:Cluster_Exp"), 
                                                                 c("ProteinGroup", "(Intercept)_SD", "Time_SD", "Cluster_Exp_SD", "Time:Cluster_Exp_SD"), skip_absent = T)
      Limma_Slopes <- Limma_Slopes |> data.table::merge.data.table(Limma_Slopes_SD, by = "ProteinGroup")
      Limma_Slopes$Fvalue = LimmaOutput$F
      
      MeanVarPlot <- data.table::data.table(Mean = LimmaOutput$Amean, Variance = sqrt(LimmaOutput$sigma)) |> 
        ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) + ggplot2::geom_point(stroke = NA) +
        ggplot2::labs(x = "Mean Log Light Protein LFQ", y = "Light Protein LFQ Variance")
      
      LightModelParameters <- limma::topTable(LimmaOutput, "Time:Cluster_Exp", number = nrow(ProtWeights_L_Imp)) |> 
        data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup") |> 
        data.table::merge.data.table(N_NAs) |> data.table::merge.data.table(Limma_Slopes) |> 
        data.table::merge.data.table(ProteinInfo, all.x = T) |> data.table::setnames("logFC", "Difference")
      LightParameters <- LightModelParameters |> data.table::copy()
      LightParameters[, `:=`(Ctl_Value = -Time, Exp_Value = -`Time:Cluster_Exp` - Time, Difference = -Difference, Parameter = "KlossL")]
      Kloss_Offset <- abs(min(LightParameters[, .(Ctl_Value, Exp_Value)]))*1.01
      LightParameters[, `:=`(Exp_Value = Exp_Value + Kloss_Offset, Ctl_Value = Ctl_Value + Kloss_Offset)]
      LightParameters[, FC := (Exp_Value + Kloss_Offset)/(Ctl_Value + Kloss_Offset)]
      LightParameters[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
      LightParameters[, Significance := data.table::fifelse(P.Value < 0.05 & Difference < 0, "Sig. Decrease", 
                                                            data.table::fifelse(P.Value < 0.05 & Difference > 0, "Sig. Increase", ""))]
      #LightParameters$Isoforms <- 1
      #for (i in 1:nrow(LightParameters)) {
      #  if (length(stringr::str_extract_all(LightParameters$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
      #    LightParameters$Isoforms[i] <- paste0(stringr::str_extract_all(LightParameters$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
      #  } else {LightParameters$Isoforms[i] <- 1}
      #}
      #LightParameters[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"),function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
      LightParameters[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      LightParameters[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
      LightParameters <- LightParameters[, .(ProteinGroup, #ProteinDescription, 
                                             Gene, P.Value, adj.P.Val, Prop_NA, Parameter, Ctl_Value, Exp_Value, FC, Log2FC, Difference)]
      if(SameInitialAbundance == T){
        LightModelledData <- ProtLFQsInput[, ProteinGroup] |> tidyr::crossing(Metadata[, .(Condition, Cluster, Time)]) |> 
          data.table::setnames(c("ProteinGroup", "Condition", "Cluster", "TimeVar")) |> data.table::data.table() |> 
          data.table::merge.data.table(LightModelParameters[, .(ProteinGroup, `(Intercept)`, Time, `Time:Cluster_Exp`)] |> data.table::setnames("Time", "TimeCoeff"))
        LightModelledData[, `:=`(ExpGroup = data.table::fifelse(grepl(ExpGroup, Cluster), 1, 0), TimeVar = as.numeric(paste(TimeVar)))]
        LightModelledData[, Abundance := exp(`(Intercept)` + TimeVar*Time + TimeVar*ExpGroup*`Time:Cluster_Exp`)] 
      } else {
        LightModelledData <- ProtLFQsInput[, ProteinGroup] |> tidyr::crossing(Metadata[, .(Condition, Cluster, Time)]) |> 
          data.table::setnames(c("ProteinGroup", "Condition", "Cluster", "TimeVar")) |> data.table::data.table() |> 
          data.table::merge.data.table(LightModelParameters[, .(ProteinGroup, `(Intercept)`, Time, Cluster_Exp, `Time:Cluster_Exp`)])
        LightModelledData[, `:=`(ExpGroup = data.table::fifelse(grepl(ExpGroup, Cluster), 1, 0), TimeVar = as.numeric(paste(TimeVar)))]
        LightModelledData[, Abundance := exp(`(Intercept)` + ExpGroup*Cluster_Exp + TimeVar*Time + TimeVar*ExpGroup*`Time:Cluster_Exp`)]
      }
    }
    if(GenerateDataPlots == T){
      message(paste0(ExpGroup, " vs. ", CtlGroup, ": Generating Modelled Light Data Plots"))
      {
        setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
        if (dir.exists("LightPlots") == TRUE) {
          unlink("LightPlots", recursive = TRUE)
        }
        dir.create("LightPlots", showWarnings = TRUE)
        setwd("LightPlots")
        
        for(POI in unique(LightModelledData$ProteinGroup)){
          POIModelData <- LightModelledData[ProteinGroup == POI]
          POIModelData$Isoforms <- 1
          for (i in 1:nrow(POIModelData)) {
            if (length(stringr::str_extract_all(POIModelData$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
              POIModelData$Isoforms[i] <- paste0(stringr::str_extract_all(POIModelData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
            } else {POIModelData$Isoforms[i] <- 1}
          }
          POIModelData[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"),function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
          
          POIActualData <- ModelProteins_L[ProteinGroup == POI]
          POIActualData$Isoforms <- 1
          for (i in 1:nrow(POIActualData)) {
            if (length(stringr::str_extract_all(POIActualData$ProteinGroup[i], "-\\d", simplify = T)) > 0) {
              POIActualData$Isoforms[i] <- paste0(stringr::str_extract_all(POIActualData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
            } else {POIActualData$Isoforms[i] <- 1}
          }
          POIActualData[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"),function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
          
          pdf(paste0(POI, "_LightPlot.pdf"), width = 16, height = 12)
          print(POIModelData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(paste(TimeVar)), y = log(Abundance), colour = Cluster)) + 
                  ggplot2::geom_line(ggplot2::aes(group = Cluster), linetype = "dashed") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::geom_point(data = POIActualData, ggplot2::aes(x = as.numeric(paste(Time)), y = log(Abundance))) +
                  ggplot2::facet_wrap(~paste0(ProteinGroup, " (", ProteinInfo[ProteinGroup == POI, Gene], ")"), scales = 'free_y') + 
                  ggplot2::annotate('text', x = quantile(TimeLevels, 0.75), y = log(max(ModelProteins_L[ProteinGroup == POI, Abundance], na.rm = T))*0.99, 
                                    label = glue::glue('{CtlGroup} kloss = {round(LightParameters[ProteinGroup == POI, Ctl_Value],3)}\n
                                                         {ExpGroup} kloss = {round(LightParameters[ProteinGroup == POI, Exp_Value], 3)}')) +
                  ggplot2::labs(x = "Time (hours)", y = "Log Light Protein Abundance") +
                  ggplot2::theme(legend.title = ggplot2::element_blank()))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Processing Modelled Light Data"))
    {
      LightModelledData$Isoforms <- 1
      for(i in 1:nrow(LightModelledData)){
        if(length(stringr::str_extract_all(LightModelledData$ProteinGroup[i], "-\\d", simplify = T)) > 0){
          LightModelledData$Isoforms[i] <- paste0(stringr::str_extract_all(LightModelledData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
        } else {LightModelledData$Isoforms[i] <- 1}
      }      
      LightModelledData[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      LightModelledData[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      LightModelledData[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
    }
    # Analyse Heavy Channel
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Isolating Heavy Channel Proteins")) 
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      ProtLFQsInput_H <- ProtLFQsInput |> data.table::dcast(ProteinGroup ~ Sample,  value.var =  "LFQ_H") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
      ProteinLFQs_H <- ProtLFQsInput_H |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
      ProteinLFQs_H <- data.table::melt.data.table(ProteinLFQs_H, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |>
        data.table::merge.data.table(Metadata)
      ProtWeights_H <- ProtWeights[Label =="Heavy"][,`:=`(Replicate = factor(Replicate), ProteinGroup = factor(ProteinGroup), 
                                                          Time = factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels)),
                                                          Cluster = factor(Cluster, levels = ClusterLevels))]
      
      ModelProteins_H <- ProteinLFQs_H[,.(ProteinGroup, Condition, Cluster, Time, Replicate, Abundance)] 
      ModelProteins_H[, Replicate := factor(Replicate)]
      ModelProteins_H[Time == 0, Abundance := 0]
      ModelProteins_H <- na.omit(ModelProteins_H)
      ModelProteins_H[, N_Quant := .N, by = .(ProteinGroup, Condition)] # Max. N_Quant is 3 (3 Replicates)
      ModelProteins_H[, N_QuantTotal := .N, by = .(ProteinGroup)] # Max. N_QuantTotal is 36 (3 Replicates, 6 Timepoints, 2 Conditions: 3x6x2 = 36)
      
      HeavyAbunBoxplot <- ModelProteins_H |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")), 
                                                                          y = Abundance, colour = Cluster)) + 
        ggplot2::geom_boxplot(outliers = F) + ggplot2::labs(x = "Time", y = "Heavy Protein Abundance") +
        ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Filtering By Heavy Channel Missingness & Variation")) 
    {
      PreFiltReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Cluster, Time, N_Quant)] |> unique()
      PreFiltTotalCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)] |> unique()
      
      ModelProteins_H <- ModelProteins_H[N_QuantTotal >= HeavyMinSamples*nrow(Metadata)] 
      ModelProteins_H[, CV := sd(Abundance)/mean(Abundance), by = .(ProteinGroup, Cluster, Time)]
      ModelProteins_H_CV <- ModelProteins_H[,.(MeanCV_H = mean(CV, na.rm =T)), by = ProteinGroup]
      ModelProteins_H <- ModelProteins_H[ProteinGroup %in% ModelProteins_H_CV[MeanCV_H <= HeavyMaxCV, ProteinGroup]]
      
      TotalReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Cluster, Time, N_Quant)][, Filtering := "Post-Filtering"] |> unique() |>
        rbind(PreFiltReplicateData[, Filtering := "Pre-Filtering"])
      TotalProteinCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)][, Filtering := "Post-Filtering"] |> unique() |>
        rbind(PreFiltTotalCountData[, Filtering := "Pre-Filtering"])
      TotalReplicateData[, Filtering := factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))]
      
      HeavyCVPlot <- ModelProteins_H_CV |> ggplot2::ggplot(ggplot2::aes(x = MeanCV_H)) + ggplot2::geom_histogram() + 
        ggplot2::labs(x = "Mean Protein CV", y = "No. Proteins") + ggplot2::geom_vline(xintercept = HeavyMaxCV, linetype = "dashed", colour = "red") +
        ggplot2::annotate("rect", xmin = HeavyMaxCV, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.3) + 
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0))
      
      HeavyCompleteness <- TotalProteinCountData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(N_QuantTotal))) + 
        ggplot2::geom_histogram() + ggplot2::labs(x = "No. Samples", y = "No. Proteins") + 
        ggplot2::geom_vline(xintercept = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, linetype = "dashed", colour = "red") +
        ggplot2::facet_wrap(~factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))) + 
        ggplot2::annotate("rect", xmin = -Inf, xmax = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, ymin = -Inf, ymax = Inf, 
                          fill = "red", alpha = 0.3) + 
        ggplot2::scale_x_continuous(limits = c(0, nrow(Metadata)), expand = c(0,0)) + ggplot2::scale_y_continuous(expand = c(0,0))
      
      HeavyCounts <- TotalReplicateData |> ggplot2::ggplot(ggplot2::aes(x = N_Quant, fill = Cluster)) + 
        ggplot2::geom_bar(stat = "count", position = ggplot2::position_dodge()) + 
        ggplot2::facet_grid(ggplot2::vars(Filtering), ggplot2::vars(factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")))) + 
        ggplot2::scale_fill_brewer(palette = "Set1", name = "Condition") + ggplot2::coord_cartesian(ylim = c(0, 5000)) +
        ggplot2::labs(x = "Replicate", y = "No. Proteins") + ggplot2::theme(strip.text.y = ggplot2::element_text(size = 26), 
                                                                            legend.title = ggplot2::element_blank(),
                                                                            legend.position = "inside", 
                                                                            legend.position.inside = c(0.2, 0.9))
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Filtering By Heavy Channel Monotonicity")) 
    {
      MonotonicityData <- ModelProteins_H[,.(Mean_Abundance = mean(Abundance, na.rm = T)), by = .(ProteinGroup, Time, Cluster)][(order(ProteinGroup, Cluster, Time))]
      MonotonicityData[, CumSum := cummax(Mean_Abundance), by = .(ProteinGroup, Cluster)]
      MonotonicityData[, Monotonic := CumSum == Mean_Abundance]
      MonotonicitySummary <- MonotonicityData[, .(N_Monotonic = sum(Monotonic)), by = "ProteinGroup"]
      MonotonicitySummary[, MonotonicityProp := (N_Monotonic/length(ConditionLevels))]
      #data.table::fwrite(MonotonicityData, file = paste0(ExpGroup, "_vs_", CtlGroup, "_HeavyMonotonicity.csv"))
      
      MonotonicityPlot <- MonotonicitySummary |> ggplot2::ggplot(ggplot2::aes(x = MonotonicityProp)) + 
        ggplot2::geom_histogram() + ggplot2::labs(x = "Prop. Monotonicity", y = "No. Proteins") + 
        ggplot2::geom_vline(xintercept = HeavyMinMonotonicity-0.033, linetype = "dashed", colour = "red") +
        ggplot2::annotate("rect", xmin = -Inf, xmax = HeavyMinMonotonicity-0.033, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.3) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(expand = c(0, 0))
      
      ModelProteins_H <- ModelProteins_H[ProteinGroup %in% MonotonicitySummary[MonotonicityProp >= HeavyMinMonotonicity, ProteinGroup]]
      ModelProteins_H[, ProteinGroup := as.factor(ProteinGroup)]
      
      HeavyMissingnessData <- ModelProteins_H[, .(ProteinGroup, N_QuantTotal)] |> data.table::copy() |> unique()
      HeavyMissingnessData[, PropNAs := N_QuantTotal/nrow(Metadata)]
      
      # Processed Heavy Data PCA
      All_PCA <- ModelProteins_H[, .(ProteinGroup, Condition, Replicate, Abundance)] |> data.table::merge.data.table(Metadata[, .(Condition, Replicate, Sample)]) |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Abundance, names_from = Sample, values_fill = NA) |> 
        tidyr::drop_na() |> data.frame(row.names = c("ProteinGroup")) |> t() |> stats::prcomp(scale. = TRUE)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
      
      HeavyProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Cluster, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_brewer(palette = "Set1") + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Processed Heavy Data")
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Running ", HeavyModel, " Modelling")) 
    {
      KlossParameters <- LightModelParameters[, .(ProteinGroup, Time, `Time:Cluster_Exp`)] |> data.table::copy()
      KlossParameters[, `:=`(Kloss_Ctl = -Time, Kloss_Exp = -`Time:Cluster_Exp` - Time)]
      Run_ProteinNLS <- function(POI){
        tryCatch(
          expr = {
            POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProtWeights_H[, .(ProteinGroup, Condition, Replicate, Weight, Sample)], all.x = T)
            POIData <- POIData[order(Cluster, Time)]
            POIData[,`:=`(Comparison = data.table::fifelse(Cluster == CtlGroup, 0, 1))]
            POIData[, Time := as.numeric(paste(Time))]
            POIData[, VAR := var(Abundance), by = .(Cluster, Time)]
            POIData[is.na(VAR), VAR := max(POIData$VAR, na.rm  = T)]
            # Define T0 Data
            T0Data <- POIData[, head(.SD,3), by = Cluster]
            T0Data[,`:=`(Time = 0, Abundance = 0, Weight = min(POIData$Weight, na.rm =T), VAR = min(POIData$VAR, na.rm = T))]
            POIData <- POIData[Time != 0] |> rbind(T0Data)
            
            T0Data <- POIData[,head(.SD,1), by = Cluster]
            T0Data[,`:=`(Time = 0, Abundance = 0)]
            POIData <- POIData |> rbind(T0Data)
            
            KsynStart <- sapply(ClusterLevels, function(COI){
              mean(POIData[order(Time)][Time != 0 & Cluster == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Cluster == COI][,head(.SD,2)]$Time)), na.rm =T)
            })
            
            AbundancePlateau <- sapply(ClusterLevels, function(COI){
              mean(POIData[order(Time)][Cluster == COI][,tail(.SD,3)]$Abundance,na.rm =T)
            })
            
            StartVals <- c(KsynStart, KsynStart / AbundancePlateau) # nls requires start values of approx. params
            # nlsLM more stable than nls
            POIFit <- minpack.lm::nlsLM(Abundance ~ (Ksyn_Ctl+(Ksyn_Exp*Comparison))/(Kloss_Ctl+(Kloss_Exp*Comparison)) * (1-exp(-(Kloss_Ctl+(Kloss_Exp*Comparison))*Time)), data = POIData,
                                        start = list(Ksyn_Ctl = StartVals[[1]], 
                                                     Ksyn_Exp = StartVals[[2]]/10, 
                                                     Kloss_Ctl = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Ctl], StartVals[[3]]), 
                                                     Kloss_Exp = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Exp], StartVals[[4]]/10)),
                                        weight = 1/POIData$Weight, control = nls.control(maxiter = 200, warnOnly = TRUE))
            POIData <- POIData[, Fitted := predict(POIFit)]
            
            FitSummary <- summary(POIFit)
            if(GenerateDataPlots == T){
              pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
              print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Cluster)) + ggplot2::geom_point() + 
                      ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), by = .(Cluster, Time, ProteinGroup)], 
                                         ggplot2::aes(y = Fitted, group = interaction(Cluster)), linetype = "dashed") + 
                      ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
                      ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", ProteinInfo[ProteinGroup == POI, Gene], ")"),
                                    subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients[3,1],3)}, ksyn = {round(FitSummary$coefficients[1,1],3)}
                                                            {ExpGroup}: kloss = {round(FitSummary$coefficients[3,1]+ FitSummary$coefficients[4,1],3)}, ksyn = {round(FitSummary$coefficients[1,1]+ FitSummary$coefficients[2,1],3)}')) + 
                      ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
              Proteopedia::Reset_Dev()
            }
            #message(POI, ": Model Applied")
            return(list(POIModel = POIFit, FittedData = POIData))
          },
          error = function(e){
            message(POI, ": Error Caught")
          })
      }
      Run_ProteinNLME <- function(POI){
        tryCatch(
          expr = {
            POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProtWeights_H[, .(ProteinGroup, Cluster, Replicate, Time, Weight)], all.x = T)
            POIData[, Acquisition := data.table::fifelse(as.numeric(paste(Time)) > max(TimeLevels)/2, "Late", "Early")]
            POIData[is.na(Acquisition), Acquisition := "Early"]
            POIData[, Time := as.numeric(paste(Time))]
            POIData[Time == 0, Weight := max(POIData$Weight, na.rm = T)]
            # Define T0 Data
            T0Data <- POIData[,head(.SD,1), by = Cluster]
            T0Data[,`:=`(Time = 0, Abundance = 0)]
            POIData <- POIData[Time != 0] |> rbind(T0Data)
            
            KsynStart <- sapply(ClusterLevels, function(COI){
              mean(POIData[order(Time)][Time != 0 & Cluster == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Cluster == COI][,head(.SD,2)]$Time)), na.rm =T)
            })
            
            AbundancePlateau <- sapply(ClusterLevels, function(COI){
              mean(POIData[order(Time)][Cluster == COI][,tail(.SD,3)]$Abundance,na.rm =T)
            })
            
            StartVals <- c(KsynStart, KsynStart / AbundancePlateau)
            
            POIFit <- nlme::nlme(Abundance ~ (Ksyn/Kloss)*(1 - exp(-Kloss*as.numeric(paste(Time)))), data = POIData, fixed = list(Ksyn ~ Cluster, Kloss ~ Cluster),
                                 random = list(Acquisition = nlme::pdDiag(Ksyn+Kloss ~1)), 
                                 start = StartVals, weights = nlme::varFixed(~Weight))
            POIData <- POIData[, Fitted := predict(POIFit)]
            
            FitSummary <- summary(POIFit)
            if(GenerateDataPlots == T){
              pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
              print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Cluster)) + ggplot2::geom_point() + 
                      ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), by = .(Cluster, Time, ProteinGroup)], 
                                         ggplot2::aes(y = Fitted, group = interaction(Cluster)), linetype = "dashed") + 
                      ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
                      ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", ProteinInfo[ProteinGroup == POI, Gene], ")"),
                                    subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients$fixed[3],3)}, ksyn = {round(FitSummary$coefficients$fixed[1],3)}
                                                              {ExpGroup}: kloss = {round(FitSummary$coefficients$fixed[3] + FitSummary$coefficients$fixed[4],3)}, ksyn = {round(FitSummary$coefficients$fixed[1] + FitSummary$coefficients$fixed[2],3)}')) + 
                      ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
              Proteopedia::Reset_Dev()
            }
            #message(paste0(POI, ": Model Applied"))
            return(list(POIModel = POIFit, FittedData = POIData))
          },
          error = function(e){
            message(paste0(POI, ": Error Caught"))
          }
        )
      }
      
      if(GenerateDataPlots == T){
        if (dir.exists("HeavyPlots") == TRUE) {
          unlink("HeavyPlots", recursive = TRUE)
        }
        dir.create("HeavyPlots", showWarnings = TRUE)
        setwd("HeavyPlots")
      }
      
      HeavyModelledData <- data.table::data.table()
      HeavyModelParameters <- data.table::data.table()
      ProteinSigmas <- data.table::data.table()
      if(HeavyModel == "NLS"){
        for(POI in levels(ModelProteins_H$ProteinGroup)){
          NLSOutput <- Run_ProteinNLS(POI)
          HeavyModelledData <- HeavyModelledData |> rbind(NLSOutput$FittedData)
          if(!is.null(NLSOutput$POIModel$convInfo$isConv)){
            if(NLSOutput$POIModel$convInfo$isConv == T){
              POIModelSummary <- summary(NLSOutput$POIModel)$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "t_Value", "P.Value"))
              ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = sigma(NLSOutput$POIModel)))
              HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
            }
          }
        }
      }
      if(HeavyModel == "NLME"){
        for(POI in levels(ModelProteins_H$ProteinGroup)){
          NLMEOutput <- Run_ProteinNLME(POI)
          HeavyModelledData <- HeavyModelledData |> rbind(NLMEOutput$FittedData)
          if(class(NLMEOutput$POIModel)[1] == "nlme"){
            POIModelSummary <- summary(NLMEOutput$POIModel)$tTable |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "DF", "t_Value", "P.Value"))
            ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = NLMEOutput$POIModel$sigma))
            HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
            HeavyModelParameters[, Effect := gsub(".Cluster.*", "_Exp", gsub("..(Intercept.)", "_Ctl", HeavyModelParameters$Effect))]
          }
        }
      }
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      HeavyModelParameters[, adj.P.Val := p.adjust(P.Value, 'BH'), by = 'Effect']
      HeavyModelParameters <- HeavyModelParameters |> data.table::merge.data.table(ProteinInfo) |> 
        data.table::merge.data.table(HeavyMissingnessData)
      
      HeavyParameters <- HeavyModelParameters |> data.table::copy()
      HeavyParameters[stringr::str_detect(Effect,'_Ctl$'), Ctl_Value := Estimate]
      HeavyParameters[, Parameter := gsub("Kloss", "KlossH", stringr::str_remove(Effect,'_.*'))]
      HeavyParameters[, Ctl_Value := mean(Ctl_Value, na.rm = T), by = .(ProteinGroup, Parameter)]
      HeavyParameters <- HeavyParameters[stringr::str_detect(Effect,'_Exp')]
      HeavyParameters[, Exp_Value := Estimate + Ctl_Value]
      KlossH_Offset <- abs(min(HeavyParameters[, .(Ctl_Value, Exp_Value)]))*1.01
      HeavyParameters[stringr::str_detect(Effect,'Kloss'), `:=`(Ctl_Value = Ctl_Value + KlossH_Offset, Exp_Value = Exp_Value + KlossH_Offset)]
      HeavyParameters[, Difference := Exp_Value - Ctl_Value]
      HeavyParameters[, FC := Exp_Value/Ctl_Value]
      HeavyParameters[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
      HeavyParameters <- HeavyParameters[, .(ProteinGroup, Gene, Parameter, Ctl_Value, Exp_Value, Difference, FC, Log2FC, P.Value, adj.P.Val, PropNAs)]
      HeavyParameters[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", 
                                                            data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", ""))]
      #HeavyParameters$Isoforms <- 1
      #for(i in 1:nrow(HeavyParameters)){
      #  if(length(stringr::str_extract_all(HeavyParameters$ProteinGroup[i], "-\\d", simplify = T)) > 0){
      #    HeavyParameters$Isoforms[i] <- paste0(stringr::str_extract_all(HeavyParameters$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
      #  } else {HeavyParameters$Isoforms[i] <- 1}
      #}
      #HeavyParameters[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      HeavyParameters[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      HeavyParameters[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
    }
    # Analyse Abundance Data 
    if(SameInitialAbundance == F){
      message(paste0(ExpGroup, " vs. ", CtlGroup, ": Compiling Abundance Data"))
      {
        AbundanceData <- merge(limma::topTable(LimmaOutput, "Cluster_Exp", number = nrow(ProtWeights_L_Imp)) |> as.data.frame() |> 
                                 tibble::rownames_to_column("ProteinGroup"), N_NAs) |> data.table::data.table() |> 
          merge(Limma_Slopes, by = "ProteinGroup") |> data.table::data.table() |> 
          data.table::merge.data.table(ProteinInfo, all.x = T) |> data.table::setnames("logFC", "logFC") |>
          data.table::merge.data.table(T0ModelledData <- LightModelledData[TimeVar == 0, .(ProteinGroup, Isoforms, Cluster, Abundance)] |> 
                                         data.table::dcast(ProteinGroup+Isoforms ~ Cluster, value.var = "Abundance") |> 
                                         data.table::setnames(c(paste(CtlGroup), paste(ExpGroup)), 
                                                              c("Ctl_Abundance", "Exp_Abundance")), all.x = T)
        AbundanceData[, FC := Exp_Abundance/Ctl_Abundance]
        AbundanceData[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
        AbundanceData[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", 
                                                            data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", ""))]
      }
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Data Files"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      data.table::fwrite(ProteinLFQs_L[, .(Sample, Condition, Cluster, Cell, Drug, Conc, Time, Replicate, Abundance, MeanAbundance, NormAbundance)], "LightInputLFQs.csv")
      data.table::fwrite(LightModelParameters, "LightModelOutput.csv")
      data.table::fwrite(LightParameters, "LightParameters.csv")
      data.table::fwrite(LightModelledData, "LightModelledData.csv")
      if(SameInitialAbundance == F){
        data.table::fwrite(AbundanceData, "AbundanceData.csv")
      }
      data.table::fwrite(ModelProteins_H, "HeavyInputLFQs.csv")
      data.table::fwrite(HeavyModelParameters, "HeavyModelOutput.csv")
      data.table::fwrite(HeavyParameters, "HeavyParameters.csv")
      data.table::fwrite(HeavyModelledData, "HeavyModelledData.csv")
      
      AnalysisSummary <- data.table::data.table(Parameter = c("Light kloss", "Heavy kloss", "ksyn", "0hr Abundance"), 
                                                N_Proteins = c(nrow(LightParameters), nrow(HeavyParameters[Parameter == "KlossH"]) , 
                                                               nrow(HeavyParameters[Parameter == "Ksyn"]), nrow(AbundanceData)),
                                                Model = c("Limma", HeavyModel, HeavyModel, "Limma"))
      AnalysisSummary[, Prop_Proteins := N_Proteins/length(ProtLFQsInput[, ProteinGroup] |> unique())]
      data.table::fwrite(AnalysisSummary, "Analysis_Summary.csv")
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Output Plots"))
    {
      pdf("LightOutputPlots.pdf", width = 16, height = 10)
      print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = Prop_NA)) + 
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() + 
              ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "black", max.overlaps = 10) + 
              Proteopedia::Add_KlossAxes())
      print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = Prop_NA)) + 
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "black", max.overlaps = 10) + 
              Proteopedia::Add_KlossAxes(scale = "Log2FC"))
      print(LightParameters[, .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |> 
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |> 
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Cluster", value.name = "KlossL") |> 
              ggplot2::ggplot(ggplot2::aes(x = KlossL, fill = Cluster)) + ggplot2::geom_density(alpha = 0.7) + 
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = "Density of Proteins") +
              ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside", legend.position.inside = c(0.8, 0.8)))
      print(LightParameters[, .(ProteinGroup, P.Value, adj.P.Val)] |> data.table::copy() |> 
              data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |> 
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Adjustment", value.name = "P") |> 
              ggplot2::ggplot(ggplot2::aes(x = P)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
              ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
              ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
      Proteopedia::Reset_Dev()
      
      if(SameInitialAbundance == F){
        pdf("AbundanceOutputPlots.pdf", width = 16, height = 10)
        print(AbundanceData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = Prop_NA)) + 
                ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() + 
                Proteopedia::Add_AbundanceAxes() + ggrepel::geom_text_repel(data = AbundanceData[P.Value < 0.05], colour = "black"))
        print(AbundanceData[, .(ProteinGroup, P.Value, adj.P.Val)] |> data.table::copy() |> 
                data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |> 
                data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Adjustment", value.name = "P") |> 
                ggplot2::ggplot(ggplot2::aes(x = P)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
                ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
                ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
        Proteopedia::Reset_Dev()
      }
      
      pdf("HeavyOutputPlots.pdf", width = 16, height = 10)
      print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) + 
              ggplot2::geom_point() + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05]) + 
              Proteopedia::Add_KlossAxes())
      print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) + 
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05], colour = "black", max.overlaps = 10) + 
              Proteopedia::Add_KlossAxes(scale = "Log2FC"))
      print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) + 
              ggplot2::geom_point() + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05]) + 
              Proteopedia::Add_KsynAxes())
      print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) + 
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05], colour = "black", max.overlaps = 10) + 
              Proteopedia::Add_KsynAxes(scale = "Log2FC"))
      print(HeavyParameters[Parameter == "KlossH", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |> 
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |> 
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Cluster", value.name = "KlossH") |> 
              ggplot2::ggplot(ggplot2::aes(x = KlossH, fill = Cluster)) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)"), y = "Density of Proteins") +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                        legend.position.inside = c(0.8, 0.8)))
      print(HeavyParameters[Parameter == "Ksyn", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |> 
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |> 
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Cluster", value.name = "Ksyn") |> 
              ggplot2::ggplot(ggplot2::aes(x = Ksyn, fill = Cluster)) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~") (Heavy Channel)"), y = "Density of Proteins") +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                        legend.position.inside = c(0.8, 0.8)))
      print(HeavyParameters[, .(ProteinGroup, P.Value, adj.P.Val, Parameter)] |> data.table::copy() |> 
              data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |> 
              data.table::melt.data.table(id.vars = c("ProteinGroup", "Parameter"), variable.name = "Adjustment", value.name = "P") |> 
              ggplot2::ggplot(ggplot2::aes(x = P, fill = Parameter)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
              ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0))+ 
              ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette) +
              ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
      Proteopedia::Reset_Dev()
    }
    message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting QC Plots & Parameter Analysis"))
    {
      pdf("QC_Plots.pdf", width = 18, height = 20)
      print(patchwork::free(LightAbunBoxplot + HeavyAbunBoxplot + patchwork::plot_layout(guides = "collect")) +
              HeavyCVPlot + MonotonicityPlot + patchwork::free(HeavyCounts) + HeavyCompleteness + 
              MeanVarPlot + patchwork::plot_layout(design = "AAAAAAAAAA\nBBCCCDDDDD\nEEEEEDDDDD\nFFFFFDDDDD") +
              patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "" ,"H", ""))))
      Proteopedia::Reset_Dev()
      
      CorrData <- HeavyParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)] |>
        rbind(LightParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)]) |>
        data.table::melt.data.table(id.vars = c("ProteinGroup", "Gene", "Parameter"), 
                                    variable.name = "Cluster", value.name = "Measure")
      CorrData[, Cluster := data.table::fifelse(grepl("Ctl", Cluster), CtlGroup, ExpGroup)]
      CorrData <- CorrData |> data.table::dcast(ProteinGroup+Gene+Cluster ~ Parameter, value.var = "Measure")
      
      pdf("CorrelationPlots.pdf", width = 16, height = 12)
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = KlossL, y = KlossH, colour = Cluster)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Rsq(T) + Proteopedia::Add_XYLine("grey") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + 
              ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
              ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossH, colour = Cluster)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Rsq(T) + Proteopedia::Add_XYLine("grey") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
              ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossL, colour = Cluster)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Rsq(T) + Proteopedia::Add_XYLine("grey") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) +
              ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Light Channel)")) +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
      Proteopedia::Reset_Dev()
    }
  }
  Proteopedia::End_Timer(start.time)
}
########### MS Analysis: TPP DIA-NN To Protein Data Functions ####################################################################################################################################
#' @export
Process_TPP_DIANN <- function(InputDirectory, CtlGroup, SILAC, ProteotypicFiltering = F){
  start.time <- Sys.time()
  set.seed(123)
  if(SILAC == F){
    message("Importing DIA-NN Report File")
    {
      setwd(InputDirectory)
      if(length(list.files(pattern = "report.tsv")) > 0){
        InputFile <- list.files(pattern = "report.tsv")[1]
        PrecursorData <- data.table::fread(InputFile)[, .(Run, Protein.Group, Protein.Ids, First.Protein.Description, Genes, 
                                                          Stripped.Sequence, Precursor.Id, Proteotypic, Precursor.Normalised,    
                                                          Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value, 
                                                          Lib.PG.Q.Value)] |> data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), 
                                                                                                   c("ProteinGroup", "Gene", "ProteinDescription"))
      } else {return(message("ERROR: No Input File Found"))}
      
      if(file.exists("Sample_Rename.csv") == TRUE){
        Sample_Metadata <- data.table::fread("Sample_Rename.csv") 
        PrecursorData$Sample <- Sample_Metadata$Renamed[match(unlist(PrecursorData$Run), Sample_Metadata$Run)]
        PrecursorData <- PrecursorData[!is.na(Sample)]
      }
      
      if(any(grepl("SILAC-", PrecursorData$Precursor.Id))){
        return(message("ERROR: Label-Free Processing on SILAC Data"))
      }
    }
    message("Defining Metadata")
    {
      PrecursorData[, Drug := gsub("(.*)_(.*)_R(\\d)","\\1", Sample)]
      PrecursorData[, Temp := gsub("(.*)_(.*)_R(\\d)","\\2", Sample)]
      PrecursorData[, Replicate := gsub("(.*)_(.*)_R(\\d)","\\3", Sample)]
      PrecursorData[, Sample := gsub("_0", "", Sample)]
      PrecursorData[, Condition := gsub("(.*)_R\\d", "\\1", Sample)]
      # Order Samples
      Sample_Order <- unique(PrecursorData[, .(Sample, Condition, Replicate)]) |> dplyr::arrange(!grepl(CtlGroup, Condition), Condition, Replicate)
      PrecursorData$Sample <- factor(PrecursorData$Sample, levels = Sample_Order$Sample)
    }
    message("Filtering Precursors")
    {
      if(ProteotypicFiltering == T){PrecursorData <- PrecursorData[Proteotypic >= 1]} else {PrecursorData <- PrecursorData}
      PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01] 
      PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
    }
    message("Compiling & Exporting Data")
    {
      data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
      LFQ <- Proteopedia::Calculate_LFQ(PrecursorData, "LFQ")
      tot_intensities <- PrecursorData[, .(Intensity = sum(Precursor.Normalised)), .(ProteinGroup, Sample)]
      PrecursorCounts <- PrecursorData[, .(N_precursors = data.table::uniqueN(Precursor.Id)), .(ProteinGroup, Sample)]
      proteotypic_counts <- PrecursorData[, .(N_precursors_proteotypic = sum(Proteotypic)), .(ProteinGroup, Sample)]
      annotations <- unique(PrecursorData[, .(ProteinGroup, Sample, Condition, Replicate, ProteinDescription, Gene)])
      ProteinData <- Reduce(Proteopedia::Merge_PrecursorData, list(LFQ, tot_intensities, PrecursorCounts, proteotypic_counts, annotations))
      data.table::fwrite(ProteinData, "LF_DIANN_Output.csv.gz")
    }
    message("Plotting Intensities")
    {
      Intensities_Data <- data.table::rbindlist(list(
        PrecursorData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Precursor.Normalised), Type = "Precursor Quantity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ), Type = "Protein MaxLFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity), Type = "Protein Intensity")]
      ))
      Intensities_Data[, Type := factor(Type, levels = c("Precursor Quantity", "Protein MaxLFQ", "Protein Intensity"))]   # Set plotting order
      
      IntensityPlot <- Intensities_Data |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = `log2 quantity`, colour = Condition)) + 
        ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + 
        ggplot2::facet_wrap("Type", scales = "free_x") + ggplot2::ylab("Log2 Value") + ggplot2::coord_flip() + 
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
    message("Plotting Precursor, Peptide & Protein Counts")
    {
      all_counts <- PrecursorData[, lapply(.SD, data.table::uniqueN), .(Sample, Condition), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")]
      facet_labels <- ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups")) # Create plot labellers for facets
      
      CountPlot <- data.table::melt.data.table(all_counts, id.vars = c("Sample","Condition"), value.name = "IDs") |> 
        ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = IDs/1000, fill = Condition, label = format(IDs, big.mark = ",", scientific = FALSE))) +
        ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.2) +
        ggplot2::facet_wrap("variable", scales = "free_x", labeller = facet_labels) + ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + 
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
    message("Plotting Data Completeness")
    {
      data_completeness <- rbind(Count_Proteins(ProteinData, "All"),
                                 Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                                 Count_Proteins(ProteinData[N_precursors_proteotypic >= 2], "≥ 2 Proteotypic"))
      
      NAsPlot <- data_completeness |> ggplot2::ggplot(ggplot2::aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors))+
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("All" = "black", "≥ 2" = "darkgrey", "≥ 2 Proteotypic" = "orange3")) +
        ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") +
        ggplot2::scale_x_continuous( breaks = seq(1, 1000, 1)) +
        ggplot2::scale_y_continuous( limits = c(0, max(data_completeness$cumulative_protein_N )/1000))+
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       legend.position = "inside", legend.position.inside = c(0.25, 0.25))
    }
    message("Calculating Missed Trypsinisation Sites")
    {
      PrecursorData <- PrecursorData[, MissedCleavage := grepl("[RK][^P]", Stripped.Sequence)]
      PrecursorCount <- PrecursorData[, N := .N, by = Sample]
      TrypsinData <- PrecursorData[MissedCleavage == TRUE, .(N = .N), by = .(Sample, Condition, MissedCleavage)]
      
      TrypsinData$PercentPrecursors <- 0
      for(i in 1:nrow(TrypsinData)){
        TrypsinData$PercentPrecursors[i] <- TrypsinData$N[i]/PrecursorCount$N[which(PrecursorCount$Sample == TrypsinData$Sample[i])]
      }
      
      MissedCleavagePlot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentPrecursors*100, fill = Condition)) + 
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::labs(x = "", y = "Precursors with Missed Tryptic Sites (%) ") + ggplot2::coord_flip() + 
        ggplot2::lims(y = c(0, max(TrypsinData$PercentPrecursors*200)))
    }
    message("Plotting Precursor & Protein Variation")
    {
      precursor_CVs <- PrecursorData[, .(CV = Calculate_CV(Precursor.Normalised), N = .N), .(Precursor.Id, Condition)]
      precursor_CVs <- precursor_CVs[, Rank := data.table::frank(CV), Condition] 
      precursor_CVs$ID <- "Precursors"
      
      protein_CVs <- ProteinData[, .(CV = Calculate_CV(LFQ), N = .N), .(ProteinGroup, Condition)]  
      protein_CVs <- protein_CVs[, Rank := data.table::frank(CV), Condition] 
      protein_CVs$ID <- "Protein Groups"
      
      all_CVs <- precursor_CVs[, Precursor.Id := NULL] |> rbind(protein_CVs[, ProteinGroup := NULL])
      
      VariationPlot <- all_CVs |> ggplot2::ggplot(ggplot2::aes(x = Rank/1000, y = CV, colour = Condition))+
        ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Coeff. of Variation [%]") +
        ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette) + ggplot2::coord_cartesian(ylim = c(0,50)) + 
        ggplot2::facet_wrap(~ID, scales = "free") +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       legend.position = "inside", legend.position.inside = c(0.6, 0.6))
    }
  } else {
    message("Importing DIA-NN Report File")
    {
      setwd(InputDirectory)
      if(length(list.files(pattern = "report.tsv")[1]) > 0){
        InputFile <- list.files(pattern = "report.tsv")[1]
        PrecursorData <- data.table::fread(InputFile)[, .(Run, Protein.Group, Protein.Ids, First.Protein.Description, Genes, 
                                                          Stripped.Sequence, Precursor.Id, Proteotypic, Precursor.Quantity,    
                                                          Precursor.Translated, Channel.Q.Value,
                                                          Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value, 
                                                          Lib.PG.Q.Value)] |> data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), 
                                                                                                   c("ProteinGroup", "Gene", "ProteinDescription"))
      } else {return(message("ERROR: No InputFile Found"))}
      
      if(file.exists("Sample_Rename.csv") == TRUE){
        PrecursorData <- PrecursorData |> data.table::merge.data.table(data.table::fread("Sample_Rename.csv") |> data.table::setnames("Renamed", "Sample"))
      }
    }
    message("Defining Metadata")
    {
      PrecursorData[, Drug := gsub("(.*)_(.*)_R(\\d)","\\1", Sample)]
      PrecursorData[, Temp := gsub("(.*)_(.*)_R(\\d)","\\2", Sample)]
      PrecursorData[, Replicate := gsub("(.*)_(.*)_R(\\d)","\\3", Sample)]
      PrecursorData[, Sample := gsub("_0", "", Sample)]
      PrecursorData[, Condition := gsub("(.*)_R\\d", "\\1", Sample)]
      
      Sample_Order <- dplyr::arrange(unique(PrecursorData[, .(Sample, Condition, Replicate)]), !grepl(CtlGroup, Condition), Condition, Replicate)
      PrecursorData[, Sample := factor(Sample, levels = Sample_Order$Sample)]
    }
    message("Filtering Precursors")
    {
      if(ProteotypicFiltering == T){PrecursorData <- PrecursorData[Proteotypic >= 1]} else {PrecursorData <- PrecursorData}
      PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.01]
      PrecursorData[Precursor.Quantity == 0, Precursor.Quantity   := NA]
      PrecursorData[Precursor.Translated == 0, Precursor.Translated := NA]
      PrecursorData <- PrecursorData[!is.na(Precursor.Quantity)]
    }
    message("Normalising Precursor Quantities")
    {
      PrecursorData[, Precursor.Quantity   := Precursor.Quantity/sum(Precursor.Quantity)*PrecursorData[, sum(Precursor.Quantity), Run][, median(V1)], by = Run]
      PrecursorData[, Precursor.Translated := Precursor.Translated/sum(Precursor.Translated, na.rm = TRUE)*PrecursorData[, sum(Precursor.Quantity), Run][, median(V1)], Run]
    }
    message("Processing SILAC Labels")
    {
      if(PrecursorData[Precursor.Id %like% "SILAC-.-L" & Precursor.Id %like% "SILAC-.-H", .N] != 0){
        message("ERROR: Multi-Label Precursors Detected")
        MultiLabelDetection = T} else {MultiLabelDetection = F}
      
      PrecursorData[data.table::like(Precursor.Id, "SILAC-.-L"), Label := "L"]
      PrecursorData[data.table::like(Precursor.Id, "SILAC-.-H"), Label := "H"]
      PrecursorData[,Precursor.Id.nolabels := gsub("SILAC-.-.", "SILAC", Precursor.Id)]
      
      FullIsotopeRatio <- PrecursorData[!is.na(Precursor.Translated) & !is.na(Label), .(LabelIntensity = sum(Precursor.Translated)), .(Sample, Condition, Replicate, Label)] |> 
        data.table::merge.data.table(PrecursorData[!is.na(Precursor.Translated) & !is.na(Label), .(TotalIntensity = sum(Precursor.Translated)), .(Sample, Condition, Replicate)])
      FullIsotopeRatio[, Prop := LabelIntensity/TotalIntensity]
      
      IsotopeIncorporation <- FullIsotopeRatio |> ggplot2::ggplot(ggplot2::aes(x = gsub("_", " ", Sample), y = Prop, fill = Condition, alpha = Label)) +
        ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::scale_alpha_manual(values = c("H" = 0.5, "L" = 1), guide = "none") + ggplot2::ylab("Label Incorporation Ratio") +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      IntensitiesRatioData <- PrecursorData |> data.table::dcast(Run+ProteinGroup+Protein.Ids+ProteinDescription+Gene+Stripped.Sequence+
                                                                   Proteotypic+Sample+Drug+Temp+Replicate+Condition+Precursor.Id.nolabels ~Label, 
                                                                 value.var = "Precursor.Translated")
      
      IntensitiesRatioData <- IntensitiesRatioData[!is.na(H) & H > 0 & !is.na(L) & L > 0]
      IntensitiesRatioData <- IntensitiesRatioData[,.(H = sum(H, na.rm = T), L = sum(L, na.rm = T), .N), .(ProteinGroup, ProteinDescription, Gene, Sample)]
      IntensitiesRatioData[, Ratio := H/L]
      data.table::fwrite(IntensitiesRatioData, "IntensitiesRatioData.csv")
    }
    message("Compiling & Exporting Data")
    {
      PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
      data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
      
      LFQ_T <- Proteopedia::Calculate_LFQ(PrecursorData, "LFQ", SILAC = T)
      Intensity_T <- PrecursorData[,.(Intensity = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_T <- PrecursorData[, .(N_precursors = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic = sum(Proteotypic)), .(ProteinGroup, Sample)]
      
      LFQ_L <- Proteopedia::Calculate_LFQ(PrecursorData[Label == "L"], "LFQ_L", SILAC = T)
      Intensity_L <- PrecursorData[Label == "L",.(Intensity_L = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_L <- PrecursorData[Label == "L" , .(N_precursors_L = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_L = sum(Proteotypic)), .(ProteinGroup, Sample)]
      
      LFQ_H <- Proteopedia::Calculate_LFQ(PrecursorData[Label == "H"], "LFQ_H", SILAC = T)
      Intensity_H <- PrecursorData[Label == "H",.(Intensity_H = sum(Precursor.Quantity)), .(ProteinGroup, Sample)]
      Counts_H <- PrecursorData[Label == "H" , .(N_precursors_H = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_H = sum(Proteotypic)), .(ProteinGroup, Sample)]
      
      Annotations <- unique(PrecursorData[, .(ProteinGroup, Run, Sample, Condition, Drug, Temp, Replicate, ProteinDescription, Gene)])
      ProteinData <- Reduce(Proteopedia::Merge_PrecursorData, list(LFQ_T, LFQ_H, LFQ_L, Intensity_T, Intensity_H, Intensity_L, Counts_T, Counts_H, Counts_L, Annotations))                               
      ProteinData[, LFQ_Ratio := LFQ_H/LFQ_L]
      data.table::fwrite(ProteinData, "SILAC_DIANN_Output.csv.gz")
    }
    message("Plotting Intensities")
    {  
      PrecursorData[, Label := gsub("L", "Light", Label)]
      PrecursorData[, Label := gsub("H", "Heavy", Label)]
      
      IntensitiesData <- data.table::rbindlist(list(
        PrecursorData[!is.na(Label), .(Sample, Condition, Replicate, Label, `log2 quantity` = log2(Precursor.Quantity), Type = "Precursor Quantity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ_L), Label = "Light", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity_L), Label = "Light", Type = "Protein Intensity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ_H), Label = "Heavy", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity_H), Label = "Heavy", Type = "Protein Intensity")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ), Label = "Total", Type = "Max. Protein LFQ")],
        ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity), Label = "Total", Type = "Protein Intensity")]
      ), use.names = TRUE)
      IntensitiesData[, Type := factor(Type, levels = c("Precursor Quantity", "Max. Protein LFQ", "Protein Intensity"))]
      
      IntensityPlot <- IntensitiesData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = `log2 quantity`, colour = Condition, alpha = Label)) + 
        ggplot2::geom_boxplot(outliers = FALSE) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1, "Total" = 1), guide = "none") + 
        ggplot2::facet_grid(cols = ggplot2::vars(Type), rows = ggplot2::vars(Label), scales = "free_x") + 
        ggplot2::ylab(expression("Log"[2]~"Value")) + ggplot2::coord_flip() + ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Plotting Precursor, Peptide & Protein Counts")
    {
      CountPlot <- data.table::melt(PrecursorData[!is.na(Label), lapply(.SD, data.table::uniqueN), .(Sample, Condition, Label), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")], 
                                    id.vars = c("Sample","Condition", "Label"), value.name = "IDs") |> 
        ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = IDs/1000, fill = Condition, alpha = Label, 
                                     label = format(IDs, big.mark = ",", scientific = FALSE))) +
        ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.1) +
        ggplot2::facet_grid(ggplot2::vars(Label), ggplot2::vars(variable), scales = "free_x", 
                            labeller = ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups", 
                                                              Light = "Light", Heavy = "Heavy"))) + 
        ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                       strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Calculating Data Completeness")
    {
      CompletenessData <- rbind(Proteopedia::Count_Proteins(ProteinData, "All"), Proteopedia::Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                                Proteopedia::Count_Proteins(ProteinData[N_precursors_proteotypic >= 2], "≥ 2 Proteotypic"))
      
      NAsPlot <- CompletenessData |> ggplot2::ggplot(ggplot2::aes(x = N_samples, y = cumulative_protein_N/1000, colour = Precursors))+
        ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_colour_manual(values = c("All" = "black", "≥ 2" = "darkgrey", "≥ 2 Proteotypic" = "orange3")) +
        ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") + ggplot2::scale_x_continuous( breaks = seq(1, 1000, 1)) +
        ggplot2::scale_y_continuous( limits = c(0, max(CompletenessData$cumulative_protein_N)/1000))+
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       legend.position = "inside", legend.position.inside = c(0.25, 0.25))
    }
    message("Plotting Missed Trypsinisation Sites")
    {
      TrypsinData <- PrecursorData |> data.table::copy()
      TrypsinData[, MissedTrypsin := grepl("[RK][^P]", Stripped.Sequence)]
      TrypsinData[, N_Trypsin := .N, by = .(Sample, MissedTrypsin, Label)]
      TrypsinData[, N_Sample := .N, by = .(Sample, Label)]
      TrypsinData <- TrypsinData[MissedTrypsin == T, .(Sample, Condition, Replicate, Label, N_Trypsin, N_Sample)] |> dplyr::distinct()
      suppressWarnings(TrypsinData[, PercentTrypsin := (N_Trypsin/N_Sample)*100])
      
      MissedCleavagePlot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentTrypsin, fill = Condition, alpha = Label)) + 
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
        ggplot2::labs(y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::coord_flip() + ggplot2::theme(axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(), 
                                               strip.text.y = ggplot2::element_text(size = 26))
    }
    message("Calculating Precursor & Protein Variation")
    {
      PrecursorCVs <- PrecursorData[, .(CV = Proteopedia::Calculate_CV(Precursor.Translated), N = .N), .(Precursor.Id, Condition, Label)]
      PrecursorCVs <- PrecursorCVs[N >= 2]  
      PrecursorCVs[, rank := data.table::frank(CV), .(Condition, Label)] 
      PrecursorCVs[, ID := "Precursors"]
      
      ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_L), N = .N), .(ProteinGroup, Condition)][, Label := "Light"] |> 
        rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_H), N = .N), .(ProteinGroup, Condition)][, Label := "Heavy"])
      ProteinCVs <- ProteinCVs[N >= 2]
      ProteinCVs <- ProteinCVs[, rank := data.table::frank(CV), .(Condition, Label)] 
      ProteinCVs[, ID := "Protein Groups"]
      
      AllCVs <- PrecursorCVs[, Precursor.Id := NULL] |> rbind(ProteinCVs[, ProteinGroup := NULL])
      
      VariationPlot <- AllCVs |> ggplot2::ggplot(ggplot2::aes(x = rank/1000, y = CV, colour = gsub("_", " ", Condition), alpha = Label)) +
        ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Variation [%]") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, name = "Condition") +
        ggplot2::coord_cartesian(ylim = c(0,50)) + ggplot2::facet_wrap(~ID, scales = "free") + 
        ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                       legend.position = "inside", legend.position.inside = c(0.7, 0.7), strip.background = ggplot2::element_blank(), 
                       strip.text.y = ggplot2::element_text(size = 26))
    }
  }
  message("Exporting QC Plots")
  {
    pdf("PrecursorQC_DIANN_Plot.pdf", width = 18, height = 20)
    print(IntensityPlot + CountPlot + patchwork::free(NAsPlot, type = "label") + 
            patchwork::free(MissedCleavagePlot, type = "label") + patchwork::free(VariationPlot, type = "label") +
            patchwork::plot_layout(design = "AAAA\nBBBB\nCCDD\nEEEE") + patchwork::plot_annotation(tag_levels = "A"))
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Analyse_TPP_Proteins <- function(InputDirectory, CtlGroup, ExpGroup, LFQChannel, PropMissingness = 0.5, Log2FCvs = "Max",
                                 ProteotypicFiltering = F, ExcludedSamples = c()){
  set.seed(123)
  start.time <- Sys.time()
  message("Loading & Formatting Data")
  {
    setwd(InputDirectory)
    ProteinData <- data.table::fread(list.files(pattern = "DIANN_Output.csv.gz")[1]) |> data.table::setorder(LFQ)
    ProteinData <- ProteinData[Sample %!in% ExcludedSamples]
    ProteinInfo <- ProteinData[, .(ProteinGroup, Gene, ProteinDescription)] |> data.table::copy() |> unique()
    BCAData <- data.table::fread("BCAData.csv")
    
    if(LFQChannel == "Heavy"){
      ProteinData[, `:=`(Replicate = gsub(".*_(R\\d+)", "\\1", Replicate), Temp = gsub(".*_(.*)C", "\\1", Condition), 
                         Drug = gsub("(.*)_.*C", "\\1", Condition), Log2LFQ = log2(LFQ_H), DetectLevel = mean(head(LFQ_H, 200), na.rm = T)), Sample]
      if(dir.exists(paste0(getwd(), "/Heavy_Output")) == TRUE){unlink(paste0(getwd(), "/Heavy_Output"), recursive = TRUE)}
      dir.create(paste0(getwd(), "/Heavy_Output"), showWarnings = TRUE)
      setwd(paste0(getwd(), "/Heavy_Output"))
    } else if(LFQChannel == "Light"){
      ProteinData[, `:=`(Replicate = gsub(".*_(R\\d+)", "\\1", Replicate), Temp = gsub(".*_(.*)C", "\\1", Condition), 
                         Drug = gsub("(.*)_.*C", "\\1", Condition), Log2LFQ = log2(LFQ_L), DetectLevel = mean(head(LFQ_L, 200), na.rm = T)), Sample]
      if(dir.exists(paste0(getwd(), "/Light_Output")) == TRUE){unlink(paste0(getwd(), "/Light_Output"), recursive = TRUE)}
      dir.create(paste0(getwd(), "/Light_Output"), showWarnings = TRUE)
      setwd(paste0(getwd(), "/Light_Output"))
    } else {
      ProteinData[, `:=`(Replicate = gsub(".*_(R\\d+)", "\\1", Replicate), Temp = gsub(".*_(.*)C", "\\1", Condition), 
                         Drug = gsub("(.*)_.*C", "\\1", Condition), Log2LFQ = log2(LFQ), DetectLevel = mean(head(LFQ, 200), na.rm = T)), Sample]
      if(dir.exists(paste0(getwd(), "/Total_Output")) == TRUE){unlink(paste0(getwd(), "/Total_Output"), recursive = TRUE)}
      dir.create(paste0(getwd(), "/Total_Output"), showWarnings = TRUE)
      setwd(paste0(getwd(), "/Total_Output"))
    }
    Metadata <- ProteinData[, .(Sample, Condition, Drug, Temp, Replicate)] |> unique() |> data.table::setorderv("Sample")
    TempLevels <- Metadata$Temp |> unique() |> as.numeric() |> sort() |> as.character()
    Metadata[, Temp := factor(Temp, levels = TempLevels)]
    data.table::fwrite(Metadata, "Sample_Metadata.csv")
    
    PCAData <- ProteinData[, .(ProteinGroup, Sample, LFQ)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "LFQ") |> tidyr::drop_na() |>
      data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC3"]*100,0), "%]", sep = ""),
                                                  y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC4"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
    pdf("RawLFQ_PCA.pdf", width = 12, height = 6)
    print(Output_PCA)
    Proteopedia::Reset_Dev()
  }
  message("Performing Median Normalisation")
  {
    SpectraMean <- median(ProteinData$Log2LFQ, na.rm = TRUE)
    ProteinData[, `:=`(NormLog2LFQ, Log2LFQ - median(Log2LFQ, na.rm = TRUE) + SpectraMean), Sample]
    ProteinDataWide <- ProteinData |> data.table::dcast(ProteinGroup ~ Sample, value.var = "NormLog2LFQ") |> tibble::column_to_rownames("ProteinGroup")
    
    ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(2^Log2LFQ), Stage = "Pre-Normalisation"), .(ProteinGroup, Condition)] |>
      rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(2^NormLog2LFQ), Stage = "Post-Normalisation"), .(ProteinGroup, Condition)])
    
    pdf("ProteinNormalisationCVs.pdf", width = 16, height = 12)
    print(ProteinCVs |> data.table::merge.data.table(Metadata[, .(Condition, Temp, Drug)] |> unique(), by = "Condition") |> 
            ggplot2::ggplot(ggplot2::aes(x = Temp, y = CV, colour = Drug)) + ggplot2::geom_boxplot(outliers = F) + 
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggplot2::facet_wrap(~forcats::fct_rev(Stage)) +
            ggplot2::labs(x = Temp.~degree~C, y = "CV (%)"))
    Proteopedia::Reset_Dev()
    
    PCAData <- ProteinData[, .(ProteinGroup, Sample, NormLog2LFQ)] |>
      data.table::dcast(ProteinGroup ~ Sample, value.var = "NormLog2LFQ") |> tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |>
      t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC3"]*100,0), "%]", sep = ""),
                                                  y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC4"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
    pdf("NormLFQ_PCA.pdf", width = 12, height = 6)
    print(Output_PCA)
    Proteopedia::Reset_Dev()
  }
  message("Checking Missingness & Calculating FC vs Median")
  {
    Missingness <- ProteinData[, .(NAs = sum(is.na(NormLog2LFQ))), .(ProteinGroup, Temp, Drug)] |> data.table::dcast(ProteinGroup+Temp~Drug, value.var = "NAs")
    Missingness[, NAsDiff := get(ExpGroup) - get(CtlGroup)]
    Missingness <- Missingness[, .(MeanNAsDiff = mean(NAsDiff)), ProteinGroup]
    
    MissingnessMap <- is.na(ProteinDataWide) |> apply(2, as.numeric)
    rownames(MissingnessMap) <- rownames(ProteinDataWide)
    N_MissingProteins <- MissingnessMap |> matrixStats::rowSums2()
    MissingnessMap <- MissingnessMap[N_MissingProteins <= nrow(Metadata)*PropMissingness,] # Pre-Imputation Missingness
    MissingnessHeatmap <- pheatmap::pheatmap(MissingnessMap, show_rownames = F, show_colnames = F, annotation_col = Metadata |> tibble::column_to_rownames("Sample"))
    ProteinDrugs <- data.table::data.table(ProteinGroup = rownames(MissingnessMap), Cluster = as.factor(cutree(MissingnessHeatmap$tree_row, k = 5)))
    pheatmap::pheatmap(MissingnessMap, show_rownames = F, show_colnames = F, color = c("#090", "#000"), legend = F,
                       annotation_col = Metadata[, .(Sample, Drug, Temp, Replicate)] |> tibble::column_to_rownames("Sample"),
                       annotation_row = ProteinDrugs |> tibble::column_to_rownames("ProteinGroup"), filename = "MissingnessHeatmap.pdf")
    Proteopedia::Reset_Dev()
    
    ProteinData[, SampleMed := median(Log2LFQ, na.rm = T), by = .(Sample)]
    ProteinData[, Log2FCvsMedian := Log2LFQ - SampleMed]
    
    PCAData <- ProteinData[, .(ProteinGroup, Sample, Log2FCvsMedian)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2FCvsMedian") |>
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC3"]*100,0), "%]", sep = ""),
                                                  y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC4"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
    pdf("Log2FCvsMed_PCA.pdf", width = 12, height = 6)
    print(Output_PCA)
    Proteopedia::Reset_Dev()
  }
  message("Perform BCA Correction")
  {
    BCAStandards <- BCAData[grepl("Standard", Sample)][, Conc := as.numeric(gsub("Standard_", "", Sample))]
    Standard_Equation <- summary(lm(Absorbance ~ Conc, data = BCAStandards))
    Standard_c <- round(Standard_Equation$coefficients[1], 3)
    Standard_m <- round(Standard_Equation$coefficients[2], 3)
    Standard_R2 <- round(Standard_Equation$r.squared, 3)
    
    BCASamples <- BCAData[!grepl("Standard", Sample)][, Conc := (Absorbance - Standard_c) / Standard_m] |> data.table::merge.data.table(Metadata)
    
    pdf("RawBCAPlot.pdf", width = 16, height = 16)
    print(BCAStandards |> ggplot2::ggplot(ggplot2::aes(x = Conc, y = Absorbance)) +
            ggplot2::geom_point(stroke = NA) +ggplot2::stat_smooth(method = "lm", se = F, colour = "#AAA", linetype = "dashed", size = 1) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.5, label = paste0("R^2 == ", Standard_R2), parse = T) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.45, label = paste0("y = ", Standard_m, "x + ", Standard_c)) +
            ggplot2::geom_rug(data = BCASamples, ggplot2::aes(colour = Drug, alpha = factor(Temp))) +
            ggplot2::scale_alpha_manual(values = c(1, 0.85, 0.7, 0.55, 0.4, 0.25), name = Temp.~degreee~C) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
            ggplot2::labs(x = "Protein Concentration (mg/mL)", y = "562 nm Absorbance (AU)"))
    Proteopedia::Reset_Dev()
    
    StandardFiltering <- c(BCAStandards$Conc |> unique(), max(BCASamples$Conc)) |> sort()
    StandardFiltering <- StandardFiltering[1:(which(StandardFiltering == max(BCASamples$Conc)) + 1)]
    
    BCAStandards <- BCAStandards[Conc %in% StandardFiltering]
    Standard_Equation <- summary(lm(Absorbance ~ Conc, data = BCAStandards))
    Standard_c <- round(Standard_Equation$coefficients[1], 3)
    Standard_m <- round(Standard_Equation$coefficients[2], 3)
    Standard_R2 <- round(Standard_Equation$r.squared, 3)
    
    BCASamples <- BCAData[!grepl("Standard", Sample)][, Conc := (Absorbance - Standard_c) / Standard_m] |> data.table::merge.data.table(Metadata)
    
    pdf("AdjBCAPlot.pdf", width = 16, height = 16)
    print(BCAStandards |> ggplot2::ggplot(ggplot2::aes(x = Conc, y = Absorbance)) +
            ggplot2::geom_point(stroke = NA) +ggplot2::stat_smooth(method = "lm", se = F, colour = "#AAA", linetype = "dashed", size = 1) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.5, label = paste0("R^2 == ", Standard_R2), parse = T) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.45, label = paste0("y = ", Standard_m, "x + ", Standard_c)) +
            ggplot2::geom_rug(data = BCASamples, ggplot2::aes(colour = Drug, alpha = factor(Temp))) +
            ggplot2::scale_alpha_manual(values = c(1, 0.85, 0.7, 0.55, 0.4, 0.25), name = Temp.~degree~C) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
            ggplot2::labs(x = "Protein Concentration (mg/mL)", y = "562 nm Absorbance (AU)"))
    Proteopedia::Reset_Dev()
    
    BCASamples[, Correction := Conc/BCASamples[, .(MedConc = median(Conc)), .(Drug, Temp)]$MedConc |> max(), Sample]
    BCASamples <- BCASamples[, .(Correction = median(Correction)), Temp]
    BCASamples[, Temp := factor(Temp)]
    
    pdf("CorrectionCoeffsPlot.pdf", width = 16, height = 12)
    print(BCASamples |> unique() |> ggplot2::ggplot(ggplot2::aes(x = Temp, y = Correction, fill = Temp)) +
            ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2(width = 0.9)) + 
            ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)], guide = "none") +
            ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::labs(x = Temp.~degree~C, y = "BCA Correction Coeff."))
    Proteopedia::Reset_Dev()
  }
  message("Imputing Non-Quantified Proteins")
  {  
    ProcessedData <- data.table::copy(ProteinData) |> data.table::dcast(ProteinGroup ~ Sample, value.var = "NormLog2LFQ") |>
      data.table::melt.data.table(id.vars = "ProteinGroup", value.name = "NormLog2LFQ", variable.name = "Sample") |>
      data.table::merge.data.table(Metadata, by = "Sample")
    
    ProcessedFCData <- data.table::copy(ProteinData) |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2FCvsMedian") |>
      data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Log2FCvsMedian") |> 
      data.table::merge.data.table(Metadata, by = "Sample")
    
    ProcessedDataMerged <- ProcessedData |> data.table::merge.data.table(ProcessedFCData, by = c("Sample", "ProteinGroup", "Condition", "Drug", "Temp", "Replicate"))
    ProcessedDataMerged[, DetectLevel := mean(head(sort(NormLog2LFQ, decreasing = F), 200), na.rm = T), Sample]
    ProcessedDataMerged[, N_Missing := sum(is.na(Log2FCvsMedian)), .(ProteinGroup, Condition)]
    
    # Calculating for same timepoint how many missing in each condition to see if missingness is biologically informative and need imputation
    DiffMissing <- unique(ProcessedDataMerged[, .(ProteinGroup, Drug, Temp, N_Missing)]) |> data.table::dcast(ProteinGroup + Temp ~ Drug, value.var = "N_Missing")
    DiffMissing[, DiffMissing := get(ExpGroup) - get(CtlGroup)]
    ProcessedDataMerged <- data.table::merge.data.table(ProcessedDataMerged, DiffMissing[, .(ProteinGroup, Temp, DiffMissing)], by = c("ProteinGroup", "Temp"))
    
    # Determine & Execute Imputation ####
    ProcessedDataMerged[, Temp := factor(Temp, levels = TempLevels)]
    ProcessedDataMerged[, NextTemp := Proteopedia::Calculate_RelativeTemp(Temp, TempLevels, Shift = 1) |> as.character(), Temp]
    
    NextTempData <- ProcessedDataMerged[, .(NextTemp_NMissing = mean(N_Missing)), .(Temp, ProteinGroup, Drug)] |> data.table::copy() |> data.table::setnames("Temp", "NextTemp")
    
    ProcessedDataMerged <- data.table::merge.data.table(ProcessedDataMerged, NextTempData, by = c("NextTemp", "ProteinGroup", "Drug"), all.x = T)
    
    # Imputed proteins quantified in one condition but not other (N_Missing == 2, Diff_detected == 2) & missing from next Temp
    ProcessedDataMerged[NextTemp_NMissing == 2 & N_Missing == 2 & abs(DiffMissing) == 2, Impute := T]
    ProcessedDataMerged[is.na(Impute), Impute := F]
    ProcessedDataMerged[, ImpLog2LFQ := data.table::fifelse(is.na(NormLog2LFQ) & Impute == T, rnorm(.N, mean = DetectLevel, 0.1) - 1, NormLog2LFQ)]
  }
  message("Plotting LFQs & Summarising")
  {
    TempLabels <- paste0(TempLevels, "~degree*C")
    names(TempLabels) <- paste0(TempLevels)
    TempLabels <- ggplot2::as_labeller(TempLabels, default = ggplot2::label_parsed)
    
    pdf("LFQHistogram.pdf", width = 16, height = 12)
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = ImpLog2LFQ, fill = is.na(Log2FCvsMedian))) + 
            ggplot2::geom_histogram() + ggplot2::scale_fill_manual(values  = c("#000", "#F10"), name = "Imputed") + 
            ggplot2::geom_vline(ggplot2::aes(xintercept = DetectLevel), colour = "#02F", linetype = "dashed") +
            ggplot2::labs(x = expression("Log"[2]~"LFQ"), y = "No. Proteins") + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
            ggplot2::facet_grid(ggplot2::vars(Drug), ggplot2::vars(Temp), labeller = ggplot2::labeller(Temp = TempLabels)))
    Proteopedia::Reset_Dev()
    
    ProcessedDataMerged <- ProcessedDataMerged |> data.table::merge.data.table(BCASamples, by = "Temp") |>
      data.table::merge.data.table(ProteinData[, .(Sample, ProteinGroup, Log2LFQ)], by = c("Sample", "ProteinGroup"), all.x = T)
    ProcessedDataMerged[, BCACorrLog2LFQ := ImpLog2LFQ*Correction]
    
    pdf("LFQProcessingBoxplots.pdf", width = 16, height = 12)
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Drug), colour = as.factor(Replicate), y = Log2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Raw Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Drug), colour = as.factor(Replicate), y = NormLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Normalised Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Drug), colour = as.factor(Replicate), y = ImpLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Post-Imputation Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Drug), colour = as.factor(Replicate), y = BCACorrLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("BCA-Corrected Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    Proteopedia::Reset_Dev()
    
    ProcessedProteinLFQs <- ProcessedDataMerged[, .(ProteinGroup, Sample, Condition, Replicate, Drug, Temp, NextTemp, DetectLevel, 
                                                    Impute, Log2LFQ, NormLog2LFQ, ImpLog2LFQ, BCACorrLog2LFQ)] |> data.table::copy()
    data.table::fwrite(ProcessedProteinLFQs, "ProcessedProteinLFQs.csv")
  }
  message("Preparing Limma Input")
  {
    ProtLFQsInput <- as.matrix(tibble::column_to_rownames(data.table::dcast(ProcessedProteinLFQs, ProteinGroup ~ Sample, value.var = "BCACorrLog2LFQ"), "ProteinGroup"))
    LimmaInputData <- data.table::data.table(tibble::rownames_to_column(as.data.frame(ProtLFQsInput), "ProteinGroup"))
    LimmaInputData <- data.table::merge.data.table(data.table::melt.data.table(LimmaInputData, id.vars = "ProteinGroup", variable.name = "Sample", 
                                                                               value.name = "Log2LFQ"), Metadata)
    LimmaInputData[, `:=`(MeanLog2LFQ, mean(Log2LFQ, na.rm = T)), by = .(ProteinGroup, Drug, Temp)]
    
    LimmaInputMatrix <- as.matrix(tibble::column_to_rownames(data.table::dcast(LimmaInputData, ProteinGroup ~ Sample, value.var = "Log2LFQ"), "ProteinGroup"))
    LimmaInputMatrix <- LimmaInputMatrix[matrixStats::rowMeans2(is.na(LimmaInputMatrix)) <= PropMissingness,]
    LimmaInputData <- LimmaInputData[ProteinGroup %in% rownames(LimmaInputMatrix)]
    
    PCAData <- LimmaInputMatrix[matrixStats::rowMeans2(is.na(LimmaInputMatrix)) == 0,] |> t() |> stats::prcomp(scale. = TRUE)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |>
      data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Drug, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC3"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC4"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()))
    pdf("LimmaInput_PCA.pdf", width = 12, height = 6)
    print(Output_PCA)
    Proteopedia::Reset_Dev()
  }
  message("Performing Linear Modelling")
  {
    Targets <- Metadata[, .(Sample, Drug, Temp)]
    Targets <- Targets[match(colnames(LimmaInputMatrix), Sample)]
    Temp <- Targets$Temp
    Drug <- factor(Targets$Drug)
    design <- model.matrix(~ Temp + Temp:Drug)
    colnames(design) <- gsub(paste0("Drug", ExpGroup), "Drug_Exp", colnames(design))
    LimmaOutput <- limma::eBayes(limma::lmFit(Biobase::ExpressionSet(assayData = LimmaInputMatrix), design))
    
    if(!is.finite(LimmaOutput$df.prior)){message("Warning: Limma df.prior is Infinite")}
    Limma_Slopes <- data.table::setnames(data.table::data.table(LimmaOutput$coefficients, keep.rownames = T), "rn", "ProteinGroup")
    Limma_Slopes_SD <- data.table::copy(data.table::data.table(LimmaOutput$stdev.unscaled, keep.rownames = T))
    Limma_Slopes_SD <- data.table::setnames(data.table::copy(data.table::data.table(LimmaOutput$stdev.unscaled, keep.rownames = T)), 
                                            c(colnames(Limma_Slopes_SD)), c("ProteinGroup", paste0(colnames(Limma_Slopes_SD)[-1], "_SD")), skip_absent = T)
    Limma_Slopes <- data.table::merge.data.table(Limma_Slopes, Limma_Slopes_SD, by = "ProteinGroup")
    Limma_Slopes$Fvalue = LimmaOutput$F
    
    pdf("MeanVariancePlot.pdf", width = 12, height = 8)
    print(ggplot2::ggplot(data.table::data.table(Mean = LimmaOutput$Amean, Variance = sqrt(LimmaOutput$sigma)), ggplot2::aes(x = Mean, y = Variance)) + 
            ggplot2::geom_point(stroke = NA) + ggplot2::labs(x = expression("Mean Log"[2]~"LFQ"), y = "Variance"))
    Proteopedia::Reset_Dev()
    
    LimmaStats <- data.table::data.table()
    for(Coefficient in colnames(LimmaOutput$coefficients)){
      LimmaCoeffStats <- limma::topTable(LimmaOutput, Coefficient, number = Inf) |> data.table::data.table(keep.rownames = T) |> 
        data.table::setnames("rn", "ProteinGroup")
      LimmaCoeffStats[, Coef := Coefficient]
      LimmaStats <- rbind(LimmaStats, LimmaCoeffStats)
    }
    LimmaStats <- LimmaStats[!is.na(logFC)] |> data.table::merge.data.table(ProteinInfo) |> data.table::setorder(adj.P.Val)
    LimmaStats[grepl(":", Coef), Drug := "Reversine"]
    LimmaStats[grepl("Temp\\d+$", Coef), Drug := "DMSO"]
    LimmaStats[grepl("Temp", Coef), Temp := gsub("Temp(\\d+).*", "\\1", Coef)]
    LimmaStats[, Significance := data.table::fifelse(adj.P.Val < 0.05, "Adj. P-Value < 0.05", 
                                                     data.table::fifelse(P.Value < 0.05, "P-Value < 0.05", "None"))] 
    LimmaStats[grepl(":", Coef), SigMark := data.table::fifelse(P.Value < 0.001, "***", data.table::fifelse(P.Value < 0.01, "**", data.table::fifelse(P.Value < 0.05, "*", "")))]
    LimmaStats[, Temp := factor(Temp, levels = TempLevels)]
    LimmaStats[, Significance := factor(Significance, levels = c("Adj. P-Value < 0.05", "P-Value < 0.05", "None"))]
    
    data.table::fwrite(LimmaStats, "LimmaStats.csv")
    
    CoefLabels <- paste0(TempLevels, "~degree*C")
    names(CoefLabels) <- paste0("Temp",TempLevels, ":Drug_Exp")
    CoefLabels <- ggplot2::as_labeller(CoefLabels, default = ggplot2::label_parsed)
    
    pdf("LimmaPValueHistograms.pdf", width = 16, height = 12)
    print(LimmaStats[grepl(":", Coef)] |> ggplot2::ggplot(ggplot2::aes(x = P.Value, fill = Significance)) + ggplot2::geom_histogram() + 
            ggplot2::facet_wrap(~Coef, labeller = ggplot2::labeller(Coef = CoefLabels), scales = "free_y") +
            ggplot2::scale_fill_manual(values = c("Adj. P-Value < 0.05" = "#FB0", "P-Value < 0.05" = "#F10", "None" = "#300")) + 
            ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
            ggplot2::labs(x = "Adjusted P-Value", y = "No. Proteins") + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2)))  
    print(LimmaStats[grepl(":", Coef)] |> ggplot2::ggplot(ggplot2::aes(x = adj.P.Val, fill = Significance)) + ggplot2::geom_histogram() + 
            ggplot2::facet_wrap(~Coef, labeller = ggplot2::labeller(Coef = CoefLabels), scales = "free_y") +
            ggplot2::scale_fill_manual(values = c("Adj. P-Value < 0.05" = "#FB0", "P-Value < 0.05" = "#F10", "None" = "#300")) + 
            ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::scale_x_continuous(expand = c(0, 0)) + 
            ggplot2::labs(x = "Adjusted P-Value", y = "No. Proteins") + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2)))  
    Proteopedia::Reset_Dev()
    
    AffectedProteins <- LimmaStats[stringr::str_detect(Coef, ":") & P.Value < 0.05][stringr::str_detect(Coef, TempLevels[1], negate = T)]
    AffectedProteinData <- LimmaInputData[ProteinGroup %in% AffectedProteins$ProteinGroup & !is.na(Log2LFQ)]
    data.table::fwrite(AffectedProteinData, "AffectedProteinsLFQs.csv")
  }
  message(paste0("Calculating Fold Changes vs ", Log2FCvs))
  {
    ProcessedDataMerged[Temp == min(TempLevels), ProteinLFQMax := mean(BCACorrLog2LFQ, na.rm = T), .(ProteinGroup, Drug)]
    ProcessedDataMerged[, ProteinLFQMax := mean(ProteinLFQMax, na.rm = T), .(ProteinGroup, Drug)]
    ProcessedDataMerged[, Log2FCvsMax := BCACorrLog2LFQ - ProteinLFQMax]
    
    ProcessedDataMerged[Temp == max(TempLevels), ProteinLFQMin := mean(BCACorrLog2LFQ, na.rm = T), .(ProteinGroup, Drug)]
    ProcessedDataMerged[, ProteinLFQMin := mean(ProteinLFQMin, na.rm = T), .(ProteinGroup, Drug)]
    ProcessedDataMerged[, Log2FCvsMin := BCACorrLog2LFQ - ProteinLFQMin]
    ProcessedProteinFCs <- ProcessedDataMerged[, .(ProteinGroup, Sample, Condition, Replicate, Drug, Temp, NextTemp, DetectLevel, 
                                                   Impute, Log2FCvsMin, Log2FCvsMedian, Log2FCvsMax)] |> data.table::copy()
    data.table::fwrite(ProcessedProteinFCs, "ProcessedProteinLog2FCs.csv")
    
    FoldChangeData <- ProcessedProteinFCs[, .(ProteinGroup, Drug, Temp, Condition, Replicate, Sample, 
                                              get(colnames(ProcessedProteinFCs)[grepl(Log2FCvs, colnames(ProcessedProteinFCs))]))] |> 
      data.table::copy() |> data.table::setnames("V7", "Log2FC")
    FoldChangeData <- FoldChangeData[, MeanLog2FC := mean(Log2FC, na.rm = T), .(ProteinGroup, Temp, Drug)]
    FoldChangeData[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
    FoldChangeDataWide <- FoldChangeData[, .(MeanLog2FC = mean(Log2FC, na.rm = T)), by = .(ProteinGroup, Temp, Drug)] |>
      data.table::dcast(ProteinGroup + Temp ~ Drug, value.var = "MeanLog2FC") |> data.table::merge.data.table(ProteinDrugs, by = "ProteinGroup")
    
    pdf(paste0("Log2FCvsTemp", Log2FCvs, "Dotplot.pdf"), width = 16, height = 12)
    print(FoldChangeDataWide[order(Cluster)] |> ggplot2::ggplot(ggplot2::aes(x = get(CtlGroup), y = get(ExpGroup), colour = factor(Cluster))) +
            ggplot2::geom_point(stroke = NA, alpha = 1) + ggplot2::geom_abline(linetype = "dashed", colour = "#AAA") +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, name = "Cluster") + ggplot2::facet_wrap("Temp") + 
            ggplot2::labs(x = Log[2]~FC~vs.~Control~Mean~at~37~degree~C, y = Log[2]~FC~vs.~37~degree~C))
    Proteopedia::Reset_Dev()
  }
  message("Identifying Most Sig. Affected Proteins")
  {
    PlotsOI = list()
    for(POI in (AffectedProteins$ProteinGroup |> unique())[1:20]){
      PlotsOI[[POI]] = LimmaInputData[ProteinGroup == POI] |> ggplot2::ggplot(ggplot2::aes(x = factor(Temp), y = Log2LFQ, colour = Drug, group = Drug)) + 
        ggplot2::geom_point() + ggplot2::geom_line(data = LimmaInputData[ProteinGroup == POI, .(MeanLog2LFQ, Drug, Temp)] |> unique(), 
                                                   ggplot2::aes(x = factor(Temp), y = MeanLog2LFQ, colour = Drug, group = Drug)) + 
        ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::labs(title = paste0(ProteinInfo[ProteinGroup == POI, Gene]  |> unique()), x = Temp.~degree~C, y = Log[2]~LFQ)
    }
    
    pdf("Top20Proteins.pdf", width = 20, height = 16)
    print(ggpubr::ggarrange(plotlist = PlotsOI, common.legend = T))
    Proteopedia::Reset_Dev()
    
    Log2FCPlotsOI = list()
    for(POI in (AffectedProteins$ProteinGroup |> unique())[1:20]){
      Log2FCPlotsOI[[POI]] = FoldChangeData[ProteinGroup == POI] |> ggplot2::ggplot(ggplot2::aes(x = factor(Temp), y = Log2FC, colour = Drug, group = Drug)) + 
        ggplot2::geom_point() + ggplot2::geom_line(data = FoldChangeData[ProteinGroup == POI, .(MeanLog2FC, Drug, Temp)] |> unique(), 
                                                   ggplot2::aes(x = factor(Temp), y = MeanLog2FC, colour = Drug, group = Drug)) + 
        ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::labs(title = paste0(ProteinInfo[ProteinGroup == POI, Gene]  |> unique()), x = Temp.~degree~C, y = expression("Log"[2]~"FC in Abundance vs ", Log2FCvs))
    }
    pdf("Top20ProteinLog2FCs.pdf", width = 20, height = 16)
    print(ggpubr::ggarrange(plotlist = Log2FCPlotsOI, common.legend = T))
    Proteopedia::Reset_Dev()
  }
  message("Exporting Annotated Protein Data")
  {
    FoldChangeData[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
    AnnotatedProteinData <- FoldChangeData |> data.table::merge.data.table(Proteopedia::Proteopedia, by.x = "GeneGroup", by.y = "ProteinGroup")
    AnnotatedProteinData[, GeneGroup.y := NULL]
    data.table::fwrite(AnnotatedProteinData, "AnnotatedProteinData.csv")
  }
  Proteopedia::End_Timer(start.time)
}
########### MS Analysis: Protein Data to Mapping Functions ####################################################################################################################################
#' @export
Map_LF_Proteins <- function(InputFile, SubsetColour = "red"){
  start.time <- Sys.time()
  set.seed(123)
  message("Loading Limma File")
  {
    setwd(gsub("(.*)/.*.csv.*", "\\1", InputFile))
    LimmaData <- data.table::fread(InputFile) 
    LimmaData |> data.table::setnames(
      c(colnames(LimmaData)[grepl("protein.*group", ignore.case = T, colnames(LimmaData))], 
        colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & !grepl("adj", ignore.case = T, colnames(LimmaData))], 
        colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & grepl("adj", ignore.case = T, colnames(LimmaData))], 
        colnames(LimmaData)[grepl("protein.*desc", ignore.case = T, colnames(LimmaData))], 
        colnames(LimmaData)[grepl("Gene", ignore.case = T, colnames(LimmaData)) & !grepl("group", ignore.case = T, colnames(LimmaData))]),
      c("ProteinGroup", "P.Value", "adj.P.Val", "ProteinDescription", "Gene"))
    
    if(length(LimmaData$ProteinGroup[grepl("\\-", LimmaData$ProteinGroup)]) > 0){
      LimmaData$Isoforms <- 1
      for(i in 1:nrow(LimmaData)){
        if(length(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T)) > 0){
          LimmaData$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
        } else {LimmaData$Isoforms[i] <- 1}
      }    
      LimmaData[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
    }
    LimmaData[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
    
    if (dir.exists(paste0(getwd(),"/Subset_Output")) == T){
      unlink(paste0(getwd(),"/Subset_Output"), recursive = T)
    }
    dir.create(paste0(getwd(),"/Subset_Output"), showWarnings = T)
    
    message("Importing Proteopedia")
    setwd(paste0(getwd(),"/Subset_Output"))
    LimmaData <- LimmaData |> merge(Proteopedia::Proteopedia, all.x = T) |> data.table::data.table()
    LimmaData$ENTREZID[is.na(LimmaData$ENTREZID)] <- "Unmapped"
    LimmaData$Strand[is.na(LimmaData$Strand)] <- "*"
    LimmaData$Chromosome[is.na(LimmaData$Chromosome)] <- "Unmapped"
    LimmaData$MedianLociStart[is.na(LimmaData$MedianLociStart)] <- "Unmapped"
    LimmaData$Deg_Profile[is.na(LimmaData$Deg_Profile)] <- "UN"
    LimmaData$Experimental_Evidence_ComplexPortal[is.na(LimmaData$Experimental_Evidence_ComplexPortal)] <- "No"
    LimmaData$Experimental_Evidence_CORUM[is.na(LimmaData$Experimental_Evidence_CORUM)] <- "No"
    LimmaData$N_ComplexPortal[is.na(LimmaData$N_ComplexPortal)] <- 0
    LimmaData$N_CORUM[is.na(LimmaData$N_CORUM)] <- 0
    LimmaData[, which(colnames(LimmaData) == "ER"):ncol(LimmaData)][is.na(LimmaData[, which(colnames(LimmaData) == "ER"):ncol(LimmaData)])] <- F
  }
  # Biochemical Trends (Numerical)
  {
    for(ColIndex in c(which(colnames(LimmaData) == "Length"), which(colnames(LimmaData) == "Length"):ncol(LimmaData))){
      if(is.numeric(LimmaData[, get(colnames(LimmaData)[ColIndex])]) & !grepl("N_", colnames(LimmaData)[ColIndex])){
        message(paste0("Analysing ", colnames(LimmaData)[ColIndex]), " Trend")
        SubsetData = LimmaData[, .(ProteinGroup, Sequence, Length, Log2FC, get(colnames(LimmaData)[ColIndex]))]
        SubsetData |> data.table::setnames("V5", "Subset")
        
        pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[2,4]
        
        TrendPlot <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Subset, y = Log2FC)) + 
          ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Add_Rsq() +
          ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Log2FC], na.rm = T)*0.93, 
                            label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) +
          ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2]~"FC in Protein Abundance")) + 
          ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
        
        pdf(paste0(colnames(LimmaData)[ColIndex],"_Trend.pdf"), width = 12, height = 10)
        print(TrendPlot)
        print(TrendPlot + ggplot2::labs(x = "", y = ""))
        Proteopedia::Reset_Dev()
      }
    }
  }
  # Cellular Subsets (TRUE/FALSE-Based)
  {
    for(ColIndex in which(colnames(LimmaData) == "ER"):ncol(LimmaData)){
      if(is.logical(LimmaData[, get(colnames(LimmaData)[ColIndex])])){
        message(paste0("Analysing ", colnames(LimmaData)[ColIndex]), " Proteins")
        SubsetData = LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, get(colnames(LimmaData)[ColIndex]))]
        SubsetData |> data.table::setnames("V5", "Subset")
        
        MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)
        
        Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
          Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes()
        
        pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[4]
        
        Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
          ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
          ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() +
          ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
          ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
          ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T))*0.93, 
                            label = paste0("P-Value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
          Proteopedia::Add_AbundanceAxes()
        pdf(paste0(colnames(LimmaData)[ColIndex],"_Volcano.pdf"))
        print(Volcano)
        print(Volcano_Rug)
        print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
        Proteopedia::Reset_Dev()
      }
    }
  }
  message("Analysing Degradation Profiles")
  {      
    DegProfileSummary <- LimmaData[, .N, by = Deg_Profile]
    
    pdf("DegradataionProfileBoxplot.pdf", width = 12, height = 10)
    print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Deg_Profile, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
            ggpubr::geom_signif(comparison = list(c("NED", "ED"), c("NED", "UN"), c("ED", "UN")), y_position = c(1, 1.25, 1.5), tip_length = 0) + 
            ggplot2::labs(x = "Degradation Profile", y = expression("Log"[2] ~ "FC in Protein Abundance")) + 
            ggplot2::geom_text(data = DegProfileSummary, ggplot2::aes(y = -1, x = Deg_Profile, label = paste0("N = ", N))))
    Proteopedia::Reset_Dev()
    
    for(Deg_Type in unique(LimmaData[!is.na(Deg_Profile), Deg_Profile])){
      MeanLog2FC <- LimmaData[Deg_Profile == Deg_Type, Log2FC] |> mean(na.rm = T)
      
      Volcano <- LimmaData[Deg_Profile == Deg_Type] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
        ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=as.character(Gene))) + 
        Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() 
      
      pval <- summary(stats::lm(Log2FC ~ Subset, data = LimmaData[, Subset := data.table::fifelse(Deg_Profile == Deg_Type, T, NA)]))$coefficients[4]
      
      Volcano_Rug <- LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") +    
        ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
        Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
        ggplot2::geom_rug(alpha = ifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
        ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(LimmaData[, P.Value], na.rm = T))*0.93, 
                          label = paste0("P-value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
        Proteopedia::Add_AbundanceAxes()      
      pdf(paste0(Deg_Type, "_DegProfile_Volcano.pdf"))  
      print(Volcano)
      print(Volcano_Rug)
      print(Volcano_Rug + ggplot2::labs(x= "", y= ""))
      Proteopedia::Reset_Dev()
    }
  }
  message("Summarising Chromosome-Based Data")
  {
    setwd(gsub("(.*)/.*.csv.*", "\\1", InputFile))
    ChromosomeSummary <- data.table::data.table("Chromosome" = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped")) |> 
      merge(LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(MeanLog2FC = mean(Log2FC), Mapped = dplyr::n()) |> 
              merge(Proteopedia::Proteopedia |> dplyr::group_by(Chromosome) |> dplyr::summarise(Total = dplyr::n())))
    ChromosomeSummary[, Coverage := Mapped/Total]
    
    pdf("Chromosome_Coverage.pdf", width = 18, height = 8)
    print(ChromosomeSummary |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), 
                                                            y = Coverage*100, fill = Coverage*100)) +
            ggplot2::geom_bar(stat= "identity") + ggplot2::scale_fill_viridis_c(guide = "none") + 
            ggplot2::labs(x = "Chromosome", y = "Protein Coverage (%)")
    )
    Proteopedia::Reset_Dev()
    
    ChromosomeSummary$Wilcoxon_p = 0
    rownumber = 1
    for(i in ChromosomeSummary$Chromosome){
      ChromosomeSummary$Wilcoxon_p[rownumber] <- wilcox.test(LimmaData[,Log2FC], LimmaData[Chromosome == i, Log2FC])$p.value
      rownumber = rownumber + 1
    }
    ChromosomeSummary[, SigSymbol := data.table::fifelse(Wilcoxon_p < 0.001, "***", data.table::fifelse(Wilcoxon_p < 0.01, "**",
                                                                                                        data.table::fifelse(Wilcoxon_p < 0.05, "*","")))]
    ChromosomeSummary[, Buffering := data.table::fifelse(MeanLog2FC < 0, paste0(">",100,"%"), 
                                                         paste0(100 - (round((MeanLog2FC/log2(3/2))*100, digits = 2)),"%"))]
    data.table::fwrite(ChromosomeSummary, "Chromosomal_Summary.csv")
    
    WhiskerTop <- LimmaData[, Log2FC] |> stats::quantile(0.75, na.rm = T) + stats::IQR(LimmaData[, Log2FC], na.rm = T)*1.5
    WhiskerBottom <- LimmaData[, Log2FC] |> stats::quantile(0.25, na.rm = T) - stats::IQR(LimmaData[, Log2FC], na.rm = T)*1.5
    
    LimmaData[, Chromosome := factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))]
    LimmaData <- LimmaData |> dplyr::arrange(Chromosome, MedianLociStart)
    LimmaData[, OrderID := as.numeric(rownames(LimmaData))]
    
    ChromosomeBorders <- LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(N_Proteins = dplyr::n()) |> 
      dplyr::mutate(Upper = cumsum(N_Proteins)) |> data.table::data.table()
    ChromosomeBorders$Midpoint <- 0
    for(i in 1:nrow(ChromosomeBorders)){
      ChromosomeBorders$Midpoint[i] <- ifelse(i > 1, ChromosomeBorders$Upper[i-1] + ChromosomeBorders$N_Proteins[i]/2, ChromosomeBorders$N_Proteins[i]/2)
    }
    ChromosomeBorders <- ChromosomeBorders[Chromosome %!in% c("Unmapped", "X/Y", "Y", "M")]
    
    ChromosomeCols <- rep(RColorBrewer::brewer.pal(n = 8, "Set2"), length.out = length(c(seq(1:22), "X")))
    names(ChromosomeCols) <- c(seq(1:22), "X")
    
    ChromosomeDotData <- LimmaData[Chromosome %!in% c("X/Y", "Y", "M", "Unmapped")] |> 
      merge(Proteopedia::Proteopedia[, .(ProteinGroup, Gene, Chromosome)], all = T, by = c("ProteinGroup", "Gene", "Chromosome"))
    ChromosomeDotData$Log2FC[is.na(ChromosomeDotData$Log2FC)] <- 0
    ChromosomeDotData <- ChromosomeDotData |> merge(data.table::data.table("Colour" = ChromosomeCols, "Chromosome" = names(ChromosomeCols)), by = "Chromosome")
    ChromosomeDotData[, Colour := data.table::fifelse(Log2FC == 0, "white", Colour)]
    
    ChromosomeDotplot <- ChromosomeDotData |> ggplot2::ggplot(ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) + 
      ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "blue") + 
      ggplot2::scale_colour_identity() +
      ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), 
                                                                label = ifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))), 
                         colour = "black", angle = 90) +
      ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + ggplot2::scale_x_continuous(expand = 0) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(x = "Chromosome", y = expression("Log"[2]~ "FC in Protein Abundance")) + ggplot2::guides(fill="none")
    
    pdf("Chromosome_Dotplot.pdf", width = 18, height = 8)
    print(ChromosomeDotplot)
    Proteopedia::Reset_Dev()
    
    ChromosomeBoxplot <- LimmaData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "Y", "M")), y = Log2FC)) + 
      ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop*1.25, label = SigSymbol), colour = "black") + 
      ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom*1.25, label = Buffering), colour = "black", size = 4) + 
      ggplot2::geom_boxplot(colour = "black", fill = NA, outliers = F, alpha = 0.7) + ggplot2::coord_cartesian(ylim = c(WhiskerTop*1.5, WhiskerBottom*1.5)) +
      ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + 
      ggplot2::labs(x = "Chromosome", y = expression("Log"[2]~ "FC in Protein Abundance")) + ggplot2::guides(fill="none")
    
    ChrGroupingData <- LimmaData |> dplyr::group_by(Significance, Chromosome) |> dplyr::summarise(N = dplyr::n()) |> 
      merge(LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(Total_N = dplyr::n())) |> dplyr::mutate(Prop = N/Total_N)
    
    ChromosomeSigBar <- ChrGroupingData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = rev(c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))), y = Prop, fill = factor(Significance, levels = c("Sig. Increase", "None", "Sig. Decrease")))) + 
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = c("Sig. Decrease" = "#02F", "None" = "#FFF", "Sig. Increase" = "#F10"), name = "Fold Change") +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~1-.)) + ggplot2::geom_hline(yintercept = seq(0.1, 0.9, by = 0.1), linetype = "dotted") + 
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::labs(x = "Chromosome", y = "Proportion of Proteins") + ggplot2::coord_flip()
    
    pdf("Chromosome_Boxplot.pdf", width = 18, height = 8)
    print(ChromosomeBoxplot)
    print(ChromosomeSigBar)
    Proteopedia::Reset_Dev()
    
    LimmaData[, Gene := gsub("\\;.*","", Gene)]
    LimmaData[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
    LimmaData[, Scien_AdjPValue := formatC(adj.P.Val, format = "e", digits = 2)]
    
    ChrInteractive <- highcharter::hchart(LimmaData, "scatter", 
                                          highcharter::hcaes(x = OrderID, y = Log2FC, group = Chromosome)) |> 
      highcharter::hc_chart(zoomType = "xy") |> 
      highcharter::hc_xAxis(title = list(text = list("Genome Position")), 
                            lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black", gridLineWidth = 0, categories = levels(LimmaData$Chromosome)) |>
      highcharter::hc_yAxis(title = list(text = paste0("Log2 Fold-Difference")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                            tickColor = "black", gridLineWidth = 0) |>
      highcharter::hc_tooltip(headerFormat = "",
                              pointFormat = "<b>{point.Gene} | {point.GeneGroup}</b><br>Log2FC: {point.Log2FC:.2f}<br>Adj. p-value: {point.Scien_AdjPValue:.2f}") |>
      highcharter::hc_plotOptions(scatter = list(jitter = list(x = 0, y = 0), states = list(hover = list(enabled = TRUE), inactive = list(enabled = FALSE)),   # Don"t dim inactive series
                                                 point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }")))))
    htmlwidgets::saveWidget(ChrInteractive, "Interactive_Chromosome_Plot.html")
  }
  # Chromosome Volcano Plots
  {
    if(dir.exists(paste0(getwd(),"/ChromosomeVolcanos")) == TRUE){
      unlink(paste0(getwd(),"/ChromosomeVolcanos"), recursive = TRUE)}
    dir.create(paste0(getwd(),"/ChromosomeVolcanos"), showWarnings = TRUE)
    setwd(paste0(getwd(),"/ChromosomeVolcanos"))
    
    for(i in levels(LimmaData$Chromosome)){
      message(paste0("Analysing Chr", i, " Proteins"))
      MeanLog2FC <- round(LimmaData[Chromosome == i, Log2FC] |> mean(), digits = 3)  
      
      pdf(paste0("Chr", gsub("/", "", i),"_Volcano.pdf"), width = 12, height = 8)
      print(
        LimmaData[Chromosome == i] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
          ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(-log10(P.Value)> -log10(0.05),as.character(Gene),""))) +
          Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black", alpha = 0.7) + 
          Proteopedia::Add_AbundanceAxes()
      )
      Proteopedia::Reset_Dev()
    }
  }
  message("Analysing Protein Complexes")
  {
    setwd(gsub("(.*)/.*.csv.*", "\\1", InputFile))
    ComplexPortalSummary <- LimmaData[, .N, by = Experimental_Evidence_ComplexPortal] 
    ComplexPortalSummary[, Database := "Complex_Portal"] |> data.table::setnames("Experimental_Evidence_ComplexPortal", "Complexed")
    CORUMSummary <- LimmaData[, .N, by = Experimental_Evidence_CORUM] 
    CORUMSummary[, Database := "CORUM"] |> data.table::setnames("Experimental_Evidence_CORUM", "Complexed")
    ComplexSummary <- ComplexPortalSummary |> rbind(CORUMSummary)
    
    pdf("ComplexBoxplot.pdf", width = 12, height = 10)
    print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_ComplexPortal, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
            ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + 
            ggplot2::labs(x = "Complex Portal Evidence", y = expression("Log"[2] ~ "FC in Protein Abundance")) + 
            ggplot2::geom_text(data = ComplexSummary[Database == "Complex_Portal"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
            LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_CORUM, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
            ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + 
            ggplot2::labs(x = "CORUM Evidence", y = expression("Log"[2] ~ "FC in Protein Abundance")) + 
            ggplot2::geom_text(data = ComplexSummary[Database == "CORUM"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
            ggplot2::theme(axis.title.y = ggplot2::element_blank()))
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_StaticSILAC_Proteins <- function(InputDirectory, SubsetColour = "red"){
  start.time <- Sys.time()
  set.seed(123)
  for(Mapping in c("Total_Analysis", "Ratio_Analysis", "Heavy_Channel_Analysis", "Light_Channel_Analysis")){
    message("Loading Limma File")
    {
      if(Mapping %in% c("Heavy_Channel_Analysis", "Light_Channel_Analysis")){
        setwd(paste0(InputDirectory, "/", gsub(".*_(Channel_Analysis)", "\\1", Mapping)))
        LimmaData <- data.table::fread(list.files(pattern = paste0(gsub("(.*)_Channel_Analysis", "\\1", Mapping), "_Limma_Output.csv")))
        
        if (dir.exists(paste0(getwd(),"/", gsub("(.*)_Channel_Analysis", "\\1", Mapping), "_Mapping")) == T){
          unlink(paste0(getwd(),"/", gsub("(.*)_Channel_Analysis", "\\1", Mapping), "_Mapping"), recursive = T)
        }
        dir.create(paste0(getwd(),"/", gsub("(.*)_Channel_Analysis", "\\1", Mapping), "_Mapping"), showWarnings = T)
        setwd(paste0(getwd(),"/", gsub("(.*)_Channel_Analysis", "\\1", Mapping), "_Mapping"))
      } else {
        setwd(paste0(InputDirectory, "/", Mapping))
        LimmaData <- data.table::fread(list.files(pattern = "Limma_Output.csv"))
      }
      MappingRootDirectory <- getwd()
      
      LimmaData |> data.table::setnames(
        c(colnames(LimmaData)[grepl("protein.*group", ignore.case = T, colnames(LimmaData))], 
          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & !grepl("adj", ignore.case = T, colnames(LimmaData))], 
          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & grepl("adj", ignore.case = T, colnames(LimmaData))], 
          colnames(LimmaData)[grepl("protein.*desc", ignore.case = T, colnames(LimmaData))], 
          colnames(LimmaData)[grepl("Gene", ignore.case = T, colnames(LimmaData)) & !grepl("group", ignore.case = T, colnames(LimmaData))]),
        c("ProteinGroup", "P.Value", "adj.P.Val", "ProteinDescription", "Gene"))
      
      if(length(LimmaData$ProteinGroup[grepl("\\-", LimmaData$ProteinGroup)]) > 0){
        LimmaData$Isoforms <- 1
        for(i in 1:nrow(LimmaData)){
          if(length(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T)) > 0){
            LimmaData$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          } else {LimmaData$Isoforms[i] <- 1}
        }    
        LimmaData[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      }
      LimmaData[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
      
      if (dir.exists(paste0(getwd(),"/Subset_Output")) == T){
        unlink(paste0(getwd(),"/Subset_Output"), recursive = T)
      }
      dir.create(paste0(getwd(),"/Subset_Output"), showWarnings = T)
      
      message("Importing Proteopedia")
      setwd(paste0(getwd(),"/Subset_Output"))
      LimmaData <- LimmaData |> merge(Proteopedia::Proteopedia, all.x = T) |> data.table::data.table()
      LimmaData$ENTREZID[is.na(LimmaData$ENTREZID)] <- "Unmapped"
      LimmaData$Strand[is.na(LimmaData$Strand)] <- "*"
      LimmaData$Chromosome[is.na(LimmaData$Chromosome)] <- "Unmapped"
      LimmaData$MedianLociStart[is.na(LimmaData$MedianLociStart)] <- "Unmapped"
      LimmaData$Deg_Profile[is.na(LimmaData$Deg_Profile)] <- "UN"
      LimmaData$Experimental_Evidence_ComplexPortal[is.na(LimmaData$Experimental_Evidence_ComplexPortal)] <- "No"
      LimmaData$Experimental_Evidence_CORUM[is.na(LimmaData$Experimental_Evidence_CORUM)] <- "No"
      LimmaData$N_ComplexPortal[is.na(LimmaData$N_ComplexPortal)] <- 0
      LimmaData$N_CORUM[is.na(LimmaData$N_CORUM)] <- 0
      LimmaData[, which(colnames(LimmaData) == "ER"):ncol(LimmaData)][is.na(LimmaData[, which(colnames(LimmaData) == "ER"):ncol(LimmaData)])] <- F
    }
    # Biochemical Trends (Numerical)
    {
      for(ColIndex in c(which(colnames(LimmaData) == "Length"), which(colnames(LimmaData) == "Length"):ncol(LimmaData))){
        if(is.numeric(LimmaData[, get(colnames(LimmaData)[ColIndex])]) & !grepl("N_", colnames(LimmaData)[ColIndex])){
          message(paste0("Analysing ", colnames(LimmaData)[ColIndex]), " Trend")
          SubsetData = LimmaData[, .(ProteinGroup, Sequence, Length, Log2FC, get(colnames(LimmaData)[ColIndex]))]
          SubsetData |> data.table::setnames("V5", "Subset")
          
          pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[2,4]
          
          TrendPlot <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Subset, y = Log2FC)) + 
            ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Add_Rsq() +
            ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Log2FC], na.rm = T)*0.93, 
                              label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) +
            ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2]~"FC in Heavy:Light Protein LFQ Ratio")) + 
            ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
          
          pdf(paste0(colnames(LimmaData)[ColIndex],"_Trend.pdf"), width = 12, height = 10)
          print(TrendPlot)
          print(TrendPlot + ggplot2::labs(x = "", y = ""))
          Proteopedia::Reset_Dev()
        }
      }
    }
    # Cellular Subsets (TRUE/FALSE-Based)
    {
      for(ColIndex in which(colnames(LimmaData) == "ER"):ncol(LimmaData)){
        if(is.logical(LimmaData[, get(colnames(LimmaData)[ColIndex])])){
          message(paste0("Analysing ", colnames(LimmaData)[ColIndex]), " Proteins")
          SubsetData = LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, get(colnames(LimmaData)[ColIndex]))]
          SubsetData |> data.table::setnames("V5", "Subset")
          
          MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)
          pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[4]
          
          if(Mapping == "Ratio_Analysis"){
            Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
              Proteopedia::Add_NotSigBox() + Proteopedia::Add_IsotopeRatioAxes()
            
            Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
              ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() +
              ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
              ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
              ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T))*0.93, 
                                label = paste0("P-Value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
              Proteopedia::Add_IsotopeRatioAxes()
          } else {
            Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
              Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes()
            
            Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
              ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() +
              ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
              ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
              ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T))*0.93, 
                                label = paste0("P-Value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
              Proteopedia::Add_AbundanceAxes()
          }
          
          pdf(paste0(colnames(LimmaData)[ColIndex],"_Volcano.pdf"))
          print(Volcano)
          print(Volcano_Rug)
          print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message("Analysing Degradation Profiles")
    {      
      DegProfileSummary <- LimmaData[, .N, by = Deg_Profile]
      
      pdf("DegradataionProfileBoxplot.pdf", width = 12, height = 10)
      print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Deg_Profile, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
              ggpubr::geom_signif(comparison = list(c("NED", "ED"), c("NED", "UN"), c("ED", "UN")), y_position = c(1, 1.25, 1.5), tip_length = 0) + 
              ggplot2::labs(x = "Degradation Profile", y = ifelse(Mapping == "Ratio_Analysis", expression("Log"[2] ~ "FC in Heavy:Light Protein LFQ Ratio"), 
                                                                  expression("Log"[2] ~ "FC in Protein Abundance"))) + 
              ggplot2::geom_text(data = DegProfileSummary, ggplot2::aes(y = -1, x = Deg_Profile, label = paste0("N = ", N))))
      Proteopedia::Reset_Dev()
      
      for(Deg_Type in unique(LimmaData[!is.na(Deg_Profile), Deg_Profile])){
        MeanLog2FC <- LimmaData[Deg_Profile == Deg_Type, Log2FC] |> mean(na.rm = T)
        
        if(Mapping == "Ratio_Analysis"){
          Volcano <- LimmaData[Deg_Profile == Deg_Type] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=as.character(Gene))) + 
            Proteopedia::Add_NotSigBox() + Proteopedia::Add_IsotopeRatioAxes() 
        } else {
          Volcano <- LimmaData[Deg_Profile == Deg_Type] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=as.character(Gene))) + 
            Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() 
        }
        
        pval <- summary(stats::lm(Log2FC ~ Subset, data = LimmaData[, Subset := data.table::fifelse(Deg_Profile == Deg_Type, T, NA)]))$coefficients[4]
        
        if(Mapping == "Ratio_Analysis"){
          Volcano_Rug <- LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") +    
            ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
            Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
            ggplot2::geom_rug(alpha = ifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
            ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(LimmaData[, P.Value], na.rm = T))*0.93, 
                              label = paste0("P-value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
            Proteopedia::Add_IsotopeRatioAxes()
        } else {
          Volcano_Rug <- LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") +    
            ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
            Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
            ggplot2::geom_rug(alpha = ifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
            ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(LimmaData[, P.Value], na.rm = T))*0.93, 
                              label = paste0("P-value: ",ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), colour = SubsetColour, size = 4) +
            Proteopedia::Add_AbundanceAxes()
        }
        
        pdf(paste0(Deg_Type, "_DegProfile_Volcano.pdf"))  
        print(Volcano)
        print(Volcano_Rug)
        print(Volcano_Rug + ggplot2::labs(x= "", y= ""))
        Proteopedia::Reset_Dev()
      }
    }
    message("Summarising Chromosome-Based Data")
    {
      setwd(MappingRootDirectory)
      ChromosomeSummary <- data.table::data.table("Chromosome" = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped")) |> 
        merge(LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(MeanLog2FC = mean(Log2FC), Mapped = dplyr::n()) |> 
                merge(Proteopedia::Proteopedia |> dplyr::group_by(Chromosome) |> dplyr::summarise(Total = dplyr::n())))
      ChromosomeSummary[, Coverage := Mapped/Total]
      
      pdf("Chromosome_Coverage.pdf", width = 18, height = 8)
      print(ChromosomeSummary |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), 
                                                              y = Coverage*100, fill = Coverage*100)) +
              ggplot2::geom_bar(stat= "identity") + ggplot2::scale_fill_viridis_c(guide = "none") + 
              ggplot2::labs(x = "Chromosome", y = "Protein Coverage (%)")
      )
      Proteopedia::Reset_Dev()
      
      ChromosomeSummary$Wilcoxon_p = 0
      rownumber = 1
      for(i in ChromosomeSummary$Chromosome){
        ChromosomeSummary$Wilcoxon_p[rownumber] <- wilcox.test(LimmaData[,Log2FC], LimmaData[Chromosome == i, Log2FC])$p.value
        rownumber = rownumber + 1
      }
      ChromosomeSummary[, SigSymbol := data.table::fifelse(Wilcoxon_p < 0.001, "***", data.table::fifelse(Wilcoxon_p < 0.01, "**",
                                                                                                          data.table::fifelse(Wilcoxon_p < 0.05, "*","")))]
      ChromosomeSummary[, Buffering := data.table::fifelse(MeanLog2FC < 0, paste0(">",100,"%"), 
                                                           paste0(100 - (round((MeanLog2FC/log2(3/2))*100, digits = 2)),"%"))]
      data.table::fwrite(ChromosomeSummary, "Chromosomal_Summary.csv")
      
      WhiskerTop <- LimmaData[, Log2FC] |> stats::quantile(0.75, na.rm = T) + stats::IQR(LimmaData[, Log2FC], na.rm = T)*1.5
      WhiskerBottom <- LimmaData[, Log2FC] |> stats::quantile(0.25, na.rm = T) - stats::IQR(LimmaData[, Log2FC], na.rm = T)*1.5
      
      LimmaData[, Chromosome := factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))]
      LimmaData <- LimmaData |> dplyr::arrange(Chromosome, MedianLociStart)
      LimmaData[, OrderID := as.numeric(rownames(LimmaData))]
      
      ChromosomeBorders <- LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(N_Proteins = dplyr::n()) |> 
        dplyr::mutate(Upper = cumsum(N_Proteins)) |> data.table::data.table()
      ChromosomeBorders$Midpoint <- 0
      for(i in 1:nrow(ChromosomeBorders)){
        ChromosomeBorders$Midpoint[i] <- ifelse(i > 1, ChromosomeBorders$Upper[i-1] + ChromosomeBorders$N_Proteins[i]/2, ChromosomeBorders$N_Proteins[i]/2)
      }
      ChromosomeBorders <- ChromosomeBorders[Chromosome %!in% c("Unmapped", "X/Y", "Y", "M")]
      
      ChromosomeCols <- rep(RColorBrewer::brewer.pal(n = 8, "Set2"), length.out = length(c(seq(1:22), "X")))
      names(ChromosomeCols) <- c(seq(1:22), "X")
      
      ChromosomeDotData <- LimmaData[Chromosome %!in% c("X/Y", "Y", "M", "Unmapped")] |> 
        merge(Proteopedia::Proteopedia[, .(ProteinGroup, Gene, Chromosome)], all = T, by = c("ProteinGroup", "Gene", "Chromosome"))
      ChromosomeDotData$Log2FC[is.na(ChromosomeDotData$Log2FC)] <- 0
      ChromosomeDotData <- ChromosomeDotData |> merge(data.table::data.table("Colour" = ChromosomeCols, "Chromosome" = names(ChromosomeCols)), by = "Chromosome")
      ChromosomeDotData[, Colour := data.table::fifelse(Log2FC == 0, "white", Colour)]
      
      ChromosomeDotplot <- ChromosomeDotData |> ggplot2::ggplot(ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) + 
        ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "blue") + 
        ggplot2::scale_colour_identity() +
        ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), 
                                                                  label = ifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))), 
                           colour = "black", angle = 90) +
        ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + ggplot2::scale_x_continuous(expand = 0) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
        ggplot2::labs(x = "Chromosome", y = ifelse(Mapping == "Ratio_Analysis", expression("Log"[2]~ "FC in Heavy:Light Protein LFQ Ratio"), 
                                                   expression("Log"[2]~ "FC in Protein Abundance"))) + ggplot2::guides(fill="none")
      
      pdf("Chromosome_Dotplot.pdf", width = 18, height = 8)
      print(ChromosomeDotplot)
      Proteopedia::Reset_Dev()
      
      ChromosomeBoxplot <- LimmaData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "Y", "M")), y = Log2FC)) + 
        ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop*1.25, label = SigSymbol), colour = "black") + 
        ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom*1.25, label = Buffering), colour = "black", size = 4) + 
        ggplot2::geom_boxplot(colour = "black", fill = NA, outliers = F, alpha = 0.7) + ggplot2::coord_cartesian(ylim = c(WhiskerTop*1.5, WhiskerBottom*1.5)) +
        ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + 
        ggplot2::labs(x = "Chromosome", y = ifelse(Mapping == "Ratio_Analysis", expression("Log"[2]~ "FC in Heavy:Light Protein LFQ Ratio"), 
                                                   expression("Log"[2]~ "FC in Protein Abundance"))) + ggplot2::guides(fill="none")
      
      ChrGroupingData <- LimmaData |> dplyr::group_by(Significance, Chromosome) |> dplyr::summarise(N = dplyr::n()) |> 
        merge(LimmaData |> dplyr::group_by(Chromosome) |> dplyr::summarise(Total_N = dplyr::n())) |> dplyr::mutate(Prop = N/Total_N)
      
      ChromosomeSigBar <- ChrGroupingData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = rev(c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))), y = Prop, fill = factor(Significance, levels = c("Sig. Increase", "None", "Sig. Decrease")))) + 
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = c("Sig. Decrease" = "#02F", "None" = "#FFF", "Sig. Increase" = "#F10"), name = "Fold Change") +
        ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~1-.)) + ggplot2::geom_hline(yintercept = seq(0.1, 0.9, by = 0.1), linetype = "dotted") + 
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::labs(x = "Chromosome", y = "Proportion of Proteins") + ggplot2::coord_flip()
      
      pdf("Chromosome_Boxplot.pdf", width = 18, height = 8)
      print(ChromosomeBoxplot)
      print(ChromosomeSigBar)
      Proteopedia::Reset_Dev()
      
      LimmaData[, Gene := gsub("\\;.*","", Gene)]
      LimmaData[, URL := paste0("https://www.uniprot.org/uniprotkb/", GeneGroup)]
      LimmaData[, Scien_AdjPValue := formatC(adj.P.Val, format = "e", digits = 2)]
      
      ChrInteractive <- highcharter::hchart(LimmaData, "scatter", 
                                            highcharter::hcaes(x = OrderID, y = Log2FC, group = Chromosome)) |> 
        highcharter::hc_chart(zoomType = "xy") |> 
        highcharter::hc_xAxis(title = list(text = list("Genome Position")), 
                              lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black", gridLineWidth = 0, categories = levels(LimmaData$Chromosome)) |>
        highcharter::hc_yAxis(title = list(text = paste0("Log2 Fold-Difference")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", 
                              tickColor = "black", gridLineWidth = 0) |>
        highcharter::hc_tooltip(headerFormat = "",
                                pointFormat = "<b>{point.Gene} | {point.GeneGroup}</b><br>Log2FC: {point.Log2FC:.2f}<br>Adj. p-value: {point.Scien_AdjPValue:.2f}") |>
        highcharter::hc_plotOptions(scatter = list(jitter = list(x = 0, y = 0), states = list(hover = list(enabled = TRUE), inactive = list(enabled = FALSE)),   # Don"t dim inactive series
                                                   point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }")))))
      htmlwidgets::saveWidget(ChrInteractive, "Interactive_Chromosome_Plot.html")
    }
    # Chromosome Volcano Plots
    {
      if(dir.exists(paste0(getwd(),"/ChromosomeVolcanos")) == TRUE){
        unlink(paste0(getwd(),"/ChromosomeVolcanos"), recursive = TRUE)}
      dir.create(paste0(getwd(),"/ChromosomeVolcanos"), showWarnings = TRUE)
      setwd(paste0(getwd(),"/ChromosomeVolcanos"))
      
      for(i in levels(LimmaData$Chromosome)){
        message(paste0("Analysing Chr", i, " Proteins"))
        MeanLog2FC <- round(LimmaData[Chromosome == i, Log2FC] |> mean(), digits = 3)  
        
        pdf(paste0("Chr", gsub("/", "", i),"_Volcano.pdf"), width = 12, height = 8)
        print(
          if(Mapping == "Ratio_Analysis"){
            LimmaData[Chromosome == i] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(-log10(P.Value)> -log10(0.05),as.character(Gene),""))) +
              Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black", alpha = 0.7) + 
              Proteopedia::Add_IsotopeRatioAxes()
          } else {
            LimmaData[Chromosome == i] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(-log10(P.Value)> -log10(0.05),as.character(Gene),""))) +
              Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black", alpha = 0.7) + 
              Proteopedia::Add_AbundanceAxes()
          }
        )
        Proteopedia::Reset_Dev()
      }
    }
    message("Analysing Protein Complexes")
    {
      setwd(MappingRootDirectory)
      ComplexPortalSummary <- LimmaData[, .N, by = Experimental_Evidence_ComplexPortal] 
      ComplexPortalSummary[, Database := "Complex_Portal"] |> data.table::setnames("Experimental_Evidence_ComplexPortal", "Complexed")
      CORUMSummary <- LimmaData[, .N, by = Experimental_Evidence_CORUM] 
      CORUMSummary[, Database := "CORUM"] |> data.table::setnames("Experimental_Evidence_CORUM", "Complexed")
      ComplexSummary <- ComplexPortalSummary |> rbind(CORUMSummary)
      
      pdf("ComplexBoxplot.pdf", width = 12, height = 10)
      print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_ComplexPortal, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
              ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + 
              ggplot2::labs(x = "Complex Portal Evidence", y = ifelse(Mapping == "Ratio_Analysis", expression("Log"[2]~ "FC in Heavy:Light Protein LFQ Ratio"), 
                                                                      expression("Log"[2]~ "FC in Protein Abundance"))) + 
              ggplot2::geom_text(data = ComplexSummary[Database == "Complex_Portal"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
              LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_CORUM, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
              ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + 
              ggplot2::labs(x = "CORUM Evidence", y = ifelse(Mapping == "Ratio_Analysis", expression("Log"[2]~ "FC in Heavy:Light Protein LFQ Ratio"), 
                                                             expression("Log"[2]~ "FC in Protein Abundance"))) + 
              ggplot2::geom_text(data = ComplexSummary[Database == "CORUM"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
              ggplot2::theme(axis.title.y = ggplot2::element_blank()))
      Proteopedia::Reset_Dev()
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_TimecourseSILAC_Proteins <- function(InputDirectory, SubsetColour = "red"){
  set.seed(123)
  start.time <- Sys.time()
  for(Parameter in c("KlossL", "KlossH", "Ksyn")){
    message(Parameter, ": Loading Limma File")
    {
      setwd(InputDirectory)
      if(Parameter == "KlossL"){LimmaData <- data.table::fread("LightParameters.csv")
      } else if(Parameter == "KlossH"){LimmaData <- data.table::fread("HeavyParameters.csv")[Parameter == "KlossH"]
      } else {LimmaData <- data.table::fread("HeavyParameters.csv")[Parameter == "Ksyn"]} 
      
      if (dir.exists(paste0(getwd(),"/",Parameter,"_SubsetOutput")) == TRUE){unlink(paste0(getwd(),"/",Parameter,"_SubsetOutput"), recursive = TRUE)}
      dir.create(paste0(getwd(),"/",Parameter,"_SubsetOutput"), showWarnings = TRUE)
      setwd(paste0(getwd(),"/",Parameter,"_SubsetOutput"))
      
      LimmaData |> data.table::setnames(c(colnames(LimmaData)[grepl("protein.*group", ignore.case = T, colnames(LimmaData))], 
                                          colnames(LimmaData)[grepl("^p.*val.*", ignore.case = T, colnames(LimmaData)) & !grepl("adj", ignore.case = T, colnames(LimmaData))], 
                                          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & grepl("adj", ignore.case = T, colnames(LimmaData))], 
                                          colnames(LimmaData)[grepl("Gene", ignore.case = T, colnames(LimmaData)) & !grepl("group", ignore.case = T, colnames(LimmaData))]), 
                                        c("ProteinGroup", "P.Value", "adj.P.Val", "Gene"), skip_absent = T)
      if(length(LimmaData$ProteinGroup[grepl("\\-", LimmaData$ProteinGroup)]) > 0){
        LimmaData$Isoforms <- 1
        for(i in 1:nrow(LimmaData)){
          if(length(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T)) > 0){
            LimmaData$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
          } else {LimmaData$Isoforms[i] <- 1}
        }
        LimmaData[, `:=`(ProteinGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
      }
      LimmaData[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
      LimmaData[, `:=`(Significance, data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", "None")))]
      message(Parameter, ": Importing Proteopedia")
      LimmaData <- data.table::data.table(data.table::merge.data.table(LimmaData, Proteopedia::Proteopedia, all.x = T))
      LimmaData$ENTREZID[is.na(LimmaData$ENTREZID)] <- "Unmapped"
      LimmaData$Strand[is.na(LimmaData$Strand)] <- "*"
      LimmaData$Chromosome[is.na(LimmaData$Chromosome)] <- "Unmapped"
      LimmaData$MedianLociStart[is.na(LimmaData$MedianLociStart)] <- "Unmapped"
      LimmaData$Deg_Profile[is.na(LimmaData$Deg_Profile)] <- "UN"
      LimmaData$Experimental_Evidence_ComplexPortal[is.na(LimmaData$Experimental_Evidence_ComplexPortal)] <- "No"
      LimmaData$Experimental_Evidence_CORUM[is.na(LimmaData$Experimental_Evidence_CORUM)] <- "No"
      LimmaData$N_ComplexPortal[is.na(LimmaData$N_ComplexPortal)] <- 0
      LimmaData$N_CORUM[is.na(LimmaData$N_CORUM)] <- 0
      LimmaData[, which(colnames(LimmaData) == "ER"):ncol(LimmaData)][is.na(LimmaData[,which(colnames(LimmaData) == "ER"):ncol(LimmaData)])] <- F
    }
    message(Parameter, ": Analysing Trends")
    {
      for(ColIndex in c(which(colnames(LimmaData) == "Length"), which(colnames(LimmaData) == "Length"):ncol(LimmaData))) {
        if(is.numeric(LimmaData[, get(colnames(LimmaData)[ColIndex])]) & !grepl("N_", colnames(LimmaData)[ColIndex])) {
          message(paste0(Parameter, ": Analysing ", colnames(LimmaData)[ColIndex]), " Trend")
          SubsetData = LimmaData[, .(ProteinGroup, Sequence, Length, Log2FC, get(colnames(LimmaData)[ColIndex]))]
          data.table::setnames(SubsetData, "V5", "Subset")
          pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[2, 4]
          
          if(Parameter == "Ksyn"){
            TrendPlot <- ggplot2::ggplot(SubsetData, ggplot2::aes(x = Subset, y = Log2FC)) + 
              ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Proteopedia::Add_Rsq() + 
              ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Log2FC], na.rm = T) * 0.93, 
                                label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) + 
              ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + 
              ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
          } else {
            TrendPlot <- ggplot2::ggplot(SubsetData, ggplot2::aes(x = Subset, y = Log2FC)) + 
              ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Proteopedia::Add_Rsq() + 
              ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Log2FC], na.rm = T) * 0.93, 
                                label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) + 
              ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + 
              ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
          }
          pdf(paste0(colnames(LimmaData)[ColIndex], "_Trend.pdf"), width = 12, height = 10)
          print(TrendPlot)
          print(TrendPlot + ggplot2::labs(x = "", y = ""))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(Parameter, ": Analysing Subsets")
    {
      if(Parameter == "Ksyn"){
        for(ColIndex in which(colnames(LimmaData) == "ER"):ncol(LimmaData)) {
          if(is.logical(LimmaData[, get(colnames(LimmaData)[ColIndex])])){
            message(paste0(Parameter, ": Analysing ", colnames(LimmaData)[ColIndex]), " Proteins")
            SubsetData <- LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, get(colnames(LimmaData)[ColIndex]))]
            data.table::setnames(SubsetData, "V5", "Subset")
            MeanLog2FC <- mean(SubsetData[Subset == T, Log2FC], na.rm = T)
            Volcano <- ggplot2::ggplot(SubsetData[Subset == T], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
              Proteopedia::Add_NotSigBox() + Proteopedia::Add_KsynAxes(scale = "Log2FC")
            pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[4]
            Volcano_Rug <- ggplot2::ggplot(SubsetData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
              ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() + 
              ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) + 
              ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
              ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T)) * 0.93, 
                                label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), 
                                colour = SubsetColour, size = 4) + Proteopedia::Add_KsynAxes(scale = "Log2FC")
            pdf(paste0(colnames(LimmaData)[ColIndex], "_Volcano.pdf"))
            print(Volcano)
            print(Volcano_Rug)
            print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
            Proteopedia::Reset_Dev()
          }
        }
      } else {
        for(ColIndex in which(colnames(LimmaData) == "ER"):ncol(LimmaData)) {
          if(is.logical(LimmaData[, get(colnames(LimmaData)[ColIndex])])){
            message(paste0(Parameter, ": Analysing ", colnames(LimmaData)[ColIndex]), " Proteins")
            SubsetData <- LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, get(colnames(LimmaData)[ColIndex]))]
            data.table::setnames(SubsetData, "V5", "Subset")
            MeanLog2FC <- mean(SubsetData[Subset == T, Log2FC], na.rm = T)
            Volcano <- ggplot2::ggplot(SubsetData[Subset == T], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
              Proteopedia::Add_NotSigBox() + Proteopedia::Add_KlossAxes(scale = "Log2FC")
            pval <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData))$coefficients[4]
            Volcano_Rug <- ggplot2::ggplot(SubsetData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
              ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() + 
              ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) + 
              ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
              ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T)) * 0.93, 
                                label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), 
                                colour = SubsetColour, size = 4) + Proteopedia::Add_KlossAxes(scale = "Log2FC")
            pdf(paste0(colnames(LimmaData)[ColIndex], "_Volcano.pdf"))
            print(Volcano)
            print(Volcano_Rug)
            print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
            Proteopedia::Reset_Dev()
          }
        }
      }
    }
    message(Parameter, ": Analysing Degradation Profiles")
    {
      DegProfileSummary <- LimmaData[, .N, by = Deg_Profile]
      
      if(Parameter == "Ksyn"){
        pdf("DegradataionProfileBoxplot.pdf", width = 12, height = 10)
        print(ggplot2::ggplot(LimmaData, ggplot2::aes(x = Deg_Profile, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
                ggpubr::geom_signif(comparison = list(c("NED", "ED"), c("NED", "UN"), c("ED", "UN")), 
                                    y_position = c(quantile(LimmaData[, Log2FC], 0.75), quantile(LimmaData[, Log2FC], 0.75)*1.2, 
                                                   quantile(LimmaData[, Log2FC], 0.75)*1.4), tip_length = 0) + 
                ggplot2::labs(x = "Degradation Profile", y = expression("Log"[2] ~ "FC in Protein k"[syn])) + 
                ggplot2::geom_text(data = DegProfileSummary, ggplot2::aes(y = quantile(LimmaData[, Log2FC], 0.01)*1.1, x = Deg_Profile, label = paste0("N = ", N))))
        Proteopedia::Reset_Dev()
        for(Deg_Type in unique(LimmaData[!is.na(Deg_Profile), Deg_Profile])) {
          MeanLog2FC <- mean(LimmaData[Deg_Profile == Deg_Type, Log2FC], na.rm = T)
          Volcano <- ggplot2::ggplot(LimmaData[Deg_Profile == Deg_Type], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
            Proteopedia::Add_NotSigBox() + Proteopedia::Add_KsynAxes(scale = "Log2FC")
          pval <- summary(stats::lm(Log2FC ~ Subset, data = LimmaData[, `:=`(Subset, data.table::fifelse(Deg_Profile == Deg_Type, T, NA))]))$coefficients[4]
          Volcano_Rug <- ggplot2::ggplot(LimmaData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) + 
            Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) + 
            ggplot2::geom_rug(alpha = data.table::fifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
            ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(LimmaData[, P.Value], na.rm = T)) * 0.93, 
                              label = paste0("P-value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), 
                              colour = SubsetColour, size = 4) + Proteopedia::Add_KsynAxes(scale = "Log2FC")
          pdf(paste0(Deg_Type, "_DegProfile_Volcano.pdf"))
          print(Volcano)
          print(Volcano_Rug)
          print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
          Proteopedia::Reset_Dev()
        }
      } else {
        pdf("DegradataionProfileBoxplot.pdf", width = 12, height = 10)
        print(ggplot2::ggplot(LimmaData, ggplot2::aes(x = Deg_Profile, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
                ggpubr::geom_signif(comparison = list(c("NED", "ED"), c("NED", "UN"), c("ED", "UN")), 
                                    y_position = c(quantile(LimmaData[, Log2FC], 0.75), quantile(LimmaData[, Log2FC], 0.75)*1.2, 
                                                   quantile(LimmaData[, Log2FC], 0.75)*1.4), tip_length = 0) + 
                ggplot2::labs(x = "Degradation Profile", y = expression("Log"[2] ~ "FC in Protein k"[loss])) + 
                ggplot2::geom_text(data = DegProfileSummary, ggplot2::aes(y = quantile(LimmaData[, Log2FC], 0.01)*1.1, x = Deg_Profile, label = paste0("N = ", N))))
        Proteopedia::Reset_Dev()
        for(Deg_Type in unique(LimmaData[!is.na(Deg_Profile), Deg_Profile])) {
          MeanLog2FC <- mean(LimmaData[Deg_Profile == Deg_Type, Log2FC], na.rm = T)
          Volcano <- ggplot2::ggplot(LimmaData[Deg_Profile == Deg_Type], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
            Proteopedia::Add_NotSigBox() + Proteopedia::Add_KlossAxes(scale = "Log2FC")
          pval <- summary(stats::lm(Log2FC ~ Subset, data = LimmaData[, `:=`(Subset, data.table::fifelse(Deg_Profile == Deg_Type, T, NA))]))$coefficients[4]
          Volcano_Rug <- ggplot2::ggplot(LimmaData, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) + 
            Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) + 
            ggplot2::geom_rug(alpha = data.table::fifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
            ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(LimmaData[, P.Value], na.rm = T)) * 0.93, 
                              label = paste0("P-value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), 
                              colour = SubsetColour, size = 4) + Proteopedia::Add_KlossAxes(scale = "Log2FC")
          pdf(paste0(Deg_Type, "_DegProfile_Volcano.pdf"))
          print(Volcano)
          print(Volcano_Rug)
          print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
          Proteopedia::Reset_Dev()
        }
      }
    }  
    message(Parameter, ": Summarising Chromosome-Based Data")
    {
      ChromosomeSummary <- data.table::data.table(Chromosome = c(seq(1:22), "X", "X/Y", "Y", "M")) |>
        data.table::merge.data.table(LimmaData[, .(MeanLog2FC = mean(Log2FC), Mapped = .N), by = Chromosome]) |>
        data.table::merge.data.table(Proteopedia::Proteopedia[, .(Total = .N), by = Chromosome])
      ChromosomeSummary[, `:=`(Coverage = Mapped/Total, Chromosome = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")))]
      
      ChromosomeSummary$Wilcoxon_p = 0
      rownumber = 1
      for(i in ChromosomeSummary$Chromosome) {
        ChromosomeSummary$Wilcoxon_p[rownumber] <- wilcox.test(LimmaData[, Log2FC], LimmaData[Chromosome == i, Log2FC])$p.value
        rownumber = rownumber + 1
      }
      ChromosomeSummary[, `:=`(SigSymbol, data.table::fifelse(Wilcoxon_p < 0.001, "***", data.table::fifelse(Wilcoxon_p < 0.01, "**", 
                                                                                                             data.table::fifelse(Wilcoxon_p < 0.05, "*", ""))))]
      ChromosomeSummary[, `:=`(Buffering, data.table::fifelse(MeanLog2FC < 0, paste0(">", 100, "%"), paste0(100 - (round((MeanLog2FC/log2(3/2)) * 100, digits = 2)), "%")))]
      data.table::fwrite(ChromosomeSummary |> data.table::setorder(Chromosome), "Chromosomal_Summary.csv")
      
      WhiskerTop <- stats::quantile(LimmaData[, Log2FC], 0.75, na.rm = T) + stats::IQR(LimmaData[, Log2FC], na.rm = T) * 1.5
      WhiskerBottom <- stats::quantile(LimmaData[, Log2FC], 0.25, na.rm = T) - stats::IQR(LimmaData[, Log2FC], na.rm = T) * 1.5
      LimmaData[, `:=`(Chromosome, factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped")))]
      LimmaData <- LimmaData |> data.table::setorder(Chromosome, MedianLociStart)
      LimmaData[, `:=`(OrderID, as.numeric(rownames(LimmaData)))]
      ChromosomeBorders <- data.table::data.table(dplyr::mutate(dplyr::summarise(dplyr::group_by(LimmaData, Chromosome), N_Proteins = dplyr::n()), 
                                                                Upper = cumsum(N_Proteins)))
      ChromosomeBorders$Midpoint <- 0
      for (i in 1:nrow(ChromosomeBorders)) {
        ChromosomeBorders$Midpoint[i] <- ifelse(i > 1, ChromosomeBorders$Upper[i - 1] + ChromosomeBorders$N_Proteins[i]/2, ChromosomeBorders$N_Proteins[i]/2)
      }
      ChromosomeBorders <- ChromosomeBorders[Chromosome %!in% c("Unmapped", "X/Y", "Y", "M")]
      ChromosomeCols <- rep(RColorBrewer::brewer.pal(n = 8, "Set2"), length.out = length(c(seq(1:22), "X")))
      names(ChromosomeCols) <- c(seq(1:22), "X")
      ChromosomeDotData <- merge(LimmaData[Chromosome %!in% c("X/Y", "Y", "M", "Unmapped")], 
                                 Proteopedia::Proteopedia[, .(ProteinGroup, Gene, Chromosome)], all = T, by = c("ProteinGroup", "Gene", "Chromosome"))
      ChromosomeDotData$Log2FC[is.na(ChromosomeDotData$Log2FC)] <- 0
      ChromosomeDotData <- merge(ChromosomeDotData, data.table::data.table(Colour = ChromosomeCols, Chromosome = names(ChromosomeCols)), by = "Chromosome")
      ChromosomeDotData[, `:=`(Colour, data.table::fifelse(Log2FC == 0, "white", Colour))]
      
      if(Parameter == "Ksyn"){
        ChromosomeDotplot <- ggplot2::ggplot(ChromosomeDotData, ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) + 
          ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "blue") + 
          ggplot2::scale_colour_identity() + ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), label = data.table::fifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))), 
                                                                colour = "black", angle = 90) + 
          ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + ggplot2::scale_x_continuous(expand = 0) + 
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + ggplot2::guides(fill = "none")
        ChromosomeBoxplot <- ggplot2::ggplot(LimmaData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Log2FC)) + 
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop * 1.25, label = SigSymbol), colour = "black") + 
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom * 1.25, label = Buffering), colour = "black", 
                             size = 4) + ggplot2::geom_boxplot(colour = "black", fill = NA, outliers = F, alpha = 0.7) + 
          ggplot2::coord_cartesian(ylim = c(WhiskerTop * 1.5, WhiskerBottom * 1.5)) + ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + 
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + ggplot2::guides(fill = "none")
      } else {
        ChromosomeDotplot <- ggplot2::ggplot(ChromosomeDotData, ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) + 
          ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "blue") + 
          ggplot2::scale_colour_identity() + ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), label = data.table::fifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))), 
                                                                colour = "black", angle = 90) + 
          ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + ggplot2::scale_x_continuous(expand = 0) + 
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + 
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + ggplot2::guides(fill = "none")
        ChromosomeBoxplot <- ggplot2::ggplot(LimmaData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Log2FC)) + 
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop * 1.25, label = SigSymbol), colour = "black") + 
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom * 1.25, label = Buffering), colour = "black", size = 4) + 
          ggplot2::geom_boxplot(colour = "black", fill = NA, outliers = F, alpha = 0.7) + 
          ggplot2::coord_cartesian(ylim = c(WhiskerTop * 1.5, WhiskerBottom * 1.5)) + ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "red") + 
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + ggplot2::guides(fill = "none")
      }
      
      ChrGroupingData <- LimmaData[, .(Sig_N = .N), by = .(Significance, Chromosome)] |> data.table::merge.data.table(LimmaData[, .(Total_N = .N), by = Chromosome])
      ChrGroupingData[, Prop := Sig_N/Total_N]
      
      ChromosomeSigBar <- ggplot2::ggplot(ChrGroupingData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = rev(c(seq(1:22), "X", "X/Y", "Y", "M"))), 
                                                                                                  y = Prop, fill = factor(Significance, levels = c("Sig. Increase", "None", "Sig. Decrease")))) + 
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = c(`Sig. Decrease` = "#02F", None = "#FFF", `Sig. Increase` = "#F10"), name = "Fold Change") + 
        ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~1 - .), expand = c(0, 0)) + ggplot2::geom_hline(yintercept = seq(0.1, 0.9, by = 0.1), linetype = "dotted") + 
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::labs(x = "Chromosome", y = "Proportion of Proteins") + ggplot2::coord_flip()
      
      pdf("ChromosomePlots.pdf", width = 18, height = 8)
      print(ggplot2::ggplot(ChromosomeSummary[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Coverage * 100, fill = Coverage * 100)) + ggplot2::geom_bar(stat = "identity") + 
              ggplot2::scale_fill_viridis_c(guide = "none") + ggplot2::labs(x = "Chromosome", y = "Protein Coverage (%)"))
      print(ChromosomeDotplot)
      print(ChromosomeBoxplot)
      print(ChromosomeSigBar)
      Proteopedia::Reset_Dev()
      
      if(dir.exists(paste0(getwd(), "/ChromosomeVolcanos")) == TRUE) {unlink(paste0(getwd(), "/ChromosomeVolcanos"), recursive = TRUE)}
      dir.create(paste0(getwd(), "/ChromosomeVolcanos"), showWarnings = TRUE)
      setwd(paste0(getwd(), "/ChromosomeVolcanos"))
      for(i in levels(LimmaData$Chromosome)){
        message(paste0(Parameter, ": Analysing Chr", i, " Proteins"))
        MeanLog2FC <- round(mean(LimmaData[Chromosome == i, Log2FC]), digits = 3)
        if(Parameter == "Ksyn"){
          pdf(paste0("Chr", gsub("/", "", i), "_Volcano.pdf"), width = 12, height = 8)
          print(ggplot2::ggplot(LimmaData[Chromosome == i], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
                  ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(-log10(P.Value) > -log10(0.05), as.character(Gene), ""))) + 
                  Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black", alpha = 0.7) + 
                  Proteopedia::Add_KsynAxes(scale = "Log2FC"))
          Proteopedia::Reset_Dev()
        } else {
          pdf(paste0("Chr", gsub("/", "", i), "_Volcano.pdf"), width = 12, height = 8)
          print(ggplot2::ggplot(LimmaData[Chromosome == i], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
                  ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(-log10(P.Value) > -log10(0.05), as.character(Gene), ""))) + 
                  Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "black", alpha = 0.7) + 
                  Proteopedia::Add_KlossAxes(scale = "Log2FC"))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(Parameter, ": Analysing Protein Complexes")
    {
      setwd(paste0(InputDirectory,"/",Parameter,"_SubsetOutput"))
      ComplexPortalSummary <- LimmaData[, .N, by = Experimental_Evidence_ComplexPortal]
      data.table::setnames(ComplexPortalSummary[, `:=`(Database, "Complex_Portal")], "Experimental_Evidence_ComplexPortal", "Complexed")
      CORUMSummary <- LimmaData[, .N, by = Experimental_Evidence_CORUM]
      data.table::setnames(CORUMSummary[, `:=`(Database, "CORUM")], "Experimental_Evidence_CORUM", "Complexed")
      ComplexSummary <- rbind(ComplexPortalSummary, CORUMSummary)
      if(Parameter == "Ksyn"){
        pdf("ComplexBoxplot.pdf", width = 12, height = 10)
        print(ggplot2::ggplot(LimmaData, ggplot2::aes(x = Experimental_Evidence_ComplexPortal, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
                ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = quantile(LimmaData[, Log2FC], 0.99), tip_length = 0) + 
                ggplot2::labs(x = "Complex Portal Evidence", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + 
                ggplot2::geom_text(data = ComplexSummary[Database == "Complex_Portal"], ggplot2::aes(y = quantile(LimmaData[, Log2FC], 0.01), x = Complexed, label = paste0("N = ", N))) + 
                ggplot2::ggplot(LimmaData, ggplot2::aes(x = Experimental_Evidence_CORUM, y = Log2FC)) + 
                ggplot2::geom_boxplot(outliers = F) + ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = quantile(LimmaData[, Log2FC], 0.75), tip_length = 0) + 
                ggplot2::labs(x = "CORUM Evidence", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~(")"))) + 
                ggplot2::geom_text(data = ComplexSummary[Database == "CORUM"], ggplot2::aes(y = quantile(LimmaData[, Log2FC], 0.01), x = Complexed, label = paste0("N = ", N))) + 
                ggplot2::theme(axis.title.y = ggplot2::element_blank()))
        Proteopedia::Reset_Dev()
      } else {
        pdf("ComplexBoxplot.pdf", width = 12, height = 10)
        print(ggplot2::ggplot(LimmaData, ggplot2::aes(x = Experimental_Evidence_ComplexPortal, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) + 
                ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = quantile(LimmaData[, Log2FC], 0.99), tip_length = 0) + 
                ggplot2::labs(x = "Complex Portal Evidence", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + 
                ggplot2::geom_text(data = ComplexSummary[Database == "Complex_Portal"], ggplot2::aes(y = quantile(LimmaData[, Log2FC], 0.01), x = Complexed, label = paste0("N = ", N))) + 
                ggplot2::ggplot(LimmaData, ggplot2::aes(x = Experimental_Evidence_CORUM, y = Log2FC)) + 
                ggplot2::geom_boxplot(outliers = F) + ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + 
                ggplot2::labs(x = "CORUM Evidence", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + 
                ggplot2::geom_text(data = ComplexSummary[Database == "CORUM"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) + 
                ggplot2::theme(axis.title.y = ggplot2::element_blank()))
        Proteopedia::Reset_Dev()
      }
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Plot_DetoxificationProteins <- function(InputFile){
  LimmaData <- data.table::fread(InputFile)
  LimmaData <- LimmaData |> data.table::merge.data.table(Proteopedia::DetoxiProt, all.x = T)
  LimmaData[is.na(Phase), Phase := "None"]
  LimmaData[is.na(Category), Category := "None"]
  
  LimmaSummary <- LimmaData[Category != "None", .(.N, Mean = mean(Log2FC, na.rm = T)), Category]
  
  pdf(paste0(gsub("(~.*/).*.csv", "\\1", InputFile), "DetoxiProtVolcano.pdf"), width = 12, height = 12)
  print(LimmaData[Category != "None"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Category)) + 
          Proteopedia::Add_NotSigBox() + ggplot2::geom_point(stroke = NA) + 
          ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(P.Value < 0.05, Gene, ""))) +
          ggplot2::labs(x = expression("Log"[2]~"FC in Protein Abundance"), y = expression("-Log"[10]~"P-Value")) + 
          ggplot2::geom_rug(data = LimmaSummary, ggplot2::aes(x = Mean, y = -0.1, colour = Category), sides = "b", size = 2) +
          ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
          ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal", legend.spacing.y = ggplot2::unit(0.1, "pt")))
  Proteopedia::Reset_Dev()
}

########### MS Analysis: Protein Data to Ontology Functions ####################################################################################################################################
#' @export
Order_Terms <- function(x){
  if(nrow(x) > 0){term_distance = as.dist(1-enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)
  if(length(term_distance) == 0){ordered_terms <- rownames(enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)}else{
    rownames(enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)[stats::hclust(term_distance, method = "ward.D")$order]}}
}
#' @export
Process_GSEAOutput <- function(GSEA_Output, GSEADatabase, PlotColour = "black"){
  GSEA_OutputData <- GSEA_Output |> data.frame() |> data.table::data.table()
  if(nrow(GSEA_OutputData)){
    GSEA_OutputData[, Ontology := GSEADatabase]
    GSEA_OutputTop30 <- GSEA_OutputData[pvalue < 0.05] |> dplyr::slice_max(order_by = rank, n = 30, with_ties = FALSE)
    
    GSEA_Dotplot <- GSEA_OutputData |> ggplot2::ggplot(ggplot2::aes(x = NES, y = -log10(pvalue), colour = PlotColour)) + 
      ggplot2::geom_point(size = 3) + ggplot2::coord_cartesian(x = c(-ceiling(max(abs(GSEA_OutputData$NES))), ceiling(max(abs(GSEA_OutputData$NES))))) + 
      ggplot2::labs(x = "Normalised Enrichment Score (NES)", y = expression("-Log"[10]~"P-Value")) + Proteopedia::Add_NotSigBox() + 
      ggplot2::scale_colour_manual(values = PlotColour, guide = "none")
    pdf(paste0("GSEA_", GSEADatabase, "_Dotplot.pdf"), width = 10, height = 8)
    print(GSEA_Dotplot)
    Proteopedia::Reset_Dev()
    
    GenesetsTop30 <- GSEA_OutputTop30 |> ggplot2::ggplot(ggplot2::aes(x = enrichmentScore, y = stringr::str_wrap(Description, 60), colour = pvalue, size = pvalue)) +
      ggplot2::geom_point() + ggplot2::geom_segment(aes(x = 0, xend = enrichmentScore, linewidth = 0.1), show.legend = FALSE) + ylab("Description") +
      ggplot2::scale_size_continuous(range = c(12,5), name = "P-Value") + ggplot2::scale_colour_viridis_c(name = "P-Value") +
      ggplot2::scale_x_continuous(limits = c(-ceiling(max(abs(GSEA_OutputTop30$enrichmentScore))), ceiling(max(abs(GSEA_OutputTop30$enrichmentScore))))) + 
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::labs(x = "Enrichment Factor", y = "")
    pdf(paste0("Top30_", GSEADatabase, "_Lollipop.pdf"), width = 20, height = 12)
    print(GenesetsTop30)
    Proteopedia::Reset_Dev()
  } else {
    GSEA_OutputData = NULL
  }
  return(GSEA_OutputData)
}
#' @export
Process_ORAOutput <- function(ORA_Output, ORADatabase, PlotColour = "black"){
  ORA_OutputData <- ORA_Output |> data.frame() |> data.table::data.table()
  if(nrow(ORA_OutputData) > 0){
    ORA_OutputData[, Ontology := ORADatabase]
    ORAOutputTop30 <- ORA_OutputData[pvalue < 0.05] |> dplyr::slice_max(order_by = p.adjust*RichFactor, n = 30, with_ties = FALSE)
    
    GenesetsTop30 <- ORAOutputTop30 |> ggplot2::ggplot(ggplot2::aes(x = RichFactor, y = stringr::str_wrap(Description, 60), colour = Count, size = pvalue)) +
      ggplot2::geom_point() + ggplot2::geom_segment(aes(x = 0, xend = RichFactor, linewidth = 0.1), show.legend = FALSE) + ylab("Description") +
      ggplot2::scale_size_continuous(range = c(12,5), name = "P-Value") + ggplot2::scale_colour_viridis_c(name = "Count") +
      ggplot2::scale_x_continuous(limits = c(0, ceiling(max(abs(ORAOutputTop30$RichFactor))))) + 
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::labs(x = "Over-Representation Score", y = "")
    pdf(paste0("Top30_", ORADatabase, "_Lollipop.pdf"), width = 20, height = 12)
    print(GenesetsTop30)
    Proteopedia::Reset_Dev()
    return(ORA_OutputData)
  } else {return(data.table("ID" = 0, "Description" = 0, "GeneRatio" = 0, "BgRatio" = 0, "RichFactor" = 0, "FoldEnrichment" = 0,
                            "zScore" = NA, "pvalue" = 0, "p.adjust" = 0, "qvalue" = 0, "geneID" = 0, "Count" = 0))}
}
#' @export
Run_GSEA <- function(InputFile, OutputName = "GSEA_Output"){
  start.time <- Sys.time()
  setwd(gsub("(~/.*)/.*.csv", "\\1", InputFile))
  InputData <- data.table::fread(InputFile)
  if(!is.null(InputData$Parameter)){
    for(Parameter in unique(InputData$Parameter)){
      setwd(gsub("(~/.*)/.*.csv", "\\1", InputFile))
      if(dir.exists(paste0(getwd(),"/", Parameter, OutputName)) == TRUE){unlink(paste0(getwd(),"/", Parameter, OutputName), recursive = TRUE)}
      dir.create(paste0(getwd(),"/", Parameter, OutputName), showWarnings = TRUE)
      setwd(paste0(getwd(),"/", Parameter, OutputName))
      message(paste0(Parameter, ": Loading & Formatting Data"))
      {
        InputData <- InputData[Parameter == Parameter, .(Gene, Log2FC)] |> data.table::merge.data.table(clusterProfiler::bitr(InputData$Gene, from = "ALIAS", to = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db) |> 
                                                             data.table::setnames("ALIAS", "Gene"), by = "Gene")
        
        GSEAData <- InputData$Log2FC[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
        names(GSEAData) <- InputData$ENTREZID[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
        GSEAData <- BiocGenerics::sort(GSEAData[!duplicated(names(GSEAData))], decreasing = TRUE)
      }
      message(paste0(Parameter, ": Performing Enrichment Analyses"))
      {
        gseaBP <- clusterProfiler::gseGO(geneList= GSEAData, ont = "BP", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
        gseaMF <- clusterProfiler::gseGO(geneList= GSEAData, ont = "MF", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
        gseaCC <- clusterProfiler::gseGO(geneList= GSEAData, ont = "CC", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
        gseaKEGG <- clusterProfiler::gseKEGG(geneList = GSEAData, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1)
        gseaMKEGG <- clusterProfiler::gseMKEGG(geneList = GSEAData, organism = "hsa", keyType  = "ncbi-geneid", pvalueCutoff = 1)
        gseaReactome <- ReactomePA::gsePathway(geneList = GSEAData, organism = "human", pvalueCutoff = 1)
        gseaDisease <- DOSE::gseDO(geneList = GSEAData, organism = "human", ont = "HDO", pvalueCutoff = 1)
        gseaPhenotype <- DOSE::gseDO(geneList = GSEAData, organism = "human", ont = "HPO", pvalueCutoff = 1)
        gseaDisGeNET <- DOSE::gseDGN(geneList = GSEAData, pvalueCutoff = 1)
        gseaNCG <- DOSE::gseNCG(geneList = GSEAData, organism = "human", pvalueCutoff = 1)
      }
      message(paste0(Parameter, ": Plotting & Compiling Analyses"))
      {
        MergedGSEA <- data.table::rbindlist(list(Proteopedia::Process_GSEAOutput(gseaBP, "Biological_Process"), Proteopedia::Process_GSEAOutput(gseaMF, "Molecular_Function"), 
                                     Proteopedia::Process_GSEAOutput(gseaCC, "Cellular_Component"), Proteopedia::Process_GSEAOutput(gseaKEGG, "KEGG"), 
                                     Proteopedia::Process_GSEAOutput(gseaMKEGG, "MKEGG"), Proteopedia::Process_GSEAOutput(gseaReactome, "Reactome"), 
                                     Proteopedia::Process_GSEAOutput(gseaDisease, "Disease"), Proteopedia::Process_GSEAOutput(gseaPhenotype, "Phenotype"), 
                                     Proteopedia::Process_GSEAOutput(gseaDisGeNET, "DisGeNET"), Proteopedia::Process_GSEAOutput(gseaNCG, "Network_of_Cancer_Genes")))
        MergedGSEA[, N_Genes := 0]
        for(GSEARow in 1:nrow(MergedGSEA)){
          MergedGSEA$N_Genes[GSEARow] <- length(str_split(MergedGSEA$core_enrichment[GSEARow], "/")[[1]])
        }
        MergedGSEA[, PropGenes := N_Genes/setSize]   
        data.table::fwrite(MergedGSEA, file = "GSEA_Results.csv")
      }
      message(paste0(Parameter, ": Formatting Data For HTML Plots"))
      {
        MergedGSEA <- MergedGSEA |> dplyr::arrange(factor(Ontology, levels = c("Cellular_Component", "Biological_Process", "Molecular_Function", "KEGG", "MKEGG", 
                                                                         "Reactome", "Disease","Phenotype", "DisGeNET", "Network_of_Cancer_Genes"))) 
        
        #MergedGSEA[, Description := factor(Description, levels = unique(c(Proteopedia::Order_Terms(gseaBP), Proteopedia::Order_Terms(gseaMF), 
        #                                                                  Proteopedia::Order_Terms(gseaCC), Proteopedia::Order_Terms(gseaKEGG), 
        #                                                                  Proteopedia::Order_Terms(gseaMKEGG), Proteopedia::Order_Terms(gseaReactome), 
        #                                                                  Proteopedia::Order_Terms(gseaDisease), Proteopedia::Order_Terms(gseaPhenotype), 
        #                                                                  Proteopedia::Order_Terms(gseaDisGeNET), Proteopedia::Order_Terms(gseaNCG))))]
        MergedGSEA[, TermIndex := 1:nrow(MergedGSEA)]
        MergedGSEA[, NegLog10PAdj := -log10(p.adjust)]
        MergedGSEA[, SciPVal := formatC(p.adjust, format = "e", digits = 2)]
      }
      message(paste0(Parameter, ": Exporting HTML Plots"))
      {
        GSEA_Significance <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = TermIndex, y = -log10(p.adjust), group = Ontology)) |>
          highcharter::hc_chart(zoomType = "xy") |> 
          highcharter::hc_xAxis(title = list(text = NULL), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", tickWidth = 0 ) |>
          highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                                gridLineWidth = 0 ) |>
          highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Normalised Enrichment Score: {point.NES:.2f}
                                  <br>P-Value: {point.SciPVal} <br>Gene count: {point.setSize}")
        htmlwidgets::saveWidget(GSEA_Significance, file = "Manhattan_Plot.html")
        
        GSEA_EnrichSig <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = NES, y = -log10(p.adjust), group = Ontology)) |>
          highcharter::hc_chart(zoomType = "xy") |>
          highcharter::hc_xAxis(title = list(text = "Normalised Enrichment Score"), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", 
                                tickWidth = 0 ) |> 
          highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                                gridLineWidth = 0 ) |> 
          highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Normalised Enrichment Score: {point.NES:.2f}
                                  <br>P-Value: {point.SciPVal} <br>Gene count: {point.setSize}")
        htmlwidgets::saveWidget(GSEA_EnrichSig, file = "Volcano_Plot.html")
      }
    }
  } else {
    if(dir.exists(paste0(getwd(),"/",OutputName)) == TRUE){unlink(paste0(getwd(),"/", OutputName), recursive = TRUE)}
    dir.create(paste0(getwd(),"/", OutputName), showWarnings = TRUE)
    setwd(paste0(getwd(),"/", OutputName))
    message("Loading & Formatting Data")
    {
      InputData <- InputData[, .(Gene, Log2FC)] |> data.table::merge.data.table(clusterProfiler::bitr(InputData$Gene, from = "ALIAS", to = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db) |> 
                                                                                  data.table::setnames("ALIAS", "Gene"), by = "Gene")
      
      GSEAData <- InputData$Log2FC[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
      names(GSEAData) <- InputData$ENTREZID[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
      GSEAData <- BiocGenerics::sort(GSEAData[!duplicated(names(GSEAData))], decreasing = TRUE)
    }
    message("Performing Enrichment Analyses")
    {
      gseaBP <- clusterProfiler::gseGO(geneList= GSEAData, ont = "BP", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
      gseaMF <- clusterProfiler::gseGO(geneList= GSEAData, ont = "MF", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
      gseaCC <- clusterProfiler::gseGO(geneList= GSEAData, ont = "CC", keyType = "ENTREZID", minGSSize = 3, maxGSSize = 800, verbose = TRUE, 
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 1)
      gseaKEGG <- clusterProfiler::gseKEGG(geneList = GSEAData, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1)
      gseaMKEGG <- clusterProfiler::gseMKEGG(geneList = GSEAData, organism = "hsa", keyType  = "ncbi-geneid", pvalueCutoff = 1)
      gseaReactome <- ReactomePA::gsePathway(geneList = GSEAData, organism = "human", pvalueCutoff = 1)
      gseaDisease <- DOSE::gseDO(geneList = GSEAData, organism = "human", ont = "HDO", pvalueCutoff = 1)
      gseaPhenotype <- DOSE::gseDO(geneList = GSEAData, organism = "human", ont = "HPO", pvalueCutoff = 1)
      gseaDisGeNET <- DOSE::gseDGN(geneList = GSEAData, pvalueCutoff = 1)
      gseaNCG <- DOSE::gseNCG(geneList = GSEAData, organism = "human", pvalueCutoff = 1)
    }
    message("Plotting & Compiling Analyses")
    {
      MergedGSEA <- data.table::rbindlist(list(Proteopedia::Process_GSEAOutput(gseaBP, "Biological_Process"), Proteopedia::Process_GSEAOutput(gseaMF, "Molecular_Function"), 
                                               Proteopedia::Process_GSEAOutput(gseaCC, "Cellular_Component"), Proteopedia::Process_GSEAOutput(gseaKEGG, "KEGG"), 
                                               Proteopedia::Process_GSEAOutput(gseaMKEGG, "MKEGG"), Proteopedia::Process_GSEAOutput(gseaReactome, "Reactome"), 
                                               Proteopedia::Process_GSEAOutput(gseaDisease, "Disease"), Proteopedia::Process_GSEAOutput(gseaPhenotype, "Phenotype"), 
                                               Proteopedia::Process_GSEAOutput(gseaDisGeNET, "DisGeNET"), Proteopedia::Process_GSEAOutput(gseaNCG, "Network_of_Cancer_Genes")))
      MergedGSEA[, N_Genes := 0]
      for(GSEARow in 1:nrow(MergedGSEA)){
        MergedGSEA$N_Genes[GSEARow] <- length(str_split(MergedGSEA$core_enrichment[GSEARow], "/")[[1]])
      }
      MergedGSEA[, PropGenes := N_Genes/setSize]   
      data.table::fwrite(MergedGSEA, file = "GSEA_Results.csv")
    }
    message("Formatting Data For HTML Plots")
    {
      MergedGSEA <- MergedGSEA |> dplyr::arrange(factor(Ontology, levels = c("Cellular_Component", "Biological_Process", "Molecular_Function", "KEGG", "MKEGG", 
                                                                             "Reactome", "Disease","Phenotype", "DisGeNET", "Network_of_Cancer_Genes"))) 
      
      #MergedGSEA[, Description := factor(Description, levels = unique(c(Proteopedia::Order_Terms(gseaBP), Proteopedia::Order_Terms(gseaMF), 
      #                                                                  Proteopedia::Order_Terms(gseaCC), Proteopedia::Order_Terms(gseaKEGG), 
      #                                                                  Proteopedia::Order_Terms(gseaMKEGG), Proteopedia::Order_Terms(gseaReactome), 
      #                                                                  Proteopedia::Order_Terms(gseaDisease), Proteopedia::Order_Terms(gseaPhenotype), 
      #                                                                  Proteopedia::Order_Terms(gseaDisGeNET), Proteopedia::Order_Terms(gseaNCG))))]
      MergedGSEA[, TermIndex := 1:nrow(MergedGSEA)]
      MergedGSEA[, NegLog10PAdj := -log10(p.adjust)]
      MergedGSEA[, SciPVal := formatC(p.adjust, format = "e", digits = 2)]
    }
    message("Exporting HTML Plots")
    {
      GSEA_Significance <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = TermIndex, y = -log10(p.adjust), group = Ontology)) |>
        highcharter::hc_chart(zoomType = "xy") |> 
        highcharter::hc_xAxis(title = list(text = NULL), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", tickWidth = 0 ) |>
        highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                              gridLineWidth = 0 ) |>
        highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Normalised Enrichment Score: {point.NES:.2f}
                                  <br>P-Value: {point.SciPVal} <br>Gene count: {point.setSize}")
      htmlwidgets::saveWidget(GSEA_Significance, file = "Manhattan_Plot.html")
      
      GSEA_EnrichSig <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = NES, y = -log10(p.adjust), group = Ontology)) |>
        highcharter::hc_chart(zoomType = "xy") |>
        highcharter::hc_xAxis(title = list(text = "Normalised Enrichment Score"), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", 
                              tickWidth = 0 ) |> 
        highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                              gridLineWidth = 0 ) |> 
        highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Normalised Enrichment Score: {point.NES:.2f}
                                  <br>P-Value: {point.SciPVal} <br>Gene count: {point.setSize}")
      htmlwidgets::saveWidget(GSEA_EnrichSig, file = "Volcano_Plot.html")
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Run_ORA <- function(SearchProteins, Proteome, OutputDirectory, OutputName = "ORA_Output"){
  start.time <- Sys.time()
  message("Loading & Formatting Data")
  {
    setwd(OutputDirectory)
    if(dir.exists(paste0(getwd(),"/", OutputName)) == TRUE){unlink(paste0(getwd(),"/",OutputName), recursive = TRUE)}
    dir.create(paste0(getwd(),"/", OutputName), showWarnings = TRUE)
    setwd(paste0(getwd(),"/", OutputName))
    
    Proteome <- clusterProfiler::bitr(Proteome, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)$ENTREZID
    Proteome <- unique(Proteome[which(!is.na(Proteome))])
    SearchProteins <- clusterProfiler::bitr(SearchProteins, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)$ENTREZID
  }
  message("Performing Over-Representation Analyses")
  {
    oraBP <- clusterProfiler::enrichGO(gene = SearchProteins, keyType = "ENTREZID", universe = Proteome, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP")
    oraMF <- clusterProfiler::enrichGO(gene = SearchProteins, keyType = "ENTREZID", universe = Proteome, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "MF") 
    oraCC <- clusterProfiler::enrichGO(gene = SearchProteins, keyType = "ENTREZID", universe = Proteome, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "CC") 
    oraKEGG <- clusterProfiler::enrichKEGG(gene = SearchProteins, organism = "hsa", universe = Proteome, pvalueCutoff = 1, keyType = "ncbi-geneid") 
    oraMKEGG <- clusterProfiler::enrichMKEGG(gene = SearchProteins, organism = "hsa", universe = Proteome, keyType = "ncbi-geneid", minGSSize = 3) 
    oraReactome <- ReactomePA::enrichPathway(gene = SearchProteins, universe = Proteome, organism = "human") 
    oraDisease <- DOSE::enrichDO(gene = SearchProteins, universe = Proteome, ont = "HDO")
    oraPhenotype <- DOSE::enrichDO(gene = SearchProteins, universe = Proteome, ont = "HPO")
    oraDisGeNET <- DOSE::enrichDGN(gene = SearchProteins, universe = Proteome)
    oraNCG <- DOSE::enrichNCG(gene = SearchProteins, universe = Proteome)
  }
  message("Plotting & Compiling Analyses")
  {
    MergedORA <- data.table::rbindlist(list(Proteopedia::Process_ORAOutput(oraBP, "Biological_Process"), Proteopedia::Process_ORAOutput(oraMF, "Molecular_Function"), 
                                            Proteopedia::Process_ORAOutput(oraCC, "Cellular_Component"), Proteopedia::Process_ORAOutput(oraKEGG, "KEGG"), 
                                            Proteopedia::Process_ORAOutput(oraMKEGG, "MKEGG"), Proteopedia::Process_ORAOutput(oraReactome, "Reactome"), 
                                            Proteopedia::Process_ORAOutput(oraDisease, "Disease"), Proteopedia::Process_ORAOutput(oraPhenotype, "Phenotype"), 
                                            Proteopedia::Process_ORAOutput(oraDisGeNET, "DisGeNET"), Proteopedia::Process_ORAOutput(oraNCG, "Network_of_Cancer_Genes")),
                                       fill = TRUE)
    MergedORA <- MergedORA[ID != 0]
    data.table::fwrite(MergedORA, file = "ORA_Results.csv")
  }
  message("Formatting Data For HTML Plots")
  {
    MergedORA <- MergedORA |> dplyr::arrange(factor(Ontology, levels = c("Cellular_Component", "Biological_Process", "Molecular Function", "KEGG", "KEGG Module", 
                                                                         "Reactome", "Disease","Phenotype","DisGeNET", "Network of Cancer Genes")))
    MergedORA[, Description := factor(Description, levels = unique(c(Proteopedia::Order_Terms(oraBP), Proteopedia::Order_Terms(oraMF), 
                                                                     Proteopedia::Order_Terms(oraCC), Proteopedia::Order_Terms(oraKEGG), 
                                                                     Proteopedia::Order_Terms(oraMKEGG), Proteopedia::Order_Terms(oraReactome), 
                                                                     Proteopedia::Order_Terms(oraDisease), Proteopedia::Order_Terms(oraPhenotype), 
                                                                     Proteopedia::Order_Terms(oraDisGeNET), Proteopedia::Order_Terms(oraNCG))))]
    MergedORA[, TermIndex := 1:nrow(MergedORA)]
    MergedORA[, NegLog10PAdj := -log10(p.adjust)]
    MergedORA[, SciPVal := formatC(p.adjust, format = "e", digits = 2)]
  }
  message("Exporting HTML Plots")
  {
    ORA_Significance <- highcharter::hchart(MergedORA, "scatter", highcharter::hcaes(x = TermIndex, y = -log10(p.adjust), group = Ontology)) |>
      highcharter::hc_chart(zoomType = "xy") |> 
      highcharter::hc_xAxis(title = list(text = NULL), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", tickWidth = 0 ) |>
      highcharter::hc_yAxis(title = list(text = "-Log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                            gridLineWidth = 0 ) |>
      highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Enrichment Factor: {point.RichFactor:.2f}
                            <br>P-Value: {point.SciPVal} <br> Gene count: {point.Count}")
    htmlwidgets::saveWidget(ORA_Significance, file = "Manhattan_Plot.html")
    
    ORA_EnrichSig <- highcharter::hchart(MergedORA, "scatter", highcharter::hcaes(x = RichFactor, y = -log10(p.adjust), group = Ontology)) |>
      highcharter::hc_chart(zoomType = "xy") |>
      highcharter::hc_xAxis(title = list(text = "Normalised Enrichment Score"), labels = list(enabled = FALSE), lineWidth = 0.5, lineColor = "black", 
                            tickWidth = 0 ) |> 
      highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "black", tickColor = "black",  
                            gridLineWidth = 0 ) |> 
      highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Enrichment Factor: {point.RichFactor:.2f}
                            <br>P-Value: {point.SciPVal} <br>Gene count: {point.Count}")
    htmlwidgets::saveWidget(ORA_EnrichSig, file = "Volcano_Plot.html")
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Plot_LFGenesetVolcano <- function(InputFile, SubsetColour = "red", Adj_PValueCutoff = 0.05){
  
  initiation.time <- Sys.time()
  setwd(gsub("(.*)/GSEA_Output/.*.csv.*", "\\1", InputFile))
  GenesetData <- data.table::fread(InputFile)[p.adjust < Adj_PValueCutoff]
  LimmaData <- data.table::fread(paste0(gsub("(.*)/GSEA_Output.*.csv.*", "\\1", InputFile), "/Limma_Output.csv")) 
  if(dir.exists(paste0(getwd(),"/GSEA_Volcanoes")) == TRUE){unlink(paste0(getwd(),"/GSEA_Volanoes"), recursive = TRUE)}
  dir.create(paste0(getwd(),"/GSEA_Volcanoes"), showWarnings = TRUE)
  setwd(paste0(getwd(),"/GSEA_Volcanoes"))
  for(GSEARow in 1:nrow(GenesetData)){
    SubsetProteins <- bitr(unlist(str_split(GenesetData$core_enrichment[GSEARow], "/")), fromType = "ENTREZID", toType = "UNIPROT", 
                           OrgDb = org.Hs.eg.db::org.Hs.eg.db)$UNIPROT
    
    SubsetData = LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, GeneGroup)]
    SubsetData[, Subset := fifelse(GeneGroup %in% SubsetProteins, T, F)]
    
    MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)
    
    Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
      ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
      Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() + 
      labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))
    
    adjpval <- GenesetData$p.adjust[GSEARow]
    
    Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
      ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() +
      ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
      ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
      ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T))*0.93, 
                        label = paste0("Adj. P-Value: ",ifelse(adjpval < 0.01, formatC(adjpval, format = "e", digits = 2), round(adjpval, digits = 2))), colour = SubsetColour, size = 4) +
      Proteopedia::Add_AbundanceAxes() + labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))
    pdf(paste0(gsub(":", "", GenesetData$ID[GSEARow]),"_Volcano.pdf"))
    print(Volcano)
    print(Volcano_Rug)
    print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
    Proteopedia::Reset_Dev()
  }
  completion.time <- Sys.time()
  duration.time <- completion.time - initiation.time
  duration.time
}
#' @export
Plot_SILACRatioGenesetVolcano <- function(InputFile, SubsetColour = "red", Adj_PValueCutoff = 0.05){
    initiation.time <- Sys.time()
    setwd(gsub("(.*)/GSEA_Output/.*.csv.*", "\\1", InputFile))
    GenesetData <- data.table::fread(InputFile)[p.adjust < Adj_PValueCutoff]
    LimmaData <- data.table::fread(paste0(gsub("(.*)/GSEA_Output.*.csv.*", "\\1", InputFile), "/Limma_Output.csv")) 
    if(dir.exists(paste0(getwd(),"/GSEA_Volcanoes")) == TRUE){unlink(paste0(getwd(),"/GSEA_Volanoes"), recursive = TRUE)}
    dir.create(paste0(getwd(),"/GSEA_Volcanoes"), showWarnings = TRUE)
    setwd(paste0(getwd(),"/GSEA_Volcanoes"))
    for(GSEARow in 1:nrow(GenesetData)){
      SubsetProteins <- bitr(unlist(str_split(GenesetData$core_enrichment[GSEARow], "/")), fromType = "ENTREZID", toType = "UNIPROT", 
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db)$UNIPROT
      
      SubsetData = LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, GeneGroup)]
      SubsetData[, Subset := fifelse(GeneGroup %in% SubsetProteins, T, F)]
      
      MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)
      
      Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + 
        ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) + 
        Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() + 
        labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))
      
      adjpval <- GenesetData$p.adjust[GSEARow]
      
      Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "black") + 
        ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) + Proteopedia::Add_NotSigBox() +
        ggplot2::geom_vline(xintercept = MeanLog2FC, colour = SubsetColour, linetype = "dashed", linewidth = 1) +
        ggplot2::geom_rug(alpha = ifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") + 
        ggplot2::annotate("label", x = MeanLog2FC, y = -log10(min(SubsetData[, P.Value], na.rm = T))*0.93, 
                          label = paste0("Adj. P-Value: ",ifelse(adjpval < 0.01, formatC(adjpval, format = "e", digits = 2), round(adjpval, digits = 2))), colour = SubsetColour, size = 4) +
        Proteopedia::Add_IsotopeRatioAxes() + labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))
      pdf(paste0(gsub(":", "", GenesetData$ID[GSEARow]),"_Volcano.pdf"))
      print(Volcano)
      print(Volcano_Rug)
      print(Volcano_Rug + ggplot2::labs(x = "", y = ""))
      Proteopedia::Reset_Dev()
    }
    completion.time <- Sys.time()
    duration.time <- completion.time - initiation.time
    duration.time
}

########### MS Analysis: Complete Functions ####################################################################################################################################
#' @export
Run_LF_CompleteAnalysis <- function(InputDirectory, ExpGroup, CtlGroup, DIANN = F, ProteotypicFiltering = F, MinPrecursors = 2, 
                                 ImputationQ = 0.01, ImputationSigma = 1, SubsetColour = "red", FunctionalAnalysis = F){
  start.Time <- Sys.time()
  if(DIANN == T){
    Proteopedia::Process_LF_DIANN(InputDirectory, CtlGroup, ProteotypicFiltering = ProteotypicFiltering)
    Proteopedia::Map_LF_PrecursorBiochemistry(InputDirectory)
  }
  Proteopedia::Analyse_LF_Proteins(InputDirectory, ExpGroup, CtlGroup, MinPrecursors = MinPrecursors, 
                                   ImputationQ = ImputationQ, ImputationSigma = ImputationSigma)
  Proteopedia::Map_LF_Proteins(InputFile = paste0(InputDirectory,"/",ExpGroup,"_vs_",CtlGroup,"_Output/Limma_Output.csv"), SubsetColour = SubsetColour)
  if(FunctionalAnalysis == T){
  Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory,"/",ExpGroup,"_vs_",CtlGroup,"_Output/Limma_Output.csv"))
  Proteopedia::Run_ORA(SearchProteins = data.table::fread(paste0(InputDirectory,"/",ExpGroup,"_vs_",CtlGroup,"_Output/Limma_Output.csv"))[Log2FC > 1 & P.Value < 0.05, ProteinGroup], 
                       Proteome = data.table::fread(paste0(InputDirectory,"/",ExpGroup,"_vs_",CtlGroup,"_Output/Limma_Output.csv"))[, ProteinGroup], 
                       OutputDirectory = gsub("(~/.*)/.*.csv", "\\1", InputFile))
  }
  Proteopedia::End_Timer(Start = start.Time)
}
#' @export
Run_StaticSILAC_CompleteAnalysis <- function(InputDirectory, ExpGroup, CtlGroup, DIANN = F, ProteotypicFiltering = F, MinPrecursors = 2, 
                                 ImputationQ = 0.01, ImputationSigma = 1, SubsetColour = "red", FunctionalAnalysis = F){
  start.time <- Sys.time()
  if(DIANN == T){
    Proteopedia::Process_StaticSILAC_DIANN(InputDirectory, CtlGroup, ProteotypicFiltering = ProteotypicFiltering)
  }
  Proteopedia::Analyse_StaticSILAC_Proteins(InputDirectory, ExpGroup, CtlGroup, MinPrecursors = MinPrecursors, 
                                             ImputationQ = ImputationQ, ImputationSigma = ImputationSigma)
  
  Proteopedia::Map_StaticSILAC_Proteins(InputDirectory, SubsetColour = SubsetColour)
  if(FunctionalAnalysis == T){
    Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/Total_Analysis/Limma_Output.csv"))
    Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/Ratio_Analysis/Limma_Output.csv"))
    Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/Channel_Analysis/Heavy_Limma_Output.csv"))
    Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/Channel_Analysis/Light_Limma_Output.csv"))
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Run_TimecourseSILAC_CompleteAnalysis <- function(InputDirectory, CtlGroup, ExpGroups, DIANN = F, ProteotypicFiltering = F, GenerateDataPlots = F, 
                                                 SameInitialAbundance = T, LightMinSamples = 0.5, HeavyModel = "NLS", HeavyMinMonotonicity = 0.5, 
                                                 HeavyMinSamples = 0.5, HeavyMaxCV = 0.3, UseLightKloss = F, ReplicatesUsed = c(1, 2, 3), SubsetColour = "red"){
  if(DIANN == T){
    Proteopedia::Process_TimecourseSILAC_DIANN(InputDirectory, CtlGroup, ProteotypicFiltering)
  }
  for(CtlGroup in ExpGroups){
    Proteopedia::Analyse_TimecourseSILAC_Proteins(InputDirectory, CtlGroup, ExpGroups, GenerateDataPlots, SameInitialAbundance, LightMinSamples, 
                                                  HeavyModel, HeavyMinMonotonicity, HeavyMinSamples, HeavyMaxCV, UseLightKloss, ReplicatesUsed)
    Proteopedia::SILAC_Timecourse_ProteinMapping(InputDirectory = paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"), 
                                                 SubsetColour = SubsetColour)
    Proteopedia::Map_TimecourseSILAC_Proteins(InputDirectory = paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
    if(FunctionalAnalysis == T){
      Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output/LightParameters.csv"))
      Proteopedia::Run_GSEA(InputFile = paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output/HeavyParameters.csv"))
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}

########### Coulter Analysis Functions ####################################################################################################################################
#' @export
Process_CoulterFiles <- function(InputDirectory, FilterBinVolumes = T){ 
  setwd(InputDirectory)
  CoulterFiles <- list.files(pattern = "*.#m4") |> purrr::map(readr::read_lines) |> purrr::set_names(basename(list.files(pattern = "*.#m4"))) 
  AnalysedData <- data.table::data.table()
  for(Run in names(CoulterFiles)){
    InputData <- unlist(CoulterFiles[Run]) |> data.table::data.table() 
    BinVol <- ((4/3)*pi)*((unlist(InputData[seq(which(InputData[,1] == "[#Bindiam]")+1,which(InputData[,1] == "[Binunits]")-1),1]) |> as.numeric()/2)^3)
    CellCount <- unlist(InputData[seq(which(InputData[,1] == "[#Binheight]")+1,which(InputData[,1] == "[SizeStats]")-1),1]) |> as.numeric()
    TotalCount <- gsub("SampleSize= ", "", unlist(InputData[which(grepl("SampleSize=", unlist(InputData)))])) |> as.numeric()
    SampleData <- data.table::data.table(BinVol, CellCount)
    SampleData[, CellCount := data.table::fifelse(FilterBinVolumes == T & (BinVol < 800 | BinVol > 20000), 0, CellCount)]
    SampleData[, TotalCount := sum(SampleData$CellCount)]
    SampleData[, Run := Run]
    SampleData[, CellProp := CellCount/TotalCount]
    AnalysedData <- AnalysedData |> rbind(SampleData)
  }
  data.table::fwrite(AnalysedData, "AnalysedCoulterData.csv")
}
#' @export
Add_CoulterLowerLimit <- function(){
  ggplot2::annotate("rect", xmin = -Inf, xmax = log10(800), ymin = -Inf, ymax = Inf, alpha = 0.1) 
}
#' @export
Add_CoulterUpperLimit <- function(){
  ggplot2::annotate("rect", xmax = Inf, xmin = log10(20000), ymin = -Inf, ymax = Inf, alpha = 0.1)
}
#' @export
Summarise_CoulterData <- function(InputDirectory, CtlGroup = NA){
  setwd(InputDirectory)
  if(length(list.files(pattern = "Meta")) == 1){
    CoulterData <- data.table::fread("AnalysedCoulterData.csv") |> data.table::merge.data.table(data.table::fread("Coulter_Meta.csv"))
    if(is.na(ControlCondition)){CtlGroup <- CoulterData$Condition[1]}
    CoulterData[, Condition:= factor(Condition, levels = c(CtlGroup, CoulterData[Condition != CtlGroup, Condition] |> unique()))]
    BinSummary <- CoulterData[, .(N = .N, MeanProp = mean(CellProp), SDProp = sd(CellProp)), by = .(Condition, BinVol)]
    BinSummary[, SEProp := SDProp/sqrt(N)]
    
    ConditionSummary <- BinSummary[, .(MeanProp = max(MeanProp)), Condition] |> data.table::merge.data.table(BinSummary[, .(Condition, BinVol, MeanProp)]) |> data.table::copy() |>
      data.table::setnames(c("MeanProp", "BinVol"), c("MaxProp", "MaxBinVol"))
    ConditionSummary$VolLabel <- paste0(signif(ConditionSummary$MaxBinVol, 4), "~fL~(µm^{3})")
  
    ConditionData <- CoulterData[, .(CellCount = mean(CellCount)), .(Condition, BinVol)]
    ExpandedConditionData <- splitstackshape::expandRows(ConditionData[CellCount > 0, .(Condition, BinVol, CellCount)], "CellCount")
    ExpandedConditionSummary <- ExpandedConditionData[, .(N = .N, Mean = mean(BinVol), SD = sd(BinVol)), .(Condition)]
    ExpandedConditionSummary[, SE := SD/sqrt(N)]
    data.table::fwrite(ExpandedConditionSummary, "ConditionSummary.csv")
    
    CoulterStats <- data.table::data.table(Condition = levels(ExpandedConditionData$Condition), P.Value = 0)
    for(ExpGroup in CoulterStats$Condition){
      CoulterStats[Condition == ExpGroup]$P.Value <- wilcox.test(ExpandedConditionData[Condition == ExpGroup, BinVol], ExpandedConditionData[Condition == CtlGroup, BinVol])$p.value
    }
    CoulterStats[, SigSymbol := data.table::fifelse(P.Value < 0.001, "***", data.table::fifelse(P.Value < 0.01, "**", data.table::fifelse(P.Value < 0.05, "*","")))]
    data.table::fwrite(CoulterStats, "Wilcoxon_Output.csv")
    
    pdf("CoulterBoxplot.pdf", width = 12, height = 8)
      print(ExpandedConditionData |> ggplot2::ggplot(ggplot2::aes(x = Condition, y = BinVol, colour = Condition)) + ggplot2::geom_boxplot(outliers = F) +
        ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::geom_text(data = CoulterStats, ggplot2::aes(y = 450, label = SigSymbol, vjust = 0.8), size = 10) + 
        ggplot2::coord_flip() + ggplot2::labs(y = Bin~Volume~(fL)~(µm^3)) + ggplot2::theme(axis.title.y = ggplot2::element_blank()))
    Proteopedia::Reset_Dev() 

    pdf("CoulterAnalysis.pdf", width = 20, height = 18)
      print(BinSummary |> ggplot2::ggplot(ggplot2::aes(x = log10(BinVol), y = MeanProp, colour = Condition)) + ggplot2::geom_line() + 
              ggplot2::labs(x = expression("Log"[10]~"Bin Volume (fL) (µm"^"3"*")"), y = "Proportion of Cells") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
              Add_CoulterLowerLimit() + Add_CoulterUpperLimit() + ggplot2::geom_errorbar(ggplot2::aes(ymax = MeanProp+SEProp, ymin = MeanProp-SEProp), width = 0.005) + 
              #ggplot2::geom_vline(data = CoulterMeans, ggplot2::aes(xintercept = log10(BinVol), colour = Condition), size = 2, linetype = "dotted") + 
              ggplot2::coord_cartesian(ylim = c(0, max(ConditionSummary$MaxProp)*1.2), clip = "off") + ggplot2::scale_x_continuous(limits = c(2.9, 4.35)) + 
              ggrepel::geom_label_repel(data = ConditionSummary, ggplot2::aes(x = log10(ConditionSummary$MaxBinVol), y = max(ConditionSummary$MaxProp)*1.1, label = VolLabel), 
                                  size = 8, parse = T, nudge_x = 0, min.segment.length = 10, force_pull = 0, max.overlaps = Inf) + ggplot2::scale_y_continuous(expand = c(0,0)) + 
              patchwork::inset_element(CoulterData[, .(Condition, Replicate, TotalCount)] |> unique() |> 
                                         ggplot2::ggplot(ggplot2::aes(x = Condition, y = TotalCount/1000, fill = Condition)) + 
                                         ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2(width = 0.9)) + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") + 
                                         ggplot2::labs(y = "Cells Analysed (x1000)") + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::theme(axis.title.x = ggplot2::element_blank()), 0.60, 0.5, 0.9, 0.99))
      Proteopedia::Reset_Dev() 
  } else {
    message("ERROR: No/Multiple Meta Files Detected")
  }
}
