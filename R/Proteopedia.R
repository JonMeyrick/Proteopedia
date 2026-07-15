########### Import Libraries ####################################################################################################################################
#' @import arrow
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
#' @import iq
#' @import Proteopedia
#' @import limma
#' @import mzR
#' @import org.Hs.eg.db
#' @import patchwork
#' @import Peptides
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
#' Update Proteopedia (Developer Only)
#'
#' For developer use, updates Proteopedia from local directory
#' @export
Update_Proteopedia <- function(LocalDirectory){
  library(devtools)
  library(roxygen2)
  setwd(LocalDirectory)
  devtools::document()
  devtools::build()
  devtools::install()
}
#' %!in%
#'
#' Extract variables not in string
#' @param x Search variable
#' @param y String variables
#' @return x variables not matching any y variables
#' @examples
#' c("a", "b", "c") %!in% c("a", "b");
#' "c"
#' @export
"%!in%" <- function(x,y)!("%in%"(x,y))
#' @export
"%like%" <- function(x,y){like(x,y)}
#' End Timer
#'
#' Report time difference since defined start
#' @param Start Start of timer called by assigning Sys.time to a variable name
#' @return Message of time difference between Start and calling of funtion
#' @export
End_Timer <- function(Start){
  time.taken <- format(difftime(Sys.time(), Start, unit= "secs"), scientific = F)
  Minutes <- floor(as.numeric(gsub("(.*)\\..*", "\\1", time.taken))/60)
  Seconds <- as.numeric(gsub("(.*)\\..*", "\\1", time.taken)) - floor(as.numeric(gsub("(.*)\\..*", "\\1", time.taken))/60)*60
  Milliseconds <- gsub("(..).*", "\\1", as.numeric(gsub(" secs", "", gsub(".*\\.(.*)", "\\1", time.taken))))
  message(paste0("Time Elapsed: ", Minutes, " min ", Seconds, " s ", Milliseconds, " ms"))
}
#' Check Aggregate Rows
#'
#' Identify aggregated data.table rows using ID columns
#' @param InputData data.table to be grouped
#' @param IDColumns Grouping variables of data.table
#' @export
Check_AggregateRows <- function(InputData, IDColumns){
  InputData[, .N, IDColumns][N > 1]
}
#' Reset Developer
#'
#' Repeat the dev.off() function until the and error is generated. Prevents remaining in developer
#' @export
Reset_Dev <- function(){
  done <- F
  while(!done){
    result <- try(dev.off(), silent = T)
    done <- class(result) == "try-error"
  }
}

########### ZenoTOF Setup Functions ####################################################################################################################################
#' Generate Metadata
#'
#' Generates metadata table using condition factor, with respect to concentration, time, and temperature variables.
#' Concentration, time, and temperature should be defined as Molar (M), hours (h), and degrees celsius (C), respectively
#' @export
Generate_Metadata <- function(Condition = 0, Conc = 0, Time = 0, Reps = 3, OrderedBy = c("Condition", "Conc", "Time", "Temp", "Rep"), Output_Directory = "~/Desktop"){

  Sample_Metadata <- tidyr::crossing(factor(Condition), Conc, Temp, Time, 1:Reps) |> data.table::data.table() |>  data.table::setnames("1:Reps", "Rep") |>
    data.table::setorderv(OrderedBy) |> tibble::rowid_to_column("ID")
  suppressWarnings(
    Sample_Metadata <- Sample_Metadata[, Sample := paste0(Condition,"_",Conc,"M_",Time,"h_", Temp, "C_","R",Rep)]
  )
  data.table::fwrite(Sample_Metadata, paste0(Output_Directory,"/Sample_Metadata.csv"))
}
#' @export
Generate_ZenoTOFBatch <- function(Input_Metadata, Run_Date = paste0("E",substr(gsub("-","",lubridate::today()), 3, 8)), GenericNames = F,
                                  Sample1Position, SampleNPosition = 0, LC_Run = "30method_Berlin_Nov_2024_CCP_depth_0.5mm",
                                  Rack_Position = 1, InjectionVolume = 2, Blank_Position = "H12", Blank_Every = 6, InjectorID = "JM"){
  Sample_Metadata <- data.table::fread(Input_Metadata)
  N_Samples <- nrow(Sample_Metadata)
  BothRacks <- N_Samples > 96
  # Define Sample Positions in Well
  PlateMap <- tidyr::crossing(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12, 1:2) |> data.table::data.table() |>
    data.table::setnames(c("RowID", "ColID", "RackID")) |> data.table::setorder(RowID, ColID, RackID)
  suppressWarnings(
    PlateMap[, Well := paste0(RowID, ColID)]
  )
  if(BothRacks){
    PlateMap <- PlateMap |> dplyr::slice_head(by = "RackID", n = N_Samples/2)
  } else {
    PlateMap <- PlateMap[RackID == Rack_Position] |> dplyr::slice_head(n = N_Samples)
  }
  PlateMap <- PlateMap |> tibble::rowid_to_column("SampleN")
  if(SampleNPosition == 0){
    Sample_Positions <- PlateMap[SampleN <= N_Samples, Well]
    Rack_Positions <- PlateMap[SampleN <= N_Samples, RackID]
  } else {
    Sample_Positions <- PlateMap[ColID <= as.numeric(gsub(".(\\d+)", "\\1", SampleNPosition)), Well][1:N_Samples]
    Rack_Positions <- PlateMap[, RackID]
  }

  N_Injections <- 1 + N_Samples + (N_Samples/Blank_Every)
  Blank_Positions <- seq(from = 1, to = N_Injections, by = Blank_Every+1)
  if(BothRacks){
    Blank_Rack <- rep(c(1, 2), each = length(Blank_Positions)/2)
  } else {
    Blank_Rack <- rep(Rack_Position, times = length(Blank_Positions))
  }

  Batch_File <- data.table::data.table("Sample_Name" = rep(0, times = N_Injections), "Sample_ID" = seq(1, to = N_Injections), "Barcode_ID" = "",
                                       "MS_Method" = 0, "Processing_Method" = "", "LC_Method" = 0, "Rack_Type" = "Sample Manager", "Rack_Position" = 0,
                                       "Plate_Type" = "Custom-96-Position", "Vial_Position" = 0, "Sample_Type" = "Unknown", "Dilution_Factor" = 1, "Weight_Volume" = 0,
                                       "Data_File" = 0, "Results_File" = "", "Comment" = "", "Injection_Volume" = 0,
                                       "Marker_Well" = F)

  Batch_File[, `:=`(MS_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_12min", "uFlow_ZenoSWATH_85VW_11ms_30T_ZENO_ON_Berlin"),
                    LC_Method = data.table::fifelse(Sample_ID %in% Blank_Positions, "Blank_5grad_12method_depth_0.5mm", LC_Run),
                    Injection_Volume = data.table::fifelse(Sample_ID %in% Blank_Positions, 5, InjectionVolume))]
  # Add Blank ID
  LoopCount <- 1
  for(BlankID in Blank_Positions){
    Batch_File$Sample_Name[BlankID] <- paste0(Run_Date,"_",Batch_File$Sample_ID[BlankID],"_Zeno_",InjectorID,"_QC_12_Blank",LoopCount)
    Batch_File$Vial_Position[BlankID] <- Blank_Position
    Batch_File$Rack_Position[BlankID] <- Blank_Rack[LoopCount]
    LoopCount <- LoopCount + 1
  }

  # Add Sample ID  & Vial Positions
  Sample_IDs <- 1:N_Injections
  Sample_IDs <- Sample_IDs[Sample_IDs %!in% Blank_Positions]

  LoopCount = 1
  if(GenericNames){
    for(SampleID in Sample_IDs){
      Batch_File$Sample_Name[SampleID] <- paste0(Run_Date,"_",Batch_File$Sample_ID[SampleID],"_Zeno_",InjectorID,"_IN_30_S",LoopCount)
      Batch_File$Vial_Position[SampleID] <- Sample_Positions[LoopCount]
      Batch_File$Rack_Position[SampleID] <- Rack_Positions[LoopCount]
      LoopCount = LoopCount + 1
    }
  } else {
    Sample_Metadata$Run <- 0
    for(SampleID in Sample_IDs){
      Batch_File$Sample_Name[SampleID] <- paste0(Run_Date,"_",Batch_File$Sample_ID[SampleID],"_Zeno_",InjectorID,"_IN_30_",Sample_Metadata$Sample[LoopCount])
      Sample_Metadata$Run[LoopCount] <- paste0(Run_Date,"_",Batch_File$Sample_ID[SampleID],"_Zeno_",InjectorID,"_IN_30_",Sample_Metadata$Sample[LoopCount])
      Batch_File$Vial_Position[SampleID] <- Sample_Positions[LoopCount]
      Batch_File$Rack_Position[SampleID] <- Rack_Positions[LoopCount]
      LoopCount = LoopCount + 1
    }
  }

  if(!GenericNames){
    setwd(gsub("/Sample_Metadata.csv","", paste0(Input_Metadata)))
    data.table::fwrite(Sample_Metadata, paste0(Input_Metadata))
  }

  Batch_File <- Batch_File[, Data_File := paste0("\\", substr(Run_Date, 1,3), "\\", substr(Run_Date, 1,5), "\\", Run_Date, "\\", Sample_Name)]
  suppressWarnings(
    Batch_File[, Sample_ID := NULL]
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
                                       "Data_File" = 0, "Results_File" = "", "Comment" = "", "Injection_Volume" = 0, "Marker_Well" = F)

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

########### MS Analysis: Data QC & Wrangling Functions ####################################################################################################################################
#' Calculate Area Under Curve (AUC)
#'
#' Calculates area under curve (AUC) of .mzML of .raw files from Sciex or ThermoFisher mass spectrometers
#' @export
Calculate_AUC <- function(InputDirectory, FileType = "mzML"){
  setwd(InputDirectory)
  AUCs = data.table::data.table()
  for(FileIndex in list.files(pattern = paste0(".", FileType, "$"))){
  TIC_data <- mzR::tic(mzR::openMSfile(paste0(InputDirectory, "/",FileIndex)))
  AUCs = rbind(AUCs, data.table::data.table(File = FileIndex, AUC = DescTools::AUC(TIC_data$rtime, TIC_data$intensity)))}
  AUCs[, Run := gsub(paste0(".", FileType, ".*"), "", File)]
  data.table::fwrite(AUCs, "AUC_Data.csv")
}
#' Simplify Geneset Data
#'
#' @export
Simplify_Data <- function(x){data.table::data.table(clusterProfiler::simplify(x, cutoff = 0.7))}
#' Calculate Label-Free Quantitation (LFQ)
#'
#' Calculates the label-free quantitation (LFQ) of proteins using the fastMaxLFQ algorithm from the iq package
#' @export
Calculate_LFQ <- function(InputData, LFQ_Col){
  tmp <- iq::fast_MaxLFQ(InputData[, .(protein_list = ProteinGroup, sample_list = Sample, id = Precursor.Id, quant = log2(Precursor.Normalised))])
  tmp <- data.table::data.table(tmp$estimate, Precursor_group = tmp$annotation, keep.rownames = "ProteinGroup")
  tmp <- data.table::melt.data.table(tmp, id.vars = c("ProteinGroup", "Precursor_group"), variable.name = "Sample", value.name = LFQ_Col)
  tmp <- tmp[ Precursor_group == "" & !is.na( get(LFQ_Col) ) ][, -"Precursor_group" ]
  tmp <- tmp[, (LFQ_Col) := 2^get(LFQ_Col)]
}
#' @export
Merge_PrecursorData <- function(x, y){merge(x, y, by = c("ProteinGroup", "Sample"), all = T)}
#' @export
Count_Proteins <- function(dt, var_name){
  tmp <- dt[, .(N_Samples = .N), ProteinGroup][, .(N_Proteins = .N), N_Samples]
  tmp[order(-N_Samples), CumulativeProteins_N := cumsum(N_Proteins)]
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
Calculate_HalfLife <- function(Kloss){log(2)/-Kloss}
#' @export
Calculate_DensityPeak <- function(x){max(density(x, na.rm = T)$y)}
#' @export
Calculate_VolcanoLog2FC <- function(x){log2(abs(x))*sign(x)}
#' @export
Add_ProteinInfo <- function(InputData, ProteinInfoFile){
  ProteinInfo <- data.table::setnames(data.table::fread(ProteinInfoFile), "Protein.Id", "ProteinGroup")[, .(ProteinGroup, Description, Gene)]
  ProteinInfo[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
  InputData[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
  InputData <- InputData |> data.table::merge.data.table(unique(ProteinInfo[, .(GeneGroup, Gene)]), by = "GeneGroup", all.x = T)
  InputData[grepl(";", GeneGroup) & is.na(Gene), Gene := sapply(strsplit(GeneGroup, ";"), function(x){paste(unique(ProteinInfo[, .(GeneGroup, Gene)])[, Gene][match(x, unique(ProteinInfo[, GeneGroup]))], collapse = ";")})]
  InputData[, URL := paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup))]
  return(InputData)
}
#' @export
Separate_Isoforms <- function(InputData){
  LimmaData <- InputData |> data.table::copy()
  if(length(grepl("Isoforms", colnames(LimmaData)) |> unique()) == 1){
    LimmaData$Isoforms <- 1
      for(i in 1:nrow(LimmaData)){
        if(length(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T)) > 0){
          LimmaData$Isoforms[i] <- paste0(stringr::str_extract_all(LimmaData$ProteinGroup[i], "-\\d", simplify = T), collapse = ", ")
        } else {LimmaData$Isoforms[i] <- 1}
      }
    LimmaData[, ProteinGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
    LimmaData[, GeneGroup := sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+","", x)), collapse = ";")})]
    LimmaData[, URL := paste0("https://www.uniprot.org/uniprotkb/", gsub(";.*", "", GeneGroup))]
  } else {
    message("Isoforms Already Separated")
  }
  return(LimmaData)
}
#' @export
Perform_Kmeans <- function(InputData, VariableIdentifier, MaxClusters = 25){
  if(length(colnames(InputData)[grepl(VariableIdentifier, colnames(InputData))]) < 2){
    message("ERROR: No. Variables columns <2\nCheck Wide Formatting")
  } else {
    set.seed(123)
    KmeansTest <- InputData |> data.table::copy()
    KmeansTest |> data.table::setnames(colnames(KmeansTest)[grepl(VariableIdentifier, colnames(KmeansTest))], c("Var1", "Var2"))
    WSS_Data <- data.table::data.table(Clusters = 1:MaxClusters, WSS = 0)
    for(ClusterIndex in 1:MaxClusters) {
      WSS_Data$WSS[ClusterIndex] <- stats::kmeans(KmeansTest[, .(Var1, Var2)], centers = ClusterIndex, nstart = 20)$tot.withinss
    }
    WSS_Start <- WSS_Data[Clusters == 1, WSS]
    WSS_End <- WSS_Data[Clusters == MaxClusters, WSS]
    WSS_Data <- WSS_Data[, Distance := abs((Clusters*(WSS_End - WSS_Start)) - (WSS*(MaxClusters - 1)) + (MaxClusters*WSS_Start) + (WSS_End*1))/sqrt((WSS_End - WSS_Start)^2 + (MaxClusters - 1)^2)]
    ElbowPoint <- WSS_Data[Distance == max(WSS_Data$Distance), Clusters]

    KmeansActual <- InputData |> data.table::copy()
    KmeansActual |> data.table::setnames(colnames(KmeansActual)[grepl(VariableIdentifier, colnames(KmeansActual))], c("Var1", "Var2"))
    InputData$Cluster <- stats::kmeans(KmeansActual[, .(Var1, Var2)], centers = ElbowPoint, nstart = 20)$cluster
    return(InputData)
  }
}
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
Add_Pearsons <- function(Subgroups = T){if(!Subgroups){ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(r.label),  group = 1), geom = "text")
} else {ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(r.label)), geom = "text")}}
#' @export
Add_R2 <- function(Subgroups = T){if(!Subgroups){ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(rr.label),  group = 1), geom = "text")
} else {ggpubr::stat_cor(ggplot2::aes(label = ggplot2::after_stat(rr.label)), geom = "text")}}
#' @export
Add_NotSigBox <- function(cutoff = 0.05){ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(cutoff), alpha = 0.1)}
#' @export
Add_XYLine <- function(colour = "#000"){ggplot2::geom_abline(linetype = "dashed", colour = colour)}
#' @export
Clean_SideDensities <- function(x = T, y = T){ggplot2::theme(ggside.axis.line = ggplot2::element_blank(), ggside.axis.ticks = ggplot2::element_blank(),
                                                             ggside.axis.text = ggplot2::element_blank())}
#' @export
Add_GSEAAxes <- function(){ggplot2::labs(x = "Normalised Enrichment Score (NES)", y = expression("-Log"[10]~"Adj. P-Value"))}
#' @export
Add_Isotope_Colour <- function(){ggplot2::scale_colour_manual(values = c("Heavy" = "#90F", "Light" = "#F09"))}
#' @export
Add_Isotope_Fill <- function(){ggplot2::scale_fill_manual(values = c("Heavy" = "#90F", "Light" = "#F09"))}
########### MS Analsis: Statistics ##########################################################################################
#' @export
Calculate_WilcoxonByVar <- function(InputData, Category, Measure){
  WilcoxSummary <- data.table::data.table("Category" = unique(InputData[, get(Category)]))
  WilcoxSummary$Wilcoxon_p = 0
  RowIndex = 1
  for(VarIndex in WilcoxSummary$Category){
    if(length(na.omit(InputData[get(Category) == VarIndex, get(Measure)])) > 0){
      WilcoxSummary$Wilcoxon_p[RowIndex] <- wilcox.test(InputData[, get(Measure)], InputData[get(Category) == VarIndex, get(Measure)], na.action = remove)$p.value
      } else {WilcoxSummary$Wilcoxon_p[RowIndex] <- NA}
    RowIndex = RowIndex + 1
    }
  WilcoxSummary[, SigSymbol := data.table::fifelse(Wilcoxon_p < 0.001, "***", data.table::fifelse(Wilcoxon_p < 0.01, "**", data.table::fifelse(Wilcoxon_p < 0.05, "*","")))]
  return(WilcoxSummary |> data.table::setnames("Category", paste(Category)))
}
########### MS Analysis: Static Analysis Functions #####################################################################################################################################
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
Process_LabelFree_DIANN <- function(InputDirectory, ProteotypicFiltering = F, DIANNVersion = 2.2){
  start.time <- Sys.time()
  set.seed(123)
  message("Importing DIA-NN Report File")
  {
    setwd(InputDirectory)
    if(DIANNVersion >= 2.2){
      if(length(list.files(pattern = "report.parquet")) > 0){
        PrecursorData <- arrow::read_parquet(list.files(pattern = "report.parquet")[1],
                                             col_select = c("Run", "Channel", "Protein.Group", "Stripped.Sequence", "Precursor.Id", "Proteotypic",
                                                            "Precursor.Normalised", "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                                                            "Global.Q.Value", "Global.PG.Q.Value", "Quantity.Quality", "Channel.Q.Value")) |>
          data.table::setDT() |> data.table::setnames("Protein.Group", "ProteinGroup")
      } else {return(message("ERROR: No Input File Found"))}
      if(nrow(PrecursorData[Channel == "H"]) > 0){return(message("ERROR: Label-Free Processing on SILAC Data"))}
    } else {
      if(length(list.files(pattern = "report.tsv")) > 0){
        PrecursorData <- data.table::fread(list.files(pattern = "report.tsv")[1])[, .(Run, Protein.Group, Protein.Ids, Stripped.Sequence,
                                                                                      Precursor.Id, Proteotypic, Precursor.Normalised,
                                                                                      Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value,
                                                                                      Lib.Q.Value, Lib.PG.Q.Value)] |> data.table::setnames("Protein.Group", "ProteinGroup")

      data.table::fwrite(unique(data.table::fread(list.files(pattern = "report.tsv")[1])[, .(Protein.Group, First.Protein.Description, Genes)]) |>
                           data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), c("Protein.Id", "Gene", "Description")),
                         "report.protein_description.tsv.gz")
      } else {return(message("ERROR: No Input File Found"))}
      if(any(grepl("SILAC-", PrecursorData$Precursor.Id))){return(message("ERROR: Label-Free Processing on SILAC Data"))}
    }
    if(file.exists("Sample_Rename.csv")){
      SampleRenaming <- data.table::fread("Sample_Rename.csv")
      PrecursorData$Sample <- SampleRenaming$Renamed[match(unlist(PrecursorData$Run), SampleRenaming$Run)]
      PrecursorData <- PrecursorData[!is.na(Sample)]
    }
  }
  message("Defining Metadata")
  {
    #PrecursorData[, Cell := gsub("(.*)_.*$","\\1", Sample)]
    #PrecursorData[, Drug := gsub("(.*)_(.*)_(.*)_(.*)_R(\\d)","\\3", Sample)]
    PrecursorData[, Conc := gsub(".*_(.*M)_.*","\\1", Sample)]
    PrecursorData[, Temp := gsub(".*_(\\d+)C_.*","\\1", Sample)]
    PrecursorData[, Time := gsub(".*_(.*h)_.*$","\\1", Sample)]
    PrecursorData[, Replicate := gsub(".*_(.*)$","\\1", Sample)]
    PrecursorData[, Sample := gsub("_0", "", Sample)]
    PrecursorData[, Condition := gsub("(.*)_.*$", "\\1", Sample)]
    Metadata <- unique(PrecursorData[, .(Sample, Conc, Time, Temp, Replicate, Condition)]) |> dplyr::arrange(Condition, Replicate)
    PrecursorData$Sample <- factor(PrecursorData$Sample, levels = Metadata$Sample)
    data.table::fwrite(Metadata, "Sample_Metadata.csv")
  }
  message("Filtering Precursors")
  {
    if(ProteotypicFiltering){PrecursorData <- PrecursorData[Proteotypic >= 1]}
    PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01]
    PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
  }
  message("Compiling & Exporting Data")
  {
    data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")
    ProteinData <- Reduce(Merge_PrecursorData, list(Calculate_LFQ(PrecursorData, "LFQ"), PrecursorData[, .(Intensity = sum(Precursor.Normalised)), .(ProteinGroup, Sample)],
                                                    PrecursorData[, .(N_precursors = data.table::uniqueN(Precursor.Id)), .(ProteinGroup, Sample)],
                                                    PrecursorData[, .(N_precursors_proteotypic = sum(Proteotypic)), .(ProteinGroup, Sample)]))
    ProteinData <- Proteopedia::Add_ProteinInfo(ProteinData, paste0(InputDirectory, "/report.protein_description.tsv.gz")) |> data.table::merge.data.table(Metadata, by = "Sample")
    data.table::fwrite(ProteinData, "LF_DIANN_Output.csv.gz")
  }
  message("Plotting Intensities")
  {
    IntensitiesData <- data.table::rbindlist(list(
      PrecursorData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Precursor.Normalised), Type = "Precursor Quantity")],
      ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(LFQ), Type = "Protein MaxLFQ")],
      ProteinData[, .(Sample, Condition, Replicate, `log2 quantity` = log2(Intensity), Type = "Protein Intensity")]
    ))
    IntensitiesData[, Type := factor(Type, levels = c("Precursor Quantity", "Protein MaxLFQ", "Protein Intensity"))]

    IntensitiesPlot <- IntensitiesData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = `log2 quantity`, colour = Condition)) +
      ggplot2::geom_boxplot(outliers = F) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::facet_wrap("Type", scales = "free_x") + ggplot2::ylab("Log2 Value") + ggplot2::coord_flip() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }
  message("Plotting Precursor, Peptide & Protein Counts")
  {
    CountsPlot <- PrecursorData[, lapply(.SD, data.table::uniqueN), .(Sample, Condition), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")] |>
      data.table::melt.data.table(id.vars = c("Sample","Condition"), value.name = "IDs") |>
      ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = IDs/1000, fill = Condition, label = format(IDs, big.mark = ",", scientific = F))) +
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.2) +
      ggplot2::scale_y_continuous(expand = 0) + ggplot2::facet_wrap("variable", scales = "free_x", labeller = ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups"))) +
      ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  }
  message("Plotting Data Completeness")
  {
    CompletenessData <- rbind(Count_Proteins(ProteinData, "All"), Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                              Count_Proteins(ProteinData[N_precursors_proteotypic >= 2], "≥ 2 Proteotypic"))
    CompletenessDataLabels <- CompletenessData[Precursors == "All" & N_Samples == 1][, Tag := "All Precursors"] |> rbind(CompletenessData[Precursors == "≥ 2" & N_Samples == round(nrow(Metadata)/2)][, Tag := "≥ 2 Precursors"]) |>
                                rbind(CompletenessData[Precursors == "≥ 2 Proteotypic" & N_Samples == nrow(Metadata)][, Tag := "≥ 2 Proteotypic\nPrecursors"])

     NAsPlot <- CompletenessData |> ggplot2::ggplot(ggplot2::aes(x = N_Samples, y = CumulativeProteins_N/1000, colour = Precursors))+
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") +
      ggplot2::scale_colour_manual(values = c("All" = "#000", "≥ 2" = "#999", "≥ 2 Proteotypic" = "#F63"), guide = "none") +
      ggplot2::scale_x_continuous(breaks = seq(1, 1000, 1)) + ggplot2::scale_y_continuous(limits = c(0, ceiling(max(CompletenessData$CumulativeProteins_N)/1000))) +
      ggrepel::geom_text_repel(data = CompletenessDataLabels, ggplot2::aes(label = Tag), nudge_y = -2) +
      ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
  }
  message("Calculating Missed Trypsinisation Sites")
  {
    TrypsinData <- PrecursorData |> data.table::copy()
    TrypsinData[, MissedTrypsin := grepl("[RK][^P]", Stripped.Sequence)]
    TrypsinData[, N_Trypsin := .N, .(Sample, MissedTrypsin)]
    TrypsinData[, N_Sample := .N, Sample]
    TrypsinData <- TrypsinData[MissedTrypsin == T, .(Sample, Condition, Replicate, N_Trypsin, N_Sample)] |> dplyr::distinct()
    suppressWarnings(TrypsinData[, PercentTrypsin := (N_Trypsin/N_Sample)*100])

    TrypsinisationPlot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentTrypsin, fill = Condition)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "", y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::coord_flip()
  }
  message("Plotting Precursor & Protein Variation")
  {
    PrecursorCVs <- PrecursorData[, .(CV = Calculate_CV(Precursor.Normalised), N = .N), .(Precursor.Id, Condition)]
    PrecursorCVs <- PrecursorCVs[, `:=`(Rank = data.table::frank(CV), ID = "Precursors"), Condition]
    ProteinCVs <- ProteinData[, .(CV = Calculate_CV(LFQ), N = .N), .(ProteinGroup, Condition)]
    ProteinCVs <- ProteinCVs[, `:=`(Rank = data.table::frank(CV), ID = "Protein Groups"), Condition]

    suppressWarnings(VariationPlot <- PrecursorCVs[, Precursor.Id := NULL] |> rbind(ProteinCVs[, ProteinGroup := NULL]) |> ggplot2::ggplot(ggplot2::aes(x = Rank/1000, y = CV, colour = Condition)) +
      ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Variation (%)") +
      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::coord_cartesian(ylim = c(0,50)) +
      ggplot2::facet_wrap(~ID, scales = "free"))
  }
  message("Exporting QC Plots")
  {
    ConditionIDPlot <- data.table::data.table(Condition = unique(Metadata$Condition)) |> tibble::rowid_to_column("ID") |>
      ggplot2::ggplot(ggplot2::aes(x = Condition, y = 1, colour = Condition)) + ggplot2::geom_point(size = 20) +
      ggplot2::geom_text(ggplot2::aes(label = Condition, y = 1), vjust = 3) +
      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())

    pdf("DIANN_QCPlot.pdf", width = 18, height = 20)
    suppressWarnings(
      print(IntensitiesPlot + CountsPlot +
            patchwork::free(NAsPlot, type = "label") + patchwork::free(TrypsinisationPlot, type = "label") +
            patchwork::free(VariationPlot, type = "label") +
            patchwork::plot_layout(design = "AAAA\nBBBB\nCCDD\nEEEE") + patchwork::plot_annotation(tag_levels = "A"))
    )
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Process_SILAC_DIANN <- function(InputDirectory, ProteotypicFiltering = F, TotalNormalisation = F, DIANNVersion = 2.2){
  set.seed(123)
  start.time <- Sys.time()
  message("Importing DIA-NN Report File")
  {
    setwd(InputDirectory)
    if(DIANNVersion >= 2.2){
      if(length(list.files(pattern = "report.parquet")) > 0){
        PrecursorData <- arrow::read_parquet(list.files(pattern = "report.parquet")[1],
                                             col_select = c("Run", "Channel", "Protein.Group", "Stripped.Sequence", "Precursor.Id", "Proteotypic",
                                                            "Precursor.Normalised", "Precursor.Quantity", "Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                                                            "Global.Q.Value", "Global.PG.Q.Value", "Quantity.Quality", "Channel.Q.Value")) |>
          data.table::setDT() |> data.table::setnames("Protein.Group", "ProteinGroup")
      } else {return(message("ERROR: No Input File Found"))}
    } else {
      if(length(list.files(pattern = "report.tsv")[1]) > 0){
        InputFile <- list.files(pattern = "report.tsv")[1]
        PrecursorData <- data.table::fread(InputFile)[, .(Run, Protein.Group, Protein.Ids, Stripped.Sequence, Precursor.Id,
                                                          Proteotypic, Precursor.Quantity, Precursor.Translated, Channel.Q.Value,
                                                          Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value, Lib.Q.Value,
                                                          Lib.PG.Q.Value)] |> data.table::setnames("Protein.Group", "ProteinGroup")
        PrecursorData[data.table::like(Precursor.Id, "SILAC-.-L"), Channel := "L"]
        PrecursorData[data.table::like(Precursor.Id, "SILAC-.-H"), Channel := "H"]
        PrecursorData[, Precursor.Id := gsub("-.-[LH]", "", Precursor.Id)]
        data.table::fwrite(unique(data.table::fread(list.files(pattern = "report.tsv")[1])[, .(Protein.Group, First.Protein.Description, Genes)]) |>
                             data.table::setnames(c("Protein.Group", "Genes", "First.Protein.Description"), c("Protein.Id", "Gene", "Description")),
                           "report.protein_description.tsv.gz")
      } else {return(message("ERROR: No InputFile Found"))}
    }
    PrecursorData[, Channel := gsub("H$", "Heavy", Channel)]
    PrecursorData[, Channel := gsub("L$", "Light", Channel)]

    if(file.exists("Sample_Rename.csv")){
      PrecursorData <- PrecursorData |> data.table::merge.data.table(data.table::fread("Sample_Rename.csv") |> data.table::setnames("Renamed", "Sample"))
    } else {
      PrecursorData[, Sample := gsub("_\\d+$","", Run)]
      PrecursorData[, Sample := gsub(".*_(.*_.*_.*_.*h_R.*$)","\\1", Sample)]
    }
    if(nrow(PrecursorData[Channel == "Heavy"]) == 0){return(message("ERROR: No Heavy Label Detected"))}
  }
  message("Defining Metadata")
  {
    PrecursorData[, Conc := data.table::fifelse(gsub(".*_(.*M)_.*","\\1", Sample) == Sample, "0", gsub(".*_(.*M)_.*","\\1", Sample))]
    PrecursorData[, Time := as.numeric(data.table::fifelse(gsub(".*_(.*)h_.*$","\\1", Sample) == Sample, "0", gsub(".*_(.*)h_.*$","\\1", Sample)))]
    PrecursorData[, Temp := data.table::fifelse(gsub(".*_(\\d+)C_.*","\\1", Sample) == Sample, "0", gsub(".*_(\\d+)C_.*","\\1", Sample))]
    PrecursorData[, Replicate := gsub(".*_(.*)$","\\1", Sample)]
    PrecursorData[, Sample := gsub("^0_", "", Sample)]
    PrecursorData[, Sample := gsub("_0h", "_Xh", Sample)]
    PrecursorData[, Sample := gsub("_0", "", Sample)]
    PrecursorData[, Sample := gsub("_Xh", "_0h", Sample)]
    PrecursorData[, Condition := gsub("(.*)_.*h_.*", "\\1", Sample)]
    PrecursorData[, ConditionTime := gsub("(.*)_.*", "\\1", Sample)]
    Metadata <- unique(PrecursorData[, .(Sample, Conc, Time, Temp, Replicate, Condition, ConditionTime)]) |> dplyr::arrange(Condition, Time, Replicate)
    PrecursorData$Sample <- factor(PrecursorData$Sample, levels = Metadata$Sample)
    data.table::fwrite(Metadata, "Sample_Metadata.csv")
  }
  message("Filtering Precursors")
  {
    if(ProteotypicFiltering){PrecursorData <- PrecursorData[Proteotypic == 1]}
    PrecursorData <- PrecursorData[Channel != ""]
    if(DIANNVersion >= 2.2){
      PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 &
                                       Quantity.Quality >= 0.2 & Channel.Q.Value <= 0.2]
      PrecursorData <- PrecursorData[!is.na(Precursor.Normalised) & Precursor.Normalised != 0]
    } else {
      PrecursorData <- PrecursorData[Q.Value <= 0.01 & PG.Q.Value <= 0.05 & Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & Channel.Q.Value <= 0.01]
      PrecursorData[Precursor.Quantity == 0, Precursor.Quantity   := NA]
      PrecursorData[Precursor.Translated == 0, Precursor.Translated := NA]
      PrecursorData <- PrecursorData[!is.na(Precursor.Quantity)]
      TotalNormalisation = T
    }
    if(nrow(PrecursorData[Time == 0 & Channel == "Heavy"]) > 0){message("Warning: Heavy Label in Post-Filtering 0hr Sample(s)")}
    if(TotalNormalisation){PrecursorData[, Precursor.Normalised := Precursor.Quantity/sum(Precursor.Quantity)*PrecursorData[, sum(Precursor.Quantity), Run][, median(V1)], Run]}
  }
  message("Calculate LFQs, Intensities & Counts")
  {
    PrecursorData[, Precursor.Length := nchar(Stripped.Sequence)]
    data.table::fwrite(PrecursorData, "Filtered_PrecursorData.csv.gz")

    LFQ_T <- Proteopedia::Calculate_LFQ(PrecursorData, "LFQ_T")
    Intensity_T <- PrecursorData[,.(Intensity = sum(Precursor.Normalised)), .(ProteinGroup, Sample)]
    Counts_T <- PrecursorData[ , .(N_precursors = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_T = sum(Proteotypic)), .(ProteinGroup, Sample)]

    LFQ_L <- Proteopedia::Calculate_LFQ(PrecursorData[Channel == "Light"], "LFQ_L")
    Intensity_L <- PrecursorData[Channel == "Light",.(Intensity_L = sum(Precursor.Normalised)), .(ProteinGroup, Sample)]
    Counts_L <- PrecursorData[Channel == "Light" , .(N_precursors_L = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_L = sum(Proteotypic)), .(ProteinGroup, Sample)]

    LFQ_H <- Proteopedia::Calculate_LFQ(PrecursorData[Channel == "Heavy"], "LFQ_H")
    Intensity_H <- PrecursorData[Channel == "Heavy",.(Intensity_H = sum(Precursor.Normalised)), .(ProteinGroup, Sample)]
    Counts_H <- PrecursorData[Channel == "Heavy" , .(N_precursors_H = data.table::uniqueN(Precursor.Id), N_precursors_proteotypic_H = sum(Proteotypic)), .(ProteinGroup, Sample)]
  }
  message("Calculate SILAC Ratios & Label Incorporation")
  {
    SILACRatios <- PrecursorData |> data.table::dcast(Precursor.Id+ProteinGroup+Sample ~ Channel, value.var = "Precursor.Normalised")
    SILACRatios[, HLRatio := Heavy/Light]
    SILACRatios[, Log2HLRatio := log2(HLRatio)]
    SILACRatios <- SILACRatios[!is.na(Log2HLRatio), .(Log2HLRatio = median(Log2HLRatio), .N), .(Sample, ProteinGroup)]
  }
  message("Merge & Output Protein-Level Data")
  {
    ProteinAnnotations <- unique(PrecursorData[, .(ProteinGroup, Run, Sample, Condition, ConditionTime, Conc, Time, Temp, Replicate)])
    ProteinData <- Reduce(Proteopedia::Merge_PrecursorData, list(LFQ_T, LFQ_H, LFQ_L, Intensity_T, Intensity_H, Intensity_L,
                                                                 Counts_T, Counts_H, Counts_L, SILACRatios, ProteinAnnotations))[!is.na(Intensity)]
    data.table::fwrite(ProteinData, "SILAC_DIANN_Output.csv.gz")
  }
  message("Plotting Intensities")
  {
    IntensitiesData <- data.table::rbindlist(list(
      PrecursorData[, .(Sample, Condition, ConditionTime, Replicate, Channel, `Log2 Quantity` = log2(Precursor.Normalised), Type = "Precursor Quantity")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(LFQ_L), Channel = "Light", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(Intensity_L), Channel = "Light", Type = "Protein Intensity")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(LFQ_H), Channel = "Heavy", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(Intensity_H), Channel = "Heavy", Type = "Protein Intensity")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(LFQ_T), Channel = "Total", Type = "Max. Protein LFQ")],
      ProteinData[, .(Sample, Condition, ConditionTime, Replicate, `Log2 Quantity` = log2(Intensity), Channel = "Total", Type = "Protein Intensity")]
    ), use.names = T)
    IntensitiesData[, Type := factor(Type, levels = c("Precursor Quantity", "Max. Protein LFQ", "Protein Intensity"))]
    IntensitiesData <- IntensitiesData |> data.table::merge.data.table(Metadata)
    suppressWarnings(
      IntensityPlot <- IntensitiesData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(Sample), y = `Log2 Quantity`, colour = Condition)) +
      ggplot2::geom_boxplot(outliers = F) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::facet_grid(cols = ggplot2::vars(Type), rows = ggplot2::vars(forcats::fct_rev(Channel)), scales = "free_x") + ggplot2::ylab(expression("Log"[2]~"Value")) +
      ggplot2::coord_flip() + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                                             strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
    )
  }
  message("Plotting Precursor, Peptide & Protein Counts")
  {
    CountPlot <- data.table::melt(PrecursorData[!is.na(Channel), lapply(.SD, data.table::uniqueN), .(Sample, Condition, Channel), .SDcols = c("Precursor.Id", "Stripped.Sequence", "ProteinGroup")],
                                  id.vars = c("Sample","Condition", "Channel"), value.name = "IDs") |>
      ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(Sample), y = IDs/1000, fill = Condition, alpha = Channel, label = format(IDs, big.mark = ",", scientific = F))) +
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity") + ggplot2::geom_text(size = 4, hjust = 1.1) +
      ggplot2::facet_grid(ggplot2::vars(forcats::fct_rev(Channel)), ggplot2::vars(variable), scales = "free_x",
                          labeller = ggplot2::as_labeller(c(Precursor.Id = "Precursors", Stripped.Sequence = "Peptides", ProteinGroup = "Protein Groups", Light = "Light", Heavy = "Heavy"))) +
      ggplot2::scale_y_continuous(expand = 0) + ggplot2::coord_flip() + ggplot2::ylab("No. IDs [x1,000]") + ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                     strip.text.y = ggplot2::element_text(size = 26))
  }
  message("Calculating Data Completeness")
  {
    CompletenessData <- rbind(Proteopedia::Count_Proteins(ProteinData, "All"), Proteopedia::Count_Proteins(ProteinData[N_precursors >= 2], "≥ 2"),
                              Proteopedia::Count_Proteins(ProteinData[N_precursors_proteotypic_T >= 2], "≥ 2 Proteotypic"))
    CompletenessDataLabels <- CompletenessData[Precursors == "All" & N_Samples == 1][, Tag := "All Precursors"] |> rbind(CompletenessData[Precursors == "≥ 2" & N_Samples == round(nrow(Metadata)/2)][, Tag := "≥ 2 Precursors"]) |>
      rbind(CompletenessData[Precursors == "≥ 2 Proteotypic" & N_Samples == nrow(Metadata)][, Tag := "≥ 2 Proteotypic\nPrecursors"])

    NAsPlot <- CompletenessData |> ggplot2::ggplot(ggplot2::aes(x = N_Samples, y = CumulativeProteins_N/1000, colour = Precursors))+
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::labs(x = "No. Samples", y = "No. Proteins [x1,000]") +
      ggplot2::scale_colour_manual(values = c("All" = "#000", "≥ 2" = "#999", "≥ 2 Proteotypic" = "#F63"), guide = "none") +
      ggplot2::scale_x_continuous(breaks = seq(1, 1000, 1)) + ggplot2::scale_y_continuous(limits = c(0, ceiling(max(CompletenessData$CumulativeProteins_N)/1000))) +
      ggrepel::geom_text_repel(data = CompletenessDataLabels, ggplot2::aes(label = Tag), nudge_y = -2) +
      ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")  }
  message("Plotting Missed Trypsinisation Sites")
  {
    TrypsinData <- PrecursorData |> data.table::copy()
    TrypsinData[, MissedTrypsin := grepl("[RK][^P]", Stripped.Sequence)]
    TrypsinData[, N_Trypsin := .N, .(Sample, MissedTrypsin)]
    TrypsinData[, N_Sample := .N, Sample]
    TrypsinData <- TrypsinData[MissedTrypsin == T, .(Sample, Condition, Replicate, N_Trypsin, N_Sample)] |> dplyr::distinct()
    suppressWarnings(TrypsinData[, PercentTrypsin := (N_Trypsin/N_Sample)*100])

    TrypsinisationPlot <- TrypsinData |> ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(gsub("_", " ", Sample)), y = PercentTrypsin, fill = Condition)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "", y = "Precursors with Missed Tryptic Sites (%)") + ggplot2::coord_flip()
  }
  message("Plotting Isotope Incorporation")
  {
    FullIsotopeRatio <- PrecursorData[!is.na(Precursor.Normalised) & !is.na(Channel), .(ChannelIntensity = sum(Precursor.Normalised)), list(Sample, Condition, Replicate, Channel)] |>
      data.table::merge.data.table(PrecursorData[!is.na(Precursor.Normalised) & !is.na(Channel), .(TotalIntensity = sum(Precursor.Normalised)), list(Sample, Condition, Replicate)])
    FullIsotopeRatio[, Prop := ChannelIntensity/TotalIntensity]
    FullIsotopeRatio <- FullIsotopeRatio |> data.table::merge.data.table(Metadata) |> data.table::setorderv(c("Conc", "Time", "Temp", "Replicate"))
    FullIsotopeRatio[, Sample := factor(Sample, levels = Metadata$Sample)]
    FullIsotopeRatio[, TimeReplicate := factor(paste0(Time, "_", Replicate), levels = unique(paste0(Metadata$Time, "_", Metadata$Replicate)))]

    LabelBar <- FullIsotopeRatio |> ggplot2::ggplot(ggplot2::aes(x = TimeReplicate, y = Prop, fill = Condition, alpha = Channel)) + ggplot2::facet_wrap(~Condition, nrow = 1) +
      ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::scale_alpha_manual(values = c("Heavy" = 0.5, "Light" = 1), guide = "none") + ggplot2::ylab("Isotope Channel Ratio") +
      ggplot2::scale_y_continuous(expand = 0) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), panel.spacing.x = ggplot2::unit(0, "lines"),)

    LabelLine <- FullIsotopeRatio[Channel == "Heavy"] |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Prop, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
      ggplot2::geom_line(data = FullIsotopeRatio[Channel == "Heavy", .(Prop = mean(Prop)), .(Condition, Time)]) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::ylab("Median Prop. Heavy\nIntensities")
  }
  message("Calculating Precursor & Protein Variation")
  {
    PrecursorCVs <- PrecursorData[, .(CV = Proteopedia::Calculate_CV(Precursor.Normalised), N = .N), .(Precursor.Id, Condition, Channel)]
    PrecursorCVs <- PrecursorCVs[N > 1]
    PrecursorCVs[, rank := data.table::frank(CV), .(Condition, Channel)]
    PrecursorCVs[, ID := "Precursors"]

    ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_L), N = .N), .(ProteinGroup, Condition)][, Channel := "Light"] |>
      rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(LFQ_H), N = .N), .(ProteinGroup, Condition)][, Channel := "Heavy"])
    ProteinCVs <- ProteinCVs[N > 1]
    ProteinCVs <- ProteinCVs[, rank := data.table::frank(CV), .(Condition, Channel)]
    ProteinCVs[, ID := "Protein Groups"]

    VariationPlot <- PrecursorCVs[, Precursor.Id := NULL] |> rbind(ProteinCVs[, ProteinGroup := NULL]) |>
      ggplot2::ggplot(ggplot2::aes(x = rank/1000, y = CV, colour = gsub("_", " ", Condition))) +
      ggplot2::geom_line() + ggplot2::labs(x = "No. IDs [x1,000]", y = "Variation (%)") +
      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::coord_cartesian(ylim = c(0,50)) + ggplot2::facet_grid(ggplot2::vars(forcats::fct_rev(Channel)), ggplot2::vars(ID), scales = "free") +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.25, colour = "grey70", linetype = "dotted"),
                     legend.position = "inside", legend.position.inside = c(0.7, 0.7), strip.background = ggplot2::element_blank())
  }
  message("Exporting QC Plots")
  {
    pdf("PrecursorQC_DIANN_Plot.pdf", width = 18, height = 20)
    suppressWarnings(print(IntensityPlot + CountPlot + patchwork::free(NAsPlot) + patchwork::free(TrypsinisationPlot) + LabelBar + patchwork::free(VariationPlot) +
            LabelLine + patchwork::plot_layout(design = "AAAAAA\nBBBBBB\nCCDDEE\nFFFFGG") + patchwork::plot_annotation(tag_levels = "A")))
    Proteopedia::Reset_Dev()
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_PrecursorBiochemistry <- function(InputDirectory){
  start.time <- Sys.time()
  set.seed(123)
  message("Loading Precursor Data")
  {
    setwd(InputDirectory)
    PrecursorData <- data.table::fread(list.files(pattern = ".*PrecursorData.csv.*"))[, .(Run, ProteinGroup, Stripped.Sequence, Precursor.Normalised,
                                                                                          Sample, Condition, Cell, Drug, Time, Replicate, Channel,
                                                                                          Precursor.Length)] |>
      data.table::setnames(c("Stripped.Sequence", "Precursor.Length"), c("Sequence", "Length"))

    PrecursorData[is.na(Channel), Channel := as.character(Channel)]

    if(is.na(unique(PrecursorData$Channel))){
      PrecursorData[, Channel := NULL]
      PrecursorData[, Channel := "Light"]
    }
  }
  message("Annotating with Biochemical Measures")
  {
    suppressWarnings(
      PrecursorData[, `:=`(Aliphatic_Score = Peptides::aIndex(Sequence), Boman_Interaction_Score = Peptides::boman(Sequence),
                           Hydrophobicity_Score = Peptides::hydrophobicity(Sequence, scale = "KyteDoolittle"),
                           Instability_Score = Peptides::instaIndex(Sequence), MW = Peptides::mw(Sequence),
                           pI = Peptides::pI(Sequence, pKscale = "Dawson"))]
    )

    for(ColIndex in which(colnames(PrecursorData) == "Length"):ncol(PrecursorData)){
      if(is.numeric(PrecursorData[, get(colnames(PrecursorData)[ColIndex])])){
        message(paste0("Analysing ", gsub("_", " ", colnames(PrecursorData)[ColIndex])), " Trend")

        for(ChannelIndex in unique(SubsetData$Channel)){
          SubsetData <- PrecursorData[Channel == ChannelIndex, .(Sequence, Length, Condition, Replicate, Precursor.Normalised, get(colnames(PrecursorData)[ColIndex]))] |> data.table::setnames("V6", "Subset")
          pvalue <- summary(stats::lm(Precursor.Normalised ~ Subset, data = SubsetData))$coefficients[2,4]

          pdf(paste0(ChannelIndex, gsub("(.*)_.*$", "\\1", colnames(PrecursorData)[ColIndex]),"_Trend.pdf"), width = 12, height = 10)
          print(SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Subset, y = Precursor.Normalised, colour = Condition, shape = factor(Replicate))) + ggplot2::geom_smooth(method = "lm", alpha = 0.1) +
                  ggplot2::annotate("label", x = mean(SubsetData$Subset, na.rm = T), y = min(SubsetData$Precursor.Normalised, na.rm = T)*0.93,
                                    label = paste0("P-Value: ", data.table::fifelse(pvalue < 0.01, formatC(pvalue, format = "e", digits = 2), as.character(round(pvalue, digits = 2))))) +
                  ggplot2::labs(x = paste0("Precursor ", gsub("_", " ", colnames(PrecursorData)[ColIndex])), y = "Precursor Intensity") +
                  ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities())
          Proteopedia::Reset_Dev()
        }
      }
    }
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Analyse_LabelFree_Proteins <- function(InputDirectory, Formula, MinPrecursors = 2, ImputationQ = 0.01, ImputationSigma = 1){
  start.time <- Sys.time()
  set.seed(123)
  message("Defining Comparison Groups")
  {
    ExpGroupsName = gsub("(.*)\\(.*", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula)))
    ExpGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula))), "\\+"))
    CtlGroupsName = gsub("(.*)\\(.*", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula)))
    CtlGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula))), "\\+"))
    if(nchar(ExpGroupsName) == 0){ExpGroupsName = "Experiment"}
    if(nchar(CtlGroupsName) == 0){CtlGroupsName = "Control"}
    ComparativeMetadata <- data.table::data.table(Condition = c(CtlGroups, ExpGroups))[, Comparative := data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName,
                                                                                                                            data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName, "None"))]
  }
  message("Loading Input File")
  {
    setwd(InputDirectory)
    SpectraRead <- data.table::fread(list.files(pattern = "DIANN_Output.csv"))[, Log2LFQ := log2(LFQ)]
    SpectraRead |> data.table::setnames(
      c(colnames(SpectraRead)[grepl("protein.*group", ignore.case = T, colnames(SpectraRead))],
        colnames(SpectraRead)[grepl("Gene", ignore.case = T, colnames(SpectraRead)) & !grepl("group", ignore.case = T, colnames(SpectraRead))]),
      c("ProteinGroup", "Gene"))

    Metadata <- SpectraRead[Condition %in% ComparativeMetadata$Condition, .(Sample, Condition, Replicate)] |> dplyr::distinct() |>
      data.table::merge.data.table(ComparativeMetadata)

    SpectraRead <- SpectraRead |> data.table::merge.data.table(Metadata)
    SpectraRead[, Comparative := factor(Comparative, levels = c(CtlGroupsName, ExpGroupsName))]

    if(dir.exists(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))){
      unlink(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"), recursive = T)
    }
    dir.create(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"), showWarnings = T)
    setwd(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
    data.table::fwrite(Metadata, file = "Sample_Metadata.csv")
  }
  message("Performing PCA")
  {
    PCAData <- SpectraRead[, .(ProteinGroup, Sample, Log2LFQ)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
    PCASummary <- summary(PCAData)$importance
    PCAData <- data.table::data.table(PCAData$x, keep.rownames = "Sample") |> data.table::merge.data.table(Metadata)
    PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]
    PCAPlot <- (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate)) +
                  ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::labs(x = paste("PC1 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC1"] * 100, 0), "%]", sep = ""),
                                y = paste("PC2 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC2"] * 100, 0), "%]", sep = "")) +
                  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                 legend.title = ggplot2::element_blank(), legend.position = "none")) +
      (ggplot2::ggplot(PCAData, ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate)) +
         ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
         ggplot2::labs(x = paste("PC3 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC3"] * 100, 0), "%]", sep = ""),
                       y = paste("PC4 [", round(PCASummary[rownames(PCASummary) == "Proportion of Variance", "PC4"] * 100, 0), "%]", sep = "")) +
         ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                        legend.title = ggplot2::element_blank()))
  }
  message("Filtering Proteins")
  {
    FilteringData <- SpectraRead[,.(N_Samples = .N, Min_Precursors = min(N_precursors), N_Conditions = data.table::uniqueN(Condition)), ProteinGroup]
    Retained1 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample) & Min_Precursors >= MinPrecursors]
    Retained2 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample)-1 & Min_Precursors >= (MinPrecursors+1)]
    Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraRead$Sample)/2) & Min_Precursors >= (MinPrecursors+1) & N_Conditions < length(unique(Metadata$Condition))]
    RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
    ExcludedProteins <- SpectraRead[ProteinGroup %!in% RetainedProteins]
    SpectraRead <- SpectraRead[ProteinGroup %in% RetainedProteins]
    data.table::fwrite(ExcludedProteins, file = "Excluded_Proteins.csv")
    data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")
    FilteringData <- SpectraRead[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")] |> rbind(ExcludedProteins[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")])
    FilteringData <- FilteringData[, .N, .(Sample, Condition, Replicate, Inclusion)] |> data.table::setorder(Condition)

    CountsBar <- FilteringData |> ggplot2::ggplot(ggplot2::aes(x = Sample, y = N, fill = Condition, alpha = Inclusion)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
      ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") +
      ggplot2::geom_text(ggplot2::aes(label = N), colour = data.table::fifelse(FilteringData$Inclusion == "Retained", "#FFF","#000"), position = ggplot2::position_stack(), vjust = 1.5) +
      ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x", nrow = 1) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) + ggplot2::labs(x = NULL, y = "Count", fill = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), strip.text.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0,"line"))

    UpsetPlot <- SpectraRead[, .(Sample = list(gsub("_", " ", Sample))), ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
      ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) +
      ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = SpectraRead[order(Condition), unique(gsub("_", " ", Sample))]) +
      ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15)))
  }
  message("Performing Median Normalisation")
  {
    SpectraRead[, Log2LFQ_Norm := Log2LFQ - median(Log2LFQ, na.rm = T) + median(SpectraRead$Log2LFQ, na.rm = T), Sample]

    suppressWarnings(
      NormPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraRead, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                  ggplot2::aes(x = Condition, fill = Condition, y = value, group = Sample))+
        ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
        ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") +
        ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::labs(x = NULL, y = expression("Log"[2]~"LFQ"))
    )
  }
  message("Imputing Undetected Condtiion Values")
  {
    SpectraRead <- SpectraRead[,.(ProteinGroup, Gene, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
    SpectraImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma),
                                         keep.rownames = "ProteinGroup")
    SpectraAll <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead, id.vars = "ProteinGroup", value.name = "Measured_LFQ", variable.name = "Sample"),
                                               data.table::melt.data.table(SpectraImp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ", variable.name = "Sample"))
    SpectraAll[, Data := data.table::fifelse(is.na(Measured_LFQ), "Imputed","Measured")]

    ImpPlot <- ggplot2::ggplot(SpectraAll, ggplot2::aes(x = Imputed_LFQ, fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_y_continuous(expand = 0) +
      ggplot2::scale_fill_manual(values = c(Measured = "#000", Imputed = "#C0C")) + ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) +
      ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))

    ImputedNAs <- data.table::melt.data.table(SpectraRead[ProteinGroup %in% Retained3$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
    ImputedNAs <- ImputedNAs[is.na(value)]
    ImputedNAs[, value := NULL]
    ImputedNAs <- ImputedNAs[, .(Vector = paste(Sample, collapse = ", ")), ProteinGroup]
    data.table::fwrite(ImputedNAs, "Imputed_LFQs.csv")
    SpectraAll[, Log2LFQ := data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ)]
    SpectraAll <- SpectraAll[, .(ProteinGroup, Sample, Log2LFQ)]
  }
  message("Performing Paired T-Testing")
  {
    SpectraAll <- SpectraAll |> data.table::merge.data.table(Metadata[, .(Comparative, Sample)], by = "Sample")
    SpectraAll[, LFQ := 2^(Log2LFQ)]

    SpectraTtest <- SpectraAll[, .(Log2MeanLFQ = log2(mean(LFQ)), CV = Proteopedia::Calculate_CV(LFQ), .N), .(Comparative, ProteinGroup)] |>
      data.table::dcast(ProteinGroup ~ Comparative, value.var = c("Log2MeanLFQ", "CV", "N"))

    CtlColIndex <- colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(CtlGroupsName, colnames(SpectraTtest)))]
    ExpColIndex <-  colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(ExpGroupsName, colnames(SpectraTtest)))]
    SpectraTtest$Log2FC <- SpectraTtest[, get(ExpColIndex)] - SpectraTtest[, get(CtlColIndex)]

    Ttest_Output <- SpectraAll[, P.Value := stats::t.test(Log2LFQ ~ Comparative)$p.value, ProteinGroup]
    SpectraTtest <- data.table::merge.data.table(data.table::data.table(SpectraTtest), Ttest_Output, by = "ProteinGroup")
    SpectraTtest <- Proteopedia::Add_ProteinInfo(SpectraTtest, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
    SpectraTtest <- SpectraTtest |> data.table::merge.data.table(ProteinInfo)
    SpectraTtest[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
    SpectraTtest[, Log2FC := as.numeric(Log2FC)]
    data.table::fwrite(SpectraTtest |> dplyr::distinct(), "Paired_T-Test_Output.csv")
  }
  message("Fitting Linear Model")
  {
    ModelDesign <- stats::model.matrix(~0 + Comparative, data = Metadata)
    colnames(ModelDesign) <- gsub("Comparative", "", colnames(ModelDesign))
    rownames(ModelDesign) <- Metadata$Sample
    ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroupsName, colnames(ModelDesign))), which(grepl(ExpGroupsName, colnames(ModelDesign)))))]

    ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
    ContrastMatrix[,1] <- c(-1,1)

    LimmaInput <- SpectraAll |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
      data.table::setcolorder(c("ProteinGroup", Metadata$Sample))

    suppressMessages(ModelFit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, ModelDesign), ContrastMatrix)))

    if(!is.finite(ModelFit$df.prior)){message("Warning: Limma Prior is Infinite")}

    MeanVarData <- data.table::data.table(ModelFit$genes, "Mean" = ModelFit$Amean, "Variance" = sqrt(ModelFit$sigma))
    MeanVarData[, Data := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured")]

    MeanVarPlot <- MeanVarData |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) +
      ggplot2::geom_point(colour = "#000") +  ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C") +
      ggplot2::labs(x = "Mean Log2LFQ", y = "Variance", colour = NULL)

    LimmaOutput <- data.table::data.table(limma::topTable(ModelFit, coef=1, adjust.method = "BH", n=Inf)) |> data.table::setnames("logFC", "Log2FC")
    LimmaOutput <- LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = T)]
    LimmaOutput <- Proteopedia::Add_ProteinInfo(LimmaOutput, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
    LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
    LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = T)
    LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
    LimmaOutput <- Proteopedia::Separate_Isoforms(LimmaOutput)
    data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
  }
  message("Generating Volcano Plots")
  {
    LimmaVolcano <- LimmaOutput |> dplyr::arrange(desc(abs(t))) |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
      ggrepel::geom_text_repel(ggplot2::aes(label= data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
      ggplot2::geom_vline(xintercept = mean(LimmaOutput$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
      Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

    TtestVolcano <- unique(SpectraTtest[, .(ProteinGroup, Log2FC, P.Value, Gene)]) |> dplyr::arrange(desc(abs(Log2FC)*-log10(P.Value))) |>
      ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
      ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
      ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
      ggplot2::geom_vline(xintercept = mean(SpectraTtest$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
      Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())

    pdf("VolcanoPlots.pdf", height = 10, width = 14)
    print(LimmaVolcano)
    print(LimmaVolcano + ggplot2::geom_point(data = LimmaOutput[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed Limma"))
    print(TtestVolcano)
    print(TtestVolcano + ggplot2::geom_point(data = SpectraTtest[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed T-Test"))
    Proteopedia::Reset_Dev()
  }
  message("Exporting QC Plot")
  {
    pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
    suppressWarnings(
      print(patchwork::free(PCAPlot) + CountsBar + patchwork::free(UpsetPlot) + NormPlot + ImpPlot + patchwork::free(MeanVarPlot, type = "label") +
              patchwork::plot_layout(design = "AAAABBBB\nCCCCDDDD\nEEEEFFFF") + patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E", "F"))))
    )
    Proteopedia::Reset_Dev()
  }
  message("Exporting HTML Volcano Plot")
  {
    InteractiveData <- data.table::data.table("Protein" = LimmaOutput$Gene, "Log2FC" = round(LimmaOutput$Log2FC, digits = 2),
                                              "PValue" = LimmaOutput$P.Value, "Significance" = factor(LimmaOutput$Significance),
                                              "Scien_PValue" = formatC(LimmaOutput$P.Value, format = "e", digits = 2),
                                              "GeneGroup" = LimmaOutput$GeneGroup, "URL" = LimmaOutput$URL)

    HighCharterVolcano <- highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC, y = -log10(PValue), group = Significance)) |>
      highcharter::hc_chart(zoomType = "xy") |>
      highcharter::hc_xAxis(title = list(text = paste0(unique(Metadata$Condition)[1]," vs ", unique(Metadata$Condition)[2]," Log2 Fold-Change")),
                            lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
      highcharter::hc_yAxis(title = list(text = "-Log10 P-Value"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
      highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}") |>
      highcharter::hc_plotOptions(scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = T), inactive = list(enabled = F)),
                                                 point = list(events = list( click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))) |>
      highcharter::hc_colors(c("#999", "#800","#03F"))
    htmlwidgets::saveWidget(HighCharterVolcano, "InteractiveVolcanoPlot.html")
  }
  message("Exporting Analysis Parameters")
  data.table::fwrite(data.table::data.table("Experimental Condition(s)" = paste(ExpGroups, collapse = ", "), "Experimental Name" = paste0(ExpGroupsName),
                                            "Control Condition(s)" = paste(CtlGroups, collapse = ", "), "Control Name" = paste0(CtlGroupsName),
                                            "Min_Precursors" = paste0(MinPrecursors), "Imputation Q-Value" = ImputationQ,
                                            "Imputation Sigma" = ImputationSigma), "Analysis_Parameters.csv")
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Analyse_StaticSILAC_Proteins <- function(InputDirectory, Formula, MinPrecursors = 2, ImputationQ = 0.01, ImputationSigma = 1, DIANNVersion = 2.2){
  set.seed(123)
  start.time <- Sys.time()
  message("Defining Comparison Groups")
  {
    ExpGroupsName = gsub("(.*)\\(.*", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula)))
    ExpGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula))), "\\+"))
    CtlGroupsName = gsub("(.*)\\(.*", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula)))
    CtlGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula))), "\\+"))
    if(nchar(ExpGroupsName) == 0){ExpGroupsName = "Experiment"}
    if(nchar(CtlGroupsName) == 0){CtlGroupsName = "Control"}
    ComparativeMetadata <- data.table::data.table(Condition = c(CtlGroups, ExpGroups))[, Comparative := data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName,
                                                                                                                            data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName, "None"))]
  }
  message("Loading & Formatting Data")
  {
    setwd(InputDirectory)
    if(length(list.files(pattern = "SILAC_DIANN_Output.csv")) > 0){
      InputFile <- list.files(pattern = "SILAC_DIANN_Output.csv")[1]

      data.table::fread("Sample_Metadata.csv")

      Metadata <- data.table::fread(InputFile)[Condition %in% ComparativeMetadata$Condition, .(Sample, Condition, ConditionTime, Replicate)] |> dplyr::distinct() |>
        data.table::merge.data.table(ComparativeMetadata)
      SpectraRead <- data.table::fread(InputFile)[, Log2LFQ_H := log2(LFQ_H)][, Log2LFQ_L := log2(LFQ_L)][, Log2LFQ := log2(LFQ_T)][Sample %in% Metadata$Sample]
    } else {return(message("ERROR: No InputFile Found"))}
    ModelDesign <- stats::model.matrix(~0 + Comparative, data = Metadata)
    colnames(ModelDesign) <- gsub("Comparative", "", colnames(ModelDesign))
    rownames(ModelDesign) <- Metadata$Sample
    ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroupsName, colnames(ModelDesign))), which(grepl(ExpGroupsName, colnames(ModelDesign)))))]
    ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
    ContrastMatrix[, 1] <- c(-1, 1)
    data.table::fwrite(Metadata, file = "Sample_Metadata.csv")

    if(dir.exists(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))){
      unlink(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"), recursive = T)
    }
    dir.create(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"), showWarnings = T)
    setwd(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
    data.table::fwrite(Metadata, file = "Sample_Metadata.csv")  }
  message("Analysing Total Intensity")
  {
    setwd(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
    if(dir.exists("Total_Analysis")){
      unlink("Total_Analysis", recursive = T)
    }
    dir.create("Total_Analysis", showWarnings = T)
    setwd("Total_Analysis")
    SpectraReadT <- SpectraRead |> data.table::copy()
    message("Total Abundance: Performing PCA")
    {
      PCAData <- stats::prcomp(t(data.frame(tidyr::drop_na(data.table::dcast(SpectraReadT, ProteinGroup ~ Sample, value.var = "Log2LFQ", values_fill = NA)), row.names = "ProteinGroup")), scale. = T)
      PCASummary <- summary(PCAData)$importance
      PCAData <- data.table::merge.data.table(data.table::data.table(PCAData$x, keep.rownames = "Sample"), Metadata)
      PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]
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
    message("Total Abundance: Filtering Proteins")
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
      FilteringPlotData <- rbind(SpectraReadT[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")],
                                 FilteredProteins[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")])[, .N , .(Sample, Condition, Replicate, Inclusion)]

      CountsBar <- ggplot2::ggplot(data.table::setorder(FilteringPlotData, Sample), ggplot2::aes(x = gsub("_", " ", Sample), y = N, fill = Condition, alpha = Inclusion)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::scale_alpha_manual(values = c(Excluded = 0.4, Retained = 1), guide = "none") +
        ggplot2::geom_text(ggplot2::aes(label = N, colour = Inclusion), size = 6, position = ggplot2::position_stack(), vjust = 1.5) +
        ggplot2::scale_colour_manual(values = c("Retained" = "#FFF", "Excluded" = "#000"), guide = "none") + ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) + ggplot2::labs(x = NULL, y = "No. Proteins", fill = NULL) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0, "line"))

      UpsetPlot <- ggplot2::ggplot(SpectraReadT[,.(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup],
                                   ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
        ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) +
        ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = gsub("_", " ", SpectraReadT[order(Condition), unique(Sample)])) +
        ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))
    }
    message("Total Abundance: Performing Median Normalisation")
    {
      SpectraReadT[, `:=`(Log2LFQ_Norm, Log2LFQ - median(Log2LFQ, na.rm = T) + median(SpectraReadT$Log2LFQ, na.rm = T)), by = Sample]
      NormPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraReadT, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                           ggplot2::aes(x = Condition, colour = Condition,  y = value, group = Sample)) +
        ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c(Log2LFQ = "Pre", Log2LFQ_Norm = "Post"))) +
        ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::labs(x = NULL, y = expression("Log"[2]~ "LFQ")) + ggplot2::theme(strip.background = ggplot2::element_blank())
    }
    message("Total Abundance: Imputing NA Values")
    {
      SpectraReadT <- SpectraReadT[, .(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup~Sample, value.var = "Log2LFQ_Norm")
      SpectraImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraReadT, rownames = "ProteinGroup"), q = ImputationQ,
                                                                      tune.sigma = ImputationSigma), keep.rownames = "ProteinGroup")
      SpectraAll <- data.table::merge.data.table(data.table::melt.data.table(SpectraReadT, id.vars = "ProteinGroup", value.name = "Measured_LFQ",
                                                                             variable.name = "Sample"),
                                                 data.table::melt.data.table(SpectraImp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ",
                                                                             variable.name = "Sample"))
      SpectraAll[, `:=`(Data, data.table::fifelse(is.na(Measured_LFQ), "Imputed", "Measured"))]

      ImpPlot <- ggplot2::ggplot(SpectraAll, ggplot2::aes(x = Imputed_LFQ, fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_y_continuous(expand = 0) +
        ggplot2::scale_fill_manual(values = c(Measured = "#000", Imputed = "#C0C")) + ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) +
        ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))

      ImputedNAs <- data.table::melt.data.table(SpectraReadT[ProteinGroup %in% Retained3$ProteinGroup],
                                                id.vars = "ProteinGroup", variable.name = "Sample")[is.na(value)]
      ImputedNAs[, `:=`(value, NULL)]
      ImputedNAs <- ImputedNAs[, .(Vector = paste(Sample, collapse = ", ")), ProteinGroup]
      data.table::fwrite(ImputedNAs, "Imputed_LFQs.csv")
      SpectraAll[, `:=`(Log2LFQ, data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ))]
      SpectraAll <- SpectraAll[, .(ProteinGroup, Sample, Log2LFQ)]
    }
    message("Total Abundance: Performing Paired T-Testing")
    {
      SpectraAll <- SpectraAll |> data.table::merge.data.table(Metadata[, .(Comparative, Sample)], by = "Sample")
      SpectraAll[, LFQ := 2^(Log2LFQ)]

      SpectraTtest <- SpectraAll[, .(Log2MeanLFQ = log2(mean(LFQ)), CV = Proteopedia::Calculate_CV(LFQ), .N), .(Comparative, ProteinGroup)] |>
        data.table::dcast(ProteinGroup ~ Comparative, value.var = c("Log2MeanLFQ", "CV", "N"))

      CtlColIndex <- colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(CtlGroupsName, colnames(SpectraTtest)))]
      ExpColIndex <-  colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(ExpGroupsName, colnames(SpectraTtest)))]
      SpectraTtest$Log2FC <- SpectraTtest[, get(ExpColIndex)] - SpectraTtest[, get(CtlColIndex)]

      Ttest_Output <- SpectraAll[, P.Value := stats::t.test(Log2LFQ ~ Comparative)$p.value, ProteinGroup]
      Ttest_Output <- Proteopedia::Add_ProteinInfo(Ttest_Output, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      SpectraTtest <- data.table::merge.data.table(data.table::data.table(SpectraTtest), Ttest_Output)
      SpectraTtest[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
      SpectraTtest[, Log2FC := as.numeric(Log2FC)]
      data.table::fwrite(SpectraTtest |> dplyr::distinct(), "Paired_T-Test_Output.csv")
    }
    message("Total Abundance: Fitting Linear Model")
    {
      ModelDesign <- stats::model.matrix(~0 + Comparative, data = Metadata)
      colnames(ModelDesign) <- gsub("Comparative", "", colnames(ModelDesign))
      rownames(ModelDesign) <- Metadata$Sample
      ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroupsName, colnames(ModelDesign))), which(grepl(ExpGroupsName, colnames(ModelDesign)))))]

      ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
      ContrastMatrix[,1] <- c(-1,1)

      LimmaInput <- SpectraAll |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
        data.table::setcolorder(c("ProteinGroup", Metadata$Sample))

      suppressMessages(ModelFit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, ModelDesign), ContrastMatrix)))

      if(!is.finite(ModelFit$df.prior)){message("Warning: Limma Prior is Infinite")}

      MeanVarData <- data.table::data.table(ModelFit$genes, "Mean" = ModelFit$Amean, "Variance" = sqrt(ModelFit$sigma))
      MeanVarData[, Data := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured")]

      MeanVarPlot <- MeanVarData |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) +
        ggplot2::geom_point(colour = "#000") +  ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C") +
        ggplot2::labs(x = expression("Mean Log"[2]~"LFQ"), y = "Variance", colour = NULL)

      LimmaOutput <- data.table::data.table(limma::topTable(ModelFit, coef=1, adjust.method = "BH", n=Inf)) |> data.table::setnames("logFC", "Log2FC")
      LimmaOutput <- LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = T)]
      LimmaOutput <- Proteopedia::Add_ProteinInfo(LimmaOutput, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
      LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = T)
      LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
      LimmaOutput <- Proteopedia::Separate_Isoforms(LimmaOutput)
      data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
    }
    message("Total Abundance: Generating Volcano Plots")
    {
      LimmaVolcano <- LimmaOutput |> dplyr::arrange(desc(abs(t))) |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
        ggrepel::geom_text_repel(ggplot2::aes(label= data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
        ggplot2::geom_vline(xintercept = mean(LimmaOutput$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      TtestVolcano <- unique(SpectraTtest[, .(ProteinGroup, Log2FC, P.Value, Gene)]) |> dplyr::arrange(desc(abs(Log2FC)*-log10(P.Value))) |>
        ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
        ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
        ggplot2::geom_vline(xintercept = mean(SpectraTtest$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())

      pdf("VolcanoPlots.pdf", height = 10, width = 14)
      print(LimmaVolcano)
      print(LimmaVolcano + ggplot2::geom_point(data = LimmaOutput[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed Limma"))
      print(TtestVolcano)
      print(TtestVolcano + ggplot2::geom_point(data = SpectraTtest[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed T-Test"))
      Proteopedia::Reset_Dev()
    }
    message("Total Abundance: Exporting QC Plot")
    {
      pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
      suppressWarnings(
        print(patchwork::free(PCAPlot) + CountsBar + patchwork::free(UpsetPlot) + NormPlot + ImpPlot + patchwork::free(MeanVarPlot, type = "label") +
                patchwork::plot_layout(design = "AAAABBBB\nCCDDDEEE\nFFFFFFFF") + patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E", "F"))))
      )
      Proteopedia::Reset_Dev()
    }
    message("Total Abundance: Exporting HTML Volcano Plot")
    {
      InteractiveData <- data.table::data.table("Protein" = LimmaOutput$Gene, "Log2FC" = round(LimmaOutput$Log2FC, digits = 2),
                                                "PValue" = LimmaOutput$P.Value, "Significance" = factor(LimmaOutput$Significance),
                                                "Scien_PValue" = formatC(LimmaOutput$P.Value, format = "e", digits = 2),
                                                "GeneGroup" = LimmaOutput$GeneGroup, "URL" = LimmaOutput$URL)

      HighCharterVolcano <- highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC, y = -log10(PValue), group = Significance)) |>
        highcharter::hc_chart(zoomType = "xy") |>
        highcharter::hc_xAxis(title = list(text = paste0(unique(Metadata$Condition)[1]," vs ", unique(Metadata$Condition)[2]," Log2 Fold-Change")),
                              lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
        highcharter::hc_yAxis(title = list(text = "-Log10 P-Value"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
        highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}") |>
        highcharter::hc_plotOptions(scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = T), inactive = list(enabled = F)),
                                                   point = list(events = list( click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))) |>
        highcharter::hc_colors(c("#999", "#800","#03F"))
      htmlwidgets::saveWidget(HighCharterVolcano, "InteractiveVolcanoPlot.html")
    }
  }
  message("Analysing Each Channel For Stability/Synthesis Measures")
  {
    setwd(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
    if(dir.exists("Channel_Analysis")){
      unlink("Channel_Analysis", recursive = T)
    }
    dir.create("Channel_Analysis", showWarnings = T)
    setwd("Channel_Analysis")
    SpectraRead_L <- SpectraRead |> data.table::copy()
    SpectraRead_H <- SpectraRead |> data.table::copy()
    message("SILAC Channels: Performing PCA")
    {
      PCAData_L <- SpectraRead_L[, .(ProteinGroup, Sample, Log2LFQ_L)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_L") |>
        tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
      PCASummary_L <- summary(PCAData_L)$importance
      PCAData_L <- data.table::data.table(PCAData_L$x, keep.rownames = "Sample") |> data.table::merge.data.table(Metadata)
      PCAData_L[, Dataset := "Light"]

      PCAData_H <- SpectraRead_H[, .(ProteinGroup, Sample, Log2LFQ_H)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_H") |>
        tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
      PCASummary_H <- summary(PCAData_H)$importance
      PCAData_H <- data.table::data.table(PCAData_H$x, keep.rownames = "Sample") |> data.table::merge.data.table(Metadata)
      PCAData_H[, Dataset := "Heavy"]

      PCAData <- PCAData_H |> rbind(PCAData_L)
      PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]
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
      FilteringData <- SpectraRead[,.(N_Samples = .N, Min_Precursors_L = min(N_precursors_L), Min_Precursors_H = min(N_precursors_H),
                                      N_Conditions = data.table::uniqueN(Condition)), ProteinGroup]

      Retained1_L <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead_L$Sample) & Min_Precursors_L >= MinPrecursors]
      Retained2_L <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead_L$Sample)-1 & Min_Precursors_L >= (MinPrecursors+1)]
      Retained3_L <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraRead_L$Sample)/2) & Min_Precursors_L >= (MinPrecursors+1) & N_Conditions < length(unique(Metadata$Condition))]
      RetainedProteins_L <- c(Retained1_L$ProteinGroup, Retained2_L$ProteinGroup, Retained3_L$ProteinGroup)
      FilteredProteins_L <- SpectraRead_L[ProteinGroup %!in% RetainedProteins_L]
      SpectraRead_L <- SpectraRead_L[ProteinGroup %in% RetainedProteins_L] |> dplyr::select(!ends_with("_H"))
      data.table::fwrite(FilteredProteins_L, file = "Light_Filtered_Proteins.csv")
      data.table::fwrite(Retained3_L, file = "Light_Imputed_Proteins.csv")

      Retained1_H <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead_H$Sample) & Min_Precursors_H >= MinPrecursors]
      Retained2_H <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead_H$Sample)-1 & Min_Precursors_H >= (MinPrecursors+1)]
      Retained3_H <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraRead_H$Sample)/2) & Min_Precursors_H >= (MinPrecursors+1) & N_Conditions < length(unique(Metadata$Condition))]
      RetainedProteins_H <- c(Retained1_H$ProteinGroup, Retained2_H$ProteinGroup, Retained3_H$ProteinGroup)
      FilteredProteins_H <- SpectraRead_H[ProteinGroup %!in% RetainedProteins_H]
      SpectraRead_H <- SpectraRead_H[ProteinGroup %in% RetainedProteins_H] |> dplyr::select(!ends_with("_L"))
      data.table::fwrite(FilteredProteins_H, file = "Heavy_Filtered_Proteins.csv")
      data.table::fwrite(Retained3_H, file = "Heavy_Imputed_Proteins.csv")

      FilteringData <- rbind(SpectraRead_L[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained", Dataset = "Light")], FilteredProteins_L[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded", Dataset = "Light")])[, .N, .(Sample, Condition, Replicate, Inclusion, Dataset)] |>
        rbind(rbind(SpectraRead_H[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained", Dataset = "Heavy")], FilteredProteins_H[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded", Dataset = "Heavy")])[, .N, .(Sample, Condition, Replicate, Inclusion, Dataset)])

      CountsBar <- FilteringData |> ggplot2::ggplot(ggplot2::aes(x = gsub("_", " ", Sample), y = N, fill = Condition, alpha = Inclusion)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") +
        ggplot2::geom_text(ggplot2::aes(label = N), size = 6, colour = ifelse(FilteringData$Inclusion == "Retained", "#FFF","#000"), position = ggplot2::position_stack(), vjust = 1.5) +
        ggplot2::facet_wrap(~Dataset, scales = "free_x") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) +
        ggplot2::labs(x = NULL, y = "No. Proteins") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = ggplot2::element_blank(),
                                                                     strip.text.y = ggplot2::element_text(size = 26))

      UpsetPlot <- (SpectraRead_H[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
                      ggplot2::geom_text(stat="count", ggplot2::aes(label =ggplot2::after_stat(count)), vjust = -0.5, size = 3) + ggplot2::ggtitle("Light") +
                      ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = gsub("_", " ", SpectraRead_H[order(Condition), unique(Sample)])) +
                      ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) +
                      ggplot2::theme(plot.title = ggplot2::element_text(size = 24))) +
        (SpectraRead_L[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
           ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) + ggplot2::ggtitle("Heavy") +
           ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = gsub("_", " ", SpectraRead_L[order(Condition), unique(Sample)])) +
           ggplot2::labs(x = NULL, y = NULL) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) +
           ggplot2::theme(plot.title = ggplot2::element_text(size = 24), axis.text.x = ggplot2::element_blank())) + patchwork::plot_layout(guides = "collect")
    }
    message("SILAC Channels: Performing Median Normalisation")
    {
      SpectraRead_L[, Log2LFQ_Norm_L := Log2LFQ_L - median(Log2LFQ_L, na.rm = T) + median(SpectraRead_L$Log2LFQ_L, na.rm = T), by = Sample]
      colnames(SpectraRead_L) <- gsub("_L", "", colnames(SpectraRead_L))
      SpectraRead_L[ , Dataset := "Light"]

      SpectraRead_H[, Log2LFQ_Norm_H := Log2LFQ_H - median(Log2LFQ_H, na.rm = T) + median(SpectraRead_H$Log2LFQ_H, na.rm = T), by = Sample]
      colnames(SpectraRead_H) <- gsub("_H", "", colnames(SpectraRead_H))
      SpectraRead_H[ , Dataset := "Heavy"]

      NormPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraRead_L |> rbind(SpectraRead_H), measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                           ggplot2::aes(x = gsub("_", " ", Condition), colour = Condition, y = value, group = Sample))+
        ggplot2::facet_grid(c("Dataset","variable"), labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
        ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::labs(x = NULL, y = expression("Log"[2]~"LFQ")) + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))
    }
    message("SILAC Channels: Imputing NA Values")
    {
      SpectraRead_L <- SpectraRead_L[, .(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
      SpectraImpL <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead_L, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma),
                                            keep.rownames = "ProteinGroup")
      SpectraAllL <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead_L, id.vars = "ProteinGroup", value.name = "MeasuredLFQ", variable.name = "Sample"),
                                                  data.table::melt.data.table(SpectraImpL, id.vars = "ProteinGroup", value.name = "ImputedLFQ", variable.name = "Sample"))
      SpectraAllL[, `:=`(Data = data.table::fifelse(is.na(MeasuredLFQ), "Imputed","Measured"), Dataset = "Light")]

      SpectraRead_H <- SpectraRead_H[, .(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
      SpectraImpH <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead_L, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma),
                                            keep.rownames = "ProteinGroup")
      SpectraAllH <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead_H, id.vars = "ProteinGroup", value.name = "MeasuredLFQ", variable.name = "Sample"),
                                                  data.table::melt.data.table(SpectraImpH, id.vars = "ProteinGroup", value.name = "ImputedLFQ", variable.name = "Sample"))
      SpectraAllH[, `:=`(Data = data.table::fifelse(is.na(MeasuredLFQ), "Imputed","Measured"), Dataset = "Heavy")]

      SpectraAll <- SpectraAllL |> rbind(SpectraAllH)

      ImpPlot <- SpectraAll |> ggplot2::ggplot(ggplot2::aes(x = ImputedLFQ, fill = Data)) + ggplot2::scale_y_continuous(expand = 0) +
        ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_fill_manual(values = c("Measured" = "#000", "Imputed" = "#C0C")) +
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
    message("SILAC Channels: Comparing Condition LFQs")
    {
      LFQComparison <- SpectraAll |> data.table::copy() |> data.table::merge.data.table(Metadata)
      LFQComparison <- LFQComparison[, .(MeanLog2LFQ = mean(Log2LFQ)), .(Condition, Comparative, Dataset, ProteinGroup)] |>
        data.table::dcast(ProteinGroup+Dataset ~ Condition, value.var = "MeanLog2LFQ")

      LFQComparison |> data.table::setnames(c(colnames(LFQComparison)[grepl(CtlGroupsName, colnames(LFQComparison))],
                                              colnames(LFQComparison)[grepl(ExpGroupsName, colnames(LFQComparison))]),
                                            c("CtlGroups", "ExpGroups"), skip_absent = T)
      LFQComparison <- Proteopedia::Add_ProteinInfo(LFQComparison, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LFQComparison[, Histone := data.table::fifelse(grepl("^H\\d.*", Gene), T, F)]

      pdf(paste0(ExpGroupsName, "_", CtlGroupsName, "_LFQComp.pdf"), height = 10, width = 14)
      print(LFQComparison |> ggplot2::ggplot(ggplot2::aes(x = CtlGroups, y = ExpGroups)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_XYLine("#999") + Proteopedia::Add_Pearsons(Subgroups = F) + ggplot2::facet_wrap(~Dataset) +
              ggplot2::annotate("text", x = max(LFQComparison$CtlGroups, na.rm = T)*0.9, y = min(LFQComparison$ExpGroups, na.rm = T)*0.9, label = "Histones", colour = "#F0F") +
              ggrepel::geom_text_repel(ggplot2::aes(label = Gene)) + ggplot2::geom_point(data = LFQComparison[Histone ==T], stroke = NA, colour = "#F0F") +
              ggplot2::labs(x = expression("Log"[2]~"LFQ"), y = expression("Log"[2]~"LFQ"), title = paste0(ExpGroupsName, " vs. ", CtlGroupsName)))
      Proteopedia::Reset_Dev()
    }
    message("SILAC Channels: Performing Paired T-Testing")
    {
      SpectraAll <- SpectraAll |> data.table::merge.data.table(Metadata)
      SpectraAll[, LFQ := 2^(Log2LFQ)]
      # Light T-Test
      {
        SpectraTtest_L <- SpectraAll[Dataset == "Light", .(.N, Log2MeanLFQ = log2(mean(LFQ)), CV = Proteopedia::Calculate_CV(LFQ)), .(Comparative, ProteinGroup)] |>
          data.table::dcast(ProteinGroup ~ Comparative, value.var = c("N", "Log2MeanLFQ", "CV"))

        CtlColIndex <- colnames(SpectraTtest_L)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_L)) & grepl(CtlGroupsName, colnames(SpectraTtest_L)))]
        ExpColIndex <-  colnames(SpectraTtest_L)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_L)) & grepl(ExpGroupsName, colnames(SpectraTtest_L)))]
        SpectraTtest_L$Log2FC <- SpectraTtest_L[, get(ExpColIndex)] - SpectraTtest_L[, get(CtlColIndex)]

        Ttest_Output_L <- SpectraAll[, P.Value := stats::t.test(Log2LFQ ~ Comparative)$p.value, ProteinGroup]
        Ttest_Output_L <- Proteopedia::Add_ProteinInfo(Ttest_Output_L, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
        SpectraTtest_L <- data.table::merge.data.table(data.table::data.table(SpectraTtest_L), Ttest_Output_L, by = "ProteinGroup")
        SpectraTtest_L[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
        SpectraTtest_L[, Log2FC := as.numeric(Log2FC)]
        data.table::fwrite(SpectraTtest_L |> dplyr::distinct(), "Paired_T-Test_Output.csv")
      }
      # Heavy T-Test
      {
        SpectraTtest_H <- SpectraAll[Dataset == "Heavy", .(.N, Log2MeanLFQ = log2(mean(LFQ)), CV = Proteopedia::Calculate_CV(LFQ)), .(Comparative, ProteinGroup)] |>
          data.table::dcast(ProteinGroup ~ Comparative, value.var = c("N", "Log2MeanLFQ", "CV"))

        CtlColIndex <- colnames(SpectraTtest_H)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_H)) & grepl(CtlGroupsName, colnames(SpectraTtest_H)))]
        ExpColIndex <-  colnames(SpectraTtest_H)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest_H)) & grepl(ExpGroupsName, colnames(SpectraTtest_H)))]
        SpectraTtest_H$Log2FC <- SpectraTtest_H[, get(ExpColIndex)] - SpectraTtest_H[, get(CtlColIndex)]

        Ttest_Output_L <- SpectraAll[, P.Value := stats::t.test(Log2LFQ ~ Comparative)$p.value, ProteinGroup]
        Ttest_Output_L <- Proteopedia::Add_ProteinInfo(Ttest_Output_L, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
        SpectraTtest_H <- data.table::merge.data.table(data.table::data.table(SpectraTtest_H), Ttest_Output_L, by = "ProteinGroup")
        SpectraTtest_H[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
        SpectraTtest_H[, Log2FC := as.numeric(Log2FC)]
        data.table::fwrite(SpectraTtest_H |> dplyr::distinct(), "Paired_T-Test_Output.csv")
      }
      SpectraTtest <- SpectraTtest_L |> rbind(SpectraTtest_H)
      data.table::fwrite(SpectraTtest, "Paired_T-Test_Output.csv")
    }
    message("SILAC Channels: Fitting to Linear Model & Exporting Data")
    {
      LimmaInput_L <- SpectraAll[Dataset == "Light"] |> data.table::copy() |> data.table::dcast(ProteinGroup+Dataset ~ Sample, value.var = "Log2LFQ") |> data.table::setcolorder(c("ProteinGroup", rownames(ModelDesign)))
      LimmaInput_H <- SpectraAll[Dataset == "Heavy"] |> data.table::copy() |> data.table::dcast(ProteinGroup+Dataset ~ Sample, value.var = "Log2LFQ") |> data.table::setcolorder(c("ProteinGroup", rownames(ModelDesign)))

      suppressMessages(ModelFit_L <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput_L[, Dataset := NULL], ModelDesign), contrasts = ContrastMatrix)))
      suppressMessages(ModelFit_H <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput_H[, Dataset := NULL], ModelDesign), contrasts = ContrastMatrix)))

      MeanVarData <- data.table::data.table(ModelFit_L$genes, "Mean" = ModelFit_L$Amean, "Variance" = sqrt(ModelFit_L$sigma), "Dataset" = "Light") |>
        rbind(data.table::data.table(ModelFit_H$genes, "Mean" = ModelFit_H$Amean, "Variance" = sqrt(ModelFit_H$sigma), "Dataset" = "Heavy"))

      MeanVarPlot <- MeanVarData[, Data := data.table::fifelse(Dataset == "Light" & ProteinGroup %in% ImputedNAs_L$ProteinGroup, "Imputed",
                                                               data.table::fifelse(Dataset == "Heavy" & ProteinGroup %in% ImputedNAs_H$ProteinGroup, "Imputed", "Measured"))] |>
        ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) + ggplot2::geom_point(stroke = NA) +
        ggplot2::scale_colour_manual(values = c(Measured = "#000", Imputed = "#C0C"), guide = "none") +
        ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C", stroke = NA) +
        ggplot2::labs(x = expression("Mean Log"[2]~"Protein LFQ"), y = "Variance") + ggplot2::facet_wrap(~Dataset, scales = "free_y", nrow = 2, strip.position = "right") +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 20), strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26),
                       panel.spacing.x = ggplot2::unit(0,"line"))

      LimmaOutput_L <- limma::topTable(ModelFit_L, coef=1, adjust.method = "BH", n=Inf) |> data.table::data.table() |> data.table::setnames("logFC", "Log2FC")
      LimmaOutput_L <- Proteopedia::Add_ProteinInfo(LimmaOutput_L, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LimmaOutput_L[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
      LimmaOutput_L <- LimmaOutput_L |> data.table::merge.data.table(LimmaInput_L, by = "ProteinGroup", sort = F) |> unique()
      LimmaOutput_L <- Proteopedia::Separate_Isoforms(LimmaOutput_L)
      LimmaOutput_L[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs_L$ProteinGroup, "Imputed", "Measured")]
      LimmaOutput_L <- LimmaOutput_L[order(abs(LimmaOutput_L$t), decreasing = T)]
      data.table::fwrite(LimmaOutput_L, file = "Light_Limma_Output.csv")

      LimmaOutput_H <- limma::topTable(ModelFit_H, coef=1, adjust.method = "BH", n=Inf) |> data.table::data.table() |> data.table::setnames("logFC", "Log2FC")
      LimmaOutput_H <- Proteopedia::Add_ProteinInfo(LimmaOutput_H, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LimmaOutput_H[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
      LimmaOutput_H <- LimmaOutput_H |> data.table::merge.data.table(LimmaInput_H, by = "ProteinGroup", sort = F) |> unique()
      LimmaOutput_H <- Proteopedia::Separate_Isoforms(LimmaOutput_H)
      LimmaOutput_H[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs_H$ProteinGroup, "Imputed", "Measured")]
      LimmaOutput_H <- LimmaOutput_H[order(abs(LimmaOutput_H$t), decreasing = T)]
      data.table::fwrite(LimmaOutput_H, file = "Heavy_Limma_Output.csv")

      data.table::fwrite(LimmaOutput_L[, .(ProteinGroup, Isoforms, Gene, Log2FC, P.Value, AveExpr, t, adj.P.Val, B, Significance, Imputed, GeneGroup)][, Dataset := "Light"] |>
                           rbind(LimmaOutput_H[, .(ProteinGroup, Isoforms, Gene, Log2FC, P.Value, AveExpr, t, adj.P.Val, B, Significance, Imputed, GeneGroup)][, Dataset := "Heavy"]), file = "Limma_Output.csv")
    }
    message("SILAC Channels: Generating Volcano Plots")
    {
      MeanLog2FC_L <- round(mean(LimmaOutput_L[, Log2FC]), digits = 3)
      MeanLog2FC_H <- round(mean(LimmaOutput_H[, Log2FC]), digits = 3)

      LimmaVolcano_L <- ggplot2::ggplot(LimmaOutput_L, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) +
        ggplot2::geom_vline(xintercept = MeanLog2FC_L, linetype = "dashed", colour = "#000") +
        ggplot2::annotate("text", x = min(LimmaOutput_L[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_L), size = 5) +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Light Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      ImputedLimmaVolcano_L <- ggplot2::ggplot(LimmaOutput_L, ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        ggplot2::scale_colour_manual(values = c(Measured = "#000", Imputed = "#C0C"), guide = "none") +
        Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) +
        ggplot2::geom_vline(xintercept = MeanLog2FC_L, linetype = "dashed", colour = "#000") +
        ggplot2::annotate("text", x = min(LimmaOutput_L[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_L), size = 5) +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Imputed Light Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      LimmaVolcano_H <- ggplot2::ggplot(LimmaOutput_H, ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) +
        ggplot2::geom_vline(xintercept = MeanLog2FC_H, linetype = "dashed", colour = "#000") +
        ggplot2::annotate("text", x = min(LimmaOutput_H[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_H), size = 5) +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Heavy Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      ImputedLimmaVolcano_H <- ggplot2::ggplot(LimmaOutput_H, ggplot2::aes(x = Log2FC, y = -log10(P.Value), colour = Imputed)) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        ggplot2::scale_colour_manual(values = c(Measured = "#000", Imputed = "#C0C"), guide = "none") +
        Proteopedia::Add_NotSigBox() + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(Gene %in% head(Gene, 250), as.character(Gene), ""))) +
        ggplot2::geom_vline(xintercept = MeanLog2FC_H, linetype = "dashed", colour = "#000") +
        ggplot2::annotate("text", x = min(LimmaOutput_H[, Log2FC]), y = 0, label = paste0("Mean Log2FC\n", MeanLog2FC_H), size = 5) +
        Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Imputed Heavy Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      pdf("Limma_Volcanoes.pdf", height = 10, width = 14)
      print(LimmaVolcano_L)
      print(ImputedLimmaVolcano_L)
      print(LimmaVolcano_H)
      print(ImputedLimmaVolcano_H)
      Proteopedia::Reset_Dev()

      # Light Interactive Plot
      InteractiveData_L <- data.table::data.table(Protein = LimmaOutput_L$Gene, Log2FC = round(LimmaOutput_L$Log2FC, digits = 2),
                                                  PValue = LimmaOutput_L$P.Value, Significance = factor(LimmaOutput_L$Significance),
                                                  Scien_PValue = formatC(LimmaOutput_L$P.Value, format = "e", digits = 2),
                                                  GeneGroup = LimmaOutput_L$GeneGroup, URL = LimmaOutput_L$URL)
      pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData_L,
                                                                                                                                                                              "scatter", highcharter::hcaes(x = Log2FC,
                                                                                                                                                                                                            y = -log10(PValue), group = Significance)),
                                                                                                                                                          zoomType = "xy"), title = list(text = paste0(ExpGroupsName, " vs ", CtlGroupsName, " Log2 Fold-Change")),
                                                                                                                                    lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000",
                                                                                                                                    tickColor = "#000", gridLineWidth = 0), title = list(text = "-Log10 P-Value"),
                                                                                                              lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000",
                                                                                                              tickColor = "#000", gridLineWidth = 0), headerFormat = "",
                                                                                        pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}"),
                                                                scatter = list(marker = list(radius = 3),
                                                                               states = list(hover = list(enabled = T),
                                                                                             inactive = list(enabled = F)), point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))),
                                    c("#999", "#800", "#03F"))
      htmlwidgets::saveWidget(pHC, "LightInteractiveVolcanoPlot.html")

      # Heavy Interactive Plot
      InteractiveData_H <- data.table::data.table(Protein = LimmaOutput_H$Gene, Log2FC = round(LimmaOutput_H$Log2FC, digits = 2),
                                                  PValue = LimmaOutput_H$P.Value, Significance = factor(LimmaOutput_H$Significance),
                                                  Scien_PValue = formatC(LimmaOutput_H$P.Value, format = "e", digits = 2),
                                                  GeneGroup = LimmaOutput_H$GeneGroup, URL = LimmaOutput_H$URL)
      pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData_H,
                                                                                                                                                                              "scatter", highcharter::hcaes(x = Log2FC,
                                                                                                                                                                                                            y = -log10(PValue), group = Significance)),
                                                                                                                                                          zoomType = "xy"), title = list(text = paste0(ExpGroupsName, " vs ", CtlGroupsName, " Log2 Fold-Change")),
                                                                                                                                    lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000",
                                                                                                                                    tickColor = "#000", gridLineWidth = 0), title = list(text = "-Log10 P-Value"),
                                                                                                              lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000",
                                                                                                              tickColor = "#000", gridLineWidth = 0), headerFormat = "",
                                                                                        pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}"),
                                                                scatter = list(marker = list(radius = 3),
                                                                               states = list(hover = list(enabled = T),
                                                                                             inactive = list(enabled = F)), point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))),
                                    c("#999", "#800", "#03F"))
      htmlwidgets::saveWidget(pHC, "HeavyInteractiveVolcanoPlot.html")
    }
    message("SILAC Channels: Exporting QC Plot")
    {
      pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
      suppressWarnings(
        print(patchwork::free(PCAPlot) + CountsBar + patchwork::free(UpsetPlot) + NormPlot + ImpPlot + patchwork::free(MeanVarPlot, type = "label") +
                patchwork::plot_layout(design = "AAAABBBB\nCCCCDDDD\nEEEFFFFF") + patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "", "D", "E", "F", "G"))))
      )
      Proteopedia::Reset_Dev()
    }
  }
  message("Analysing Channel Ratio For Turnover Measures")
  {
    setwd(InputDirectory)
    RatioData <- data.table::fread(InputFile)[Sample %in% Metadata$Sample]
    setwd(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
    RatioData[, Ratio := 2^Log2HLRatio]
    RatioData <- RatioData[!is.na(Ratio)]
    if(dir.exists("Ratio_Analysis")){unlink("Ratio_Analysis", recursive = T)}
    dir.create("Ratio_Analysis", showWarnings = T)
    setwd("Ratio_Analysis")
    message("Channel Ratio: Performing PCA")
    {
      PCAData <- RatioData |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Ratio") |> tidyr::drop_na() |>
        data.frame(row.names = "ProteinGroup") |> t()
      PCAData <- PCAData[ , which(apply(PCAData, 2, var) != 0)] |> stats::prcomp(scale. = T)

      PCASummary <- summary(PCAData)$importance
      PCAData <- data.table::merge.data.table(data.table::data.table(PCAData$x, keep.rownames = "Sample"), Metadata)
      PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]

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
      FilteringData <- RatioData[, .(N_Samples = .N, Min_Precursors_L = min(N_precursors_L), Min_Precursors_H = min(N_precursors_H), N_conditions = data.table::uniqueN(Condition)), ProteinGroup]
      Retained1 <- FilteringData[N_Samples == data.table::uniqueN(RatioData$Sample) & Min_Precursors_L >= MinPrecursors & Min_Precursors_H >= MinPrecursors]
      Retained2 <- FilteringData[N_Samples == data.table::uniqueN(RatioData$Sample) - 1 & Min_Precursors_L >= (MinPrecursors + 1) & Min_Precursors_H >= (MinPrecursors + 1)]
      Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(RatioData$Sample)/2) & Min_Precursors_L >= (MinPrecursors + 1) & Min_Precursors_H >= (MinPrecursors + 1) & N_conditions >= 1]
      RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
      FilteredProteins <- RatioData[ProteinGroup %!in% RetainedProteins]
      RatioData <- RatioData[ProteinGroup %in% RetainedProteins]
      data.table::fwrite(FilteredProteins, file = "Filtered_Proteins.csv")
      data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")

      FilteringData <- rbind(RatioData[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")], FilteredProteins[, .(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")])[, .N, .(Sample, Condition, Replicate, Inclusion)]

      CountsBar <- ggplot2::ggplot(data.table::setorder(FilteringData, Sample), ggplot2::aes(x = gsub("_", " ", Sample), y = N, fill = Condition, alpha = Inclusion)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::scale_alpha_manual(values = c(Excluded = 0.4, Retained = 1), guide = "none") +
        ggplot2::geom_text(ggplot2::aes(label = N, colour = Inclusion), size = 6, position = ggplot2::position_stack(), vjust = 1.5) +
        ggplot2::scale_colour_manual(values = c("Retained" = "#FFF", "Excluded" = "#000"), guide = "none") + ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) + ggplot2::labs(x = NULL, y = "No. Proteins", fill = NULL) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = ggplot2::element_blank(),
                       strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0, "line"))

      UpsetPlot <- ggplot2::ggplot(RatioData[, .(Sample = list(gsub("_", " ", Sample))), by = ProteinGroup], ggplot2::aes(x = Sample)) +
        ggplot2::geom_bar() + ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) +
        ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = gsub("_", " ", RatioData[order(Condition), unique(Sample)])) +
        ggplot2::labs(x = NULL, y = "Post-Filtering Count") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))

      data.table::fwrite(RatioData, "FilteredRatioData.csv")
    }
    message("Channel Ratio: Plotting Experimental vs Control Ratio Data")
    {
      RatioDataWide <- RatioData[, .(MeanRatio = mean(Ratio)), .(ProteinGroup, Condition)] |>
        data.table::dcast(ProteinGroup ~ Condition, value.var = "MeanRatio") |>
        Proteopedia::Add_ProteinInfo(ProteinInfoFile = paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      RatioDataWide |> data.table::setnames(c(colnames(RatioDataWide)[grepl(CtlGroupsName, colnames(RatioDataWide))],
                                              colnames(RatioDataWide)[grepl(ExpGroupsName, colnames(RatioDataWide))]), c("CtlGroups", "ExpGroups"))
      RatioDataWide[, Histone := data.table::fifelse(grepl("^H\\d.*", Gene), T, F)]

      RatioDataLong <- RatioDataWide |> data.table::melt.data.table(id.vars = c("ProteinGroup", "Gene", "Histone"),
                                                                    measure.vars = c("CtlGroups", "ExpGroups"), value.name = "SILACRatio",
                                                                    variable.name = "Condition")
      RatioDataLong[, Comparative := data.table::fifelse(Condition == "CtlGroups", CtlGroupsName, ExpGroupsName)]

      DistributionPval <- wilcox.test(RatioDataLong[Comparative == CtlGroupsName, SILACRatio], RatioDataLong[Comparative == ExpGroupsName, SILACRatio])$p.value
      DistributionPval <- data.table::fifelse(DistributionPval < 0.001, "p < 0.001 ***", data.table::fifelse(DistributionPval < 0.01, "p < 0.01 **",
                                                                                                             data.table::fifelse(DistributionPval < 0.05, "p < 0.05 *", "p > 0.05")))

      pdf(paste0(ExpGroupsName, "_", CtlGroupsName, "_RatioComp.pdf"), height = 10, width = 14)
      print(RatioDataWide |> ggplot2::ggplot(ggplot2::aes(x = CtlGroups, y = ExpGroups)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_XYLine("#999") + Proteopedia::Add_Pearsons(Subgroups = F) + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
              ggplot2::annotate("text", x = max(RatioDataWide$CtlGroups, na.rm = T)*0.9, y = min(RatioDataWide$ExpGroups, na.rm = T)*0.9, label = "Histones", colour = "#F0F") +
              ggrepel::geom_text_repel(ggplot2::aes(label = Gene)) + ggplot2::geom_point(data = RatioDataWide[Histone ==T], stroke = NA, colour = "#F0F") +
              ggplot2::labs(x = paste0("Heavy:Light Protein LFQ Ratio (",  CtlGroupsName, ")"), y = paste0("Heavy:Light Protein LFQ Ratio (",  ExpGroupsName, ")")))
      print(RatioDataLong |> ggplot2::ggplot(ggplot2::aes(x = SILACRatio, fill = factor(Comparative, levels = c(ExpGroupsName, CtlGroupsName)))) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, name = "") + ggplot2::labs(x = "Heavy:Light Protein LFQ Ratio", y = "Protein Frequency") +
              ggplot2::annotate("text", x = median(RatioDataLong$SILACRatio, na.rm = T), y = Proteopedia::Calculate_DensityPeak(RatioDataLong$SILACRatio)*1.1, label = DistributionPval) +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)))
      Proteopedia::Reset_Dev()

      InteractiveData <- data.table::data.table(Protein = RatioDataWide$Gene, Log2FC_CtlGroup = round(RatioDataWide$CtlGroup, digits = 3),
                                                Log2FC_ExpGroup = round(RatioDataWide$ExpGroup, digits = 3), GeneGroup = RatioDataWide$GeneGroup, URL = RatioDataWide$URL)
      pHC <- highcharter::hc_colors(highcharter::hc_plotOptions(highcharter::hc_tooltip(highcharter::hc_yAxis(
        highcharter::hc_xAxis(highcharter::hc_chart(highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC_CtlGroup, y = Log2FC_ExpGroup)), zoomType = "xy"),
                              title = list(text = paste0("Log2 Heavy:Light Protein LFQ Ratio (", CtlGroupsName, ")")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0),
        title = list(text = paste0("Log2 Heavy:Light Protein LFQ Ratio (", ExpGroupsName, ")")), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0), headerFormat = "",
        pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC (Ctl): {point.Log2FC_CtlGroup:.3f}<br>Log2FC (Exp): {point.Log2FC_ExpGroup:.3f}"),
        scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = T), inactive = list(enabled = F)),
                       point = list(events = list(click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))),
        c("#999", "#800", "#03F")) |> highcharter::hc_title(text = paste0(ExpGroupsName, " vs. ", CtlGroupsName))
      htmlwidgets::saveWidget(pHC, "InteractiveCorrPlot.html")
    }
    message("Channel Ratio: Performing Median Normalisation")
    {
      RatioData[, `:=`(Raw = Ratio, Normalised = Ratio - median(Ratio, na.rm = T) + median(RatioData$Ratio, na.rm = T)), Sample]
      RatioData <- RatioData[Normalised > 0]
      NormPlot <- ggplot2::ggplot(data.table::melt.data.table(RatioData, measure.vars = c("Raw", "Normalised")),
                                           ggplot2::aes(x = Condition, colour = Condition, y = value, group = Sample)) +
        ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c(LFQ_Ratio = "Pre", SILACRatio_Norm = "Post"))) +
        ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
        ggplot2::labs(x = NULL, y = "Heavy:Light Protein LFQ Ratio") + ggplot2::theme(strip.background = ggplot2::element_blank())
    }
    message("Channel Ratio: Imputing NA Values")
    {
      RatioData <- RatioData[, .(ProteinGroup, Sample, Normalised)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Normalised")
      RatioImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(RatioData, rownames = "ProteinGroup"), q = ImputationQ,
                                                                    tune.sigma = ImputationSigma), keep.rownames = "ProteinGroup")
      RatioAll <- data.table::merge.data.table(data.table::melt.data.table(RatioData, id.vars = "ProteinGroup", value.name = "MeasuredRatio", variable.name = "Sample"),
                                               data.table::melt.data.table(RatioImp, id.vars = "ProteinGroup", value.name = "ImputedRatio", variable.name = "Sample"))
      RatioAll[, `:=`(Data, data.table::fifelse(is.na(MeasuredRatio), "Imputed", "Measured"))]

      ImpPlot <- ggplot2::ggplot(RatioAll, ggplot2::aes(x = log2(ImputedRatio), fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_y_continuous(expand = 0) +
        ggplot2::scale_fill_manual(values = c(Measured = "#000", Imputed = "#C0C")) + ggplot2::labs(x = expression("Log"[2]~"Heavy:Light Ratio"), y = "Density", fill = NULL) +
        ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))

      ImputedNAs <- data.table::melt.data.table(RatioData[ProteinGroup %in% Retained3$ProteinGroup], id.vars = "ProteinGroup",
                                                variable.name = "Sample")[is.na(value)][, `:=`(value, NULL)]
      ImputedNAs <- dplyr::summarise(dplyr::group_by(ImputedNAs, ProteinGroup), vector = paste(Sample, collapse = ", "))
      data.table::fwrite(ImputedNAs, "Imputed_Ratios.csv")
      RatioAll[, `:=`(Ratio, data.table::fifelse(is.na(MeasuredRatio), ImputedRatio, MeasuredRatio))]
      RatioAll <- RatioAll[Ratio > 0, .(ProteinGroup, Sample, Ratio)]
    }
    message("Channel Ratio: Performing Paired T-Testing")
    {
      RatioAll <- RatioAll |> data.table::merge.data.table(Metadata[, .(Comparative, Sample)], by = "Sample")

      RatioTtest <- RatioAll[, .(Log2MeanRatio = log2(mean(Ratio, na.rm = T)), CV = Proteopedia::Calculate_CV(Ratio), .N), .(Comparative, ProteinGroup)] |>
        data.table::dcast(ProteinGroup ~ Comparative, value.var = c("Log2MeanRatio", "CV", "N"))

      CtlColIndex <- colnames(RatioTtest)[which(grepl("Log2MeanRatio_", colnames(RatioTtest)) & grepl(CtlGroupsName, colnames(RatioTtest)))]
      ExpColIndex <-  colnames(RatioTtest)[which(grepl("Log2MeanRatio_", colnames(RatioTtest)) & grepl(ExpGroupsName, colnames(RatioTtest)))]
      RatioTtest$Log2FC <- RatioTtest[, get(ExpColIndex)] - RatioTtest[, get(CtlColIndex)]

      Ttest_Output <- RatioAll[, P.Value := stats::t.test(Ratio ~ Comparative)$p.value, ProteinGroup]
      Ttest_Output <- Proteopedia::Add_ProteinInfo(Ttest_Output, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      RatioTtest <- data.table::merge.data.table(data.table::data.table(RatioTtest), Ttest_Output)
      RatioTtest[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
      data.table::fwrite(RatioTtest |> dplyr::distinct(), "Paired_T-Test_Output.csv")
    }
    message("Channel Ratio: Fitting Linear Model")
    {
      ModelDesign <- stats::model.matrix(~0 + Comparative, data = Metadata)
      colnames(ModelDesign) <- gsub("Comparative", "", colnames(ModelDesign))
      rownames(ModelDesign) <- Metadata$Sample
      ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroupsName, colnames(ModelDesign))), which(grepl(ExpGroupsName, colnames(ModelDesign)))))]

      ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
      ContrastMatrix[,1] <- c(-1,1)

      RatioAll[, Log2Ratio := log2(Ratio)]

      LimmaInput <- RatioAll |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2Ratio") |>
        data.table::setcolorder(c("ProteinGroup", Metadata$Sample))

      suppressMessages(ModelFit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, ModelDesign), ContrastMatrix)))

      if(!is.finite(ModelFit$df.prior)){message("Warning: Limma Prior is Infinite")}

      MeanVarData <- data.table::data.table(ModelFit$genes, "Mean" = ModelFit$Amean, "Variance" = sqrt(ModelFit$sigma))
      MeanVarData[, Data := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured")]

      MeanVarPlot <- MeanVarData |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) +
        ggplot2::geom_point(colour = "#000") +  ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C") +
        ggplot2::labs(x = expression("Mean Log"[2]~"Ratio"), y = "Variance", colour = NULL)

      LimmaOutput <- data.table::data.table(limma::topTable(ModelFit, coef=1, adjust.method = "BH", n=Inf)) |> data.table::setnames("logFC", "Log2FC")
      LimmaOutput <- LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = T)]
      LimmaOutput <- Proteopedia::Add_ProteinInfo(LimmaOutput, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
      LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = T)
      LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
      LimmaOutput <- Proteopedia::Separate_Isoforms(LimmaOutput)
      data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
    }
    message("Channel Ratio: Generating Volcano Plots")
    {
      LimmaVolcano <- LimmaOutput |> dplyr::arrange(desc(abs(t))) |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
        ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
        ggrepel::geom_text_repel(ggplot2::aes(label= data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
        ggplot2::geom_vline(xintercept = mean(LimmaOutput$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
        ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Ratio"), y = expression("-Log"[10]~"P-Value")) + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

      TtestVolcano <- unique(SpectraTtest[, .(ProteinGroup, Log2FC, P.Value, Gene)]) |> dplyr::arrange(desc(abs(Log2FC)*-log10(P.Value))) |>
        ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
        ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
        ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
        ggplot2::geom_vline(xintercept = mean(SpectraTtest$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
        ggplot2::labs(x = expression("Log"[2]~"FC in Heavy:Light Ratio"), y = expression("-Log"[10]~"P-Value")) + ggplot2::ggtitle("T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())

      pdf("VolcanoPlots.pdf", height = 10, width = 14)
      print(LimmaVolcano)
      print(LimmaVolcano + ggplot2::geom_point(data = LimmaOutput[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed Limma"))
      print(TtestVolcano)
      print(TtestVolcano + ggplot2::geom_point(data = RatioTtest[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed T-Test"))
      Proteopedia::Reset_Dev()
    }
    message("Channel Ratio: Exporting QC Plot")
    {
      pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
      suppressWarnings(
        print(patchwork::free(PCAPlot) + CountsBar + patchwork::free(UpsetPlot) + NormPlot + ImpPlot + patchwork::free(MeanVarPlot, type = "label") +
                patchwork::plot_layout(design = "AAAABBBB\nCCCCDDDD\nEEEEFFFF") + patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E", "F"))))
      )
      Proteopedia::Reset_Dev()
    }
    message("Channel Ratio: Exporting HTML Volcano Plot")
    {
      InteractiveData <- data.table::data.table("Protein" = LimmaOutput$Gene, "Log2FC" = round(LimmaOutput$Log2FC, digits = 2),
                                                "PValue" = LimmaOutput$P.Value, "Significance" = factor(LimmaOutput$Significance),
                                                "Scien_PValue" = formatC(LimmaOutput$P.Value, format = "e", digits = 2),
                                                "GeneGroup" = LimmaOutput$GeneGroup, "URL" = LimmaOutput$URL)

      HighCharterVolcano <- highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC, y = -log10(PValue), group = Significance)) |>
        highcharter::hc_chart(zoomType = "xy") |>
        highcharter::hc_xAxis(title = list(text = paste0(unique(Metadata$Condition)[1]," vs ", unique(Metadata$Condition)[2]," Log2 Fold-Change in Ratio")),
                              lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
        highcharter::hc_yAxis(title = list(text = "-Log10 P-Value"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
        highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}") |>
        highcharter::hc_plotOptions(scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = T), inactive = list(enabled = F)),
                                                   point = list(events = list( click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))) |>
        highcharter::hc_colors(c("#999", "#800","#03F"))
      htmlwidgets::saveWidget(HighCharterVolcano, "InteractiveVolcanoPlot.html")
    }
  }
  setwd(paste0(InputDirectory,"/",ExpGroupsName,"_vs_",CtlGroupsName,"_Output"))
  data.table::fwrite(data.table::data.table("Experimental Condition(s)" = paste(ExpGroups, collapse = ", "), "Experimental Name" = paste0(ExpGroupsName),
                                            "Control Condition(s)" = paste(CtlGroups, collapse = ", "), "Control Name" = paste0(CtlGroupsName),
                                            "Min_Precursors" = paste0(MinPrecursors), "Imputation Q-Value" = ImputationQ,
                                            "Imputation Sigma" = ImputationSigma), "Analysis_Parameters.csv")
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Map_ProteinSubsets <- function(InputDirectory, SILAC = F, SubsetColour = "#F00"){
  start.time <- Sys.time()
  set.seed(123)
  for(Mapping in c("Total_Analysis", "Ratio_Analysis", "Heavy_Channel_Analysis", "Light_Channel_Analysis")){
    message(gsub("_", " ", Mapping), ": Loading Limma File")
    {
      if(SILAC){
        setwd(paste0(InputDirectory, "/", gsub(".*_(Channel_Analysis)", "\\1", Mapping)))
      } else {
        setwd(InputDirectory)
      }
      LimmaData <- data.table::fread("Limma_Output.csv")
      MappingDirectory <- getwd()
      LimmaData |> data.table::setnames(
        c(colnames(LimmaData)[grepl("protein.*group", ignore.case = T, colnames(LimmaData))],
          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & !grepl("adj", ignore.case = T, colnames(LimmaData))],
          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & grepl("adj", ignore.case = T, colnames(LimmaData))],
          colnames(LimmaData)[grepl("Gene", ignore.case = T, colnames(LimmaData)) & !grepl("group", ignore.case = T, colnames(LimmaData))]),
        c("ProteinGroup", "P.Value", "adj.P.Val", "Gene"))

      if(length(LimmaData$ProteinGroup[grepl("\\-", LimmaData$ProteinGroup)]) > 0){LimmaData <- Proteopedia::Separate_Isoforms(LimmaData)}

      if(length(colnames(LimmaData)[grepl("Channel", colnames(LimmaData))]) < 1){
        LimmaData$Channel <- gsub("_Analysis", "", Mapping)
      }
      LimmaData <- LimmaData |> data.table::merge.data.table(Proteopedia::Proteopedia, all.x = T)
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
    message(gsub("_", " ", Mapping), ": Biochemical Trends")
    {
      for(ColIndex in c(which(colnames(LimmaData) == "Length"), which(colnames(LimmaData) == "Length"):ncol(LimmaData))){
        if(is.numeric(LimmaData[, get(colnames(LimmaData)[ColIndex])]) & !grepl("N_", colnames(LimmaData)[ColIndex])){
          message(paste0(gsub("_", " ", Mapping), ": ", colnames(LimmaData)[ColIndex]), " Trend")
          SubsetData <- LimmaData[, .(ProteinGroup, Channel, Sequence, Length, Log2FC, get(colnames(LimmaData)[ColIndex]))]
          SubsetData |> data.table::setnames("V6", "Subset")

          SubsetStats <- data.table::data.table(Channel = unique(LimmaData$Channel), MeanSubsetVal = 0, MinLog2FC = 0, P.Value = 0)
          for(RowIndex in nrow(SubsetStats)){
            SubsetStats$Subset[RowIndex] <- quantile(SubsetData[Channel == SubsetStats$Channel[RowIndex], Subset], 0.75, na.rm = T)*1.2
            SubsetStats$Log2FC[RowIndex] <- min(SubsetData[Channel == SubsetStats$Channel[RowIndex], Log2FC], na.rm = T)*0.93
            SubsetStats$P.Value[RowIndex] <- summary(stats::lm(Log2FC ~ Subset, data = SubsetData[Channel == SubsetStats$Channel[RowIndex]]))$coefficients[2,4]
          }
          SubsetStats[, P.Value := data.table::fifelse(P.Value < 0.01, formatC(P.Value, format = "e", digits = 2), as.character(round(P.Value, digits = 2)))]

          pdf(paste0(colnames(LimmaData)[ColIndex],"_Trend.pdf"), width = 12, height = 10)
          print(SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Subset, y = Log2FC)) +
                  ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Add_Pearsons() + ggplot2::facet_wrap(~Channel) +
                  ggplot2::geom_text(data = SubsetStats, ggplot2::aes(label = paste0("p-value = ", P.Value)), size = 6) +
                  ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2]~"FC in Protein Abundance")) +
                  ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities())
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(gsub("_", " ", Mapping), ": Cellular Trends")
    {
      for(ColIndex in which(colnames(LimmaData) == "ER"):ncol(LimmaData)){
        if(is.logical(LimmaData[, get(colnames(LimmaData)[ColIndex])])){
          message(paste0(gsub("_", " ", Mapping), ": ", colnames(LimmaData)[ColIndex]), " Proteins")
          SubsetData <- LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, Channel, get(colnames(LimmaData)[ColIndex]))]
          SubsetData |> data.table::setnames("V6", "Subset")
          if(nrow(SubsetData[Subset == T]) > 0){
            SubsetStats <- data.table::data.table(Channel = unique(LimmaData$Channel), MeanLog2FC = 0, P.Value = 0)
            for(RowIndex in nrow(SubsetStats)){
              SubsetStats$Log2FC[RowIndex] <- mean(SubsetData[Channel == SubsetStats$Channel[RowIndex] & Subset == T, Log2FC], na.rm = T)
              SubsetStats$P.Value[RowIndex] <- 0.00002
              SubsetStats$StatLabel[RowIndex] <- wilcox.test(SubsetData[Channel == SubsetStats$Channel[RowIndex] & Subset == T, Log2FC], SubsetData[Channel == SubsetStats$Channel[RowIndex] & Subset == F, Log2FC])$p.value
            }
            SubsetStats[, StatLabel := data.table::fifelse(StatLabel < 0.01, formatC(StatLabel, format = "e", digits = 2), as.character(round(StatLabel, digits = 2)))]

            pdf(paste0(colnames(LimmaData)[ColIndex],"_Volcano.pdf"))
            print(SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
                    ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) +
                    Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() + ggplot2::facet_wrap(~Channel))
            print(SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::facet_wrap(~Channel) +
                    ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3) + ggplot2::geom_point(data = SubsetData[Subset == T], colour = SubsetColour) +
                    Proteopedia::Add_NotSigBox() + ggplot2::geom_rug(alpha = data.table::fifelse(SubsetData[, Gene] %in% SubsetData[Subset == T, Gene], 1, 0), colour = SubsetColour, sides = "tr") +
                    ggplot2::geom_label(data = SubsetStats, ggplot2::aes(label = paste0("p = ", StatLabel)), colour = SubsetColour, size = 6) +
                    Proteopedia::Add_AbundanceAxes())
            Proteopedia::Reset_Dev()
          }
        }
      }
    }
    message(gsub("_", " ", Mapping), ": Degradation Profiles")
    {
      DegProfileSummary <- LimmaData[, .N, .(Deg_Profile, Channel)]

      pdf("DegradataionProfileBoxplot.pdf", width = 12, height = 10)
      print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Deg_Profile, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) +
              ggpubr::geom_signif(comparison = list(c("NED", "ED"), c("NED", "UN"), c("ED", "UN")), y_position = c(1, 1.25, 1.5), tip_length = 0) +
              ggplot2::labs(x = "Degradation Profile", y = expression("Log"[2]~"FC in Protein Abundance")) + ggplot2::facet_wrap(~Channel) +
              ggplot2::geom_text(data = DegProfileSummary, y = -1, ggplot2::aes(label = paste0("N = ", N))))
      Proteopedia::Reset_Dev()

      for(Deg_Type in unique(LimmaData[!is.na(Deg_Profile), Deg_Profile])){
        SubsetStats <- data.table::data.table(Channel = unique(LimmaData$Channel), MeanLog2FC = 0, P.Value = 0)
        for(RowIndex in nrow(SubsetStats)){
          SubsetStats$Log2FC[RowIndex] <- mean(LimmaData[Channel == SubsetStats$Channel[RowIndex] & Deg_Profile == Deg_Type, Log2FC], na.rm = T)
          SubsetStats$P.Value[RowIndex] <- 0.00002
          SubsetStats$StatLabel[RowIndex] <- wilcox.test(LimmaData[Channel == SubsetStats$Channel[RowIndex] & Deg_Profile == Deg_Type, Log2FC], SubsetData[Channel == SubsetStats$Channel[RowIndex] & Subset == F, Log2FC])$p.value
        }
        SubsetStats[, StatLabel := data.table::fifelse(StatLabel < 0.01, formatC(StatLabel, format = "e", digits = 2), as.character(round(StatLabel, digits = 2)))]

        pdf(paste0(Deg_Type, "_DegProfile_Volcano.pdf"))
        print(LimmaData[Deg_Profile == Deg_Type] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
                ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) +
                Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() + ggplot2::facet_wrap(~Channel))
        print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::facet_wrap(~Channel) +
                ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3) + ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
                Proteopedia::Add_NotSigBox() + ggplot2::geom_rug(alpha = data.table::fifelse(LimmaData[, Gene] %in% LimmaData[Deg_Profile == Deg_Type, Gene], 1, 0), colour = SubsetColour, sides = "tr") +
                ggplot2::geom_label(data = SubsetStats, ggplot2::aes(label = paste0("p = ", StatLabel)), colour = SubsetColour, size = 6) +
                Proteopedia::Add_AbundanceAxes())
        Proteopedia::Reset_Dev()
      }
    }
    message(gsub("_", " ", Mapping), ": Chromosome-Based Data")
    {
      ChromosomeSummary <- data.table::data.table()
      for(ChannelIndex in unique(LimmaData$Channel)){
        ChromosomeSummary <- ChromosomeSummary <- rbind(Calculate_WilcoxonByVar(InputData = LimmaData[Channel == ChannelIndex], Category = "Chromosome", Measure = "Log2FC"))
      }
      ChromosomeSummary <- ChromosomeSummary |> data.table::merge.data.table(LimmaData[, .(MeanLog2FC = mean(Log2FC, na.rm = T), Mapped = .N,
                                                                                           WhiskerTop = (stats::quantile(Log2FC, 0.75, na.rm = T) + stats::IQR(Log2FC, na.rm = T))*1.1,
                                                                                           WhiskerBottom = (stats::quantile(Log2FC, 0.25, na.rm = T) - stats::IQR(Log2FC, na.rm = T))*1.1), .(Chromosome, Channel)]) |>
        data.table::merge.data.table(Proteopedia::Proteopedia[, .(Total = .N), Chromosome])
      ChromosomeSummary[, `:=`(Coverage = Mapped/Total, Chromosome = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")))]
      data.table::fwrite(ChromosomeSummary |> data.table::setorder(Chromosome), "Chromosomal_Summary.csv")

      LimmaData[, Chromosome := factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))]
      LimmaData <- LimmaData |> dplyr::arrange(Chromosome, MedianLociStart)
      suppressWarnings(LimmaData[, OrderID := as.numeric(rownames(LimmaData))])

      ChromosomeDotData <- unique(LimmaData[Chromosome != "Unmapped", .(ProteinGroup, Log2FC, Gene, Chromosome, MedianLociStart, Channel)])
      ChromosomeDotData[, MedianLociStart := as.numeric(MedianLociStart)]
      ChromosomeDotData <- ChromosomeDotData |> data.table::setorder(Chromosome, MedianLociStart) |>
        data.table::merge.data.table(ChromosomeDotData[, .(MaxChrPos = .N), .(Chromosome, Channel)])
      ChromosomeDotData[, ChrPos := seq_len(.N), .(Chromosome, Channel)]
      ChromosomeDotData[, RelChrPos := ChrPos/MaxChrPos]

      ChromosomeBorders <- ChromosomeDotData[, .(N_Proteins = .N), .(Chromosome, Channel)] |> tibble::rowid_to_column("ChrIndex")
      suppressWarnings(ChromosomeBorders[, Lower := cumsum(N_Proteins)-(N_Proteins-1)])
      ChromosomeBorders[, Midpoint := (cumsum(N_Proteins)+Lower)/2]
      ChromosomeBorders[, Upper := cumsum(N_Proteins)]
      ChromosomeBorders[, Lower := (Lower/max(ChromosomeBorders$Upper)) + ChrIndex]
      ChromosomeBorders[, Midpoint := (Midpoint/max(ChromosomeBorders$Upper)) + ChrIndex]
      ChromosomeBorders[, Upper := (Upper/max(ChromosomeBorders$Upper)) + ChrIndex]
      ChromosomeDotData <- ChromosomeDotData |> data.table::merge.data.table(ChromosomeBorders)
      ChromosomeDotData[, RelPos := RelChrPos + ChrIndex]

      ChrGroupingData <- LimmaData[Chromosome != "Unmapped", .N, .(Significance, Chromosome, Channel)] |> data.table::merge.data.table(LimmaData[Chromosome != "Unmapped", .(Total_N = .N), Chromosome])
      ChrGroupingData[, `:=`(Prop = N/Total_N, Significance = data.table::fifelse(grepl("Sig.", Significance), Significance, "None"))]

      pdf("ChromosomePlots.pdf", width = 18, height = 8)
      print(ChromosomeDotData |> ggplot2::ggplot(ggplot2::aes(x = RelPos, y = Log2FC, colour = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")))) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_manual(values = rep(Proteopedia::NiceColourPalette, 2), guide = "none") + ggplot2::facet_wrap(~Channel, ncol = 1) +
              ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) + ggplot2::scale_y_continuous(limits = c(-max(ceiling(abs(ChromosomeDotData$Log2FC))), max(ceiling(abs(ChromosomeDotData$Log2FC))))) +
              ggplot2::labs(x = "Chromosome", y = expression("Log"[2]~ "FC in Protein Abundance")) + ggplot2::geom_vline(xintercept = ChromosomeBorders[Chromosome != 1, ChrIndex], colour = "#000") +
              ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = ChrIndex+0.5, y = I(0.05), label = data.table::fifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))), angle = 90) +
              ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()))
      print(ChromosomeSummary[!is.na(Chromosome)] |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Coverage*100, fill = Coverage*100)) +
              ggplot2::geom_bar(stat = "identity") + ggplot2::scale_fill_viridis_c(guide = "none") + ggplot2::labs(x = "Chromosome", y = "Protein Coverage (%)") + ggplot2::facet_wrap(~Channel, ncol = 1))
      print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "Y", "M")), y = Log2FC)) + ggplot2::facet_wrap(~Channel, ncol = 1) +
              ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop*1.1, label = SigSymbol)) +
              ggplot2::geom_boxplot(fill = NA, outliers = F) + ggplot2::labs(x = "Chromosome", y = expression("Log"[2]~ "FC in Protein Abundance")))
      print(ChrGroupingData |> ggplot2::ggplot(ggplot2::aes(x = factor(Chromosome, levels = rev(c(seq(1:22), "X", "X/Y", "Y", "M", "Unmapped"))), y = Prop, fill = factor(Significance, levels = c("Sig. Increase", "None", "Sig. Decrease")))) +
              ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = c("Sig. Decrease" = "#02F", "None" = "#FFF", "Sig. Increase" = "#F10"), name = "Fold Change") +
              ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~1-.)) + ggplot2::geom_hline(yintercept = seq(0.1, 0.9, by = 0.1), linetype = "dotted") + ggplot2::facet_wrap(~Channel, ncol = 1) +
              ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::labs(x = "Chromosome", y = "Proportion of Proteins") + ggplot2::coord_flip())
      Proteopedia::Reset_Dev()

      if(dir.exists(paste0(MappingDirectory,"/ChromosomeVolcanos"))){unlink(paste0(MappingDirectory,"/ChromosomeVolcanos"), recursive = T)}
      dir.create(paste0(MappingDirectory,"/ChromosomeVolcanos"), showWarnings = T)
      setwd(paste0(MappingDirectory,"/ChromosomeVolcanos"))
      for(Chr in levels(LimmaData$Chromosome)){
        if(nrow(LimmaData[Chromosome == Chr]) > 0){
          message(paste0(gsub("_", " ", Mapping), ": Chr", Chr, " Proteins"))
          pdf(paste0("Chr", gsub("/", "", Chr),"_Volcano.pdf"), width = 12, height = 8)
          print(LimmaData[Chromosome == Chr] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA) +
                  ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(-log10(P.Value)> -log10(0.05),as.character(Gene),""))) +
                  Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = mean(LimmaData[Chromosome == Chr, Log2FC]), linetype = "dashed", alpha = 0.7) +
                  ggplot2::facet_wrap(~Channel) + Proteopedia::Add_AbundanceAxes())
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(gsub("_", " ", Mapping), ": Protein Complexes")
    {
      setwd(MappingDirectory)
      ComplexPortalSummary <- LimmaData[, .N, .(Experimental_Evidence_ComplexPortal, Channel)]
      ComplexPortalSummary[, Database := "Complex_Portal"] |> data.table::setnames("Experimental_Evidence_ComplexPortal", "Complexed")
      CORUMSummary <- LimmaData[, .N, .(Experimental_Evidence_CORUM, Channel)]
      CORUMSummary[, Database := "CORUM"] |> data.table::setnames("Experimental_Evidence_CORUM", "Complexed")
      ComplexSummary <- ComplexPortalSummary |> rbind(CORUMSummary)

      pdf("ComplexBoxplot.pdf", width = 12, height = 10)
      print(LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_ComplexPortal, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) +
              ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + ggplot2::facet_wrap(~Channel) +
              ggplot2::labs(x = "Complex Portal Evidence", y = expression("Log"[2] ~ "FC in Protein Abundance")) +
              ggplot2::geom_text(data = ComplexSummary[Database == "Complex_Portal"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
              LimmaData |> ggplot2::ggplot(ggplot2::aes(x = Experimental_Evidence_CORUM, y = Log2FC)) + ggplot2::geom_boxplot(outliers = F) +
              ggpubr::geom_signif(comparison = list(c("Yes", "No")), y_position = 1, tip_length = 0) + ggplot2::facet_wrap(~Channel) +
              ggplot2::labs(x = "CORUM Evidence", y = expression("Log"[2] ~ "FC in Protein Abundance")) +
              ggplot2::geom_text(data = ComplexSummary[Database == "CORUM"], ggplot2::aes(y = -1.5, x = Complexed, label = paste0("N = ", N))) +
              ggplot2::theme(axis.title.y = ggplot2::element_blank()))
      Proteopedia::Reset_Dev()
    }
    if(Mapping == "Light_Channel_Analysis"){SILAC = F}
    if(!SILAC){return(Proteopedia::End_Timer(Start = start.time))}
  }
}
########### MS Analysis: Timecourse Analysis Functions ########################################################################################################################################
#' @export
Analyse_TimecourseSILAC_Proteins <- function(InputDirectory, Formula, GenerateDataPlots = F, SameInitialAbundance = T, LightMinSamples = 0.45,
                                             HeavyModel = "NLS", HeavyMinMonotonicity = 0.5, HeavyMinSamples = 0.6, MaxCV = 0.3, MeanCentring = F,
                                             OffsetKloss = F, UseLightKloss = F, ExcludeTimepoints = NULL, ReplicatesUsed = c(1, 2, 3)){
  set.seed(123)
  start.time <- Sys.time()
  Relative_Time <- function(CurrentTime, Times = TimeLevels, Shift = 1){
    Times <- gsub("h", "", Times)
    NextPos = which(Times == unique(CurrentTime)) + Shift
    if(data.table::between(NextPos,1, length(TimeLevels))){return(Times[NextPos])}else{return(NA_character_)}
  }
  message("Defining Comparison Groups")
  {
    ExpGroupsName = gsub("(.*)\\(.*", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula)))
    ExpGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub("(.*)-.*", "\\1", gsub(" ", "", Formula))), "\\+"))
    CtlGroupsName = gsub("(.*)\\(.*", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula)))
    CtlGroups = unlist(stringr::str_split(gsub(".*\\((.*))", "\\1", gsub(".*-(.*)", "\\1", gsub(" ", "", Formula))), "\\+"))
    if(nchar(ExpGroupsName) == 0){ExpGroupsName = "Experiment"}
    if(nchar(CtlGroupsName) == 0){CtlGroupsName = "Control"}
    ComparativeMetadata <- data.table::data.table(Condition = c(CtlGroups, ExpGroups))[, Comparative := data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName,
                                                                                                                            data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName, NA))]
  }
  if(length(c(CtlGroups,ExpGroups)) > 2){
    setwd(InputDirectory)
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Loading Data"))
    {
      ProtLFQsInput <- data.table::fread(list.files(pattern = "DIANN_Output"))[Condition %in% c(CtlGroups, ExpGroups)] |> data.table::setnames("Protein_group", "ProteinGroup", skip_absent = T)
      ProtLFQsInput |> data.table::setnames(colnames(ProtLFQsInput)[grepl("Protein.*group", colnames(ProtLFQsInput), ignore.case = T)], "ProteinGroup")
      ProtLFQsInput <- ProtLFQsInput[gsub(".*(\\d+)$", "\\1", Replicate) %in% ReplicatesUsed]

      ProtLFQsInput <- ProtLFQsInput |> data.table::merge.data.table(ProtLFQsInput[Time == 0, .(BaselineMean = mean(LFQ_L, na.rm = T)), Condition])
      ProtLFQsInput[, RelLFQ_L := LFQ_L/BaselineMean]

      for(ColumnID in c("Conc")){
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

      Metadata <- data.table::fread("Sample_Metadata.csv")[Condition %in% c(CtlGroups, ExpGroups)] |> data.table::setorderv(c("Conc", "Time", "Replicate"))
      SampleLevels <- unique(Metadata$Sample)
      ConcLevels <- unique(Metadata$Conc)
      TimeLevels <- unique(Metadata$Time)
      ConditionLevels <- unique(Metadata$Condition)
      Metadata[, Sample := factor(Sample, levels = SampleLevels)]
      Metadata[, Conc := factor(Conc, levels = ConcLevels)]
      Metadata[, Time := factor(Time, levels = TimeLevels)]
      Metadata[, Condition := factor(Condition, levels = ConditionLevels)]
      Metadata[, Replicate := factor(Replicate)]
      Metadata[, Grouping := data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName,
                                                 data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName, "Exclude"))]
      Metadata[, Grouping := factor(Grouping, levels = c(CtlGroupsName, ExpGroupsName))]

      ProteinWeights <- data.table::fread("Filtered_PrecursorData.csv.gz")[,.(Sum_Intensities  = sum(log2(Precursor.Normalised), na.rm = T)), .(Sample, Channel, ProteinGroup)] |>
        data.table::merge.data.table(Metadata)
      ProteinWeights[, MeanSum := mean(Sum_Intensities, na.rm =T), .(ProteinGroup, Channel)]
      ProteinWeights[, RowID := paste0(ProteinGroup, "_", Channel)]
      ProteinWeights[, Weight := Sum_Intensities/MeanSum]

      # Protein Weights PCA
      All_PCA <- ProteinWeights |> data.table::dcast(RowID ~ Sample, value.var = "Weight") |> tidyr::drop_na() |> data.frame(row.names = "RowID") |> t() |> stats::prcomp(scale. = T)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

      ProteinWeights_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Protein Weights")

      ProteinWeight_Density <- ProteinWeights |> ggplot2::ggplot(ggplot2::aes(x = Weight, fill = Time)) +
        ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap("Channel", nrow = 2, strip.position = "right", scales = "free_x") +
        ggplot2::scale_x_log10() + ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(1, length(Proteopedia::ThermalPalette), length.out = length(TimeLevels))]) +
        ggplot2::labs(title = "Protein Weights for Limma", x = "Protein Weight", y = "No. Proteins", fill = "Time") +
        ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))

      ProteinWeights |> data.table::setnames(c("Channel", colnames(ProteinWeights)[grepl("Protein.*group", colnames(ProteinWeights), ignore.case = T)], colnames(ProteinWeights)[grepl("weight", colnames(ProteinWeights), ignore.case = T)],
                                               colnames(ProteinWeights)[grepl("sum.*int", colnames(ProteinWeights), ignore.case = T)], colnames(ProteinWeights)[grepl("mean.*sum", colnames(ProteinWeights), ignore.case = T)]),
                                             c("Channel", "ProteinGroup", "Weight", "SumIntensities", "MeanSum"), skip_absent = T)
      ProteinWeights <- ProteinWeights[gsub(".*(\\d+)$", "\\1", Replicate) %in% ReplicatesUsed]

      if(length(unique(grepl("Cluster", colnames(ProteinWeights)))) == 1){
        ProteinWeights[, Cluster := gsub("(.*)_\\d+h", "\\1", Condition)]
      }
      if(length(unique(grepl("Time", colnames(ProteinWeights)))) == 1){
        ProteinWeights[, Time := as.numeric(gsub(".*_(\\d+)h", "\\1", Condition))]
      }

      ProteinWeights[, Channel := gsub("L$", "Light", Channel)]
      ProteinWeights[, Channel := gsub("H$", "Heavy", Channel)]

      if (dir.exists(paste0(getwd(), "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"))) {
        unlink(paste0(getwd(), "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"), recursive = T)
      }
      dir.create(paste0(getwd(), "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"), showWarnings = T)
      setwd(paste0(getwd(), "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"))
    }
    # Analyse Light Channel
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Isolating Light Channel Proteins"))
    {
      ProtLFQsInput_L <- ProtLFQsInput |> data.table::copy()
      ProtLFQsInput_L <- ProtLFQsInput_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = "RelLFQ_L") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
      ProtWeights_L <- ProteinWeights[Channel == "Light"] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Weight") |>
        tibble::column_to_rownames("ProteinGroup") |> as.matrix()

      ProteinLFQs_L <- ProtLFQsInput_L |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
      ProteinLFQs_L <- data.table::melt.data.table(ProteinLFQs_L, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |>
        data.table::merge.data.table(Metadata)

      ProteinLFQs_L[, N_Values := sum(is.finite(Abundance)), .(ProteinGroup, Condition, Time)]
      ProteinLFQs_L[, MeanAbundance := mean(Abundance, na.rm = T), ProteinGroup]
      ProteinLFQs_L[, NormAbundance := Abundance - MeanAbundance]
      ProteinLFQs_L[, Sum_Values := sum(is.finite(Abundance)), .(ProteinGroup, Time)]
      ProteinLFQs_L[, Diff_Detected := Sum_Values - N_Values]
      ProteinLFQs_L[, Time := factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels))]

      LightAbunBoxplot <- ProteinLFQs_L |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")),
                                                                        y = Abundance, colour = Condition)) +
        ggplot2::geom_boxplot(outliers = F) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::labs(x = "Time", y = "Light Protein Abundance") + ggplot2::theme(legend.title = ggplot2::element_blank())

      ProteinLFQs_L[, NextTime := Relative_Time(Time, Times = TimeLevels), Time]
      NextTime <- ProteinLFQs_L[,.(NextTimeSamples = mean(N_Values)), .(Time, ProteinGroup, Condition)] |> data.table::setnames("Time","NextTime")
      ProteinLFQs_L <- ProteinLFQs_L |> data.table::merge.data.table(NextTime, by = c("NextTime","ProteinGroup","Condition"), all.x = T)
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Filtering Light Channel Proteins By Variation"))
    {
      ProteinLFQs_L[, N_Quant := .N, .(ProteinGroup, Condition, Time)]
      ProteinLFQs_L[, N_QuantTotal := .N, ProteinGroup]
      ProteinLFQs_L <- ProteinLFQs_L[N_QuantTotal >= LightMinSamples*nrow(Metadata)]
      ProteinLFQs_L[, CV := sd(Abundance)/mean(Abundance), .(ProteinGroup, Condition, Time)]
      ProteinLFQs_L_CV <- ProteinLFQs_L[,.(MeanCV = mean(CV, na.rm =T)), ProteinGroup]
      ProteinLFQs_L <- ProteinLFQs_L[ProteinGroup %in% ProteinLFQs_L_CV[MeanCV <= MaxCV, ProteinGroup]]
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Imputing Light Abundances & Weights"))
    {
      ProteinLFQs_L[NextTimeSamples == 0 & N_Values == 0 & Diff_Detected == 3, Impute :=  T]
      ProteinLFQs_L[is.na(Impute), Impute := F]

      ModelProteins_L <-  imputeLCMD::impute.MinProb(ProtLFQsInput_L, q = 0.01) |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |>
        data.table::data.table() |> data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Sample", value.name = "ImpAbundance")
      ModelProteins_L <- ProteinLFQs_L |> data.table::merge.data.table(ModelProteins_L, by = c("ProteinGroup", "Sample"))
      ModelProteins_L[Impute == T, Abundance := ImpAbundance]
      ModelProteins_L[, ImpAbundance := NULL]
      ModelProteins_L |> data.table::setorder(Condition)

      AbundanceMatrix_L <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Abundance') |>
        tibble::column_to_rownames('ProteinGroup') |> as.matrix()
      AbundanceMatrix_L <- AbundanceMatrix_L[ , SampleLevels]
      AbundanceMatrix_L <- AbundanceMatrix_L[matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) <= LightMinSamples,]

      ProtWeights_L_Imp <- ProtWeights_L[rownames(AbundanceMatrix_L), colnames(AbundanceMatrix_L)]

      ImputationMatrix <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Impute') |>
        tibble::column_to_rownames('ProteinGroup') |> as.matrix()
      ImputationMatrix <- ImputationMatrix[rownames(ProtWeights_L_Imp), colnames(ProtWeights_L_Imp)]

      ProtWeights_L_Imp[ImputationMatrix==T] <- (ProtWeights_L_Imp |> min(na.rm = T))*0.1

      N_NAs <- matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) |> tibble::enframe(name = "ProteinGroup", value = "Prop_NA") |> data.table::data.table()

      # Plot PCA on Processed Light Data
      All_PCA <- AbundanceMatrix_L |> data.frame() |> tidyr::drop_na() |> t() |> stats::prcomp(scale. = T)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

      LightProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Processed Light Data")
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Modelling Light Data"))
    {
      Targets <- Metadata[, .(Sample, Grouping, Time)][ , Time := as.numeric(gsub("h", "", Time))]
      Targets <- Targets[match(ProtWeights_L_Imp |> colnames(), Sample)]
      Time <- Targets$Time
      Grouping <- factor(Targets$Grouping)
      ModelDesign <- model.matrix(~Time + Time:Grouping)
      colnames(ModelDesign) <- gsub(paste0("Grouping", ExpGroupsName), "Grouping_Exp", colnames(ModelDesign))
      colnames(ModelDesign) <- gsub(paste0("Grouping", CtlGroupsName), "Grouping_Ctl", colnames(ModelDesign))

      LimmaOutput <- limma::eBayes(limma::lmFit(Biobase::ExpressionSet(assayData = log(AbundanceMatrix_L)), ModelDesign, method = "robust", weights = ProtWeights_L_Imp))
      Limma_Slopes <- LimmaOutput$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup")
      Limma_Slopes_SD <- LimmaOutput$stdev.unscaled |> data.table::data.table(keep.rownames = T) |> data.table::copy()
      Limma_Slopes_SD <- Limma_Slopes_SD |> data.table::setnames(c("rn", "(Intercept)", "Time", "Grouping_Exp", "Time:Grouping_Exp"),
                                                                 c("ProteinGroup", "(Intercept)_SD", "Time_SD", "Grouping_Exp_SD", "Time:Grouping_Exp_SD"), skip_absent = T)
      Limma_Slopes <- Limma_Slopes |> data.table::merge.data.table(Limma_Slopes_SD, by = "ProteinGroup")
      Limma_Slopes$Fvalue = LimmaOutput$F

      LightMeanVarPlot <- data.table::data.table(Mean = LimmaOutput$Amean, Variance = sqrt(LimmaOutput$sigma)) |>
        ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) + ggplot2::geom_point(stroke = NA) +
        ggplot2::labs(x = "Mean Light Protein Log LFQ Ratio vs. T0 Baseline", y = "Light Protein LFQ Ratio Variance")

      LightModelParameters <- limma::topTable(LimmaOutput, "Time:Grouping_Exp", number = nrow(ProtWeights_L_Imp)) |>
        data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup") |>
        data.table::merge.data.table(N_NAs) |> data.table::merge.data.table(Limma_Slopes) |> data.table::setnames("logFC", "Difference")

      LightModelParameters <- Proteopedia::Add_ProteinInfo(LightModelParameters, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
      LightParameters <- LightModelParameters |> data.table::copy()
      LightParameters[, `:=`(Ctl_Value = -Time, Exp_Value = -`Time:Grouping_Exp` - Time, Difference = -Difference, Parameter = "KlossL")]

      if(OffsetKloss){
        Kloss_Offset <- abs(min(LightParameters[, .(Ctl_Value, Exp_Value)]))*1.01
        LightParameters[, `:=`(Exp_Value = Exp_Value + Kloss_Offset, Ctl_Value = Ctl_Value + Kloss_Offset)]
      } else {
        LightParameters <- LightParameters[Ctl_Value > 0 & Exp_Value > 0]
      }
      LightParameters[, FC := Exp_Value/Ctl_Value]
      LightParameters[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
      LightParameters[, Significance := data.table::fifelse(P.Value < 0.05 & Difference < 0, "Sig. Decrease",
                                                            data.table::fifelse(P.Value < 0.05 & Difference > 0, "Sig. Increase", ""))]

      LightParameters <- LightParameters[, .(ProteinGroup, Gene, P.Value, adj.P.Val, Prop_NA, Parameter, Ctl_Value, Exp_Value, FC, Log2FC, Difference)]
      LightModelledData <- ProtLFQsInput[, ProteinGroup] |> tidyr::crossing(Metadata[, .(Condition, Time)]) |>
        data.table::setnames(c("ProteinGroup", "Condition", "TimeVar")) |> data.table::data.table() |>
        data.table::merge.data.table(LightModelParameters[, .(ProteinGroup, `(Intercept)`, Time, `Time:Grouping_Exp`)] |> data.table::setnames("Time", "TimeCoeff"))
      LightModelledData[, `:=`(ExpGroup = data.table::fifelse(Condition %in% ExpGroups, 1, 0), TimeVar = as.numeric(paste(TimeVar)))]
      LightModelledData[, Abundance := exp(`(Intercept)` + TimeVar*TimeCoeff + TimeVar*ExpGroup*`Time:Grouping_Exp`)]
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Calculating Light Mean Absolute Percentage Errors (MAPEs)"))
    {
      suppressWarnings(LightMAPEData <- unique(ModelProteins_L[, .(ProteinGroup, Grouping, MeanAbundance)]) |> data.table::merge.data.table(LightModelledData[, .(ProteinGroup, Grouping, Abundance)]))
      LightMAPEData <- LightMAPEData[, .(MAPE = MetricsWeighted::mape(MeanAbundance, Abundance)), .(ProteinGroup, Grouping)]
      LightMAPEData[, MAPEBin := data.table::fifelse(MAPE > 50, "MAPE > 50", data.table::fifelse(MAPE > 25, "MAPE > 25",
                                                                                                 data.table::fifelse(MAPE > 10, "MAPE > 10", "MAPE ≤ 10")))]
      LightMAPEData[, Channel := "Light"]
    }
    if(GenerateDataPlots){
      message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Generating Modelled Light Data Plots"))
      {
        setwd(paste0(InputDirectory, "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"))
        if(dir.exists("LightPlots")){unlink("LightPlots", recursive = T)}
        dir.create("LightPlots", showWarnings = T)
        setwd("LightPlots")

        for(POI in unique(LightModelledData$ProteinGroup)){
          POIModelData <- LightModelledData[ProteinGroup == POI]
          POIModelData <- Proteopedia::Separate_Isoforms(POIModelData)
          POIModelData <- Proteopedia::Add_ProteinInfo(POIModelData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
          POIModelData[, Grouping := data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName,
                                                         data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName, "Exclude"))]

          POIActualData <- ModelProteins_L[ProteinGroup == POI]
          POIActualData <- Proteopedia::Separate_Isoforms(POIActualData)
          POIActualData <- Proteopedia::Add_ProteinInfo(POIActualData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
          POIActualData[, Grouping := data.table::fifelse(Condition %in% CtlGroups, CtlGroupsName,
                                                          data.table::fifelse(Condition %in% ExpGroups, ExpGroupsName, "Exclude"))]

          pdf(paste0(POI, "_LightPlot.pdf"), width = 16, height = 12)
          print(POIModelData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(paste(TimeVar)), y = Abundance, colour = Grouping)) +
                  ggplot2::geom_line(ggplot2::aes(group = Grouping), linetype = "dashed") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::geom_point(data = POIActualData, ggplot2::aes(x = as.numeric(paste(Time)), y = Abundance)) +
                  ggplot2::facet_wrap(~paste0(ProteinGroup, " (", unique(POIActualData[ProteinGroup == POI, Gene]), ")"), scales = 'free_y') +
                  ggplot2::annotate('text', x = quantile(TimeLevels, 0.75), y = max(ModelProteins_L[ProteinGroup == POI, Abundance], na.rm = T)*0.99,
                                    label = glue::glue('{CtlGroupsName} kloss = {round(LightParameters[ProteinGroup == POI, Ctl_Value],3)}\n
                                                           {ExpGroupsName} kloss = {round(LightParameters[ProteinGroup == POI, Exp_Value], 3)}')) +
                  ggplot2::labs(x = "Time (hours)", y = "Log Light Protein Abundance") +
                  ggplot2::theme(legend.title = ggplot2::element_blank()))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting Light Data Files"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      data.table::fwrite(ProteinLFQs_L[, .(ProteinGroup, Sample, Condition, Conc, Time, Replicate, Abundance, MeanAbundance, NormAbundance)], "LightInputLFQs.csv")
      data.table::fwrite(LightModelParameters, "LightModelOutput.csv")
      data.table::fwrite(Proteopedia::Separate_Isoforms(LightParameters), "LightParameters.csv")
      data.table::fwrite(Proteopedia::Separate_Isoforms(LightModelledData), "LightModelledData.csv")
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting Light Output Plots"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      pdf("LightOutputPlots.pdf", width = 16, height = 10)
      print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = Prop_NA)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "#000", max.overlaps = 10) +
              Proteopedia::Add_KlossAxes())
      print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = Prop_NA)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "#000", max.overlaps = 10) +
              Proteopedia::Add_KlossAxes(scale = "Log2FC"))
      print(LightParameters[, .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "KlossL") |>
              ggplot2::ggplot(ggplot2::aes(x = KlossL, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = "Density of Proteins") +
              ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside", legend.position.inside = c(0.8, 0.8)))
      print(LightParameters[, .(ProteinGroup, P.Value, adj.P.Val)] |> data.table::copy() |>
              data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |>
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Adjustment", value.name = "P") |>
              ggplot2::ggplot(ggplot2::aes(x = P)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
              ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = 0) +
              ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
      Proteopedia::Reset_Dev()
    }
    # Analyse Heavy Channel
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Isolating Heavy Channel Proteins"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"))
      ProtLFQsInput <- ProtLFQsInput |> data.table::merge.data.table(ProtLFQsInput[Time == max(TimeLevels), .(SummitLineMean = mean(LFQ_H, na.rm = T)), Condition])
      ProtLFQsInput[, RelLFQ_H := LFQ_H/SummitLineMean]
      ProtLFQsInput_H <- ProtLFQsInput |> data.table::dcast(ProteinGroup ~ Sample,  value.var =  "RelLFQ_H") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
      ProteinLFQs_H <- ProtLFQsInput_H |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
      ProteinLFQs_H <- data.table::melt.data.table(ProteinLFQs_H, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |>
        data.table::merge.data.table(Metadata)
      ProtWeights_H <- ProteinWeights[Channel =="Heavy"][,`:=`(Replicate = factor(Replicate), ProteinGroup = factor(ProteinGroup),
                                                               Time = factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels)),
                                                               Condition = factor(Condition, levels = ConditionLevels))]

      ModelProteins_H <- ProteinLFQs_H[,.(ProteinGroup, Condition, Time, Replicate, Abundance)]
      ModelProteins_H[, Replicate := factor(Replicate)]
      ModelProteins_H[Time == 0, Abundance := 0]
      ModelProteins_H <- na.omit(ModelProteins_H)
      ModelProteins_H[, N_Quant := .N, .(ProteinGroup, Condition, Time)] # Max. N_Quant is 3 (3 Replicates)
      ModelProteins_H[, N_QuantTotal := .N, ProteinGroup] # Max. N_QuantTotal is 36 (3 Replicates, 6 Timepoints, 2 Conditions: 3x6x2 = 36)

      HeavyAbunBoxplot <- ModelProteins_H |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")),
                                                                          y = Abundance, colour = Condition)) +
        ggplot2::geom_boxplot(outliers = F) + ggplot2::labs(x = "Time", y = "Heavy Protein Abundance") +
        ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Filtering By Heavy Channel Missingness & Variation"))
    {
      PreFiltReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Condition, Time, N_Quant)] |> unique()
      PreFiltTotalCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)] |> unique()

      ModelProteins_H <- ModelProteins_H[N_QuantTotal >= HeavyMinSamples*nrow(Metadata)]
      ModelProteins_H[, CV := sd(Abundance)/mean(Abundance), .(ProteinGroup, Condition, Time)]
      ModelProteins_H_CV <- ModelProteins_H[,.(MeanCV = mean(CV, na.rm =T)), ProteinGroup]
      ModelProteins_H <- ModelProteins_H[ProteinGroup %in% ModelProteins_H_CV[MeanCV <= MaxCV, ProteinGroup]]

      TotalReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Condition, Time, N_Quant)][, Filtering := "Post-Filtering"] |> unique() |>
        rbind(PreFiltReplicateData[, Filtering := "Pre-Filtering"])
      TotalProteinCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)][, Filtering := "Post-Filtering"] |> unique() |>
        rbind(PreFiltTotalCountData[, Filtering := "Pre-Filtering"])
      TotalReplicateData[, Filtering := factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))]

      HeavyCompleteness <- TotalProteinCountData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(N_QuantTotal))) +
        ggplot2::geom_histogram() + ggplot2::labs(x = "No. Samples", y = "No. Proteins") +
        ggplot2::geom_vline(xintercept = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, linetype = "dashed", colour = "#F00") +
        ggplot2::facet_wrap(~factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, ymin = -Inf, ymax = Inf,
                          fill = "#F00", alpha = 0.3) +
        ggplot2::scale_x_continuous(limits = c(0, nrow(Metadata)), expand = 0) + ggplot2::scale_y_continuous(expand = 0)

      HeavyCounts <- TotalReplicateData |> ggplot2::ggplot(ggplot2::aes(x = N_Quant, fill = Condition)) +
        ggplot2::geom_bar(stat = "count", position = ggplot2::position_dodge()) +
        ggplot2::facet_grid(ggplot2::vars(Filtering), ggplot2::vars(factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")))) +
        ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, name = "Condition") + ggplot2::labs(x = "Replicates", y = "No. Proteins") +
        ggplot2::theme(strip.text.y = ggplot2::element_text(size = 26),legend.title = ggplot2::element_blank(),
                       legend.position = "inside", legend.position.inside = c(0.2, 0.9))
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Filtering By Heavy Channel Monotonicity"))
    {
      MonotonicityData <- ModelProteins_H[,.(Mean_Abundance = mean(Abundance, na.rm = T)), .(ProteinGroup, Time, Condition)][(order(ProteinGroup, Condition, Time))]
      MonotonicityData[, CumSum := cummax(Mean_Abundance), .(ProteinGroup, Condition)]
      MonotonicityData[, Monotonic := CumSum == Mean_Abundance]
      MonotonicitySummary <- MonotonicityData[, .(N_Monotonic = sum(Monotonic)), "ProteinGroup"]
      MonotonicitySummary[, MonotonicityProp := N_Monotonic/(length(TimeLevels)*length(ConditionLevels))]

      MonotonicityPlot <- MonotonicitySummary |> ggplot2::ggplot(ggplot2::aes(x = MonotonicityProp)) +
        ggplot2::geom_histogram() + ggplot2::labs(x = "Prop. Monotonicity", y = "No. Proteins") +
        ggplot2::geom_vline(xintercept = HeavyMinMonotonicity-0.033, linetype = "dashed", colour = "#F00") +
        ggplot2::annotate("rect", xmin = -Inf, xmax = HeavyMinMonotonicity-0.033, ymin = -Inf, ymax = Inf, fill = "#F00", alpha = 0.3) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(expand = 0)

      ModelProteins_H <- ModelProteins_H[ProteinGroup %in% MonotonicitySummary[MonotonicityProp >= HeavyMinMonotonicity, ProteinGroup]]
      ModelProteins_H[, ProteinGroup := as.factor(ProteinGroup)]

      HeavyMissingnessData <- ModelProteins_H[, .(ProteinGroup, N_QuantTotal)] |> data.table::copy() |> unique()
      HeavyMissingnessData[, PropNAs := N_QuantTotal/nrow(Metadata)]

      # Processed Heavy Data PCA
      All_PCA <- ModelProteins_H[, .(ProteinGroup, Condition, Time, Replicate, Abundance)] |> data.table::merge.data.table(Metadata[, .(Condition, Time, Replicate, Sample)]) |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Abundance, names_from = Sample, values_fill = NA) |>
        tidyr::drop_na() |> data.frame(row.names = c("ProteinGroup")) |> t() |> stats::prcomp(scale. = T)
      SummaryPCA <- summary(All_PCA)$importance
      All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

      HeavyProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                      y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
        All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
        ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
        ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                      y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
        patchwork::plot_annotation(title = "Processed Heavy Data")
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Running ", HeavyModel, " Modelling"))
    {
      KlossParameters <- LightModelParameters[, .(ProteinGroup, Time, `Time:Grouping_Exp`)] |> data.table::copy()
      KlossParameters[, `:=`(Kloss_Ctl = -Time, Kloss_Exp = -`Time:Grouping_Exp` - Time)]
      Run_ProteinNLS <- function(POI, Progress = NULL){
        tryCatch(
          expr = {
            POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProteinWeights[Channel == "Heavy", .(ProteinGroup, Condition, Replicate, Weight)], all.x = T)
            POIData <- POIData[order(Condition, Time)]
            POIData[,`:=`(Comparison = data.table::fifelse(Condition == CtlGroup, 0, 1))]
            POIData[, Time := as.numeric(paste(Time))]
            POIData[, VAR := var(Abundance), .(Condition, Time)]
            POIData[is.na(VAR), VAR := max(POIData$VAR, na.rm  = T)]
            # Define T0 Data
            T0Data <- POIData[, head(.SD,3), Condition]
            T0Data[,`:=`(Time = 0, Abundance = 0, Weight = min(POIData$Weight, na.rm =T), VAR = min(POIData$VAR, na.rm = T))]
            POIData <- POIData[Time != 0] |> rbind(T0Data)

            T0Data <- POIData[,head(.SD,1), by = Condition]
            T0Data[,`:=`(Time = 0, Abundance = 0)]
            POIData <- POIData |> rbind(T0Data)

            KsynStart <- sapply(ConditionLevels, function(COI){
              mean(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Time)), na.rm =T)
            })

            AbundancePlateau <- sapply(ConditionLevels, function(COI){
              mean(POIData[order(Time)][Condition == COI][,tail(.SD,3)]$Abundance,na.rm =T)
            })

            StartVals <- c(KsynStart, KsynStart / AbundancePlateau) # nls requires start values of approx. params
            # nlsLM more stable than nls
            POIFit <- minpack.lm::nlsLM(Abundance ~ (Ksyn_Ctl+(Ksyn_Exp*Comparison))/(Kloss_Ctl+(Kloss_Exp*Comparison)) * (1-exp(-(Kloss_Ctl+(Kloss_Exp*Comparison))*Time)), data = POIData,
                                        start = list(Ksyn_Ctl = StartVals[[1]],
                                                     Ksyn_Exp = StartVals[[2]]/10,
                                                     Kloss_Ctl = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Ctl], StartVals[[3]]),
                                                     Kloss_Exp = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Exp], StartVals[[4]]/10)),
                                        weight = 1/POIData$Weight, control = nls.control(maxiter = 200, warnOnly = T))
            POIData <- POIData[, Fitted := predict(POIFit)]
            POIData <- Proteopedia::Add_ProteinInfo(POIData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))

            FitSummary <- summary(POIFit)
            if(GenerateDataPlots){
              pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
              print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Condition)) + ggplot2::geom_point() +
                      ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), .(Condition, Time, ProteinGroup)],
                                         ggplot2::aes(y = Fitted, group = interaction(Condition)), linetype = "dashed") +
                      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                      ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", POIData[ProteinGroup == POI, Gene], ")"),
                                    subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients[3,1],3)}, ksyn = {round(FitSummary$coefficients[1,1],3)}
                                                                {ExpGroup}: kloss = {round(FitSummary$coefficients[3,1]+ FitSummary$coefficients[4,1],3)}, ksyn = {round(FitSummary$coefficients[1,1]+ FitSummary$coefficients[2,1],3)}')) +
                      ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
              Proteopedia::Reset_Dev()
            }
            message(POI, ": Model Applied ", Progress)
            return(list(POIModel = POIFit, FittedData = POIData))
          },
          error = function(e){
            message(POI, ": Error Caught ", Progress)
          })
      }
      Run_ProteinNLME <- function(POI, Progress = NULL){
        tryCatch(
          expr = {
            POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProtWeights_H[, .(ProteinGroup, Condition, Replicate, Time, Weight)], all.x = T)
            POIData[, Acquisition := data.table::fifelse(as.numeric(paste(Time)) > max(TimeLevels)/2, "Late", "Early")]
            POIData[is.na(Acquisition), Acquisition := "Early"]
            POIData[, Time := as.numeric(paste(Time))]
            POIData[Time == 0, Weight := max(POIData$Weight, na.rm = T)]
            # Define T0 Data
            T0Data <- POIData[,head(.SD,1), Condition]
            T0Data[,`:=`(Time = 0, Abundance = 0)]
            POIData <- POIData[Time != 0] |> rbind(T0Data)

            KsynStart <- sapply(ConditionLevels, function(COI){
              mean(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Time)), na.rm =T)
            })

            AbundancePlateau <- sapply(ConditionLevels, function(COI){
              mean(POIData[order(Time)][Condition == COI][,tail(.SD,3)]$Abundance,na.rm =T)
            })

            StartVals <- c(KsynStart, KsynStart / AbundancePlateau)

            POIFit <- nlme::nlme(Abundance ~ (Ksyn/Kloss)*(1 - exp(-Kloss*as.numeric(paste(Time)))), data = POIData, fixed = list(Ksyn ~ Condition, Kloss ~ Condition),
                                 random = list(Acquisition = nlme::pdDiag(Ksyn+Kloss ~1)),
                                 start = StartVals, weights = nlme::varFixed(~Weight))
            POIData <- POIData[, Fitted := predict(POIFit)]

            FitSummary <- summary(POIFit)
            if(GenerateDataPlots){
              pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
              print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Condition)) + ggplot2::geom_point() +
                      ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), .(Condition, Time, ProteinGroup)],
                                         ggplot2::aes(y = Fitted, group = interaction(Condition)), linetype = "dashed") +
                      ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                      ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", ProteinInfo[ProteinGroup == POI, Gene], ")"),
                                    subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients$fixed[3],3)}, ksyn = {round(FitSummary$coefficients$fixed[1],3)}
                                                                  {ExpGroup}: kloss = {round(FitSummary$coefficients$fixed[3] + FitSummary$coefficients$fixed[4],3)}, ksyn = {round(FitSummary$coefficients$fixed[1] + FitSummary$coefficients$fixed[2],3)}')) +
                      ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
              Proteopedia::Reset_Dev()
            }
            message(paste0(POI, ": Model Applied ", Progress))
            return(list(POIModel = POIFit, FittedData = POIData))
          },
          error = function(e){
            message(paste0(POI, ": Error Caught ", Progress))
          }
        )
      }

      if(GenerateDataPlots){
        if (dir.exists("HeavyPlots")) {
          unlink("HeavyPlots", recursive = T)
        }
        dir.create("HeavyPlots", showWarnings = T)
        setwd("HeavyPlots")
      }

      HeavyModelledData <- data.table::data.table()
      HeavyModelParameters <- data.table::data.table()
      ProteinSigmas <- data.table::data.table()
      if(HeavyModel == "NLS"){
        for(POI in levels(ModelProteins_H$ProteinGroup)){
          NLSOutput <- Run_ProteinNLS(POI, paste0("(", round(which(levels(ModelProteins_H$ProteinGroup) == POI)/length(levels(ModelProteins_H$ProteinGroup)), digits = 2)*100, "%)"))
          HeavyModelledData <- HeavyModelledData |> rbind(NLSOutput$FittedData)
          if(!is.null(NLSOutput$POIModel$convInfo$isConv)){
            if(NLSOutput$POIModel$convInfo$isConv){
              POIModelSummary <- summary(NLSOutput$POIModel)$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "t_Value", "P.Value"))
              ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = sigma(NLSOutput$POIModel)))
              HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
            }
          }
        }
      }
      if(HeavyModel == "NLME"){
        for(POI in levels(ModelProteins_H$ProteinGroup)){
          NLMEOutput <- Run_ProteinNLME(POI, paste0("(", round(which(levels(ModelProteins_H$ProteinGroup) == POI)/length(levels(ModelProteins_H$ProteinGroup)), digits = 2)*100, "%)"))
          HeavyModelledData <- HeavyModelledData |> rbind(NLMEOutput$FittedData)
          if(class(NLMEOutput$POIModel)[1] == "nlme"){
            POIModelSummary <- summary(NLMEOutput$POIModel)$tTable |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "DF", "t_Value", "P.Value"))
            ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = NLMEOutput$POIModel$sigma))
            HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
            HeavyModelParameters[, Effect := gsub(".Grouping.*", "_Exp", gsub("..(Intercept.)", "_Ctl", HeavyModelParameters$Effect))]
          }
        }
      }
      setwd(paste0(InputDirectory, "/", ExpGroupsName, "_vs_", CtlGroupsName, "_Output"))
      HeavyModelParameters[, adj.P.Val := p.adjust(P.Value, 'BH'), 'Effect']
      HeavyModelParameters <- Proteopedia::Add_ProteinInfo(HeavyModelParameters, paste0(InputDirectory, "/report.protein_description.tsv.gz")) |>
        data.table::merge.data.table(HeavyMissingnessData, by = "ProteinGroup")

      HeavyParameters <- HeavyModelParameters |> data.table::copy()
      HeavyParameters[stringr::str_detect(Effect,'_Ctl$'), Ctl_Value := Estimate]
      HeavyParameters[, Parameter := gsub("Kloss", "KlossH", stringr::str_remove(Effect,'_.*'))]
      HeavyParameters[, Ctl_Value := mean(Ctl_Value, na.rm = T), .(ProteinGroup, Parameter)]
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
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Calculating Heavy Mean Absolute Percentage Errors (MAPEs)"))
    {
      HeavyMAPEData <- HeavyModelledData[Time != 0, .(MeanAbundance = mean(Abundance, na.rm = T), MeanFitted = mean(Fitted, na.rm = T)), .(ProteinGroup, Grouping)]
      HeavyMAPEData <- HeavyMAPEData[, .(MAPE = MetricsWeighted::mape(MeanAbundance, MeanFitted)), .(ProteinGroup, Grouping)]
      HeavyMAPEData[, MAPEBin := data.table::fifelse(MAPE > 50, "MAPE > 50", data.table::fifelse(MAPE > 25, "MAPE > 25", data.table::fifelse(MAPE > 10, "MAPE > 10", "MAPE ≤ 10")))]
      HeavyMAPEData[, Channel := "Heavy"]
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting Heavy Data Files"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      data.table::fwrite(ModelProteins_H, "HeavyInputLFQs.csv")
      data.table::fwrite(HeavyModelParameters, "HeavyModelOutput.csv")
      data.table::fwrite(Proteopedia::Separate_Isoforms(HeavyParameters), "HeavyParameters.csv")
      data.table::fwrite(HeavyModelledData, "HeavyModelledData.csv")
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting Heavy Output Plots"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      pdf("HeavyOutputPlots.pdf", width = 16, height = 10)
      print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05]) +
              Proteopedia::Add_KlossAxes())
      print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05], colour = "#000", max.overlaps = 10) +
              Proteopedia::Add_KlossAxes(scale = "Log2FC"))
      print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05]) +
              Proteopedia::Add_KsynAxes())
      print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
              ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
              ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05], colour = "#000", max.overlaps = 10) +
              Proteopedia::Add_KsynAxes(scale = "Log2FC"))
      print(HeavyParameters[Parameter == "KlossH", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "KlossH") |>
              ggplot2::ggplot(ggplot2::aes(x = KlossH, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)"), y = "Density of Proteins") +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                        legend.position.inside = c(0.8, 0.8)))
      print(HeavyParameters[Parameter == "Ksyn", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
              data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
              data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "Ksyn") |>
              ggplot2::ggplot(ggplot2::aes(x = Ksyn, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~") (Heavy Channel)"), y = "Density of Proteins") +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                        legend.position.inside = c(0.8, 0.8)))
      print(HeavyParameters[, .(ProteinGroup, P.Value, adj.P.Val, Parameter)] |> data.table::copy() |>
              data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |>
              data.table::melt.data.table(id.vars = c("ProteinGroup", "Parameter"), variable.name = "Adjustment", value.name = "P") |>
              ggplot2::ggplot(ggplot2::aes(x = P, fill = Parameter)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
              ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0))+
              ggplot2::scale_y_continuous(expand = 0) + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette) +
              ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
      Proteopedia::Reset_Dev()
    }
    # Analyse Abundance Data
    if(!SameInitialAbundance){
      message(paste0("Running ", ExpGroupsName, " vs. ", CtlGroupsName, ": 0hr Comparison"))
      {
        MinPrecursors = 2; ImputationQ = 0.01; ImputationSigma = 1
        message("Loading Input File")
        {
          setwd(InputDirectory)
          SpectraRead <- data.table::fread(list.files(pattern = "DIANN_Output.csv"))[Time == 0][, Log2LFQ := log2(LFQ_L)]
          SpectraRead |> data.table::setnames(colnames(SpectraRead)[grepl("protein.*group", ignore.case = T, colnames(SpectraRead))], "ProteinGroup")

          AbundanceMetadata <- SpectraRead[Condition %in% ComparativeMetadata$Condition, .(Sample, Condition, Replicate)] |> dplyr::distinct() |>
            data.table::merge.data.table(ComparativeMetadata)

          SpectraRead <- SpectraRead |> data.table::merge.data.table(Metadata)
          SpectraRead[, Comparative := factor(Comparative, levels = c(CtlGroupsName, ExpGroupsName))]

          if(dir.exists(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_0hrOutput"))){
            unlink(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_0hrOutput"), recursive = T)
          }
          dir.create(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_0hrOutput"), showWarnings = T)
          setwd(paste0(getwd(),"/",ExpGroupsName,"_vs_",CtlGroupsName,"_0hrOutput"))
          data.table::fwrite(Metadata, file = "Sample_Metadata.csv")
        }
        message("Performing PCA")
        {
          PCAData <- SpectraRead[, .(ProteinGroup, Sample, Log2LFQ)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
            tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
          PCASummary <- summary(PCAData)$importance
          PCAData <- data.table::data.table(PCAData$x, keep.rownames = "Sample") |> data.table::merge.data.table(AbundanceMetadata)
          PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]
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
        message("Filtering Proteins")
        {
          FilteringData <- SpectraRead[,.(N_Samples = .N, Min_Precursors = min(N_precursors), N_Conditions = data.table::uniqueN(Condition)), ProteinGroup]
          Retained1 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample) & Min_Precursors >= MinPrecursors]
          Retained2 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample)-1 & Min_Precursors >= (MinPrecursors+1)]
          Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraRead$Sample)/2) & Min_Precursors >= (MinPrecursors+1) & N_Conditions < length(unique(Metadata$Condition))]
          RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
          ExcludedProteins <- SpectraRead[ProteinGroup %!in% RetainedProteins]
          SpectraRead <- SpectraRead[ProteinGroup %in% RetainedProteins]
          data.table::fwrite(ExcludedProteins, file = "Excluded_Proteins.csv")
          data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")
          FilteringData <- SpectraRead[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")] |> rbind(ExcludedProteins[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")])
          FilteringData <- FilteringData[, .N, .(Sample, Condition, Replicate, Inclusion)] |> data.table::setorder(Condition)

          CountsBar <- FilteringData |> ggplot2::ggplot(ggplot2::aes(x = Sample, y = N, fill = Condition, alpha = Inclusion)) +
            ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") +
            ggplot2::geom_text(ggplot2::aes(label = N), colour = data.table::fifelse(FilteringData$Inclusion == "Retained", "#FFF","#000"), position = ggplot2::position_stack(), vjust = 1.5) +
            ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x", nrow = 1) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) + ggplot2::labs(x = NULL, y = "Count", fill = NULL) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), strip.text.x = ggplot2::element_blank(),
                           strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0,"line"))

          UpsetPlot <- SpectraRead[, .(Sample = list(gsub("_", " ", Sample))), ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
            ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) +
            ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = SpectraRead[order(Condition), unique(gsub("_", " ", Sample))]) +
            ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15)))
        }
        message("Performing Median Normalisation")
        {
          SpectraRead[, Log2LFQ_Norm := Log2LFQ - median(Log2LFQ, na.rm = T) + median(SpectraRead$Log2LFQ, na.rm = T), Sample]

          suppressWarnings(
            NormPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraRead, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                        ggplot2::aes(x = Condition, fill = Condition, y = value, group = Sample))+
              ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
              ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") +
              ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::labs(x = NULL, y = expression("Log"[2]~"LFQ"))
          )
        }
        message("Imputing Undetected Condtiion Values")
        {
          SpectraRead <- SpectraRead[,.(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2LFQ_Norm")
          SpectraImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma),
                                               keep.rownames = "ProteinGroup")
          SpectraAll <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead, id.vars = "ProteinGroup", value.name = "Measured_LFQ", variable.name = "Sample"),
                                                     data.table::melt.data.table(SpectraImp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ", variable.name = "Sample"))
          SpectraAll[, Data := data.table::fifelse(is.na(Measured_LFQ), "Imputed","Measured")]

          ImpPlot <- ggplot2::ggplot(SpectraAll, ggplot2::aes(x = Imputed_LFQ, fill = Data)) + ggplot2::geom_density(adjust = 2, alpha = 0.8) + ggplot2::scale_y_continuous(expand = 0) +
            ggplot2::scale_fill_manual(values = c(Measured = "#000", Imputed = "#C0C")) + ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) +
            ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))

          ImputedNAs <- data.table::melt.data.table(SpectraRead[ProteinGroup %in% Retained3$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
          ImputedNAs <- ImputedNAs[is.na(value)]
          ImputedNAs[, value := NULL]
          ImputedNAs <- ImputedNAs[, .(Vector = paste(Sample, collapse = ", ")), ProteinGroup]
          data.table::fwrite(ImputedNAs, "Imputed_LFQs.csv")
          SpectraAll[, Log2LFQ := data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ)]
          SpectraAll <- SpectraAll[, .(ProteinGroup, Sample, Log2LFQ)]
        }
        message("Performing Paired T-Testing")
        {
          SpectraAll <- SpectraAll |> data.table::merge.data.table(Metadata[, .(Comparative, Sample)], by = "Sample")
          SpectraAll[, LFQ := 2^(Log2LFQ)]

          SpectraTtest <- SpectraAll[, .(Log2MeanLFQ = log2(mean(LFQ)), CV = Proteopedia::Calculate_CV(LFQ), .N), .(Comparative, ProteinGroup)] |>
            data.table::dcast(ProteinGroup ~ Comparative, value.var = c("Log2MeanLFQ", "CV", "N"))

          CtlColIndex <- colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(CtlGroupsName, colnames(SpectraTtest)))]
          ExpColIndex <-  colnames(SpectraTtest)[which(grepl("Log2MeanLFQ_", colnames(SpectraTtest)) & grepl(ExpGroupsName, colnames(SpectraTtest)))]
          SpectraTtest$Log2FC <- SpectraTtest[, get(ExpColIndex)] - SpectraTtest[, get(CtlColIndex)]

          Ttest_Output <- SpectraAll[, P.Value := stats::t.test(Log2LFQ ~ Comparative)$p.value, ProteinGroup]
          SpectraTtest <- data.table::merge.data.table(data.table::data.table(SpectraTtest), Ttest_Output, by = "ProteinGroup")
          SpectraTtest <- Proteopedia::Add_ProteinInfo(SpectraTtest, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
          SpectraTtest[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
          SpectraTtest[, Log2FC := as.numeric(Log2FC)]
          data.table::fwrite(SpectraTtest |> dplyr::distinct(), "Paired_T-Test_Output.csv")
        }
        message("Fitting Linear Model")
        {
          ModelDesign <- stats::model.matrix(~0 + Comparative, data = AbundanceMetadata)
          colnames(ModelDesign) <- gsub("Comparative", "", colnames(ModelDesign))
          rownames(ModelDesign) <- AbundanceMetadata$Sample
          ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroupsName, colnames(ModelDesign))), which(grepl(ExpGroupsName, colnames(ModelDesign)))))]

          ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
          ContrastMatrix[,1] <- c(-1,1)

          LimmaInput <- SpectraAll |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
            data.table::setcolorder(c("ProteinGroup", Metadata$Sample))

          suppressMessages(ModelFit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, ModelDesign), ContrastMatrix)))

          if(!is.finite(ModelFit$df.prior)){message("Warning: Limma Prior is Infinite")}

          MeanVarData <- data.table::data.table(ModelFit$genes, "Mean" = ModelFit$Amean, "Variance" = sqrt(ModelFit$sigma))
          MeanVarData[, Data := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured")]

          MeanVarPlot <- MeanVarData |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) +
            ggplot2::geom_point(colour = "#000") +  ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C") +
            ggplot2::labs(x = "Mean Log2LFQ", y = "Variance", colour = NULL)

          LimmaOutput <- data.table::data.table(limma::topTable(ModelFit, coef=1, adjust.method = "BH", n=Inf)) |> data.table::setnames("logFC", "Log2FC")
          LimmaOutput <- LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = T)]
          LimmaOutput <- Proteopedia::Add_ProteinInfo(LimmaOutput, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
          LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
          LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = T)
          LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")]
          LimmaOutput <- Proteopedia::Separate_Isoforms(LimmaOutput)
          data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
        }
        message("Generating Volcano Plots")
        {
          LimmaVolcano <- LimmaOutput |> dplyr::arrange(desc(abs(t))) |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
            ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) + ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
            ggrepel::geom_text_repel(ggplot2::aes(label= data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
            ggplot2::geom_vline(xintercept = mean(LimmaOutput$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
            Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

          TtestVolcano <- unique(SpectraTtest[, .(ProteinGroup, Log2FC, P.Value, Gene)]) |> dplyr::arrange(desc(abs(Log2FC)*-log10(P.Value))) |>
            ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
            ggplot2::scale_colour_manual("#000") + Proteopedia::Add_NotSigBox() +
            ggrepel::geom_text_repel(ggplot2::aes(label = data.table::fifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
            ggplot2::geom_vline(xintercept = mean(SpectraTtest$Log2FC, na.rm = T), linetype = "dashed", colour = "#000") +
            Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("T-Test") + ggplot2::theme(legend.title = ggplot2::element_blank())

          pdf("VolcanoPlots.pdf", height = 10, width = 14)
          print(LimmaVolcano)
          print(LimmaVolcano + ggplot2::geom_point(data = LimmaOutput[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed Limma"))
          print(TtestVolcano)
          print(TtestVolcano + ggplot2::geom_point(data = SpectraTtest[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed T-Test"))
          Proteopedia::Reset_Dev()
        }
        message("Exporting QC Plot")
        {
          pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
          suppressWarnings(
            print(patchwork::free(PCAPlot, type= "label") + CountsBar + patchwork::free(UpsetPlot, type = "label") +
                    patchwork::free(NormPlot) + patchwork::free(ImpPlot) + patchwork::free(MeanVarPlot, type = "label") +
                    patchwork::plot_layout(design = "AAAABBBB\nCCDDEEEE\nFFFFFFFF") + patchwork::plot_annotation(tag_levels = "A"))
          )
          Proteopedia::Reset_Dev()
        }
        message("Exporting HTML Volcano Plot")
        {
          InteractiveData <- data.table::data.table("Protein" = LimmaOutput$Gene, "Log2FC" = round(LimmaOutput$Log2FC, digits = 2),
                                                    "PValue" = LimmaOutput$P.Value, "Significance" = factor(LimmaOutput$Significance),
                                                    "Scien_PValue" = formatC(LimmaOutput$P.Value, format = "e", digits = 2),
                                                    "GeneGroup" = LimmaOutput$GeneGroup, "URL" = LimmaOutput$URL)

          HighCharterVolcano <- highcharter::hchart(InteractiveData, "scatter", highcharter::hcaes(x = Log2FC, y = -log10(PValue), group = Significance)) |>
            highcharter::hc_chart(zoomType = "xy") |>
            highcharter::hc_xAxis(title = list(text = paste0(unique(Metadata$Condition)[1]," vs ", unique(Metadata$Condition)[2]," Log2 Fold-Change")),
                                  lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
            highcharter::hc_yAxis(title = list(text = "-Log10 P-Value"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000", gridLineWidth = 0 ) |>
            highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Protein} | {point.GeneGroup} </b> <br>Log2FC: {point.Log2FC:.2f}<br>p-value: {point.Scien_PValue:.2f}") |>
            highcharter::hc_plotOptions(scatter = list(marker = list(radius = 3), states = list(hover = list(enabled = T), inactive = list(enabled = F)),
                                                       point = list(events = list( click = htmlwidgets::JS("function() { window.open(this.URL, '_blank'); }"))))) |>
            highcharter::hc_colors(c("#999", "#800","#03F"))
          htmlwidgets::saveWidget(HighCharterVolcano, "InteractiveVolcanoPlot.html")
        }
        message("Exporting Analysis Parameters")
        data.table::fwrite(data.table::data.table("Experimental Condition(s)" = paste(ExpGroups, collapse = ", "), "Experimental Name" = paste0(ExpGroupsName),
                                                  "Control Condition(s)" = paste(CtlGroups, collapse = ", "), "Control Name" = paste0(CtlGroupsName),
                                                  "Min_Precursors" = paste0(MinPrecursors), "Imputation Q-Value" = ImputationQ,
                                                  "Imputation Sigma" = ImputationSigma), "Analysis_Parameters.csv")
      }
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting Analysis Summary"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      AnalysisSummary <- data.table::data.table(Parameter = c("Light kloss", "Heavy kloss", "ksyn", "0hr Abundance"),
                                                N_Proteins = c(nrow(LightParameters), nrow(HeavyParameters[Parameter == "KlossH"]) ,
                                                               nrow(HeavyParameters[Parameter == "Ksyn"]), nrow(AbundanceData)),
                                                Model = c("Limma", HeavyModel, HeavyModel, "Limma"))
      AnalysisSummary[, Prop_Proteins := N_Proteins/length(ProtLFQsInput[, ProteinGroup] |> unique())]
      data.table::fwrite(AnalysisSummary, "Analysis_Summary.csv")
    }
    message(paste0(ExpGroupsName, " vs. ", CtlGroupsName, ": Exporting QC Data, Plots & Comparisons"))
    {
      setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
      MAPEData <- Proteopedia::Separate_Isoforms(LightMAPEData |> rbind(HeavyMAPEData))
      MAPEData[, MAPEBin := factor(MAPEBin, levels = c("MAPE ≤ 10", "MAPE > 10", "MAPE > 25", "MAPE > 50"))]
      data.table::fwrite(MAPEData, "MAPEData.csv")
      WideMAPEData <- MAPEData |> data.table::dcast(ProteinGroup+Isoforms+Channel+Time ~ Condition, value.var = "MAPE")
      WideMAPEData[, DiffMAPE := get(ExpGroup) - get(CtlGroup)]

      pdf("MAPEPlots.pdf", width = 16, height = 12)
      print(MAPEData[Channel == "Light"] |> ggplot2::ggplot(ggplot2::aes(x = MAPE, fill = MAPEBin)) + ggplot2::geom_histogram() +
              ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(16, length(Proteopedia::ThermalPalette), length.out = 4)]) +
              ggplot2::facet_grid(ggplot2::vars(Condition), ggplot2::vars(factor(paste(Time, "h"), levels = paste(TimeLevels, "h"))), scales = "free") +
              ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "Mean Absolute Percentage Error (MAPE)", y = "No. Light Proteins") +
              ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(MAPEData[Channel == "Heavy"] |> ggplot2::ggplot(ggplot2::aes(x = MAPE, fill = MAPEBin)) + ggplot2::geom_histogram() +
              ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(16, length(Proteopedia::ThermalPalette), length.out = 4)]) +
              ggplot2::facet_grid(ggplot2::vars(Condition), ggplot2::vars(paste(Time, "h")), scales = "free") +
              ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "Mean Absolute Percentage Error (MAPE)", y = "No. Light Proteins") +
              ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(WideMAPEData |> ggplot2::ggplot(ggplot2::aes(x = DiffMAPE, fill = Channel)) + ggplot2::geom_density(alpha = 0.3) +
              Proteopedia::Add_Isotope_Fill() + ggplot2::facet_wrap(~paste(Time, "h"), scales = "free", nrow = 1) +
              ggplot2::labs(x = paste0("Difference in MAPE (", ExpGroup, " - ", CtlGroup, ")"), y = "No. Proteins") +
              ggplot2::scale_y_continuous(expand = 0))
      Proteopedia::Reset_Dev()

      CVPlot <- ProteinLFQs_L_CV[, Channel := "Light"] |> rbind(ModelProteins_H_CV[, Channel := "Heavy"]) |>
        ggplot2::ggplot(ggplot2::aes(x = MeanCV)) + ggplot2::geom_histogram() +
        ggplot2::facet_wrap(~forcats::fct_rev(Channel), nrow = 3, strip.position = "right") +
        ggplot2::labs(x = "Mean Protein Variation", y = "No. Proteins") + ggplot2::geom_vline(xintercept = MaxCV, linetype = "dashed", colour = "#F00") +
        ggplot2::annotate("rect", xmin = MaxCV, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F00", alpha = 0.3) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = 0) + ggplot2::scale_y_continuous(expand = 0)

      pdf("QC_Plots.pdf", width = 18, height = 20)
      print(patchwork::free(LightAbunBoxplot + HeavyAbunBoxplot + patchwork::plot_layout(guides = "collect")) +
              ProteinWeights_PCA + ProteinWeight_Density + CVPlot +
              LightMeanVarPlot + MonotonicityPlot + patchwork::free(HeavyCounts) + HeavyCompleteness) +
        patchwork::plot_layout(design = "AAAAAAAAAA\nBBBBBBCCCC\nEEEEEDDDDD\nFFFFFDDDDD") +
        patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "" ,"H", "")))
      Proteopedia::Reset_Dev()

      CorrData <- HeavyParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)] |>
        rbind(LightParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)]) |>
        data.table::melt.data.table(id.vars = c("ProteinGroup", "Gene", "Parameter"),
                                    variable.name = "Condition", value.name = "Measure")
      CorrData[, Condition := data.table::fifelse(grepl("Ctl", Condition), CtlGroup, ExpGroup)]
      CorrData <- CorrData |> data.table::dcast(ProteinGroup+Gene+Condition ~ Parameter, value.var = "Measure")

      pdf("CorrelationPlots.pdf", width = 16, height = 12)
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = KlossL, y = KlossH, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
              ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
              ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossH, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
              ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
      print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossL, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
              Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
              ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Light Channel)")) +
              ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
      Proteopedia::Reset_Dev()
    }
  } else {
    for(ExpGroup in ExpGroups){
      for(CtlGroup in CtlGroups){
        setwd(InputDirectory)
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Loading Data"))
        {
          ProtLFQsInput <- data.table::fread(list.files(pattern = "DIANN_Output"))[Time %!in% ExcludeTimepoints & (grepl(CtlGroup, Sample)|grepl(ExpGroup, Sample))] |> data.table::setnames("Protein_group", "ProteinGroup", skip_absent = T)
          ProtLFQsInput |> data.table::setnames(colnames(ProtLFQsInput)[grepl("Protein.*group", colnames(ProtLFQsInput), ignore.case = T)], "ProteinGroup")
          ProtLFQsInput <- ProtLFQsInput[gsub(".*(\\d+)$", "\\1", Replicate) %in% ReplicatesUsed]

          for(ColumnID in c("Conc")){
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

          Metadata <- data.table::fread("Sample_Metadata.csv")[Time %!in% ExcludeTimepoints][grepl(CtlGroup, Sample)|grepl(ExpGroup, Sample)] |> data.table::setorderv(c("Conc", "Time", "Replicate"))
          Metadata[, Cluster := gsub("(.*)_\\d+h", "\\1", Condition)]
          Metadata[, Cluster := factor(Cluster, levels = c(CtlGroup, ExpGroup))]
          Metadata |> data.table::setorder(Cluster)
          SampleLevels <- unique(Metadata$Sample)
          ConcLevels <- unique(Metadata$Conc)
          TimeLevels <- unique(Metadata$Time)
          ConditionLevels <- unique(Metadata$Condition)
          ClusterLevels <- unique(Metadata$Cluster)
          Metadata[, Sample := factor(Sample, levels = SampleLevels)]
          Metadata[, Conc := factor(Conc, levels = ConcLevels)]
          Metadata[, Time := factor(Time, levels = TimeLevels)]
          Metadata[, Condition := factor(Condition, levels = ConditionLevels)]
          Metadata[, Replicate := factor(Replicate)]
          Metadata[, Cluster := factor(Cluster, levels = ClusterLevels)]

          ProteinWeights <- data.table::fread("Filtered_PrecursorData.csv.gz")[Time %!in% ExcludeTimepoints,.(Sum_Intensities  = sum(log2(Precursor.Normalised), na.rm = T)), .(Sample, Channel, ProteinGroup)] |>
            data.table::merge.data.table(Metadata)
          ProteinWeights[, MeanSum := mean(Sum_Intensities, na.rm =T ), .(ProteinGroup, Channel)]
          ProteinWeights[, RowID := paste0(ProteinGroup, "_", Channel)]
          ProteinWeights[, Weight := Sum_Intensities/MeanSum]

          # Protein Weights PCA
          All_PCA <- ProteinWeights |> data.table::dcast(RowID~Sample, value.var = "Weight") |> tidyr::drop_na() |> data.frame(row.names = "RowID") |> t() |> stats::prcomp(scale. = T)
          SummaryPCA <- summary(All_PCA)$importance
          All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

          ProteinWeights_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
            ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
            ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                          y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
            All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
            ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
            ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                          y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
            patchwork::plot_annotation(title = "Protein Weights")

          ProteinWeight_Density <- ProteinWeights |> ggplot2::ggplot(ggplot2::aes(x = Weight, fill = Time)) +
            ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap("Channel", nrow = 2, strip.position = "right", scales = "free_x") +
            ggplot2::scale_x_log10() + ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(1, length(Proteopedia::ThermalPalette), length.out = length(TimeLevels))]) +
            ggplot2::labs(title = "Protein Weights for Limma", x = "Protein Weight", y = "No. Proteins", fill = "Time") +
            ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 26))

          ProteinWeights |> data.table::setnames(c("Channel", colnames(ProteinWeights)[grepl("Protein.*group", colnames(ProteinWeights), ignore.case = T)], colnames(ProteinWeights)[grepl("weight", colnames(ProteinWeights), ignore.case = T)],
                                                colnames(ProteinWeights)[grepl("sum.*int", colnames(ProteinWeights), ignore.case = T)], colnames(ProteinWeights)[grepl("mean.*sum", colnames(ProteinWeights), ignore.case = T)]),
                                              c("Channel", "ProteinGroup", "Weight", "SumIntensities", "MeanSum"), skip_absent = T)
          ProteinWeights <- ProteinWeights[gsub(".*(\\d+)$", "\\1", Replicate) %in% ReplicatesUsed]

          if(length(unique(grepl("Cluster", colnames(ProteinWeights)))) == 1){
            ProteinWeights[, Cluster := gsub("(.*)_\\d+h", "\\1", Condition)]
          }
          if(length(unique(grepl("Time", colnames(ProteinWeights)))) == 1){
            ProteinWeights[, Time := as.numeric(gsub(".*_(\\d+)h", "\\1", Condition))]
          }

          ProteinWeights[, Channel := gsub("L$", "Light", Channel)]
          ProteinWeights[, Channel := gsub("H$", "Heavy", Channel)]

          if (dir.exists(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"))) {
            unlink(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"), recursive = T)
          }
          dir.create(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"), showWarnings = T)
          setwd(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
        }
        # Analyse Light Channel
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Isolating Light Channel Proteins"))
        {
          ProtLFQsInput_L <- ProtLFQsInput |> data.table::copy()

          if(MeanCentring){
            ProtLFQsInput_L[grepl(CtlGroup, Condition), LFQ_Mean := mean(log(LFQ_L),na.rm = T), .(ProteinGroup, Condition, Time)]
            ProtLFQsInput_L[, LFQ_Mean := mean(LFQ_Mean,na.rm = T), .(ProteinGroup, Time)]
            ProtLFQsInput_L[, Diff := log(LFQ_L) - LFQ_Mean]
            ProtLFQsInput_L[, MeanDiff := mean(Diff, na.rm = T), .(Time, Replicate, Condition)]
            ProtLFQsInput_L[, LFQ_H:= exp(log(LFQ_L) - MeanDiff)]
          }

          ProtLFQsInput_L <- ProtLFQsInput_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = "LFQ_L") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
          ProtWeights_L <- ProteinWeights[Channel == "Light"] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Weight") |>
            tibble::column_to_rownames("ProteinGroup") |> as.matrix()

          ProteinLFQs_L <- ProtLFQsInput_L |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
          ProteinLFQs_L <- data.table::melt.data.table(ProteinLFQs_L, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |>
            data.table::merge.data.table(Metadata)

          ProteinLFQs_L[, N_Values := sum(is.finite(Abundance)), .(ProteinGroup, Condition, Time)]
          ProteinLFQs_L[, MeanAbundance := mean(Abundance, na.rm = T), ProteinGroup]
          ProteinLFQs_L[, NormAbundance := Abundance - MeanAbundance]
          ProteinLFQs_L[, Sum_Values := sum(is.finite(Abundance)), .(ProteinGroup, Time)]
          ProteinLFQs_L[, Diff_Detected := Sum_Values - N_Values]
          ProteinLFQs_L[, Time := factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels))]

          LightAbunBoxplot <- ProteinLFQs_L |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")),
                                                                            y = Abundance, colour = Condition)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
            ggplot2::labs(x = "Time", y = "Light Protein Abundance") +
            ggplot2::theme(legend.title = ggplot2::element_blank())

          ProteinLFQs_L[, NextTime := Relative_Time(Time, Times = TimeLevels), Time]
          NextTime <- ProteinLFQs_L[,.(NextTimeSamples = mean(N_Values)), .(Time, ProteinGroup, Condition)] |> data.table::setnames("Time","NextTime")
          ProteinLFQs_L <- ProteinLFQs_L |> data.table::merge.data.table(NextTime, by = c("NextTime","ProteinGroup","Condition"), all.x = T)
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Filtering Light Channel Proteins By Variation"))
        {
          ProteinLFQs_L[, N_Quant := .N, .(ProteinGroup, Condition, Time)]
          ProteinLFQs_L[, N_QuantTotal := .N, ProteinGroup]
          ProteinLFQs_L <- ProteinLFQs_L[N_QuantTotal >= LightMinSamples*nrow(Metadata)]
          ProteinLFQs_L[, CV := Proteopedia::Calculate_CV(Abundance), .(ProteinGroup, Condition, Time)]
          ProteinLFQs_L_CV <- ProteinLFQs_L[,.(MeanCV = mean(CV, na.rm =T)), ProteinGroup]
          ProteinLFQs_L <- ProteinLFQs_L[ProteinGroup %in% ProteinLFQs_L_CV[MeanCV <= (MaxCV*100), ProteinGroup]]
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
          ModelProteins_L |> data.table::setorder(Condition)

          AbundanceMatrix_L <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Abundance') |>
            tibble::column_to_rownames('ProteinGroup') |> as.matrix()
          AbundanceMatrix_L <- AbundanceMatrix_L[ , SampleLevels]
          AbundanceMatrix_L <- AbundanceMatrix_L[matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) <= LightMinSamples,]

          ProtWeights_L_Imp <- ProtWeights_L[rownames(AbundanceMatrix_L), colnames(AbundanceMatrix_L)]

          ImputationMatrix <- ModelProteins_L |> data.table::dcast(ProteinGroup ~ Sample, value.var = 'Impute') |>
            tibble::column_to_rownames('ProteinGroup') |> as.matrix()
          ImputationMatrix <- ImputationMatrix[rownames(ProtWeights_L_Imp), colnames(ProtWeights_L_Imp)]

          ProtWeights_L_Imp[ImputationMatrix==T] <- (ProtWeights_L_Imp |> min(na.rm = T))*0.1

          N_NAs <- matrixStats::rowMeans2(is.na(AbundanceMatrix_L)) |> tibble::enframe(name = "ProteinGroup", value = "Prop_NA") |> data.table::data.table()

          # Plot PCA on Processed Light Data
          All_PCA <- log(AbundanceMatrix_L) |> data.frame() |> tidyr::drop_na() |> t() |> stats::prcomp(scale. = T)
          SummaryPCA <- summary(All_PCA)$importance
          All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

          LightProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
            ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
            ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
                          y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
            All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
            ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
            ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
                          y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
            patchwork::plot_annotation(title = "Processed Light Data")
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Modelling Light Data"))
        {
          Targets <- Metadata[, .(Sample, Condition, Time)][ , Time := as.numeric(gsub("h", "", Time))]
          Targets <- Targets[match(ProtWeights_L_Imp |> colnames(), Sample)]
          Time <- Targets$Time
          Condition <- factor(Targets$Condition)
          if(SameInitialAbundance){ModelDesign <- model.matrix(~Time + Time:Condition)} else {ModelDesign <- model.matrix(~Time + Condition + Time:Condition)}
          colnames(ModelDesign) <- gsub(paste0("Condition", ExpGroup), "Condition_Exp", colnames(ModelDesign))
          colnames(ModelDesign) <- gsub(paste0("Condition", CtlGroup), "Condition_Ctl", colnames(ModelDesign))

          LimmaOutput <- limma::eBayes(limma::lmFit(Biobase::ExpressionSet(assayData = log(AbundanceMatrix_L)), ModelDesign, method = "robust",
                                                    weights = ProtWeights_L_Imp))
          Limma_Slopes <- LimmaOutput$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup")
          Limma_Slopes_SD <- LimmaOutput$stdev.unscaled |> data.table::data.table(keep.rownames = T) |> data.table::copy()
          Limma_Slopes_SD <- Limma_Slopes_SD |> data.table::setnames(c("rn", "(Intercept)", "Time", "Condition_Exp", "Time:Condition_Exp"),
                                                                     c("ProteinGroup", "(Intercept)_SD", "Time_SD", "Condition_Exp_SD", "Time:Condition_Exp_SD"), skip_absent = T)
          Limma_Slopes <- Limma_Slopes |> data.table::merge.data.table(Limma_Slopes_SD, by = "ProteinGroup")
          Limma_Slopes$Fvalue = LimmaOutput$F

          LightMeanVarPlot <- data.table::data.table(Mean = LimmaOutput$Amean, Variance = sqrt(LimmaOutput$sigma)) |>
            ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) + ggplot2::geom_point(stroke = NA) +
            ggplot2::labs(x = "Mean Log Light Protein LFQ", y = "Light Protein LFQ Variance")

          LightModelParameters <- limma::topTable(LimmaOutput, "Time:Condition_Exp", number = nrow(ProtWeights_L_Imp)) |>
            data.table::data.table(keep.rownames = T) |> data.table::setnames("rn", "ProteinGroup") |>
            data.table::merge.data.table(N_NAs) |> data.table::merge.data.table(Limma_Slopes) |> data.table::setnames("logFC", "Difference")

          LightModelParameters <- Proteopedia::Add_ProteinInfo(LightModelParameters, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
          LightParameters <- LightModelParameters |> data.table::copy()
          LightParameters[, `:=`(Ctl_Value = -Time, Exp_Value = -`Time:Condition_Exp` - Time, Difference = -Difference, Parameter = "KlossL")]

          if(OffsetKloss){
            Kloss_Offset <- abs(min(LightParameters[, .(Ctl_Value, Exp_Value)]))*1.01
            LightParameters[, `:=`(Exp_Value = Exp_Value + Kloss_Offset, Ctl_Value = Ctl_Value + Kloss_Offset)]
          } else {
            LightParameters <- LightParameters[Ctl_Value > 0 & Exp_Value > 0]
          }
          LightParameters[, FC := Exp_Value/Ctl_Value]
          LightParameters[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
          LightParameters[, Significance := data.table::fifelse(P.Value < 0.05 & Difference < 0, "Sig. Decrease",
                                                                data.table::fifelse(P.Value < 0.05 & Difference > 0, "Sig. Increase", ""))]

          LightParameters <- LightParameters[, .(ProteinGroup, Gene, P.Value, adj.P.Val, Prop_NA, Parameter, Ctl_Value, Exp_Value, FC, Log2FC, Difference)]
          if(SameInitialAbundance){
            LightModelledData <- ProtLFQsInput[, ProteinGroup] |> tidyr::crossing(Metadata[, .(Condition, Time)]) |>
              data.table::setnames(c("ProteinGroup", "Condition", "TimeVar")) |> data.table::data.table() |>
              data.table::merge.data.table(LightModelParameters[, .(ProteinGroup, `(Intercept)`, Time, `Time:Condition_Exp`)] |> data.table::setnames("Time", "TimeCoeff"))
            LightModelledData[, `:=`(ExpGroup = data.table::fifelse(grepl(ExpGroup, Condition), 1, 0), TimeVar = as.numeric(paste(TimeVar)))]
            LightModelledData[, Abundance := exp(`(Intercept)` + TimeVar*TimeCoeff + TimeVar*ExpGroup*`Time:Condition_Exp`)]
          } else {
            LightModelledData <- ProtLFQsInput[, ProteinGroup] |> tidyr::crossing(Metadata[, .(Condition, Time)]) |>
              data.table::setnames(c("ProteinGroup", "Condition", "TimeVar")) |> data.table::data.table() |>
              data.table::merge.data.table(LightModelParameters[, .(ProteinGroup, `(Intercept)`, Time, Condition_Exp, `Time:Condition_Exp`)])
            LightModelledData[, `:=`(ExpGroup = data.table::fifelse(grepl(ExpGroup, Condition), 1, 0), TimeVar = as.numeric(paste(TimeVar)))]
            LightModelledData[, Abundance := exp(`(Intercept)` + ExpGroup*Condition_Exp + TimeVar*Time + TimeVar*ExpGroup*`Time:Condition_Exp`)]
          }
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Calculating Light Mean Absolute Percentage Errors (MAPEs)"))
        {
          LightMAPEData <- unique(ModelProteins_L[, .(ProteinGroup, Condition, MeanAbundance, Time)]) |> data.table::merge.data.table(LightModelledData[, .(ProteinGroup, Condition, Abundance, as.factor(TimeVar))] |> data.table::copy() |> data.table::setnames("V4", "Time"))
          LightMAPEData <- LightMAPEData[, .(MAPE = MetricsWeighted::mape(MeanAbundance, Abundance)), .(ProteinGroup, Condition, Time)]
          LightMAPEData[, MAPEBin := data.table::fifelse(MAPE > 50, "MAPE > 50", data.table::fifelse(MAPE > 25, "MAPE > 25",
                                                                                                     data.table::fifelse(MAPE > 10, "MAPE > 10", "MAPE ≤ 10")))]
          LightMAPEData[, Channel := "Light"]
        }
        if(GenerateDataPlots){
          message(paste0(ExpGroup, " vs. ", CtlGroup, ": Generating Modelled Light Data Plots"))
          {
            setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
            if (dir.exists("LightPlots")) {
              unlink("LightPlots", recursive = T)
            }
            dir.create("LightPlots", showWarnings = T)
            setwd("LightPlots")

            for(POI in unique(LightModelledData$ProteinGroup)){
              POIModelData <- LightModelledData[ProteinGroup == POI]
              POIModelData <- Proteopedia::Separate_Isoforms(POIModelData)
              POIModelData <- Proteopedia::Add_ProteinInfo(POIModelData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))

              POIActualData <- ModelProteins_L[ProteinGroup == POI]
              POIActualData <- Proteopedia::Separate_Isoforms(POIActualData)
              POIActualData <- Proteopedia::Add_ProteinInfo(POIActualData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))

              pdf(paste0(POI, "_LightPlot.pdf"), width = 16, height = 12)
              print(POIModelData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(paste(TimeVar)), y = log(Abundance), colour = Condition)) +
                      ggplot2::geom_line(ggplot2::aes(group = Condition), linetype = "dashed") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                      ggplot2::geom_point(data = POIActualData, ggplot2::aes(x = as.numeric(paste(Time)), y = log(Abundance))) +
                      ggplot2::facet_wrap(~paste0(ProteinGroup, " (", unique(POIActualData[ProteinGroup == POI, Gene]), ")"), scales = 'free_y') +
                      ggplot2::annotate('text', x = quantile(TimeLevels, 0.75), y = log(max(ModelProteins_L[ProteinGroup == POI, Abundance], na.rm = T))*0.99,
                                        label = glue::glue('{CtlGroup} kloss = {round(LightParameters[ProteinGroup == POI, Ctl_Value],3)}\n
                                                             {ExpGroup} kloss = {round(LightParameters[ProteinGroup == POI, Exp_Value], 3)}')) +
                      ggplot2::labs(x = "Time (hours)", y = "Log Light Protein Abundance") +
                      ggplot2::theme(legend.title = ggplot2::element_blank()))
              Proteopedia::Reset_Dev()
            }
          }
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Light Data Files"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          data.table::fwrite(ProteinLFQs_L[, .(ProteinGroup, Sample, Condition, Conc, Time, Replicate, Abundance, MeanAbundance, NormAbundance)], "LightInputLFQs.csv")
          data.table::fwrite(LightModelParameters, "LightModelOutput.csv")
          data.table::fwrite(Proteopedia::Separate_Isoforms(LightParameters), "LightParameters.csv")
          data.table::fwrite(Proteopedia::Separate_Isoforms(LightModelledData), "LightModelledData.csv")
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Light Output Plots"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          pdf("LightOutputPlots.pdf", width = 16, height = 10)
          print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = Prop_NA)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "#000", max.overlaps = 10) +
                  Proteopedia::Add_KlossAxes())
          print(LightParameters |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = Prop_NA)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = LightParameters[P.Value < 0.05], colour = "#000", max.overlaps = 10) +
                  Proteopedia::Add_KlossAxes(scale = "Log2FC"))
          print(LightParameters[, .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
                  data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
                  data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "KlossL") |>
                  ggplot2::ggplot(ggplot2::aes(x = KlossL, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
                  ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = "Density of Proteins") +
                  ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside", legend.position.inside = c(0.8, 0.8)))
          print(LightParameters[, .(ProteinGroup, P.Value, adj.P.Val)] |> data.table::copy() |>
                  data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |>
                  data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Adjustment", value.name = "P") |>
                  ggplot2::ggplot(ggplot2::aes(x = P)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
                  ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = 0) +
                  ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
          Proteopedia::Reset_Dev()
        }
        # Analyse Abundance Data
        if(!SameInitialAbundance){
          message(paste0(ExpGroup, " vs. ", CtlGroup, ": Compiling Modelled Abundance Data"))
          {
            setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
            AbundanceData <- merge(limma::topTable(LimmaOutput, "Condition_Exp", number = nrow(ProtWeights_L_Imp)) |> as.data.frame() |>
                                     tibble::rownames_to_column("ProteinGroup"), N_NAs) |> data.table::data.table() |>
              merge(Limma_Slopes, by = "ProteinGroup") |> data.table::data.table() |>
              data.table::merge.data.table(LightModelledData[TimeVar == 0, .(ProteinGroup, Condition, Abundance)] |>
                                             data.table::dcast(ProteinGroup ~ Condition, value.var = "Abundance") |>
                                             data.table::setnames(c(paste(CtlGroup), paste(ExpGroup)),
                                                                  c("Ctl_Abundance", "Exp_Abundance")), all.x = T)
            AbundanceData <- Proteopedia::Separate_Isoforms(AbundanceData)
            AbundanceData <- Proteopedia::Add_ProteinInfo(AbundanceData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))
            AbundanceData[, FC := Exp_Abundance/Ctl_Abundance]
            AbundanceData[, Log2FC := Proteopedia::Calculate_VolcanoLog2FC(FC)]
            AbundanceData[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease",
                                                                data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", ""))]
            setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
            data.table::fwrite(AbundanceData, "AbundanceData.csv")
          }
          message(paste0("Running ", ExpGroup, " vs. ", CtlGroup, ": 0hr Comparison"))
          {
            MinPrecursors = 2; ImputationQ = 0.01; ImputationSigma = 1
            setwd(InputDirectory)
            SpectraRead <- data.table::fread(list.files(pattern = "DIANN_Output"))[Time == 0][, Log2LFQ := log2(LFQ_L)]
            # Loading
            {
              SpectraRead |> data.table::setnames(colnames(SpectraRead)[grepl("protein.*group", ignore.case = T, colnames(SpectraRead))], "ProteinGroup")
              AbundanceMetadata <- SpectraRead[grepl(ExpGroup, Condition)|grepl(CtlGroup, Condition), .(Sample, Condition, Replicate)] |> dplyr::distinct() |>
                dplyr::arrange(c(which(grepl(CtlGroup, Condition)), which(grepl(ExpGroup, Condition))))

              SpectraRead <- SpectraRead[Sample %in% AbundanceMetadata$Sample]
              SpectraRead[, Condition := factor(Condition, levels = c(unique(SpectraRead[grepl(CtlGroup, Condition), Condition]),
                                                                      unique(SpectraRead[grepl(ExpGroup, Condition), Condition])))]
              if (dir.exists(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_0hrOutput"))) {
                unlink(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_0hrOutput"), recursive = T)
              }
              dir.create(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_0hrOutput"), showWarnings = T)
              setwd(paste0(getwd(), "/", ExpGroup, "_vs_", CtlGroup, "_0hrOutput"))
              data.table::fwrite(Metadata, file = "Sample_Metadata.csv")

              PCAData <- SpectraRead[, .(ProteinGroup, Sample, Log2LFQ)] |>
                tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Log2LFQ, names_from = Sample, values_fill = NA) |>
                tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
              PCASummary <- summary(PCAData)$importance
              PCAData <- data.table::data.table(PCAData$x, keep.rownames = "Sample") |> data.table::merge.data.table(AbundanceMetadata)
              PCAData[, Replicate := paste0("Rep. ", gsub("R", "", Replicate))]

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
            # Filtering
            {
              FilteringData <- SpectraRead[,.(N_Samples = .N, Min_Precursors = min(N_precursors), N_Conditions = data.table::uniqueN(Condition)), ProteinGroup]
              Retained1 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample) & Min_Precursors >= MinPrecursors]
              Retained2 <- FilteringData[N_Samples == data.table::uniqueN(SpectraRead$Sample)-1 & Min_Precursors >= (MinPrecursors+1)]
              Retained3 <- FilteringData[N_Samples == floor(data.table::uniqueN(SpectraRead$Sample)/2) & Min_Precursors >= (MinPrecursors+1) & N_Conditions < length(unique(Metadata$Condition))]
              RetainedProteins <- c(Retained1$ProteinGroup, Retained2$ProteinGroup, Retained3$ProteinGroup)
              ExcludedProteins <- SpectraRead[ProteinGroup %!in% RetainedProteins]
              SpectraRead <- SpectraRead[ProteinGroup %in% RetainedProteins]
              data.table::fwrite(ExcludedProteins, file = "Excluded_Proteins.csv")
              data.table::fwrite(Retained3, file = "Imputed_Proteins.csv")
              FilteringData <- SpectraRead[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Retained")] |> rbind(ExcludedProteins[,.(ProteinGroup, Sample, Condition, Replicate, Inclusion = "Excluded")])
              FilteringData <- FilteringData[, .N, .(Sample, Condition, Replicate, Inclusion)] |> data.table::setorder(Condition)

              CountsBar <- FilteringData |> ggplot2::ggplot(ggplot2::aes(x = Sample, y = N, fill = Condition, alpha = Inclusion)) +
                ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
                ggplot2::scale_alpha_manual(values = c("Excluded" = 0.4, "Retained" = 1), guide = "none") +
                ggplot2::geom_text(ggplot2::aes(label = N), colour = data.table::fifelse(FilteringData$Inclusion == "Retained", "#FFF","#000"), position = ggplot2::position_stack(), vjust = 1.5) +
                ggplot2::facet_wrap(~Condition, strip.position = "bottom", scales = "free_x", nrow = 1) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15))) + ggplot2::labs(x = NULL, y = "Count", fill = NULL) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), strip.text.x = ggplot2::element_blank(),
                               strip.background = ggplot2::element_blank(), panel.spacing.x = grid::unit(0,"line"))

              UpsetPlot <- SpectraRead[, .(Sample = list(gsub("_", " ", Sample))), ProteinGroup] |> ggplot2::ggplot(ggplot2::aes(x = Sample)) + ggplot2::geom_bar() +
                ggplot2::geom_text(stat="count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5, size = 3) +
                ggupset::scale_x_upset(order_by = "degree", reverse = T, sets = SpectraRead[order(Condition), unique(gsub("_", " ", Sample))]) +
                ggplot2::labs(x = NULL, y = stringr::str_wrap("Post-Filtering Count", 10)) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.15)))
            }
            # Normalisation
            {
              SpectraRead[, Log2LFQ_Norm := Log2LFQ - median(Log2LFQ, na.rm = T) + median(SpectraRead$Log2LFQ, na.rm = T), by = Sample]

              NormPlot <- ggplot2::ggplot(data.table::melt.data.table(SpectraRead, measure.vars = c("Log2LFQ", "Log2LFQ_Norm")),
                                          ggplot2::aes(x = Condition, fill = Condition, y = value, group = Sample)) +
                ggplot2::facet_wrap("variable", labeller = ggplot2::labeller(variable = c("Log2LFQ" = "Pre", "Log2LFQ_Norm" = "Post")))+
                ggplot2::geom_boxplot(outliers = F) + ggplot2::ggtitle("Normalisation") +
                ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, guide = "none") + ggplot2::labs(x = NULL, y = expression("Log"[2]~"LFQ"))
            }
            # Imputation
            {
              SpectraRead <- SpectraRead[,.(ProteinGroup, Sample, Log2LFQ_Norm)] |> data.table::dcast(ProteinGroup~Sample, value.var = "Log2LFQ_Norm", fill = NA)
              SpectraImp <- data.table::data.table(imputeLCMD::impute.MinProb(as.matrix(SpectraRead, rownames = "ProteinGroup"), q = ImputationQ, tune.sigma = ImputationSigma),
                                                   keep.rownames = "ProteinGroup")
              SpectraAll <- data.table::merge.data.table(data.table::melt.data.table(SpectraRead, id.vars = "ProteinGroup", value.name = "Measured_LFQ", variable.name = "Sample"),
                                                         data.table::melt.data.table(SpectraImp, id.vars = "ProteinGroup", value.name = "Imputed_LFQ", variable.name = "Sample"))
              SpectraAll[, Data := data.table::fifelse(is.na(Measured_LFQ), "Imputed","Measured")]

              ImpPlot <- SpectraAll |> ggplot2::ggplot(ggplot2::aes(x = Imputed_LFQ, fill = Data)) + ggplot2::scale_y_continuous(expand = 0) +
                ggplot2::geom_density(adjust = 2, alpha = 0.5) + ggplot2::scale_fill_manual(values = c("Measured" = "#000", "Imputed" = "#C0C")) +
                ggplot2::labs(x = expression("Log"[2]~"LFQ Intensity"), y = "Density", fill = NULL) + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8))
              ImputedNAs <- data.table::melt.data.table(SpectraRead[ProteinGroup %in% Retained3$ProteinGroup], id.vars = "ProteinGroup", variable.name = "Sample")
              ImputedNAs <- ImputedNAs[is.na(value)]
              ImputedNAs[, value := NULL]
              ImputedNAs <- ImputedNAs[, .(vector = paste(Sample, collapse = ", ")), ProteinGroup]
              data.table::fwrite(ImputedNAs, "Imputed_LFQs.csv")
              SpectraAll[, Log2LFQ := data.table::fifelse(is.na(Measured_LFQ), Imputed_LFQ, Measured_LFQ)]
              SpectraAll <- SpectraAll[, .(ProteinGroup, Sample, Log2LFQ)]
            }
            # Modelling
            {
              ModelDesign <- stats::model.matrix(~0 + Condition, data = AbundanceMetadata)
              colnames(ModelDesign) <- gsub("Condition", "", colnames(ModelDesign))
              rownames(ModelDesign) <- AbundanceMetadata$Sample
              ModelDesign <- ModelDesign[, (c(which(grepl(CtlGroup, colnames(ModelDesign))), which(grepl(ExpGroup, colnames(ModelDesign)))))]

              ContrastMatrix <- matrix(nrow = 2, ncol = 1, dimnames = list("Levels" = colnames(ModelDesign), "Contrasts" = "comp"))
              ContrastMatrix[,1] <- c(-1,1)

              LimmaInput <- SpectraAll |> data.table::dcast(formula = ProteinGroup ~ Sample, value.var = "Log2LFQ") |>
                data.table::setcolorder(neworder = c("ProteinGroup", AbundanceMetadata$Sample))

              suppressMessages(ModelFit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(LimmaInput, ModelDesign), ContrastMatrix)))
              if(!is.finite(ModelFit$df.prior)){message("Warning: Limma Prior is Infinite")}
              MeanVarData <- data.table::data.table(ModelFit$genes, "Mean" = ModelFit$Amean, "Variance" = sqrt(ModelFit$sigma))
              MeanVarData[, Data := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Imputed", "Measured")]
              MeanVarPlot <- MeanVarData |> ggplot2::ggplot(ggplot2::aes(x = Mean, y = Variance)) +
                ggplot2::geom_point(colour = "#000") +  ggplot2::geom_point(data = MeanVarData[Data == "Imputed"], colour = "#C0C") +
                ggplot2::labs(x = "Mean Log2LFQ", y = "Variance", colour = NULL)

              LimmaOutput <- limma::topTable(ModelFit, coef=1, adjust.method = "BH", n=Inf) |>  data.table::data.table() |> data.table::setnames("logFC", "Log2FC")
              LimmaOutput <- Proteopedia::Add_ProteinInfo(LimmaOutput[order(abs(LimmaOutput$Log2FC), decreasing = T),], paste0(InputDirectory, "/report.protein_description.tsv.gz"))
              LimmaOutput[, Significance := data.table::fifelse(P.Value < 0.05 & Log2FC < 0, "Sig. Decrease", data.table::fifelse(P.Value < 0.05 & Log2FC > 0, "Sig. Increase", "None"))]
              LimmaOutput <- data.table::merge.data.table(LimmaOutput, LimmaInput, by = "ProteinGroup", all.x = T)
              LimmaOutput <- Proteopedia::Separate_Isoforms(LimmaOutput)
              suppressWarnings(LimmaOutput[, Imputed := data.table::fifelse(ProteinGroup %in% ImputedNAs$ProteinGroup, "Yes","No")])
              data.table::fwrite(LimmaOutput, file = "Limma_Output.csv")
            }
            # Plotting
            {
              LimmaVolcano <- LimmaOutput |> dplyr::arrange(desc(abs(t))) |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) + ggplot2::geom_point(alpha = 0.7, stroke = NA, size = 2) +
                ggplot2::scale_colour_manual("#000") + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), alpha = 0.1) +
                ggrepel::geom_text_repel(ggplot2::aes(label=ifelse(Gene %in% head(Gene,250), as.character(Gene),""))) +
                ggplot2::geom_vline(xintercept = round(mean(LimmaOutput$Log2FC, na.rm = T), digits = 3), linetype = "dashed", colour = "#000") +
                Proteopedia::Add_AbundanceAxes() + ggplot2::ggtitle("Limma") + ggplot2::theme(legend.title = ggplot2::element_blank())

              pdf("LimmaVolcanoes.pdf", height = 10, width = 14)
              print(LimmaVolcano)
              print(LimmaVolcano + ggplot2::geom_point(data = LimmaOutput[Imputed == "Yes"], colour = "#F00", alpha = 0.7, stroke = NA, size = 2) + ggplot2::ggtitle("Imputed Limma"))
              Proteopedia::Reset_Dev()

              AbundanceCorr <- LimmaOutput[, .(ProteinGroup, Isoforms, Log2FC)][, Measure := "0hr"] |> rbind(AbundanceData[, .(ProteinGroup, Isoforms, Log2FC)][, Measure := "Modelled"]) |>
                data.table::dcast(ProteinGroup+Isoforms ~ Measure, value.var = "Log2FC") |> ggplot2::ggplot(ggplot2::aes(x = `0hr`, y = Modelled)) + ggplot2::geom_point(stroke = NA) +
                Proteopedia::Add_XYLine() + Proteopedia::Add_Pearsons() + ggplot2::labs(x = expression("Log"[2]~"FC in Protein Abundance (at 0 hr)"),
                                                                                        y = expression(atop("Log"[2]~"FC in Protein", "Abundance (Modelled)")))

              pdf("LimmaQC_Plot.pdf", width = 18, height = 20)
              print(patchwork::free(PCAPlot) + CountsBar + patchwork::free(UpsetPlot) +
                      NormPlot + ImpPlot + patchwork::free(MeanVarPlot, type = "label") +
                      AbundanceCorr + patchwork::plot_layout(design = "AAAABBBB\nCCDDEEEE\nFFFFGGGG") +
                      patchwork::plot_annotation(tag_levels = list(c("A", "", "B", "C", "D", "E", "F", "G"))))
              Proteopedia::Reset_Dev()

            }
            data.table::fwrite(data.table::data.table("Experimental_Condition" = paste0(ExpGroup), "Control_Condition" = paste0(CtlGroup),
                                                      "Min_Precursors" = paste0(MinPrecursors), "Imputation_Q-value" = ImputationQ,
                                                      "Imputation_Sigma" = ImputationSigma), "Analysis_Parameters.csv")
          }
          message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Abundance Output Plots"))
          {
            setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
            pdf("AbundanceOutputPlots.pdf", width = 16, height = 10)
            print(AbundanceData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = Prop_NA)) +
                    ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                    Proteopedia::Add_AbundanceAxes() + ggrepel::geom_text_repel(data = AbundanceData[P.Value < 0.05], colour = "#000"))
            print(AbundanceData[, .(ProteinGroup, P.Value, adj.P.Val)] |> data.table::copy() |>
                    data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |>
                    data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Adjustment", value.name = "P") |>
                    ggplot2::ggplot(ggplot2::aes(x = P)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
                    ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = 0) +
                    ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
            Proteopedia::Reset_Dev()
          }
        }
        # Analyse Heavy Channel
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Isolating Heavy Channel Proteins"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          ProtLFQsInput_H <- ProtLFQsInput |> data.table::copy()
          if(MeanCentring){
            ProtLFQsInput_H[grepl(CtlGroup, Condition), LFQ_Mean := mean(log(LFQ_H),na.rm = T), .(ProteinGroup, Condition, Time)]
            ProtLFQsInput_H[, LFQ_Mean := mean(LFQ_Mean,na.rm = T), .(ProteinGroup, Time)]
            ProtLFQsInput_H[, Diff := log(LFQ_H) - LFQ_Mean]
            ProtLFQsInput_H[, MeanDiff := mean(Diff, na.rm = T), .(Time, Replicate, Condition)]
            ProtLFQsInput_H[, LFQ_H:= exp(log(LFQ_H) - MeanDiff)]
          }
          ProtLFQsInput_H <- ProtLFQsInput_H |> data.table::dcast(ProteinGroup ~ Sample,  value.var =  "LFQ_H") |> tibble::column_to_rownames("ProteinGroup") |> as.matrix()
          ProteinLFQs_H <- ProtLFQsInput_H |> as.data.frame() |> tibble::rownames_to_column("ProteinGroup") |> data.table::data.table()
          ProteinLFQs_H <- data.table::melt.data.table(ProteinLFQs_H, id.vars = "ProteinGroup", variable.name = "Sample", value.name = "Abundance") |>
            data.table::merge.data.table(Metadata)
          ProtWeights_H <- ProteinWeights[Channel =="Heavy"][,`:=`(Replicate = factor(Replicate), ProteinGroup = factor(ProteinGroup),
                                                              Time = factor(gsub("h", "", Time), levels = gsub("h", "", TimeLevels)),
                                                              Condition = factor(Condition, levels = ConditionLevels))]

          ModelProteins_H <- ProteinLFQs_H[,.(ProteinGroup, Condition, Time, Replicate, Abundance)]
          ModelProteins_H[, Replicate := factor(Replicate)]
          ModelProteins_H[Time == 0, Abundance := 0]
          ModelProteins_H <- na.omit(ModelProteins_H)
          ModelProteins_H[, N_Quant := .N, .(ProteinGroup, Condition, Time)] # Max. N_Quant is 3 (3 Replicates)
          ModelProteins_H[, N_QuantTotal := .N, ProteinGroup] # Max. N_QuantTotal is 36 (3 Replicates, 6 Timepoints, 2 Conditions: 3x6x2 = 36)

          HeavyAbunBoxplot <- ModelProteins_H |> ggplot2::ggplot(ggplot2::aes(x = factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")),
                                                                              y = Abundance, colour = Condition)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::labs(x = "Time", y = "Heavy Protein Abundance") +
            ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
            ggplot2::theme(legend.title = ggplot2::element_blank())
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Filtering By Heavy Channel Missingness & Variation"))
        {
          PreFiltReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Condition, Time, N_Quant)] |> unique()
          PreFiltTotalCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)] |> unique()

          ModelProteins_H <- ModelProteins_H[N_QuantTotal >= HeavyMinSamples*nrow(Metadata)]
          ModelProteins_H[, CV := sd(Abundance)/mean(Abundance), .(ProteinGroup, Condition, Time)]
          ModelProteins_H_CV <- ModelProteins_H[,.(MeanCV = mean(CV, na.rm =T)), ProteinGroup]
          ModelProteins_H <- ModelProteins_H[ProteinGroup %in% ModelProteins_H_CV[MeanCV <= MaxCV, ProteinGroup]]

          TotalReplicateData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, Condition, Time, N_Quant)][, Filtering := "Post-Filtering"] |> unique() |>
            rbind(PreFiltReplicateData[, Filtering := "Pre-Filtering"])
          TotalProteinCountData <- ModelProteins_H[Abundance != 0, .(ProteinGroup, N_QuantTotal)][, Filtering := "Post-Filtering"] |> unique() |>
            rbind(PreFiltTotalCountData[, Filtering := "Pre-Filtering"])
          TotalReplicateData[, Filtering := factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))]

          HeavyCompleteness <- TotalProteinCountData |> ggplot2::ggplot(ggplot2::aes(x = as.numeric(N_QuantTotal))) +
            ggplot2::geom_histogram() + ggplot2::labs(x = "No. Samples", y = "No. Proteins") +
            ggplot2::geom_vline(xintercept = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, linetype = "dashed", colour = "#F00") +
            ggplot2::facet_wrap(~factor(Filtering, levels = c("Pre-Filtering", "Post-Filtering"))) +
            ggplot2::annotate("rect", xmin = -Inf, xmax = ceiling(HeavyMinSamples*nrow(Metadata))-0.5, ymin = -Inf, ymax = Inf,
                              fill = "#F00", alpha = 0.3) +
            ggplot2::scale_x_continuous(limits = c(0, nrow(Metadata)), expand = 0) + ggplot2::scale_y_continuous(expand = 0)

          HeavyCounts <- TotalReplicateData |> ggplot2::ggplot(ggplot2::aes(x = N_Quant, fill = Condition)) +
            ggplot2::geom_bar(stat = "count", position = ggplot2::position_dodge()) +
            ggplot2::facet_grid(ggplot2::vars(Filtering), ggplot2::vars(factor(paste0(Time, "h"), levels = paste0(TimeLevels, "h")))) +
            ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette, name = "Condition") + ggplot2::labs(x = "Replicates", y = "No. Proteins") +
            ggplot2::theme(strip.text.y = ggplot2::element_text(size = 26),legend.title = ggplot2::element_blank(),
                           legend.position = "inside", legend.position.inside = c(0.2, 0.9))
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Filtering By Heavy Channel Monotonicity"))
        {
          MonotonicityData <- ModelProteins_H[,.(Mean_Abundance = mean(Abundance, na.rm = T)), .(ProteinGroup, Time, Condition)][(order(ProteinGroup, Condition, Time))]
          MonotonicityData[, CumSum := cummax(Mean_Abundance), .(ProteinGroup, Condition)]
          MonotonicityData[, Monotonic := CumSum == Mean_Abundance]
          MonotonicitySummary <- MonotonicityData[, .(N_Monotonic = sum(Monotonic)), "ProteinGroup"]
          MonotonicitySummary[, MonotonicityProp := N_Monotonic/(length(TimeLevels)*length(ConditionLevels))]

          MonotonicityPlot <- MonotonicitySummary |> ggplot2::ggplot(ggplot2::aes(x = MonotonicityProp)) +
            ggplot2::geom_histogram() + ggplot2::labs(x = "Prop. Monotonicity", y = "No. Proteins") +
            ggplot2::geom_vline(xintercept = HeavyMinMonotonicity-0.033, linetype = "dashed", colour = "#F00") +
            ggplot2::annotate("rect", xmin = -Inf, xmax = HeavyMinMonotonicity-0.033, ymin = -Inf, ymax = Inf, fill = "#F00", alpha = 0.3) +
            ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(expand = 0)

          ModelProteins_H <- ModelProteins_H[ProteinGroup %in% MonotonicitySummary[MonotonicityProp >= HeavyMinMonotonicity, ProteinGroup]]
          ModelProteins_H[, ProteinGroup := as.factor(ProteinGroup)]

          HeavyMissingnessData <- ModelProteins_H[, .(ProteinGroup, N_QuantTotal)] |> data.table::copy() |> unique()
          HeavyMissingnessData[, PropNAs := 1-(N_QuantTotal/nrow(Metadata))]

          # Processed Heavy Data PCA
          #All_PCA <- ModelProteins_H[, .(ProteinGroup, Condition, Time, Replicate, Abundance)] |> data.table::merge.data.table(Metadata[, .(Condition, Time, Replicate, Sample)]) |> tidyr::pivot_wider(id_cols = ProteinGroup, values_from = Abundance, names_from = Sample, values_fill = NA) |>
          #  tidyr::drop_na() |> data.frame(row.names = c("ProteinGroup")) |> t() |> stats::prcomp(scale. = T)
          #SummaryPCA <- summary(All_PCA)$importance
          #All_PCA <- data.table::data.table(All_PCA$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

          #HeavyProcessed_PCA <- All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
          #  ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
          #  ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC1"] * 100, 0), "%]", sep = ""),
          #                y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC2"] * 100, 0), "%]", sep = "")) +
          #  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
          #  All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, colour = Condition, shape = Replicate, label = paste0(Time, "h"))) +
          #  ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggrepel::geom_text_repel() +
          #  ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC3"] * 100, 0), "%]", sep = ""),
          #                y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance" , "PC4"] * 100, 0), "%]", sep = "")) +
          #  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank()) +
          #  patchwork::plot_annotation(title = "Processed Heavy Data")
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Running ", HeavyModel, " Modelling"))
        {
          KlossParameters <- LightModelParameters[, .(ProteinGroup, Time, `Time:Condition_Exp`)] |> data.table::copy()
          KlossParameters[, `:=`(Kloss_Ctl = -Time, Kloss_Exp = -`Time:Condition_Exp` - Time)]
          Run_ProteinNLS <- function(POI, Progress = NULL){
            tryCatch(
              expr = {
                POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProteinWeights[Channel == "Heavy", .(ProteinGroup, Condition, Replicate, Weight)], all.x = T)
                POIData <- POIData[order(Condition, Time)]
                POIData[,`:=`(Comparison = data.table::fifelse(Condition == CtlGroup, 0, 1))]
                POIData[, Time := as.numeric(paste(Time))]
                POIData[, VAR := var(Abundance), .(Condition, Time)]
                POIData[is.na(VAR), VAR := max(POIData$VAR, na.rm  = T)]
                # Define T0 Data
                T0Data <- POIData[, head(.SD,3), Condition]
                T0Data[,`:=`(Time = 0, Abundance = 0, Weight = min(POIData$Weight, na.rm =T), VAR = min(POIData$VAR, na.rm = T))]
                POIData <- POIData[Time != 0] |> rbind(T0Data)

                T0Data <- POIData[,head(.SD,1), by = Condition]
                T0Data[,`:=`(Time = 0, Abundance = 0)]
                POIData <- POIData |> rbind(T0Data)

                KsynStart <- sapply(ConditionLevels, function(COI){
                  mean(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Time)), na.rm =T)
                })

                AbundancePlateau <- sapply(ConditionLevels, function(COI){
                  mean(POIData[order(Time)][Condition == COI][,tail(.SD,3)]$Abundance,na.rm =T)
                })

                StartVals <- c(KsynStart, KsynStart / AbundancePlateau) # nls requires start values of approx. params
                # nlsLM more stable than nls
                POIFit <- minpack.lm::nlsLM(Abundance ~ (Ksyn_Ctl+(Ksyn_Exp*Comparison))/(Kloss_Ctl+(Kloss_Exp*Comparison)) * (1-exp(-(Kloss_Ctl+(Kloss_Exp*Comparison))*Time)), data = POIData,
                                            start = list(Ksyn_Ctl = StartVals[[1]],
                                                         Ksyn_Exp = StartVals[[2]]/10,
                                                         Kloss_Ctl = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Ctl], StartVals[[3]]),
                                                         Kloss_Exp = data.table::fifelse(UseLightKloss == T, KlossParameters[ProteinGroup == POI, Kloss_Exp], StartVals[[4]]/10)),
                                            weight = 1/POIData$Weight, control = nls.control(maxiter = 200, warnOnly = T))
                POIData <- POIData[, Fitted := predict(POIFit)]
                POIData <- Proteopedia::Add_ProteinInfo(POIData, paste0(InputDirectory, "/report.protein_description.tsv.gz"))

                FitSummary <- summary(POIFit)
                if(GenerateDataPlots){
                  pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
                  print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Condition)) + ggplot2::geom_point() +
                          ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), .(Condition, Time, ProteinGroup)],
                                             ggplot2::aes(y = Fitted, group = interaction(Condition)), linetype = "dashed") +
                          ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                          ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", POIData[ProteinGroup == POI, Gene], ")"),
                                        subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients[3,1],3)}, ksyn = {round(FitSummary$coefficients[1,1],3)}
                                                                {ExpGroup}: kloss = {round(FitSummary$coefficients[3,1]+ FitSummary$coefficients[4,1],3)}, ksyn = {round(FitSummary$coefficients[1,1]+ FitSummary$coefficients[2,1],3)}')) +
                          ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
                  Proteopedia::Reset_Dev()
                }
                message(POI, ": Model Applied ", Progress)
                return(list(POIModel = POIFit, FittedData = POIData))
              },
              error = function(e){
                message(POI, ": Error Caught ", Progress)
              })
          }
          Run_ProteinNLME <- function(POI, Progress = NULL){
            tryCatch(
              expr = {
                POIData <- ModelProteins_H[ProteinGroup == POI] |> data.table::copy() |> data.table::merge.data.table(ProtWeights_H[, .(ProteinGroup, Condition, Replicate, Time, Weight)], all.x = T)
                POIData[, Acquisition := data.table::fifelse(as.numeric(paste(Time)) > max(TimeLevels)/2, "Late", "Early")]
                POIData[is.na(Acquisition), Acquisition := "Early"]
                POIData[, Time := as.numeric(paste(Time))]
                POIData[Time == 0, Weight := max(POIData$Weight, na.rm = T)]
                # Define T0 Data
                T0Data <- POIData[,head(.SD,1), Condition]
                T0Data[,`:=`(Time = 0, Abundance = 0)]
                POIData <- POIData[Time != 0] |> rbind(T0Data)

                KsynStart <- sapply(ConditionLevels, function(COI){
                  mean(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Abundance/as.numeric(paste0(POIData[order(Time)][Time != 0 & Condition == COI][,head(.SD,2)]$Time)), na.rm =T)
                })

                AbundancePlateau <- sapply(ConditionLevels, function(COI){
                  mean(POIData[order(Time)][Condition == COI][,tail(.SD,3)]$Abundance,na.rm =T)
                })

                StartVals <- c(KsynStart, KsynStart / AbundancePlateau)

                POIFit <- nlme::nlme(Abundance ~ (Ksyn/Kloss)*(1 - exp(-Kloss*as.numeric(paste(Time)))), data = POIData, fixed = list(Ksyn ~ Condition, Kloss ~ Condition),
                                     random = list(Acquisition = nlme::pdDiag(Ksyn+Kloss ~1)),
                                     start = StartVals, weights = nlme::varFixed(~Weight))
                POIData <- POIData[, Fitted := predict(POIFit)]

                FitSummary <- summary(POIFit)
                if(GenerateDataPlots){
                  pdf(paste0(POI, "_HeavyPlot.pdf"), width = 16, height = 12)
                  print(POIData |> ggplot2::ggplot(ggplot2::aes(x = Time, y = Abundance, colour = Condition)) + ggplot2::geom_point() +
                          ggplot2::geom_line(data = POIData[,.(Fitted = (mean(Fitted,na.rm =T))), .(Condition, Time, ProteinGroup)],
                                             ggplot2::aes(y = Fitted, group = interaction(Condition)), linetype = "dashed") +
                          ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                          ggplot2::labs(x = "Time (hours)", y = "Heavy Protein Abundance", title = paste0(POI, " (", ProteinInfo[ProteinGroup == POI, Gene], ")"),
                                        subtitle = glue::glue('{CtlGroup}: kloss = {round(FitSummary$coefficients$fixed[3],3)}, ksyn = {round(FitSummary$coefficients$fixed[1],3)}
                                                                  {ExpGroup}: kloss = {round(FitSummary$coefficients$fixed[3] + FitSummary$coefficients$fixed[4],3)}, ksyn = {round(FitSummary$coefficients$fixed[1] + FitSummary$coefficients$fixed[2],3)}')) +
                          ggplot2::theme(plot.title = ggplot2::element_text(size = 26), legend.title = ggplot2::element_blank()))
                  Proteopedia::Reset_Dev()
                }
                message(paste0(POI, ": Model Applied ", Progress))
                return(list(POIModel = POIFit, FittedData = POIData))
              },
              error = function(e){
                message(paste0(POI, ": Error Caught ", Progress))
              }
            )
          }

          if(GenerateDataPlots){
            if (dir.exists("HeavyPlots")) {
              unlink("HeavyPlots", recursive = T)
            }
            dir.create("HeavyPlots", showWarnings = T)
            setwd("HeavyPlots")
          }

          HeavyModelledData <- data.table::data.table()
          HeavyModelParameters <- data.table::data.table()
          ProteinSigmas <- data.table::data.table()
          if(HeavyModel == "NLS"){
            for(POI in levels(ModelProteins_H$ProteinGroup)){
              NLSOutput <- Run_ProteinNLS(POI, paste0("(", round(which(levels(ModelProteins_H$ProteinGroup) == POI)/length(levels(ModelProteins_H$ProteinGroup)), digits = 2)*100, "%)"))
              HeavyModelledData <- HeavyModelledData |> rbind(NLSOutput$FittedData)
              if(!is.null(NLSOutput$POIModel$convInfo$isConv)){
                if(NLSOutput$POIModel$convInfo$isConv){
                  POIModelSummary <- summary(NLSOutput$POIModel)$coefficients |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "t_Value", "P.Value"))
                  ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = sigma(NLSOutput$POIModel)))
                  HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
                }
              }
            }
          } else if(HeavyModel == "NLME"){
            for(POI in levels(ModelProteins_H$ProteinGroup)){
              NLMEOutput <- Run_ProteinNLME(POI, paste0("(", round(which(levels(ModelProteins_H$ProteinGroup) == POI)/length(levels(ModelProteins_H$ProteinGroup)), digits = 2)*100, "%)"))
              HeavyModelledData <- HeavyModelledData |> rbind(NLMEOutput$FittedData)
              if(class(NLMEOutput$POIModel)[1] == "nlme"){
                POIModelSummary <- summary(NLMEOutput$POIModel)$tTable |> data.table::data.table(keep.rownames = T) |> data.table::setnames(c("Effect", "Estimate", "SE", "DF", "t_Value", "P.Value"))
                ProteinSigmas <- ProteinSigmas |> rbind(data.table::data.table(ProteinGroup = POI, Sigma = NLMEOutput$POIModel$sigma))
                HeavyModelParameters <- HeavyModelParameters |> rbind(POIModelSummary[,`:=`(ProteinGroup = POI)])
                HeavyModelParameters[, Effect := gsub(".Condition.*", "_Exp", gsub("..(Intercept.)", "_Ctl", HeavyModelParameters$Effect))]
              }
            }
          }
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          HeavyModelParameters[, adj.P.Val := p.adjust(P.Value, 'BH'), by = 'Effect']
          HeavyModelParameters <- Proteopedia::Add_ProteinInfo(HeavyModelParameters, paste0(InputDirectory, "/report.protein_description.tsv.gz")) |>
            data.table::merge.data.table(HeavyMissingnessData, by = "ProteinGroup")

          HeavyParameters <- HeavyModelParameters |> data.table::copy()
          HeavyParameters[stringr::str_detect(Effect,'_Ctl$'), Ctl_Value := Estimate]
          HeavyParameters[, Parameter := gsub("Kloss", "KlossH", stringr::str_remove(Effect,'_.*'))]
          HeavyParameters[, Ctl_Value := mean(Ctl_Value, na.rm = T), .(ProteinGroup, Parameter)]
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
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Calculating Heavy Mean Absolute Percentage Errors (MAPEs)"))
        {
          HeavyMAPEData <- HeavyModelledData[Time != 0, .(MeanAbundance = mean(Abundance, na.rm = T), MeanFitted = mean(Fitted, na.rm = T)), .(ProteinGroup, Condition, Time)]
          HeavyMAPEData <- HeavyMAPEData[, .(MAPE = MetricsWeighted::mape(MeanAbundance, MeanFitted)), .(ProteinGroup, Condition, Time)]
          HeavyMAPEData[, MAPEBin := data.table::fifelse(MAPE > 50, "MAPE > 50", data.table::fifelse(MAPE > 25, "MAPE > 25",
                                                                                                     data.table::fifelse(MAPE > 10, "MAPE > 10", "MAPE ≤ 10")))]
          HeavyMAPEData[, Channel := "Heavy"]
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Heavy Data Files"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          data.table::fwrite(ModelProteins_H, "HeavyInputLFQs.csv")
          data.table::fwrite(HeavyModelParameters, "HeavyModelOutput.csv")
          data.table::fwrite(Proteopedia::Separate_Isoforms(HeavyParameters), "HeavyParameters.csv")
          data.table::fwrite(HeavyModelledData, "HeavyModelledData.csv")
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Heavy Output Plots"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          pdf("HeavyOutputPlots.pdf", width = 16, height = 10)
          print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05]) +
                  Proteopedia::Add_KlossAxes())
          print(HeavyParameters[Parameter == "KlossH"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "KlossH" & P.Value < 0.05], colour = "#000", max.overlaps = 10) +
                  Proteopedia::Add_KlossAxes(scale = "Log2FC"))
          print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Difference, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05]) +
                  Proteopedia::Add_KsynAxes())
          print(HeavyParameters[Parameter == "Ksyn"] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value), label = Gene, colour = PropNAs)) +
                  ggplot2::geom_point(stroke = NA) + ggplot2::scale_colour_viridis_c(name = "Prop. NAs") + Proteopedia::Add_NotSigBox() +
                  ggrepel::geom_text_repel(data = HeavyParameters[Parameter == "Ksyn" & P.Value < 0.05], colour = "#000", max.overlaps = 10) +
                  Proteopedia::Add_KsynAxes(scale = "Log2FC"))
          print(HeavyParameters[Parameter == "KlossH", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
                  data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
                  data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "KlossH") |>
                  ggplot2::ggplot(ggplot2::aes(x = KlossH, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
                  ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)"), y = "Density of Proteins") +
                  ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                            legend.position.inside = c(0.8, 0.8)))
          print(HeavyParameters[Parameter == "Ksyn", .(ProteinGroup, Ctl_Value, Exp_Value)] |> data.table::copy() |>
                  data.table::setnames(c("Ctl_Value", "Exp_Value"), c(CtlGroup, ExpGroup)) |>
                  data.table::melt.data.table(id.vars = "ProteinGroup", variable.name = "Condition", value.name = "Ksyn") |>
                  ggplot2::ggplot(ggplot2::aes(x = Ksyn, fill = Condition)) + ggplot2::geom_density(alpha = 0.7) +
                  ggplot2::scale_fill_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~") (Heavy Channel)"), y = "Density of Proteins") +
                  ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "inside",
                                                            legend.position.inside = c(0.8, 0.8)))
          print(HeavyParameters[, .(ProteinGroup, P.Value, adj.P.Val, Parameter)] |> data.table::copy() |>
                  data.table::setnames(c("P.Value", "adj.P.Val"), c("Raw", "Adjusted")) |>
                  data.table::melt.data.table(id.vars = c("ProteinGroup", "Parameter"), variable.name = "Adjustment", value.name = "P") |>
                  ggplot2::ggplot(ggplot2::aes(x = P, fill = Parameter)) + ggplot2::geom_histogram() + ggplot2::facet_wrap(~Adjustment, scales = "free_y") +
                  ggplot2::labs(y = "No. Proteins", x = "P-Value") + ggplot2::scale_x_continuous(expand = c(0, 0))+
                  ggplot2::scale_y_continuous(expand = 0) + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette) +
                  ggplot2::annotate("rect", xmin = -Inf, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0F0", alpha = 0.3))
          Proteopedia::Reset_Dev()
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting Analysis Summary"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          AnalysisSummary <- data.table::data.table(Parameter = c("Light kloss", "Heavy kloss", "ksyn", "0hr Abundance"),
                                                    N_Proteins = c(nrow(LightParameters), nrow(HeavyParameters[Parameter == "KlossH"]) ,
                                                                   nrow(HeavyParameters[Parameter == "Ksyn"]), nrow(AbundanceData)),
                                                    Model = c("Limma", HeavyModel, HeavyModel, "Limma"))
          AnalysisSummary[, Prop_Proteins := N_Proteins/length(ProtLFQsInput[, ProteinGroup] |> unique())]
          data.table::fwrite(AnalysisSummary, "Analysis_Summary.csv")
        }
        message(paste0(ExpGroup, " vs. ", CtlGroup, ": Exporting QC Data, Plots & Comparisons"))
        {
          setwd(paste0(InputDirectory, "/", ExpGroup, "_vs_", CtlGroup, "_Output"))
          MAPEData <- Proteopedia::Separate_Isoforms(LightMAPEData |> rbind(HeavyMAPEData))
          MAPEData[, MAPEBin := factor(MAPEBin, levels = c("MAPE ≤ 10", "MAPE > 10", "MAPE > 25", "MAPE > 50"))]
          data.table::fwrite(MAPEData, "MAPEData.csv")
          WideMAPEData <- MAPEData |> data.table::dcast(ProteinGroup+Isoforms+Channel+Time ~ Condition, value.var = "MAPE")
          WideMAPEData[, DiffMAPE := get(ExpGroup) - get(CtlGroup)]

          pdf("MAPEPlots.pdf", width = 16, height = 12)
          print(MAPEData[Channel == "Light"] |> ggplot2::ggplot(ggplot2::aes(x = MAPE, fill = MAPEBin)) + ggplot2::geom_histogram() +
                  ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(16, length(Proteopedia::ThermalPalette), length.out = 4)]) +
                  ggplot2::facet_grid(ggplot2::vars(Condition), ggplot2::vars(factor(paste(Time, "h"), levels = paste(TimeLevels, "h"))), scales = "free") +
                  ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "Mean Absolute Percentage Error (MAPE)", y = "No. Light Proteins") +
                  ggplot2::theme(legend.title = ggplot2::element_blank()))
          print(MAPEData[Channel == "Heavy"] |> ggplot2::ggplot(ggplot2::aes(x = MAPE, fill = MAPEBin)) + ggplot2::geom_histogram() +
                  ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(16, length(Proteopedia::ThermalPalette), length.out = 4)]) +
                  ggplot2::facet_grid(ggplot2::vars(Condition), ggplot2::vars(paste(Time, "h")), scales = "free") +
                  ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = "Mean Absolute Percentage Error (MAPE)", y = "No. Light Proteins") +
                  ggplot2::theme(legend.title = ggplot2::element_blank()))
          print(WideMAPEData |> ggplot2::ggplot(ggplot2::aes(x = DiffMAPE, fill = Channel)) + ggplot2::geom_density(alpha = 0.3) +
                  Proteopedia::Add_Isotope_Fill() + ggplot2::facet_wrap(~paste(Time, "h"), scales = "free", nrow = 1) +
                  ggplot2::labs(x = paste0("Difference in MAPE (", ExpGroup, " - ", CtlGroup, ")"), y = "No. Proteins") +
                  ggplot2::scale_y_continuous(expand = 0))
          Proteopedia::Reset_Dev()

          CVPlot <- ProteinLFQs_L_CV[, Channel := "Light"] |> rbind(ModelProteins_H_CV[, Channel := "Heavy"]) |>
            ggplot2::ggplot(ggplot2::aes(x = MeanCV)) + ggplot2::geom_histogram() +
            ggplot2::facet_wrap(~forcats::fct_rev(Channel), nrow = 3, strip.position = "right") +
            ggplot2::labs(x = "Mean Protein Variation", y = "No. Proteins") + ggplot2::geom_vline(xintercept = MaxCV, linetype = "dashed", colour = "#F00") +
            ggplot2::annotate("rect", xmin = MaxCV, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F00", alpha = 0.3) +
            ggplot2::scale_x_continuous(limits = c(0, 1), expand = 0) + ggplot2::scale_y_continuous(expand = 0)

          pdf("QC_Plots.pdf", width = 18, height = 20)
          print(patchwork::free(LightAbunBoxplot + HeavyAbunBoxplot + patchwork::plot_layout(guides = "collect")) +
                  ProteinWeights_PCA + ProteinWeight_Density + CVPlot +
                  LightMeanVarPlot + MonotonicityPlot + patchwork::free(HeavyCounts) + HeavyCompleteness) +
            patchwork::plot_layout(design = "AAAAAAAAAA\nBBBBBBCCCC\nEEEEEDDDDD\nFFFFFDDDDD") +
                  patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "" ,"H", "")))
          Proteopedia::Reset_Dev()

          CorrData <- HeavyParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)] |>
            rbind(LightParameters[, .(ProteinGroup, Gene, Ctl_Value, Exp_Value, Parameter)]) |>
            data.table::melt.data.table(id.vars = c("ProteinGroup", "Gene", "Parameter"),
                                        variable.name = "Condition", value.name = "Measure")
          CorrData[, Condition := data.table::fifelse(grepl("Ctl", Condition), CtlGroup, ExpGroup)]
          CorrData <- CorrData |> data.table::dcast(ProteinGroup+Gene+Condition ~ Parameter, value.var = "Measure")

          pdf("CorrelationPlots.pdf", width = 16, height = 12)
          print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = KlossL, y = KlossH, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
                  Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::labs(x = expression("Rate of Turnover (k"[loss]~") (Light Channel)"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
                  ggplot2::theme(legend.title = ggplot2::element_blank()))
          print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossH, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
                  Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Heavy Channel)")) +
                  ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
          print(CorrData |> ggplot2::ggplot(ggplot2::aes(x = Ksyn, y = KlossL, colour = Condition)) + ggplot2::geom_point(stroke = NA) +
                  Proteopedia::Add_Pearsons(T) + Proteopedia::Add_XYLine("#999") + ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) +
                  ggplot2::labs(x = expression("Rate of Synthesis (k"[syn]~")"), y = expression("Rate of Turnover (k"[loss]~") (Light Channel)")) +
                  ggplot2::scale_x_log10() + ggplot2::theme(legend.title = ggplot2::element_blank()))
          Proteopedia::Reset_Dev()
        }
      }
    }
  }
  Proteopedia::End_Timer(start.time)
}
#' @export
Map_TimecourseSILAC_Proteins <- function(InputDirectory, SubsetColour = "#F00"){
  set.seed(123)
  start.time <- Sys.time()
  for(Parameter in c("KlossL", "KlossH", "Ksyn")){
    message(Parameter, ": Loading Limma File")
    {
      setwd(InputDirectory)
      if(Parameter == "KlossL"){LimmaData <- data.table::fread("LightParameters.csv")
      } else if(Parameter == "KlossH"){LimmaData <- data.table::fread("HeavyParameters.csv")[Parameter == "KlossH"]
      } else {LimmaData <- data.table::fread("HeavyParameters.csv")[Parameter == "Ksyn"]}

      if (dir.exists(paste0(getwd(),"/",Parameter,"_MappingOutput"))){unlink(paste0(getwd(),"/",Parameter,"_SubsetOutput"), recursive = T)}
      dir.create(paste0(getwd(),"/",Parameter,"_MappingOutput"), showWarnings = T)
      setwd(paste0(getwd(),"/",Parameter,"_MappingOutput"))

      LimmaData |> data.table::setnames(c(colnames(LimmaData)[grepl("protein.*group", ignore.case = T, colnames(LimmaData))],
                                          colnames(LimmaData)[grepl("^p.*val.*", ignore.case = T, colnames(LimmaData)) & !grepl("adj", ignore.case = T, colnames(LimmaData))],
                                          colnames(LimmaData)[grepl("p.*val.*", ignore.case = T, colnames(LimmaData)) & grepl("adj", ignore.case = T, colnames(LimmaData))],
                                          colnames(LimmaData)[grepl("Gene", ignore.case = T, colnames(LimmaData)) & !grepl("group", ignore.case = T, colnames(LimmaData))]),
                                        c("ProteinGroup", "P.Value", "adj.P.Val", "Gene"), skip_absent = T)
      if(length(LimmaData$ProteinGroup[grepl("\\-", LimmaData$ProteinGroup)]) > 0){
        LimmaData <- Proteopedia::Separate_Isoforms(LimmaData)
      }
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
              ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Proteopedia::Add_Pearsons() +
              ggplot2::annotate("label", x = mean(SubsetData[, Subset], na.rm = T), y = min(SubsetData[, Log2FC], na.rm = T) * 0.93,
                                label = paste0("P-Value: ", ifelse(pval < 0.01, formatC(pval, format = "e", digits = 2), round(pval, digits = 2))), size = 6) +
              ggplot2::labs(x = paste0("Protein ", gsub("_", " ", colnames(LimmaData)[ColIndex])), y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) +
              ggside::geom_xsidedensity() + Proteopedia::Clean_SideDensities()
          } else {
            TrendPlot <- ggplot2::ggplot(SubsetData, ggplot2::aes(x = Subset, y = Log2FC)) +
              ggplot2::geom_smooth(method = "lm", alpha = 0.1) + Proteopedia::Add_Pearsons() +
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
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") +
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
              ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") +
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
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") + ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
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
            ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") + ggplot2::geom_point(data = LimmaData[Deg_Profile == Deg_Type], colour = SubsetColour) +
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
      ChromosomeSummary <- Calculate_WilcoxonByVar(InputData = LimmaData, Category = "Chromosome", Measure = "Log2FC") |>
        data.table::merge.data.table(LimmaData[, .(MeanLog2FC = mean(Log2FC, na.rm = T), Mapped = .N), by = Chromosome]) |>
        data.table::merge.data.table(Proteopedia::Proteopedia[, .(Total = .N), by = Chromosome])
      ChromosomeSummary[, `:=`(Coverage = Mapped/Total, Chromosome = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")))]
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
      ChromosomeCols <- rep(Proteopedia::NiceColourPalette, length.out = length(c(seq(1:22), "X")))
      names(ChromosomeCols) <- c(seq(1:22), "X")
      ChromosomeDotData <- merge(LimmaData[Chromosome %!in% c("X/Y", "Y", "M", "Unmapped")],
                                 Proteopedia::Proteopedia[, .(ProteinGroup, Gene, Chromosome)], all = T, by = c("ProteinGroup", "Gene", "Chromosome"))
      ChromosomeDotData$Log2FC[is.na(ChromosomeDotData$Log2FC)] <- 0
      ChromosomeDotData <- merge(ChromosomeDotData, data.table::data.table(Colour = ChromosomeCols, Chromosome = names(ChromosomeCols)), by = "Chromosome")
      ChromosomeDotData[, `:=`(Colour, data.table::fifelse(Log2FC == 0, "#FFF", Colour))]

      if(Parameter == "Ksyn"){
        ChromosomeDotplot <- ggplot2::ggplot(ChromosomeDotData, ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) +
          ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "#00F") +
          ggplot2::scale_colour_identity() + ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), label = data.table::fifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))),
                                                                colour = "#000", angle = 90) +
          ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "#F00") + ggplot2::scale_x_continuous(expand = 0) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + ggplot2::guides(fill = "none")
        ChromosomeBoxplot <- ggplot2::ggplot(LimmaData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Log2FC)) +
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop * 1.25, label = SigSymbol), colour = "#000") +
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom * 1.25, label = Buffering), colour = "#000",
                             size = 4) + ggplot2::geom_boxplot(colour = "#000", fill = NA, outliers = F, alpha = 0.7) +
          ggplot2::coord_cartesian(ylim = c(WhiskerTop * 1.5, WhiskerBottom * 1.5)) + ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "#F00") +
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Synthesis Rate (k"[syn]~")")) + ggplot2::guides(fill = "none")
      } else {
        ChromosomeDotplot <- ggplot2::ggplot(ChromosomeDotData, ggplot2::aes(x = OrderID, y = Log2FC, colour = Colour)) +
          ggplot2::geom_point(stroke = NA, alpha = 0.7) + ggplot2::geom_vline(xintercept = ChromosomeBorders$Upper, colour = "#00F") +
          ggplot2::scale_colour_identity() + ggplot2::geom_text(data = ChromosomeBorders, ggplot2::aes(x = Midpoint, y = I(0.05), label = data.table::fifelse(Chromosome != "Unmapped", paste0("Chr", Chromosome), paste0(Chromosome))),
                                                                colour = "#000", angle = 90) +
          ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "#F00") + ggplot2::scale_x_continuous(expand = 0) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + ggplot2::guides(fill = "none")
        ChromosomeBoxplot <- ggplot2::ggplot(LimmaData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Log2FC)) +
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerTop * 1.25, label = SigSymbol), colour = "#000") +
          ggplot2::geom_text(data = ChromosomeSummary, ggplot2::aes(x = Chromosome, y = WhiskerBottom * 1.25, label = Buffering), colour = "#000", size = 4) +
          ggplot2::geom_boxplot(colour = "#000", fill = NA, outliers = F, alpha = 0.7) +
          ggplot2::coord_cartesian(ylim = c(WhiskerTop * 1.5, WhiskerBottom * 1.5)) + ggplot2::geom_hline(yintercept = log2(3/2), linetype = "dashed", colour = "#F00") +
          ggplot2::labs(x = "Chromosome", y = expression("Log"[2] ~ "FC in Protein Turnover Rate (k"[loss]~")")) + ggplot2::guides(fill = "none")
      }

      ChrGroupingData <- LimmaData[, .(Sig_N = .N), .(Significance, Chromosome)] |> data.table::merge.data.table(LimmaData[, .(Total_N = .N), by = Chromosome])
      ChrGroupingData[, Prop := Sig_N/Total_N]

      ChromosomeSigBar <- ggplot2::ggplot(ChrGroupingData[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = rev(c(seq(1:22), "X", "X/Y", "Y", "M"))),
                                                                                                  y = Prop, fill = factor(Significance, levels = c("Sig. Increase", "None", "Sig. Decrease")))) +
        ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::scale_fill_manual(values = c(`Sig. Decrease` = "#02F", None = "#FFF", `Sig. Increase` = "#F10"), name = "Fold Change") +
        ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~1 - .), expand = c(0, 0)) + ggplot2::geom_hline(yintercept = seq(0.1, 0.9, by = 0.1), linetype = "dotted") +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") + ggplot2::labs(x = "Chromosome", y = "Proportion of Proteins") + ggplot2::coord_flip()

      pdf(paste0(Parameter, "_ChromosomePlots.pdf"), width = 18, height = 8)
      print(ggplot2::ggplot(ChromosomeSummary[Chromosome != "Unmapped"], ggplot2::aes(x = factor(Chromosome, levels = c(seq(1:22), "X", "X/Y", "Y", "M")), y = Coverage * 100, fill = Coverage * 100)) + ggplot2::geom_bar(stat = "identity") +
              ggplot2::scale_fill_viridis_c(guide = "none") + ggplot2::labs(x = "Chromosome", y = "Protein Coverage (%)"))
      print(ChromosomeDotplot)
      print(ChromosomeBoxplot)
      print(ChromosomeSigBar)
      Proteopedia::Reset_Dev()

      if(dir.exists(paste0(getwd(), "/ChromosomeVolcanos"))) {unlink(paste0(getwd(), "/ChromosomeVolcanos"), recursive = T)}
      dir.create(paste0(getwd(), "/ChromosomeVolcanos"), showWarnings = T)
      setwd(paste0(getwd(), "/ChromosomeVolcanos"))
      for(i in levels(LimmaData$Chromosome)){
        message(paste0(Parameter, ": Analysing Chr", i, " Proteins"))
        MeanLog2FC <- round(mean(LimmaData[Chromosome == i, Log2FC]), digits = 3)
        if(Parameter == "Ksyn"){
          pdf(paste0("Chr", gsub("/", "", i), "_Volcano.pdf"), width = 12, height = 8)
          print(ggplot2::ggplot(LimmaData[Chromosome == i], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
                  ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(-log10(P.Value) > -log10(0.05), as.character(Gene), ""))) +
                  Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "#000", alpha = 0.7) +
                  Proteopedia::Add_KsynAxes(scale = "Log2FC"))
          Proteopedia::Reset_Dev()
        } else {
          pdf(paste0("Chr", gsub("/", "", i), "_Volcano.pdf"), width = 12, height = 8)
          print(ggplot2::ggplot(LimmaData[Chromosome == i], ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
                  ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = ifelse(-log10(P.Value) > -log10(0.05), as.character(Gene), ""))) +
                  Proteopedia::Add_NotSigBox() + ggplot2::geom_vline(xintercept = MeanLog2FC, linetype = "dashed", colour = "#000", alpha = 0.7) +
                  Proteopedia::Add_KlossAxes(scale = "Log2FC"))
          Proteopedia::Reset_Dev()
        }
      }
    }
    message(Parameter, ": Analysing Protein Complexes")
    {
      setwd(paste0(InputDirectory,"/",Parameter,"_MappingOutput"))
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
########### MS Analysis: TPP Analysis Functions ###############################################################################################################################################
#' @export
Analyse_TPP_Proteins <- function(InputDirectory, CtlGroup, ExpGroup, LFQChannel, PropMissingness = 0.5, Log2FCvs = "Max",
                                 ProteotypicFiltering = F, ExcludedSamples = NULL){
  set.seed(123)
  start.time <- Sys.time()
  message("Loading & Formatting Data")
  {
    setwd(InputDirectory)
    Metadata <- data.table::fread("Sample_Metadata.csv")
    ProteinData <- data.table::fread(list.files(pattern = "DIANN_Output.csv.gz")[1])[Sample %!in% ExcludedSamples] |> data.table::merge.data.table(Metadata)
    BCAData <- data.table::fread("BCAData.csv")

    if(LFQChannel == "Heavy"){
      if(dir.exists(paste0(getwd(), "/Heavy_Output"))){unlink(paste0(getwd(), "/Heavy_Output"), recursive = T)}
      dir.create(paste0(getwd(), "/Heavy_Output"), showWarnings = T)
      setwd(paste0(getwd(), "/Heavy_Output"))
      ProteinData |> data.table::setnames("LFQ_H", "LFQ", skip_absent = T)
    } else if(LFQChannel == "Light"){
      if(dir.exists(paste0(getwd(), "/Light_Output"))){unlink(paste0(getwd(), "/Light_Output"), recursive = T)}
      dir.create(paste0(getwd(), "/Light_Output"), showWarnings = T)
      setwd(paste0(getwd(), "/Light_Output"))
      ProteinData |> data.table::setnames("LFQ_L", "LFQ", skip_absent = T)
    } else {
      if(dir.exists(paste0(getwd(), "/Total_Output"))){unlink(paste0(getwd(), "/Total_Output"), recursive = T)}
      dir.create(paste0(getwd(), "/Total_Output"), showWarnings = T)
      setwd(paste0(getwd(), "/Total_Output"))
      ProteinData |> data.table::setnames("LFQ_T", "LFQ", skip_absent = T)
    }
    ProteinData[, DetectLevel := mean(head(LFQ, 200), na.rm = T), Sample]
    TempLevels <- as.character(sort(as.numeric(unique(Metadata$Temp))))
    Metadata[, Temp := factor(Temp, levels = TempLevels)]

    Metadata[, Group := gsub("_\\d+C", "", Condition)]

    PCAData <- ProteinData[, .(ProteinGroup, Sample, LFQ)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "LFQ") |> tidyr::drop_na() |>
      data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)

    pdf("RawLFQ_PCA.pdf", width = 12, height = 6)
    print(patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Group, colour = Temp)) +
                            ggplot2::geom_point(size = 6, stroke = NA) +
                            ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                            ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                            ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                          y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                            All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Group, colour = Temp)) +
                            ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                            ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                            ggplot2::labs(x = paste("PC3 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC3"]*100,0), "%]", sep = ""),
                                          y = paste("PC4 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC4"]*100, 0), "%]", sep = "")) +
                            ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.title = ggplot2::element_blank())))
    Proteopedia::Reset_Dev()
  }
  message("Performing Median Normalisation")
  {
    ProteinData[, Log2LFQ := log2(LFQ)]
    SpectraMean <- median(ProteinData$Log2LFQ, na.rm = T)
    ProteinData[, `:=`(NormLog2LFQ, Log2LFQ - median(Log2LFQ, na.rm = T) + SpectraMean), Sample]
    ProteinDataWide <- ProteinData |> data.table::dcast(ProteinGroup ~ Sample, value.var = "NormLog2LFQ") |> tibble::column_to_rownames("ProteinGroup")

    ProteinCVs <- ProteinData[, .(CV = Proteopedia::Calculate_CV(2^Log2LFQ), Stage = "Pre-Normalisation"), .(ProteinGroup, Condition)] |>
      rbind(ProteinData[, .(CV = Proteopedia::Calculate_CV(2^NormLog2LFQ), Stage = "Post-Normalisation"), .(ProteinGroup, Condition)])

    pdf("ProteinNormalisationCVs.pdf", width = 16, height = 12)
    print(ProteinCVs |> data.table::merge.data.table(Metadata[, .(Condition, Temp, Group)] |> unique(), by = "Condition") |>
            ggplot2::ggplot(ggplot2::aes(x = Temp, y = CV, colour = Group)) + ggplot2::geom_boxplot(outliers = F) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggplot2::facet_wrap(~forcats::fct_rev(Stage)) +
            ggplot2::labs(x = Temp.~degree~C, y = "Variation (%)"))
    Proteopedia::Reset_Dev()

    PCAData <- ProteinData[, .(ProteinGroup, Sample, NormLog2LFQ)] |>
      data.table::dcast(ProteinGroup ~ Sample, value.var = "NormLog2LFQ") |> tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |>
      t() |> stats::prcomp(scale. = T)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Group, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Group, colour = Temp)) +
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
    ProteinData[, Group := gsub("_\\d+C", "", Condition)]
    Missingness <- ProteinData[, .(NAs = sum(is.na(NormLog2LFQ))), .(ProteinGroup, Temp, Group)] |> data.table::dcast(ProteinGroup+Temp~Group, value.var = "NAs")
    Missingness[, NAsDiff := get(ExpGroup) - get(CtlGroup)]
    Missingness <- Missingness[, .(MeanNAsDiff = mean(NAsDiff)), ProteinGroup]

    MissingnessMap <- is.na(ProteinDataWide) |> apply(2, as.numeric)
    rownames(MissingnessMap) <- rownames(ProteinDataWide)
    N_MissingProteins <- MissingnessMap |> matrixStats::rowSums2()
    MissingnessMap <- MissingnessMap[N_MissingProteins <= nrow(Metadata)*PropMissingness,] # Pre-Imputation Missingness
    MissingnessHeatmap <- pheatmap::pheatmap(MissingnessMap, show_rownames = F, show_colnames = F, annotation_col = Metadata |> tibble::column_to_rownames("Sample"))
    ProteinGroups <- data.table::data.table(ProteinGroup = rownames(MissingnessMap), Cluster = as.factor(cutree(MissingnessHeatmap$tree_row, k = 5)))
    pheatmap::pheatmap(MissingnessMap, show_rownames = F, show_colnames = F, color = c("#090", "#000"), legend = F,
                       annotation_col = Metadata[, .(Sample, Group, Temp, Replicate)] |> tibble::column_to_rownames("Sample"),
                       annotation_row = ProteinGroups |> tibble::column_to_rownames("ProteinGroup"), filename = "MissingnessHeatmap.pdf")
    Proteopedia::Reset_Dev()

    ProteinData[, SampleMed := median(Log2LFQ, na.rm = T), .(Sample)]
    ProteinData[, Log2FCvsMedian := Log2LFQ - SampleMed]

    PCAData <- ProteinData[, .(ProteinGroup, Sample, Log2FCvsMedian)] |> data.table::dcast(ProteinGroup ~ Sample, value.var = "Log2FCvsMedian") |>
      tidyr::drop_na() |> data.frame(row.names = "ProteinGroup") |> t() |> stats::prcomp(scale. = T)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |> data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Group, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) +
                                    ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Group, colour = Temp)) +
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

    BCASamples <- BCAData[!grepl("Standard", Sample)][, Conc := (Absorbance - Standard_c) / Standard_m] |> data.table::merge.data.table(Metadata[, .(Sample, Group, Temp)])

    pdf("RawBCAPlot.pdf", width = 16, height = 16)
    print(BCAStandards |> ggplot2::ggplot(ggplot2::aes(x = Conc, y = Absorbance)) +
            ggplot2::geom_point(stroke = NA) +ggplot2::stat_smooth(method = "lm", se = F, colour = "#AAA", linetype = "dashed", size = 1) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.5, label = paste0("R^2 == ", Standard_R2), parse = T) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.4, label = paste0("y = ", Standard_m, "x + ", Standard_c)) +
            ggplot2::geom_rug(data = BCASamples, ggplot2::aes(colour = Group, alpha = factor(Temp))) +
            ggplot2::scale_alpha_manual(values = c(1, 0.85, 0.7, 0.55, 0.4, 0.25), name = Temp.~degree~C) +
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

    BCASamples <- BCAData[!grepl("Standard", Sample)][, Conc := (Absorbance - Standard_c) / Standard_m] |> data.table::merge.data.table(Metadata[, .(Sample, Group, Temp)])

    pdf("AdjBCAPlot.pdf", width = 16, height = 16)
    print(BCAStandards |> ggplot2::ggplot(ggplot2::aes(x = Conc, y = Absorbance)) +
            ggplot2::geom_point(stroke = NA) +ggplot2::stat_smooth(method = "lm", se = F, colour = "#AAA", linetype = "dashed", size = 1) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.5, label = paste0("R^2 == ", Standard_R2), parse = T) +
            ggplot2::annotate("text", x = mean(BCAStandards$Conc, na.rm = T), y = mean(BCAStandards$Absorbance, na.rm = T)*0.4, label = paste0("y = ", Standard_m, "x + ", Standard_c)) +
            ggplot2::geom_rug(data = BCASamples, ggplot2::aes(colour = Group, alpha = factor(Temp))) + ggplot2::scale_alpha_manual(values = c(1, 0.85, 0.7, 0.55, 0.4, 0.25), name = Temp.~degree~C) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette) + ggplot2::labs(x = "Protein Concentration (mg/mL)", y = "562 nm Absorbance (AU)"))
    Proteopedia::Reset_Dev()

    BCASamples[, Correction := Conc/BCASamples[, .(MedConc = median(Conc)), .(Group, Temp)]$MedConc |> max(), Sample]
    BCASamples <- BCASamples[, .(Correction = median(Correction)), Temp]
    BCASamples[, Temp := factor(Temp)]

    pdf("CorrectionCoeffsPlot.pdf", width = 16, height = 12)
    print(BCASamples |> unique() |> ggplot2::ggplot(ggplot2::aes(x = Temp, y = Correction, fill = Temp)) +
            ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2(width = 0.9)) +
            ggplot2::scale_fill_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)], guide = "none") +
            ggplot2::scale_y_continuous(expand = 0) + ggplot2::labs(x = Temp.~degree~C, y = "BCA Correction Coeff."))
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

    ProcessedDataMerged <- ProcessedData |> data.table::merge.data.table(ProcessedFCData, by = c("Sample", "ProteinGroup", "Condition", "Group", "Temp", "Replicate"))
    ProcessedDataMerged[, DetectLevel := mean(head(sort(NormLog2LFQ, decreasing = F), 200), na.rm = T), Sample]
    ProcessedDataMerged[, N_Missing := sum(is.na(Log2FCvsMedian)), .(ProteinGroup, Condition)]

    # Calculating for same timepoint how many missing in each condition to see if missingness is biologically informative and need imputation
    DiffMissing <- unique(ProcessedDataMerged[, .(ProteinGroup, Group, Temp, N_Missing)]) |> data.table::dcast(ProteinGroup + Temp ~ Group, value.var = "N_Missing")
    DiffMissing[, DiffMissing := get(ExpGroup) - get(CtlGroup)]
    ProcessedDataMerged <- data.table::merge.data.table(ProcessedDataMerged, DiffMissing[, .(ProteinGroup, Temp, DiffMissing)], by = c("ProteinGroup", "Temp"))

    # Determine & Execute Imputation ####
    ProcessedDataMerged[, Temp := factor(Temp, levels = TempLevels)]
    ProcessedDataMerged[, NextTemp := Proteopedia::Calculate_RelativeTemp(Temp, TempLevels, Shift = 1) |> as.character(), Temp]

    NextTempData <- ProcessedDataMerged[, .(NextTemp_NMissing = mean(N_Missing)), .(Temp, ProteinGroup, Group)] |> data.table::copy() |> data.table::setnames("Temp", "NextTemp")

    ProcessedDataMerged <- data.table::merge.data.table(ProcessedDataMerged, NextTempData, by = c("NextTemp", "ProteinGroup", "Group"), all.x = T)

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
            ggplot2::labs(x = expression("Log"[2]~"LFQ"), y = "No. Proteins") + ggplot2::scale_y_continuous(expand = 0) +
            ggplot2::facet_grid(ggplot2::vars(Group), ggplot2::vars(Temp), labeller = ggplot2::labeller(Temp = TempLabels)))
    Proteopedia::Reset_Dev()

    ProcessedDataMerged <- ProcessedDataMerged |> data.table::merge.data.table(BCASamples, by = "Temp") |>
      data.table::merge.data.table(ProteinData[, .(Sample, ProteinGroup, Log2LFQ)], by = c("Sample", "ProteinGroup"), all.x = T)
    ProcessedDataMerged[, BCACorrLog2LFQ := ImpLog2LFQ*Correction]

    pdf("LFQProcessingBoxplots.pdf", width = 16, height = 12)
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Group), colour = as.factor(Replicate), y = Log2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Raw Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Group), colour = as.factor(Replicate), y = NormLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Normalised Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Group), colour = as.factor(Replicate), y = ImpLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("Post-Imputation Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    print(ProcessedDataMerged |> ggplot2::ggplot(ggplot2::aes(x = as.factor(Group), colour = as.factor(Replicate), y = BCACorrLog2LFQ)) +
            ggplot2::geom_boxplot(outliers = F) + ggplot2::facet_wrap(~Temp, nrow = 1, labeller = TempLabels) +
            ggplot2::scale_colour_manual(values = Proteopedia::NiceColourPalette, guide = "none") +
            ggplot2::labs(y = expression("BCA-Corrected Log"[2] ~ "LFQ")) + ggplot2::theme(axis.title.x = ggplot2::element_blank()))
    Proteopedia::Reset_Dev()

    ProcessedProteinLFQs <- ProcessedDataMerged[, .(ProteinGroup, Sample, Condition, Replicate, Group, Temp, NextTemp, DetectLevel,
                                                    Impute, Log2LFQ, NormLog2LFQ, ImpLog2LFQ, BCACorrLog2LFQ)] |> data.table::copy()
    data.table::fwrite(ProcessedProteinLFQs, "ProcessedProteinLFQs.csv")
  }
  message("Preparing Limma Input")
  {
    ProtLFQsInput <- as.matrix(tibble::column_to_rownames(data.table::dcast(ProcessedProteinLFQs, ProteinGroup ~ Sample, value.var = "BCACorrLog2LFQ"), "ProteinGroup"))
    LimmaInputData <- data.table::data.table(tibble::rownames_to_column(as.data.frame(ProtLFQsInput), "ProteinGroup"))
    LimmaInputData <- data.table::merge.data.table(data.table::melt.data.table(LimmaInputData, id.vars = "ProteinGroup", variable.name = "Sample",
                                                                               value.name = "Log2LFQ"), Metadata)
    LimmaInputData[, `:=`(MeanLog2LFQ, mean(Log2LFQ, na.rm = T)), .(ProteinGroup, Group, Temp)]

    LimmaInputMatrix <- as.matrix(tibble::column_to_rownames(data.table::dcast(LimmaInputData, ProteinGroup ~ Sample, value.var = "Log2LFQ"), "ProteinGroup"))
    LimmaInputMatrix <- LimmaInputMatrix[matrixStats::rowMeans2(is.na(LimmaInputMatrix)) <= PropMissingness,]
    LimmaInputData <- LimmaInputData[ProteinGroup %in% rownames(LimmaInputMatrix)]

    PCAData <- LimmaInputMatrix[matrixStats::rowMeans2(is.na(LimmaInputMatrix)) == 0,] |> t() |> stats::prcomp(scale. = T)
    SummaryPCA <- summary(PCAData)$importance
    All_PCA <- data.table::data.table(PCAData$x, keep.rownames = "Sample")[, .(Sample, PC1, PC2, PC3, PC4)] |>
      data.table::merge.data.table(Metadata)
    Output_PCA <- patchwork::free(All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, shape = Group, colour = Temp)) +
                                    ggplot2::geom_point(size = 6, stroke = NA) + ggplot2::scale_colour_manual(values = Proteopedia::ThermalPalette[seq(20, 46, length.out = 5)]) +
                                    ggrepel::geom_text_repel(ggplot2::aes(label = Temp)) +
                                    ggplot2::labs(x = paste("PC1 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC1"]*100, 0), "%]", sep = ""),
                                                  y = paste("PC2 [", round(SummaryPCA[rownames(SummaryPCA) == "Proportion of Variance", "PC2"]*100, 0), "%]", sep = "")) +
                                    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), legend.position = "none") +
                                    All_PCA |> ggplot2::ggplot(ggplot2::aes(x = PC3, y = PC4, shape = Group, colour = Temp)) +
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
    Targets <- Metadata[, .(Sample, Group, Temp)]
    Targets <- Targets[match(colnames(LimmaInputMatrix), Sample)]
    Temp <- Targets$Temp
    Group <- factor(Targets$Group)
    ModelDesign <- model.matrix(~ Temp + Temp:Group)
    colnames(ModelDesign) <- gsub(paste0("Group", ExpGroup), "Group_Exp", colnames(ModelDesign))
    LimmaOutput <- limma::eBayes(limma::lmFit(Biobase::ExpressionSet(assayData = LimmaInputMatrix), ModelDesign))

    if(!is.finite(LimmaOutput$df.prior)){message("Warning: Limma Prior is Infinite")}
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
    LimmaStats[grepl(":", Coef), Group := "Reversine"]
    LimmaStats[grepl("Temp\\d+$", Coef), Group := "DMSO"]
    LimmaStats[grepl("Temp", Coef), Temp := gsub("Temp(\\d+).*", "\\1", Coef)]
    LimmaStats[, Significance := data.table::fifelse(adj.P.Val < 0.05, "Adj. P-Value < 0.05",
                                                     data.table::fifelse(P.Value < 0.05, "P-Value < 0.05", "None"))]
    LimmaStats[grepl(":", Coef), SigMark := data.table::fifelse(P.Value < 0.001, "***", data.table::fifelse(P.Value < 0.01, "**", data.table::fifelse(P.Value < 0.05, "*", "")))]
    LimmaStats[, Temp := factor(Temp, levels = TempLevels)]
    LimmaStats[, Significance := factor(Significance, levels = c("Adj. P-Value < 0.05", "P-Value < 0.05", "None"))]

    data.table::fwrite(LimmaStats, "LimmaStats.csv")

    CoefLabels <- paste0(TempLevels, "~degree*C")
    names(CoefLabels) <- paste0("Temp",TempLevels, ":Group_Exp")
    CoefLabels <- ggplot2::as_labeller(CoefLabels, default = ggplot2::label_parsed)

    pdf("LimmaPValueHistograms.pdf", width = 16, height = 12)
    print(LimmaStats[grepl(":", Coef)] |> ggplot2::ggplot(ggplot2::aes(x = P.Value, fill = Significance)) + ggplot2::geom_histogram() +
            ggplot2::facet_wrap(~Coef, labeller = ggplot2::labeller(Coef = CoefLabels), scales = "free_y") +
            ggplot2::scale_fill_manual(values = c("Adj. P-Value < 0.05" = "#FB0", "P-Value < 0.05" = "#F10", "None" = "#300")) +
            ggplot2::scale_y_continuous(expand = 0) + ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::labs(x = "Adjusted P-Value", y = "No. Proteins") + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2)))
    print(LimmaStats[grepl(":", Coef)] |> ggplot2::ggplot(ggplot2::aes(x = adj.P.Val, fill = Significance)) + ggplot2::geom_histogram() +
            ggplot2::facet_wrap(~Coef, labeller = ggplot2::labeller(Coef = CoefLabels), scales = "free_y") +
            ggplot2::scale_fill_manual(values = c("Adj. P-Value < 0.05" = "#FB0", "P-Value < 0.05" = "#F10", "None" = "#300")) +
            ggplot2::scale_y_continuous(expand = 0) + ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::labs(x = "Adjusted P-Value", y = "No. Proteins") + ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2)))
    Proteopedia::Reset_Dev()

    AffectedProteins <- LimmaStats[stringr::str_detect(Coef, ":") & P.Value < 0.05][stringr::str_detect(Coef, TempLevels[1], negate = T)]
    AffectedProteinData <- LimmaInputData[ProteinGroup %in% AffectedProteins$ProteinGroup & !is.na(Log2LFQ)]
    data.table::fwrite(AffectedProteinData, "AffectedProteinsLFQs.csv")
  }
  message(paste0("Calculating Fold Changes vs ", Log2FCvs))
  {
    ProcessedDataMerged[Temp == min(TempLevels), ProteinLFQMax := mean(BCACorrLog2LFQ, na.rm = T), .(ProteinGroup, Group)]
    ProcessedDataMerged[, ProteinLFQMax := mean(ProteinLFQMax, na.rm = T), .(ProteinGroup, Group)]
    ProcessedDataMerged[, Log2FCvsMax := BCACorrLog2LFQ - ProteinLFQMax]

    ProcessedDataMerged[Temp == max(TempLevels), ProteinLFQMin := mean(BCACorrLog2LFQ, na.rm = T), .(ProteinGroup, Group)]
    ProcessedDataMerged[, ProteinLFQMin := mean(ProteinLFQMin, na.rm = T), .(ProteinGroup, Group)]
    ProcessedDataMerged[, Log2FCvsMin := BCACorrLog2LFQ - ProteinLFQMin]
    ProcessedProteinFCs <- ProcessedDataMerged[, .(ProteinGroup, Sample, Condition, Replicate, Group, Temp, NextTemp, DetectLevel,
                                                   Impute, Log2FCvsMin, Log2FCvsMedian, Log2FCvsMax)] |> data.table::copy()
    data.table::fwrite(ProcessedProteinFCs, "ProcessedProteinLog2FCs.csv")

    FoldChangeData <- ProcessedProteinFCs[, .(ProteinGroup, Group, Temp, Condition, Replicate, Sample,
                                              get(colnames(ProcessedProteinFCs)[grepl(Log2FCvs, colnames(ProcessedProteinFCs))]))] |>
      data.table::copy() |> data.table::setnames("V7", "Log2FC")
    FoldChangeData <- FoldChangeData[, MeanLog2FC := mean(Log2FC, na.rm = T), .(ProteinGroup, Temp, Group)]
    FoldChangeData[, `:=`(GeneGroup, sapply(strsplit(ProteinGroup, ";"), function(x){paste(unique(gsub("-.+", "", x)), collapse = ";")}))]
    FoldChangeDataWide <- FoldChangeData[, .(MeanLog2FC = mean(Log2FC, na.rm = T)), .(ProteinGroup, Temp, Group)] |>
      data.table::dcast(ProteinGroup + Temp ~ Group, value.var = "MeanLog2FC") |> data.table::merge.data.table(ProteinGroups, by = "ProteinGroup")

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
      PlotsOI[[POI]] = LimmaInputData[ProteinGroup == POI] |> ggplot2::ggplot(ggplot2::aes(x = factor(Temp), y = Log2LFQ, colour = Group, group = Group)) +
        ggplot2::geom_point() + ggplot2::geom_line(data = LimmaInputData[ProteinGroup == POI, .(MeanLog2LFQ, Group, Temp)] |> unique(),
                                                   ggplot2::aes(x = factor(Temp), y = MeanLog2LFQ, colour = Group, group = Group)) +
        ggplot2::scale_color_manual(values = Proteopedia::NiceColourPalette) +
        ggplot2::labs(title = paste0(ProteinInfo[ProteinGroup == POI, Gene]  |> unique()), x = Temp.~degree~C, y = Log[2]~LFQ)
    }

    pdf("Top20Proteins.pdf", width = 20, height = 16)
    print(ggpubr::ggarrange(plotlist = PlotsOI, common.legend = T))
    Proteopedia::Reset_Dev()

    Log2FCPlotsOI = list()
    for(POI in (AffectedProteins$ProteinGroup |> unique())[1:20]){
      Log2FCPlotsOI[[POI]] = FoldChangeData[ProteinGroup == POI] |> ggplot2::ggplot(ggplot2::aes(x = factor(Temp), y = Log2FC, colour = Group, group = Group)) +
        ggplot2::geom_point() + ggplot2::geom_line(data = FoldChangeData[ProteinGroup == POI, .(MeanLog2FC, Group, Temp)] |> unique(),
                                                   ggplot2::aes(x = factor(Temp), y = MeanLog2FC, colour = Group, group = Group)) +
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
########### MS Analysis: Protein Data to Ontology Functions ###################################################################################################################################
#' @export
Order_Terms <- function(x){
  if(!is.null(x)){
    if(nrow(x) > 0){term_distance = as.dist(1-enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)
      if(length(term_distance) == 0){ordered_terms <- rownames(enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)
      }else{rownames(enrichplot::pairwise_termsim(x, showCategory = Inf)@termsim)[stats::hclust(term_distance, method = "ward.D")$order]}
    }
  } else {NULL}
}
#' @export
Process_GSEAOutput <- function(GSEA_Output, GSEADatabase, PlotColour = "#000"){

  GSEA_OutputData <- GSEA_Output |> data.frame() |> data.table::data.table()
  if(nrow(GSEA_OutputData)){
    GSEA_OutputData[, Ontology := GSEADatabase]
    GSEA_OutputTop30 <- GSEA_OutputData[pvalue < 0.05] |> dplyr::slice_max(order_by = rank, n = 30, with_ties = F)

    GSEA_Dotplot <- GSEA_OutputData |> ggplot2::ggplot(ggplot2::aes(x = NES, y = -log10(pvalue), colour = PlotColour)) +
      ggplot2::geom_point(size = 3) + ggplot2::coord_cartesian(x = c(-ceiling(max(abs(GSEA_OutputData$NES))), ceiling(max(abs(GSEA_OutputData$NES))))) +
      ggplot2::labs(x = "Normalised Enrichment Score (NES)", y = expression("-Log"[10]~"P-Value")) + Proteopedia::Add_NotSigBox() +
      ggplot2::scale_colour_manual(values = PlotColour, guide = "none")
    pdf(paste0("GSEA_", GSEADatabase, "_Dotplot.pdf"), width = 10, height = 8)
    print(GSEA_Dotplot)
    Proteopedia::Reset_Dev()

    GenesetsTop30 <- GSEA_OutputTop30 |> ggplot2::ggplot(ggplot2::aes(x = enrichmentScore, y = stringr::str_wrap(Description, 60), colour = pvalue, size = pvalue)) +
      ggplot2::geom_point() + ggplot2::geom_segment(aes(x = 0, xend = enrichmentScore, linewidth = 0.1), show.legend = F) + ylab("Description") +
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
Process_ORAOutput <- function(ORA_Output, ORADatabase, PlotColour = "#000"){
  ORA_OutputData <- ORA_Output |> data.frame() |> data.table::data.table()
  if(nrow(ORA_OutputData) > 0){
    ORA_OutputData[, Ontology := ORADatabase]
    ORAOutputTop30 <- ORA_OutputData[pvalue < 0.05] |> dplyr::slice_max(order_by = p.adjust*RichFactor, n = 30, with_ties = F)

    GenesetsTop30 <- ORAOutputTop30 |> ggplot2::ggplot(ggplot2::aes(x = RichFactor, y = stringr::str_wrap(Description, 60), colour = Count, size = pvalue)) +
      ggplot2::geom_point() + ggplot2::geom_segment(aes(x = 0, xend = RichFactor, linewidth = 0.1), show.legend = F) + ylab("Description") +
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
Run_GSEA <- function(InputFile, OutputName = "GSEA_Output", FastGSEA = T){
  start.time <- Sys.time()
  setwd(gsub("(.*)/.*.csv", "\\1", InputFile))
  InputData <- data.table::fread(InputFile)
  if(is.null(InputData$Parameter)){
    Parameters <- "Abundance"
    InputData[, Parameter := "Abundance"]
  } else {Parameters <- unique(InputData$Parameter)}
  for(Parameter in Parameters){
    setwd(gsub("(.*)/.*.csv", "\\1", InputFile))
    if(dir.exists(paste0(getwd(),"/", Parameter, OutputName))){unlink(paste0(getwd(),"/", Parameter, OutputName), recursive = T)}
    dir.create(paste0(getwd(),"/", Parameter, OutputName), showWarnings = T)
    setwd(paste0(getwd(),"/", Parameter, OutputName))
    message(paste0(Parameter, ": Loading & Formatting Data"))
    {
      InputData <- InputData[Parameter == Parameter, .(Gene, Log2FC)] |> data.table::merge.data.table(clusterProfiler::bitr(InputData$Gene, from = "ALIAS", to = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db) |>
                                                                                                        data.table::setnames("ALIAS", "Gene"), by = "Gene")

      GSEAData <- InputData$Log2FC[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
      names(GSEAData) <- InputData$ENTREZID[which(InputData$ENTREZID %in% unique(InputData$ENTREZID))]
      GSEAData <- BiocGenerics::sort(GSEAData[!duplicated(names(GSEAData))], decreasing = T)
    }
    message(paste0(Parameter, ": Performing Enrichment Analyses"))
    {
      if(FastGSEA){
        gseaBP <- clusterProfiler::gseGO(geneList = GSEAData, ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaMF <- clusterProfiler::gseGO(geneList = GSEAData, ont = "MF", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaCC <- clusterProfiler::gseGO(geneList = GSEAData, ont = "CC", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaKEGG <- clusterProfiler::gseKEGG(geneList = GSEAData, organism = "hsa", keyType = "ncbi-geneid",
                                             minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaMKEGG <- clusterProfiler::gseMKEGG(geneList = GSEAData, organism = "hsa", keyType  = "ncbi-geneid",
                                               minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaReactome <- ReactomePA::gsePathway(geneList = GSEAData, organism = "human",
                                               minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaDisease <- DOSE::gseDO(geneList = GSEAData, ont = "HDO", organism = "human",
                                   minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaPhenotype <- DOSE::gseDO(geneList = GSEAData, ont = "HPO", organism = "human",
                                     minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaDisGeNET <- DOSE::gseDGN(geneList = GSEAData, minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaNCG <- DOSE::gseNCG(geneList = GSEAData, organism = "human", minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
        gseaHallmark <- clusterProfiler::GSEA(geneList = GSEAData, TERM2GENE = msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                                              minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH")
      } else {
        gseaBP <- clusterProfiler::gseGO(geneList = GSEAData, ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaMF <- clusterProfiler::gseGO(geneList = GSEAData, ont = "MF", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaCC <- clusterProfiler::gseGO(geneList = GSEAData, ont = "CC", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                         minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaKEGG <- clusterProfiler::gseKEGG(geneList = GSEAData, organism = "hsa", keyType = "ncbi-geneid",
                                             minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaMKEGG <- clusterProfiler::gseMKEGG(geneList = GSEAData, organism = "hsa", keyType  = "ncbi-geneid",
                                               minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaReactome <- ReactomePA::gsePathway(geneList = GSEAData, organism = "human",
                                               minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaDisease <- DOSE::gseDO(geneList = GSEAData, ont = "HDO", organism = "human",
                                   minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaPhenotype <- DOSE::gseDO(geneList = GSEAData, ont = "HPO", organism = "human",
                                     minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaDisGeNET <- DOSE::gseDGN(geneList = GSEAData, minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaNCG <- DOSE::gseNCG(geneList = GSEAData, organism = "human", minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
        gseaHallmark <- clusterProfiler::GSEA(geneList = GSEAData, TERM2GENE = msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                                              minGSSize = 3, maxGSSize = 800, verbose = F, pvalueCutoff = 1, pAdjustMethod = "BH", by = "DOSE")
      }
    }
    message(paste0(Parameter, ": Exporting Plots & Results"))
    {
      MergedGSEA <- data.table::rbindlist(list(Proteopedia::Process_GSEAOutput(gseaBP, "Biological_Process"), Proteopedia::Process_GSEAOutput(gseaMF, "Molecular_Function"),
                                               Proteopedia::Process_GSEAOutput(gseaCC, "Cellular_Component"), Proteopedia::Process_GSEAOutput(gseaKEGG, "KEGG"),
                                               Proteopedia::Process_GSEAOutput(gseaMKEGG, "MKEGG"), Proteopedia::Process_GSEAOutput(gseaReactome, "Reactome"),
                                               Proteopedia::Process_GSEAOutput(gseaDisease, "Disease"), Proteopedia::Process_GSEAOutput(gseaPhenotype, "Phenotype"),
                                               Proteopedia::Process_GSEAOutput(gseaDisGeNET, "DisGeNET"), Proteopedia::Process_GSEAOutput(gseaNCG, "Network_of_Cancer_Genes"),
                                               Proteopedia::Process_GSEAOutput(gseaHallmark, "Hallmark")))
      MergedGSEA[, N_Genes := 0]
      for(GSEARow in 1:nrow(MergedGSEA)){
        MergedGSEA$N_Genes[GSEARow] <- length(stringr::str_split(MergedGSEA$core_enrichment[GSEARow], "/")[[1]])
      }
      MergedGSEA[, PropGenes := N_Genes/setSize]
      data.table::fwrite(MergedGSEA, file = "GSEA_Results.csv")
    }
    message(paste0(Parameter, ": Exporting HTML Plots"))
    {
      MergedGSEA <- MergedGSEA |> dplyr::arrange(factor(Ontology, levels = c("Cellular_Component", "Biological_Process", "Molecular_Function", "KEGG", "MKEGG",
                                                                             "Reactome", "Disease","Phenotype", "DisGeNET", "Network_of_Cancer_Genes", "Hallmark")))
      MergedGSEA[, `:=`(TermIndex = 1:nrow(MergedGSEA), NegLog10PAdj = -log10(p.adjust), SciPVal = formatC(p.adjust, format = "e", digits = 2))]

      GSEA_Significance <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = TermIndex, y = -log10(p.adjust), group = Ontology)) |>
        highcharter::hc_chart(zoomType = "xy") |>
        highcharter::hc_xAxis(title = list(text = NULL), labels = list(enabled = F), lineWidth = 0.5, lineColor = "#000", tickWidth = 0 ) |>
        highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000",
                              gridLineWidth = 0 ) |>
        highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Normalised Enrichment Score: {point.NES:.2f}
                                <br>P-Value: {point.SciPVal} <br>Gene count: {point.setSize}")
      htmlwidgets::saveWidget(GSEA_Significance, file = "Manhattan_Plot.html")

      GSEA_EnrichSig <- highcharter::hchart(MergedGSEA, "scatter", highcharter::hcaes(x = NES, y = -log10(p.adjust), group = Ontology)) |>
        highcharter::hc_chart(zoomType = "xy") |>
        highcharter::hc_xAxis(title = list(text = "Normalised Enrichment Score"), labels = list(enabled = F), lineWidth = 0.5, lineColor = "#000",
                              tickWidth = 0 ) |>
        highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000",
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
    if(dir.exists(paste0(getwd(),"/", OutputName))){unlink(paste0(getwd(),"/",OutputName), recursive = T)}
    dir.create(paste0(getwd(),"/", OutputName), showWarnings = T)
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
    oraHallmark <- clusterProfiler::enricher(SearchProteins, TERM2GENE = msigdbr::msigdbr(db_species = "HS", collection = "H"))
  }
  message("Plotting & Compiling Analyses")
  {
    MergedORA <- data.table::rbindlist(list(Proteopedia::Process_ORAOutput(oraBP, "Biological_Process"), Proteopedia::Process_ORAOutput(oraMF, "Molecular_Function"),
                                            Proteopedia::Process_ORAOutput(oraCC, "Cellular_Component"), Proteopedia::Process_ORAOutput(oraKEGG, "KEGG"),
                                            Proteopedia::Process_ORAOutput(oraMKEGG, "MKEGG"), Proteopedia::Process_ORAOutput(oraReactome, "Reactome"),
                                            Proteopedia::Process_ORAOutput(oraDisease, "Disease"), Proteopedia::Process_ORAOutput(oraPhenotype, "Phenotype"),
                                            Proteopedia::Process_ORAOutput(oraDisGeNET, "DisGeNET"), Proteopedia::Process_ORAOutput(oraNCG, "Network_of_Cancer_Genes"),
                                            Proteopedia::Process_ORAOutput(oraHallmark, "Hallmark")), fill = T)
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
                                                                     Proteopedia::Order_Terms(oraDisGeNET), Proteopedia::Order_Terms(oraNCG),
                                                                     Proteopedia::Order_Terms(oraHallmark))))]
    suppressWarnings(MergedORA[, TermIndex := 1:nrow(MergedORA)])
    MergedORA[, NegLog10PAdj := -log10(p.adjust)]
    MergedORA[, SciPVal := formatC(p.adjust, format = "e", digits = 2)]
  }
  message("Exporting HTML Plots")
  {
    ORA_Significance <- highcharter::hchart(MergedORA, "scatter", highcharter::hcaes(x = TermIndex, y = -log10(p.adjust), group = Ontology)) |>
      highcharter::hc_chart(zoomType = "xy") |>
      highcharter::hc_xAxis(title = list(text = NULL), labels = list(enabled = F), lineWidth = 0.5, lineColor = "#000", tickWidth = 0 ) |>
      highcharter::hc_yAxis(title = list(text = "-Log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000",
                            gridLineWidth = 0 ) |>
      highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Enrichment Factor: {point.RichFactor:.2f}
                            <br>P-Value: {point.SciPVal} <br> Gene count: {point.Count}")
    htmlwidgets::saveWidget(ORA_Significance, file = "Manhattan_Plot.html")

    ORA_EnrichSig <- highcharter::hchart(MergedORA, "scatter", highcharter::hcaes(x = RichFactor, y = -log10(p.adjust), group = Ontology)) |>
      highcharter::hc_chart(zoomType = "xy") |>
      highcharter::hc_xAxis(title = list(text = "Normalised Enrichment Score"), labels = list(enabled = F), lineWidth = 0.5, lineColor = "#000",
                            tickWidth = 0 ) |>
      highcharter::hc_yAxis(title = list(text = "-log10(Adj. P-Value)"), lineWidth = 0.5, tickWidth = 0.5, lineColor = "#000", tickColor = "#000",
                            gridLineWidth = 0 ) |>
      highcharter::hc_tooltip(headerFormat = "", pointFormat = "<b>{point.Description}</b> <br>Enrichment Factor: {point.RichFactor:.2f}
                            <br>P-Value: {point.SciPVal} <br>Gene count: {point.Count}")
    htmlwidgets::saveWidget(ORA_EnrichSig, file = "Volcano_Plot.html")
  }
  Proteopedia::End_Timer(Start = start.time)
}
#' @export
Plot_LabelFreeGenesetVolcano <- function(InputFile, SubsetColour = "#F00", Adj_PValueCutoff = 0.05){

  initiation.time <- Sys.time()
  setwd(gsub("(.*)/GSEA_Output/.*.csv.*", "\\1", InputFile))
  GenesetData <- data.table::fread(InputFile)[p.adjust < Adj_PValueCutoff]
  LimmaData <- data.table::fread(paste0(gsub("(.*)/GSEA_Output.*.csv.*", "\\1", InputFile), "/Limma_Output.csv"))
  if(dir.exists(paste0(getwd(),"/GSEA_Volcanoes"))){unlink(paste0(getwd(),"/GSEA_Volanoes"), recursive = T)}
  dir.create(paste0(getwd(),"/GSEA_Volcanoes"), showWarnings = T)
  setwd(paste0(getwd(),"/GSEA_Volcanoes"))
  for(GSEARow in 1:nrow(GenesetData)){
    SubsetProteins <- bitr(unlist(str_split(GenesetData$core_enrichment[GSEARow], "/")), fromType = "ENTREZID", toType = "UNIPROT",
                           OrgDb = org.Hs.eg.db::org.Hs.eg.db)$UNIPROT

    SubsetData <- LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, GeneGroup)]
    SubsetData[, Subset := data.table::fifelse(GeneGroup %in% SubsetProteins, T, F)]

    MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)

    Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) +
      Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() +
      labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))

    adjpval <- GenesetData$p.adjust[GSEARow]

    Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") +
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
Plot_SILACRatioGenesetVolcano <- function(InputFile, SubsetColour = "#F00", Adj_PValueCutoff = 0.05){
  initiation.time <- Sys.time()
  setwd(gsub("(.*)/GSEA_Output/.*.csv.*", "\\1", InputFile))
  GenesetData <- data.table::fread(InputFile)[p.adjust < Adj_PValueCutoff]
  LimmaData <- data.table::fread(paste0(gsub("(.*)/GSEA_Output.*.csv.*", "\\1", InputFile), "/Limma_Output.csv"))
  if(dir.exists(paste0(getwd(),"/GSEA_Volcanoes"))){unlink(paste0(getwd(),"/GSEA_Volanoes"), recursive = T)}
  dir.create(paste0(getwd(),"/GSEA_Volcanoes"), showWarnings = T)
  setwd(paste0(getwd(),"/GSEA_Volcanoes"))
  for(GSEARow in 1:nrow(GenesetData)){
    SubsetProteins <- bitr(unlist(str_split(GenesetData$core_enrichment[GSEARow], "/")), fromType = "ENTREZID", toType = "UNIPROT",
                           OrgDb = org.Hs.eg.db::org.Hs.eg.db)$UNIPROT

    SubsetData = LimmaData[, .(ProteinGroup, Gene, Log2FC, P.Value, GeneGroup)]
    SubsetData[, Subset := data.table::fifelse(GeneGroup %in% SubsetProteins, T, F)]

    MeanLog2FC <- SubsetData[Subset == T, Log2FC] |> mean(na.rm = T)

    Volcano <- SubsetData[Subset == T] |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(alpha = 0.7, stroke = NA) + ggrepel::geom_text_repel(ggplot2::aes(label = as.character(Gene))) +
      Proteopedia::Add_NotSigBox() + Proteopedia::Add_AbundanceAxes() +
      labs(caption = paste0(GenesetData$ID[GSEARow], "\n", GenesetData$Description[GSEARow]))

    adjpval <- GenesetData$p.adjust[GSEARow]

    Volcano_Rug <- SubsetData |> ggplot2::ggplot(ggplot2::aes(x = Log2FC, y = -log10(P.Value))) +
      ggplot2::geom_point(stroke = NA, alpha = 0.3, size = 3, colour = "#000") +
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
########### Coulter Analysis Functions ########################################################################################################################################################
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
    BinSummary <- CoulterData[, .(N = .N, MeanProp = mean(CellProp), SDProp = sd(CellProp)), .(Condition, BinVol)]
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
            ggplot2::coord_cartesian(ylim = c(0, max(ConditionSummary$MaxProp)*1.2), clip = "off") + ggplot2::scale_x_continuous(limits = c(2.9, 4.35)) +
            ggrepel::geom_label_repel(data = ConditionSummary, ggplot2::aes(x = log10(ConditionSummary$MaxBinVol), y = max(ConditionSummary$MaxProp)*1.1, label = VolLabel),
                                      size = 8, parse = T, nudge_x = 0, min.segment.length = 10, force_pull = 0, max.overlaps = Inf) + ggplot2::scale_y_continuous(expand = 0) +
            patchwork::inset_element(CoulterData[, .(Condition, Replicate, TotalCount)] |> unique() |>
                                       ggplot2::ggplot(ggplot2::aes(x = Condition, y = TotalCount/1000, fill = Condition)) +
                                       ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge2(width = 0.9)) + ggplot2::scale_fill_manual(values =Proteopedia::NiceColourPalette, guide = "none") +
                                       ggplot2::labs(y = "Cells Analysed (x1000)") + ggplot2::scale_y_continuous(expand = 0) + ggplot2::theme(axis.title.x = ggplot2::element_blank()), 0.60, 0.5, 0.9, 0.99))
    Proteopedia::Reset_Dev()
  } else {
    message("ERROR: No/Multiple Meta Files Detected")
  }
}
