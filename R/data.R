#' Data for the Standardised Analysis of Proteomic Data (Human)
#'
#' Chromosome Lengths
#'
#' ChromosomeLengths is a data.table with the base-pair lengths of each human Chromosome
#'
#' @format A data table with 25 rows and 2 columns.
#' describe{
#' item{Chromosome}{Chromosome ID}
#' item{Length}{Length of Chromosome in base-pairs}
#' }
#' @source Created by Jon Meyrick, by extracting lengths from the Human Reference Genome GENCODE V49 (GRCh38.p14)
#' @examples
#' data(ChromosomeLengths)
"ChromosomeLengths"
#'
#' Complexed Proteins
#'
#' ComplexedProteins contains the names of protein complexes compiled from both Complex Portal and CORUM
#' @source Created by Jon Meyrick
#' @examples
#' data(Complexed_Proteins)
"Complexed_Proteins"
#'
#' Complex Portal Complexes
#'
#' Proteins assigned to complexes within the Complex Portal database
#'
#' @source Created by Jon Meyrick
#' @examples
#' data(ComplexPortal_Complexes)
"ComplexPortal_Complexes"
#'
#' CORUM Complexes
#'
#' Proteins assigned to complexes within the CORUM database
#'
#' @source Created by Jon Meyrick
#' @examples
#' data(CORUM_Complexes)
"CORUM_Complexes"
#'
#' Degradation Profiles
#' @source Created by Jon Meyrick
#' @examples
#' data(DegradationProfiles)
"DegradationProfiles"
#'
#' Genomic Mapping
#' @source Created by Jon Meyrick
#' @examples
#' data(GenomicMapping)
"GenomicMapping"
#'
#' NiceColourPalette
#'
#' NiceColourPalette is a string of hexademical colour codes which allow for plotting of 24 distinct
#' variable groups using clear, distinct colours
#'
#' @format A character string of 24 length,
#' @source Created by Jon Meyrick
#' @examples
#' data(NiceColourPalette)
"NiceColourPalette"
#'
#' Proteopedia
#'
#' Proteopedia is a comprehensive dataset of the Human proteome, with various annotations of interest in proteomic analysis
#'
#' @format A data table with 42434 rows and 50 columns.
#' describe{
#' item{ProteinGroup}{UniProt identifer for protein groups}
#' item{Gene}{Encoding NCBI gene name for protein}
#' item{GeneGroup}{UniProt identifer for gene groups (identical to ProteinGroup, without isoforms)}
#' item{Entry}{UniProt identifer)}
#' item{Alt_Names}{Alternative names for protein-encoding gene)}
#' item{ProteinDescription}{Description of protein function}
#' item{Sequence}{Amino acid sequence of protein}
#' item{Length}{Length of protein in amino acids}
#' item{ENTREZID}{ENTREZ identifier for encoding gene}
#' item{Chromosome}{Encoding chromosome of protein}
#' item{MedianLociStart}{Median base-pair position of encoding gene on chromosome}
#' item{Deg_Profile}{Degradation profile of protein as determined by McShane et al. (2014)}
#' item{Experimental_Evidence_ComplexPortal}{Experimental evidence for incorporation of protein in complex (derived from Complex Portal)}
#' item{Experimental_Evidence_CORUM}{Experimental evidence for incorporation of protein in complex (derived from CORUM)}
#' item{N_ComplexPortal}{Number of entries for protein incorporation within complexes (derived from Complex Portal)}
#' item{N_CORUM}{Number of entries for protein incorporation within complexes (derived from CORUM)}
#' item{Aliphatic_Score}{Protein aliphatic score (calculated from protein sequence using the "peptides" package)}
#' item{Boman_Interaction_Score}{Protein interaction score (calculated from protein sequence using the "peptides" package)}
#' item{Hydrophobicity_Score}{Protein hydrophobicity score (calculated from protein sequence using the "peptides" package)}
#' item{MW}{Protein molecular weight (calculated from protein sequence using the "peptides" package)}
#' item{pI}{Protein isoelectric point at pH7 (calculated from protein sequence using the "peptides" package)}
#' item{Instability_Score}{Protein instability score (calculated from protein sequence using the "peptides" package)}
#' item{ER}{Localisation of the protein to the endoplasmic reticulum (ER), derived from UniProt}
#' item{Nucleus}{Localisation of the protein to the nucleus, derived from UniProt}
#' item{Kinetochore}{Formation of the protein in the kinetochore, derived from UniProt}
#' item{Ribosome}{Formation of the protein in the cytoplasmic ribosome, manually annotated}
#' item{RQC}{Involvement of the protein in cytoplasmic ribosome quality control, manually annotated}
#' item{Proteasome}{Formation of the protein in the proteasome, manually annotated}
#' item{E3Ligase}{Role of the protein as an E3 ubiquitin ligase, manually annotated}
#' item{E3LigaseAccesories}{Role of the protein as a contributor to the function of E3 ubiquitin ligases, manually annotated}
#' item{Mitochondria}{Localisation of the protein to the mitochondria, based on inclusion in MitoCarta 3.0}
#' item{OMM}{Localisation of the protein to the outer mitochondrial membrane (OMM), based on annotations in MitoCarta 3.0}
#' item{IMS}{Localisation of the protein to the intermembrane space (IMS) in mitochondria, based on annotations in MitoCarta 3.0}
#' item{IMM}{Localisation of the protein to the inner mitochondrial membrane (IMM), based on annotations in MitoCarta 3.0}
#' item{Matrix}{Localisation of the protein to the mitochondrial matrix, based on annotations in MitoCarta 3.0}
#' item{Mitoribosome}{Formation of the protein in the mitochondrial ribosome, based on annotations in MitoCarta 3.0}
#' item{MitoribosomeAssembly}{Involvement of the protein in assmebly of the mitochondrial ribosomes, based on annotations in MitoCarta 3.0}
#' item{OXPHOS}{Formation of the protein in any of the respiratory chain proteins or ATP synthase, based on annotations in MitoCarta 3.0}
#' item{MitochondriaDetox}{Involvement of the protein in detoxification processes in mitochondria, based on annotations in MitoCarta 3.0}
#' item{Lysosome}{Localisation of the protein to the lysosome, manually annotated}
#' item{Autophagy}{Involvement of the protein in autophagy, manually annotated}
#' item{Kinase}{Function of the protein as a kinase, manually annotated}
#' item{Apoptosis}{Involvement of the protein in apoptosis, manually annotated}
#' item{MCS_ER_Mito}{Involvement of the protein in endoplasmic reticulum (ER)-mitochondria membrane contact sites (MCS), manually annotated}
#' item{MCS_ER_PlasmaMembrane}{Involvement of the protein in endoplasmic reticulum (ER)-plasma membrane membrane contact sites (MCS), manually annotated}
#' item{MCS_ER_Golgi}{Involvement of the protein in endoplasmic reticulum (ER)-Gogli apparatus membrane contact sites (MCS), manually annotated}
#' item{MCS_ER_Phagosome}{Involvement of the protein in endoplasmic reticulum (ER)-phagosome membrane contact sites (MCS), manually annotated}
#' item{MCS_ER_Autophagosome}{Involvement of the protein in endoplasmic reticulum (ER)-autophagosome membrane contact sites (MCS), manually annotated}
#' }
#' @source Created by Jon Meyrick, compiling data from various sources
#' @examples
#' data(Proteopedia)
"Proteopedia"
#'
#' ThermalPalette
#'
#' ThermalPalette is a string of hexademical colour codes which allow for plotting of 46 continuous
#' variable groups with the colours white, yellow, orange, red, and black
#' @format A character string of 46 length.
#' @source Created by Jon Meyrick
#' @examples
#' data(ThermalPalette)
"ThermalPalette"
#'
#' UnmappedProteins
#'
#' Contains proteins which cannot be mapped to Chromosomes using GenomicMapping. In most cases, this is because the Gene identifer has been changed.
#'
#' @format A data table with 114 rows and 9 columns.
#' @source Created by Jon Meyrick
#' @examples
#' data(UnmappedProteins)
"UnmappedProteins"
