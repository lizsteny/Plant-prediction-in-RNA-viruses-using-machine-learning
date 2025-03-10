# ==============================
# step 1: load metadata and fasta sequences
# ==============================

# load required libraries
library(Biostrings)
library(stringr)
library(tidyr)
library(dplyr)


# define file paths
fasta_protein_path <- "data/raw/raw_protein_seq.fasta"
fasta_nuc_path <- "data/raw/raw_nuc_seq.fasta"
metadata_path <- "data/metadata/metadata.csv"
output_metadata_path <- "data/metadata/metadata_with_sequences.csv"
output_nuc_fasta_path <- "data/sequences/nuc_sequences.fasta"
output_protein_fasta_path <- "data/sequences/protein_sequences.fasta"


# ensure input files exist before proceeding
if (!file.exists(fasta_protein_path)) stop("error: protein fasta file not found!")
if (!file.exists(fasta_nuc_path)) stop("error: nucleotide fasta file not found!")
if (!file.exists(metadata_path)) stop("error: metadata file not found!")

# load metadata
message("loading metadata")
metadata <- read.csv(metadata_path)

# load protein sequences from fasta file
message("loading protein sequences")
protein_sequences <- readAAStringSet(fasta_protein_path)

# load nucleotide sequences from fasta file
message("loading nucleotide sequences")
nuc_sequences <- readDNAStringSet(fasta_nuc_path)

# ==============================
# step 2: process and format sequence data
# ==============================

# convert protein sequences to a dataframe with sequence name and sequence
protein_df <- data.frame(
  Sequence_Name = names(protein_sequences),
  Sequence = as.character(protein_sequences),
  stringsAsFactors = FALSE
)

# convert nucleotide sequences to a dataframe with sequence name and sequence
nuc_df <- data.frame(
  Sequence_Name = names(nuc_sequences),
  Sequence = as.character(nuc_sequences),
  stringsAsFactors = FALSE
)

# split sequence_name column into separate columns for id and assembly
protein_df <- transform(
  protein_df,
  ID = sapply(strsplit(Sequence_Name, "\\|"), `[`, 1),
  Assembly = sapply(strsplit(Sequence_Name, "\\|"), `[`, 2)
)

nuc_df <- transform(
  nuc_df,
  ID = sapply(strsplit(Sequence_Name, "\\|"), `[`, 1),
  Assembly = sapply(strsplit(Sequence_Name, "\\|"), `[`, 2)
)

# remove row names
rownames(protein_df) <- NULL
rownames(nuc_df) <- NULL

# drop the sequence_name column
protein_df <- protein_df[-1]
nuc_df <- nuc_df[-1]

# ==============================
# step 3: filter sequences based on metadata assemblies
# ==============================

# retain only assemblies present in the metadata
nuc_df <- nuc_df %>%
  filter(Assembly %in% metadata$assembly)

protein_df <- protein_df %>%
  filter(Assembly %in% metadata$assembly)

# ==============================
# step 4: concatenate sequences grouped by assembly
# ==============================

# define a function to concatenate sequences grouped by assembly
# sequences separated by "nnnnnnnnnn" to indicate distinct regions
concatenate_sequences <- function(df) {
  df %>%
    group_by(Assembly) %>%
    summarise(
      Sequence = paste(Sequence, collapse = "nnnnnnnnnn"),
      .groups = "drop"
    )
}

# apply the concatenation function to nucleotide and protein dataframes
nuc_df_final <- concatenate_sequences(nuc_df)
protein_df_final <- concatenate_sequences(protein_df)

# rename columns in the concatenated nucleotide dataframe
colnames(nuc_df_final)[colnames(nuc_df_final) == "Sequence"] <- "nuc_sequence"
colnames(nuc_df_final)[colnames(nuc_df_final) == "Assembly"] <- "assembly"

# rename columns in the concatenated protein dataframe
colnames(protein_df_final)[colnames(protein_df_final) == "Sequence"] <- "pro_sequence"
colnames(protein_df_final)[colnames(protein_df_final) == "Assembly"] <- "assembly"

# ==============================
# step 5: merge sequence data with metadata
# ==============================

# merge nucleotide sequences with metadata
merged_data <- merge(metadata, nuc_df_final, by = "assembly", all.x = TRUE)

# merge protein sequences with the previous result
metadata_with_sequences <- merge(merged_data, protein_df_final, by = "assembly", all.x = TRUE)

# ==============================
# step 6: export nucleotide and protein sequences to fasta
# ==============================

# nucleotide sequences processing
nuc_fasta <- DNAStringSet(metadata_with_sequences$nuc_sequence)
names(nuc_fasta) <- metadata_with_sequences$sequence

# save nucleotide sequences to a fasta file
writeXStringSet(nuc_fasta, filepath = output_nuc_fasta_path)
message("nucleotide fasta saved in ", output_nuc_fasta_path)

# protein sequences processing
pro_fasta <- AAStringSet(metadata_with_sequences$pro_sequence)
names(pro_fasta) <- metadata_with_sequences$sequence

# save protein sequences to a fasta file
writeXStringSet(pro_fasta, filepath = output_protein_fasta_path)
message("protein fasta saved in ", output_protein_fasta_path)

# ==============================
# step 7: 6-frame translation of nucleotide sequences
# ==============================
# if you already have known protein sequences in 'protein_sequences', this step may be skipped.

message("translating nucleotide sequences in 6 frames")

# forward frames
fwd_frame1 <- nuc_fasta                  # frame 1: start at position 1
fwd_frame2 <- subseq(nuc_fasta, start=2) # frame 2: start at position 2
fwd_frame3 <- subseq(nuc_fasta, start=3) # frame 3: start at position 3

# reverse frames
rev_comp_fasta <- reverseComplement(nuc_fasta) # reverse complement of all
rev_frame1 <- rev_comp_fasta
rev_frame2 <- subseq(rev_comp_fasta, start=2)
rev_frame3 <- subseq(rev_comp_fasta, start=3)

# translate each frame set
# 'if.fuzzy.codon="X"' marks incomplete codons with 'x'
fwd_protein1 <- translate(fwd_frame1, if.fuzzy.codon="X")
fwd_protein2 <- translate(fwd_frame2, if.fuzzy.codon="X")
fwd_protein3 <- translate(fwd_frame3, if.fuzzy.codon="X")
rev_protein1 <- translate(rev_frame1, if.fuzzy.codon="X")
rev_protein2 <- translate(rev_frame2, if.fuzzy.codon="X")
rev_protein3 <- translate(rev_frame3, if.fuzzy.codon="X")

# rename sequences to indicate frame
original_names <- names(nuc_fasta)
names(fwd_protein1) <- paste0(original_names, "_frame1")
names(fwd_protein2) <- paste0(original_names, "_frame2")
names(fwd_protein3) <- paste0(original_names, "_frame3")
names(rev_protein1) <- paste0(original_names, "_frame4")
names(rev_protein2) <- paste0(original_names, "_frame5")
names(rev_protein3) <- paste0(original_names, "_frame6")

# combine and save all 6-frame translations
all_six_frames <- c(
  fwd_protein1, fwd_protein2, fwd_protein3,
  rev_protein1, rev_protein2, rev_protein3
)

writeXStringSet(all_six_frames, filepath = six_frame_fasta_path)
message("6-frame translation fasta saved in ", six_frame_fasta_path)

# ==============================
# step 8: extract the longest orf directly from nucleotide sequences
# ==============================

# load original nucleotide sequences
message("loading original nucleotide sequences")
nuc_fasta <- readDNAStringSet(output_nuc_fasta_path)

# generate the six-frame nucleotide sequences
message("generating six reading frames")
fwd_frame1 <- nuc_fasta  # frame 1 (starts at position 1)
fwd_frame2 <- subseq(nuc_fasta, start = 2)  # frame 2 (starts at position 2)
fwd_frame3 <- subseq(nuc_fasta, start = 3)  # frame 3 (starts at position 3)
rev_comp_fasta <- reverseComplement(nuc_fasta)
rev_frame1 <- rev_comp_fasta  # frame 4 (starts at position 1, reverse)
rev_frame2 <- subseq(rev_comp_fasta, start = 2)  # frame 5 (starts at position 2, reverse)
rev_frame3 <- subseq(rev_comp_fasta, start = 3)  # frame 6 (starts at position 3, reverse)

# combine all frames into a list
frames_nuc <- list(fwd_frame1, fwd_frame2, fwd_frame3, rev_frame1, rev_frame2, rev_frame3)

# function to extract the longest orf from a given nucleotide frame
extract_longest_orf_nuc <- function(nuc_seq) {
  nuc_seq <- as.character(nuc_seq)
  
  # identify start codons (atg) and stop codons (taa, tag, tga)
  starts <- unlist(gregexpr("ATG", nuc_seq, perl = TRUE))
  stops <- sort(c(
    unlist(gregexpr("TAA", nuc_seq, perl = TRUE)),
    unlist(gregexpr("TAG", nuc_seq, perl = TRUE)),
    unlist(gregexpr("TGA", nuc_seq, perl = TRUE))
  ))
  
  # check if there are any start or stop codons
  if (length(starts) == 0 || length(stops) == 0) {
    return(NA_character_)  # no valid orf found
  }
  
  # initialize variables to track the longest orf
  longest_orf_nuc <- ""
  max_length <- 0
  
  # iterate through each start codon to find the longest orf
  for (start in starts) {
    valid_stops <- stops[stops > start]  # find stop codons after the start
    if (length(valid_stops) > 0) {
      stop_pos <- min(valid_stops) + 2  # include the full stop codon
      orf_nuc <- substr(nuc_seq, start, stop_pos)
      if (nchar(orf_nuc) > max_length) {
        longest_orf_nuc <- orf_nuc
        max_length <- nchar(orf_nuc)
      }
    }
  }
  
  return(longest_orf_nuc)
}

# extract the longest orf per sequence from all six frames
message("extracting longest orf from all six frames")
longest_orfs_nuc <- lapply(frames_nuc, function(frame) {
  sapply(as.character(frame), extract_longest_orf_nuc)
})

# format results into a dataframe
longest_orf_df <- do.call(cbind, longest_orfs_nuc)
longest_orf_df <- data.frame(longest_orf_df, stringsAsFactors = FALSE)

# add sequence ids
longest_orf_df$sequence_id <- names(nuc_fasta)

# find the longest orf across all frames for each sequence
longest_orf_df$longest_orf <- apply(longest_orf_df[,1:6], 1, function(row) {
  row[which.max(nchar(row))]
})

# extract orf lengths
longest_orf_df$orf_length <- nchar(longest_orf_df$longest_orf)

# keep only relevant columns
longest_orf_df <- longest_orf_df[, c("sequence_id", "longest_orf", "orf_length")]

# save longest orfs as fasta
longest_orf_nuc_set <- DNAStringSet(longest_orf_df$longest_orf)
names(longest_orf_nuc_set) <- longest_orf_df$sequence_id

# create a dataframe for additional information if needed
longest_orf_df <- data.frame(
  sequence_id = names(longest_orf_nuc_set),
  orf_sequence = as.character(longest_orf_nuc_set),
  orf_length = nchar(as.character(longest_orf_nuc_set)),
  stringsAsFactors = FALSE
)

# save the longest nucleotide orfs to a fasta file
writeXStringSet(longest_orf_nuc_set, filepath = longest_orf_nuc_fasta_path)
message("longest nucleotide orf fasta saved in ", longest_orf_nuc_fasta_path)

# ==============================
# step 9: translate longest orfs into protein sequences
# ==============================

# translate the longest orfs
longest_orf_protein_set <- translate(longest_orf_nuc_set, if.fuzzy.codon="X")

# save translated protein sequences to fasta file
writeXStringSet(longest_orf_protein_set, filepath = longest_orf_protein_fasta_path)
message("longest orf protein fasta saved in ", longest_orf_protein_fasta_path)
