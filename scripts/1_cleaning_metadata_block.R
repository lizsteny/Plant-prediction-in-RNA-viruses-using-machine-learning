# ==============================
# load and process input metadata
# ==============================

# load required libraries
library(Biostrings)
library(tidyr)
library(dplyr)
library(conflicted)
library(stringr)

conflict_prefer("rename", "dplyr")
conflict_prefer("first", "dplyr")
conflicts_prefer(dplyr::filter)

# define file paths
virusdb_path <- "data/raw/virushostdb.csv"
metadata_path <- "data/raw/raw_metadata.csv"
protein_fasta_path <- "data/raw/raw_protein_seq.fasta"
output_path <- "data/metadata/metadata.csv"
y_path <-"data/y.csv"

# ==============================
# step 1: load data
# ==============================

# ensure input files exist before proceeding
if (!file.exists(virusdb_path)) stop("error: virus-host database file not found!")
if (!file.exists(metadata_path)) stop("error: metadata file not found!")

# Load virus-host association database
message("loading virus-host association")
virusdb <- read.csv(virusdb_path)


# Load metadata containing viral genome information
message("loading viral genome metadata")
metadataNCBI <- read.csv(metadata_path)

# ==============================
# step 2: process virus tax IDs
# ==============================

# - select relevant columns:
#     - `virus.tax.id` (column 4)
#     - `associated.IDs` (column 1) → Corresponds to `refseq.id`
# - expand rows where multiple associated IDs exist (separated by commas)
# - rename the column for clarity
taxids <- virusdb %>%
  select(virus.tax.id, associated.IDs) %>%
  separate_rows(associated.IDs, sep = ",", convert = TRUE) %>%
  rename(refseq.id = associated.IDs)

selected_df <- select(virusdb, virus.tax.id, associated.IDs)
print(head(selected_df))

separated_df <- separate_rows(selected_df, associated.IDs, sep = ",")
print(head(separated_df))


# ==============================
# step 3: merge metadata with taxonomic information
# ==============================

# - perform a left join to retain all `metadata_ncbi` rows
# - fill missing `virus.tax.id` values using `Assembly` and `Organism_Name` group-wise
merged_df <- metadataNCBI %>%
  left_join(taxids, by = "refseq.id") %>%
  group_by(Assembly) %>%
  fill(virus.tax.id, .direction = "downup") %>%
  ungroup() %>%
  group_by(Organism_Name) %>%
  fill(virus.tax.id, .direction = "downup")

# ==============================
# step 4: clean dataset
# ==============================

# - remove duplicate rows
# - remove rows where `Assembly` is missing, empty, or not marked as "complete"
# - drop the 4th column (`Nuc_Completeness`) as it is no longer needed
# - remove entries where `Assembly` starts with "set:"
# - aggregate data:
#     - combine multiple `refseq.id` values per (`virus.tax.id`, `Assembly`)
cleaned_df <- merged_df %>%
  distinct() %>%
  filter(!is.na(Assembly), Assembly != "", Nuc_Completeness == "complete") %>%
  select(-4) %>%  # drop nuc completeness col
  filter(!grepl("^set:", Assembly)) %>%
  group_by(virus.tax.id, Assembly) %>%
  summarise(
    refseq.id = paste(unique(refseq.id), collapse = ", "),  
    .groups = "drop"
  )

# ==============================
# step 5: process virus-host associations
# ==============================

# - select relevant columns:
#     - `virus.tax.id` (1)
#     - `host.tax.id` (5)
#     - `host.name` (6)
#     - `host.lineage` (7)
#     - `evidence` (9)
# - remove rows where `host.lineage` is empty or missing
# - create a binary column `infection`:
#     - 1 if "Viridiplantae" appears in `host.lineage` (indicating a plant host)
#     - 0 otherwise
# - group by `virus.tax.id` and `infection` to separate plant and non-plant hosts
# - concatenate **unique** `host.tax.id` values into a single string per `virus.tax.id`
# - concatenate **unique** `evidence` values into a single string per `virus.tax.id`
virusdb_hosts <- virusdb %>%
  select(1, 5, 6, 7, 9) %>%  # Select key columns
  filter(host.name != "", !is.na(host.lineage), host.lineage != "") %>%  
  mutate(infection = ifelse(grepl("Viridiplantae", host.lineage, ignore.case = TRUE), 1, 0)) %>%
  group_by(virus.tax.id, infection) %>%
  summarise(
    host.tax.id = paste(unique(host.tax.id), collapse = ", "),  
    evidence = paste(unique(evidence), collapse = ", "),  
    .groups = "drop"
  )

# ==============================
# step 6: merge metadata with virus-host associations
# ==============================

metadata <- merge(cleaned_df, virusdb_hosts, by = "virus.tax.id") %>%
  rename_with(tolower)

# identify `virus.tax.id` values that have both `infection = 0` and `infection = 1`
dual_infection_ids <- metadata %>%
  group_by(virus.tax.id) %>%
  summarize(infection_types = n_distinct(infection)) %>%
  filter(infection_types > 1) %>%
  pull(virus.tax.id)

# - remove rows where `virus.tax.id` has both `infection = 0` and `infection = 1`
metadata <- metadata %>%
  filter(!virus.tax.id %in% dual_infection_ids)

# ==============================
# step 7: process and integrate protein sequence data
# ==============================

# ensure protein FASTA file exists before proceeding
if (!file.exists(protein_fasta_path)) stop("Error: Protein FASTA file not found!")

# read protein FASTA file
message("Loading protein sequences...")
protein_sequences <- readAAStringSet(protein_fasta_path)

# extract sequence headers (names) from FASTA
protein_headers <- names(protein_sequences)

# split headers to extract protein ID and assembly
# - Assumes the format: "ProteinID | AssemblyID"
split_headers <- strsplit(protein_headers, " \\|")

# - create a dataframe with `protein_id` and `assembly_id`
protein_assembly_df <- do.call(rbind, lapply(split_headers, function(x) { 
  data.frame(protein_id = x[1], assembly_id = x[2], stringsAsFactors = FALSE)
}))

# - filter to keep only protein assemblies that match those in `metadata`
filtered_protein_df <- protein_assembly_df %>%
  filter(assembly_id %in% metadata$assembly)

# - group by `assembly_id` and concatenate multiple protein IDs into a single string
grouped_protein_df <- filtered_protein_df %>%
  group_by(assembly_id) %>%
  summarise(protein_id = str_c(protein_id, collapse = ","), .groups = "drop")

# merge grouped protein data with `metadata`, retaining all `metadata` rows
metadata_with_protein <- merge(metadata, grouped_protein_df, by.x = "assembly", by.y = "assembly_id", all.x = TRUE)

# remove rows where `protein_id` is NA (no associated proteins)
cleaned_metadata <- metadata_with_protein %>%
  filter(!is.na(protein_id))

# order the dataframe by `assembly`
ordered_metadata <- cleaned_metadata %>%
  arrange(assembly)

# add a sequential column `sequence_id` with "seq1", "seq2", ...
ordered_metadata$sequence <- paste0("seq", seq_len(nrow(ordered_metadata)))

# reorder columns for consistency
ordered_metadata <- ordered_metadata %>%
  select(sequence, assembly, infection, virus.tax.id, host.tax.id, evidence, refseq.id, protein_id)

# ==============================
# step 8: create dataframe with sequences and infections labels
# ==============================

# save y as CSV
yvariable <- ordered_metadata %>%
  select(1, 3) 

# ==============================
# step 9: save processed data
# ==============================

# save final processed metadata as CSV
write.csv(ordered_metadata, output_path, row.names = FALSE)
message("metadata saved in ", output_path)

write.csv(yvariable, y_path, row.names = FALSE)
message("y saved in ", y_path)