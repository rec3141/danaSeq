library(ShortRead)
library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

# Define a new transformation for the fourth root
fourth_root_trans <- trans_new(
  "fourth_root",
  transform = function(x) x^(1/4),
  inverse = function(x) x^4
)

# Set your working directory to where your FASTQ files are located
setwd("/data/CMO-20241208/no_sample_id/20241208_1217_X1_FAY90288_1d314608/")
setwd("/data/CMO-20241204")

# Suppose your FASTQ files are named something like "barcode01.fastq.gz", "barcode02.fastq.gz", etc.
fastq_files <- list.files(pattern = "fastq\\.gz$")

extract_info <- function(fname) {
  lines <- readLines(gzfile(fname), warn = FALSE)
  
  # Extract sequence lines (2nd line of each 4-line group)
  seq_lines <- lines[seq(2, length(lines), 4)]
  # Extract quality lines (4th line of each 4-line group)
  qual_lines <- lines[seq(4, length(lines), 4)]
  
  # Calculate read lengths
  read_lengths <- nchar(seq_lines)
  
  # Convert ASCII quality chars to Q-scores and compute mean per read
  avg_qscores <- sapply(qual_lines, function(qual_string) {
    qvals <- utf8ToInt(qual_string) - 33
    mean(qvals)
  })
  
  # Return a data frame with both AvgQ and Length
  data.frame(AvgQ = avg_qscores, Length = read_lengths, stringsAsFactors = FALSE)
}

qscore_data <- data.frame(Barcode = character(), AvgQ = numeric(), Length = numeric(), stringsAsFactors = FALSE)

for (f in fastq_files) {

  print(f)
  # Extract the barcode name from filename
  barcode_name <- sub("\\.fastq\\.gz$", "", f)
  
  info_df <- extract_info(f)
  info_df$Barcode <- barcode_name
  qscore_data <- rbind(qscore_data, info_df)
}


qscore_data_filtered <- qscore_data %>%
  group_by(Barcode) %>%
  filter(n() >= 10000) %>%
  ungroup()


gg = ggplot(qscore_data_filtered, aes(x = Length, y = AvgQ)) +
  facet_wrap(~ Barcode, ncol = 4) +
  geom_hex(bins=100) +
  scale_x_log10() + 
  scale_fill_viridis_c(trans=fourth_root_trans) +  # Optional, for a nice color map
  theme_minimal()

ggsave(filename = "length-qscore3.pdf",device = pdf(), plot = gg, width=12, height=6)




# Create a boxplot of Q-score by Barcode, with jittered points colored by read length
# ggplot(qscore_data, aes(x = Barcode, y = AvgQ)) +
#   geom_boxplot(outlier.shape = NA) +
#  # geom_jitter(aes(color = Length), width = 0.2, alpha = 0.7) +
#   scale_color_viridis_c(option = "plasma") + # Requires "viridis" package if you want a nice color scale
#   theme_minimal() +
#   labs(title = "Q-score Distribution by Barcode",
#        x = "Barcode",
#        y = "Average Q-score per Read") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# gg = ggplot(qscore_data, aes(x = Length, y = AvgQ)) +
#   geom_point(alpha = 0.1, size=0.1, color="blue") +
#   facet_wrap(~ Barcode, nrow = 2) +
#   theme_minimal() +
#   scale_x_log10() +
#   labs(title = "Read Length vs. Q-score by Barcode",
#        x = "Read Length (bases)",
#        y = "Average Q-score")
# 
# ggsave(filename = "length-qscore2.pdf",device = pdf(), plot = gg, width=12, height=6)

