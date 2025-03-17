library(stringr)
start_time <- Sys.time()

codons <- list(
  GCT = 'Ala', GCC = 'Ala', GCA = 'Ala', GCG = 'Ala', CGT = 'Arg', CGC = 'Arg', CGA = 'Arg', CGG = 'Arg', AGA = 'Arg', AGG = 'Arg',
  AAT = 'Asn', AAC = 'Asn', GAT = 'Asp', GAC = 'Asp',
  TGT = 'Cys', TGC = 'Cys', CAA = 'Gln', CAG = 'Gln',
  GAA = 'Glu', GAG = 'Glu', GGT = 'Gly', GGC = 'Gly', GGA = 'Gly', GGG = 'Gly',
  CAT = 'His', CAC = 'His', ATT = 'Ile', ATC = 'Ile', ATA = 'Ile',
  TTA = 'Leu', TTG = 'Leu', CTT = 'Leu', CTC = 'Leu', CTA = 'Leu', CTG = 'Leu',
  AAA = 'Lys', AAG = 'Lys', ATG = 'Met', TTT = 'Phe', TTC = 'Phe',
  CCT = 'Pro', CCC = 'Pro', CCA = 'Pro', CCG = 'Pro',
  TCT = 'Ser', TCC = 'Ser', TCA = 'Ser', TCG = 'Ser', AGT = 'Ser', AGC = 'Ser',
  ACT = 'Thr', ACC = 'Thr', ACA = 'Thr', ACG = 'Thr',
  TGG = 'Trp', TAT = 'Tyr', TAC = 'Tyr', GTT = 'Val', GTC = 'Val', GTA = 'Val', GTG = 'Val',
  TAA = '', TAG = '', TGA = ''
)
# filename <- "sequence1.fasta"
# file_str <- paste(readLines(filename), collapse = "")
# result <- gsub("\\s+", "", file_str)

results <- 'gatgttcgggggcgcctgtgcgaagccgaatgtccttaccctggttaattgtctaagtcagttatttaacatcgtcacgggtgacacgatcgctttgagctattataagtgtcagactgccacttaccgttgatagtggtcaatgtattt'
result <- toupper(results)

cat(paste("The original sequence is\n", result))
newlist <- c()
seq <- ''
for (lines in strsplit(result, "\n")[[1]]) {
    count <- 0
    if (grepl("\n", lines)) {
        newlist <- append(newlist, lines)
    }
}
for (i in 2:(length(newlist) - 1)) {
    seq <- paste0(seq, newlist[i])
}

translation <- function(final_string) {
  library(stringr)
  points <- 0
  translated <- ''
  letters <- ''
  final <- str_replace_all(final_string, ' ', '')
  for (i in 1:nchar(final)) {
    letters <- paste0(letters, substr(final, i, i))
    points <- points + 1
    if (points == 3) {
      translated <- paste0(translated, codons[[letters]], sep= ' ')
      points <- 0
      letters <- ''
    }
  }
  if (final_string == '') {
    return(NULL)
  } else {
    return(translated)
  }
}


#function for adding spaces between codons
spaces <- function(f, result) {
  seq_pos <- ""
  count <- 0
  
  for (i in (f+1):nchar(result)) {
    seq_pos <- paste0(seq_pos, substr(result, i, i))
    count <- count + 1
    if (count == 3) {
      count <- 0
      seq_pos <- paste0(seq_pos, " ")
    }
  }
  return(seq_pos)
}

del_spaces <- function(seq) {
  seqnew <- gsub(" ", "", seq)
  return(seqnew)
}

#function for making 6 reading frames
driver_main <- function(count, result, number) {
  for (f in 0:2) {
    seq_pos <- spaces(f, result)
    if (count == 0) {
      cat(paste("\n\n\nThe +", f + 1, "reading frame is:", seq_pos, "\n"))
    } else if (count == 1) {
      cat(paste("\n\n\nThe -", f + 1, "reading frame is:", seq_pos, "\n"))
    }
    orfs <- c()
    point <- 0
    found <- 1
    final_string <- ''
    for (i in 1:nchar(seq_pos) - 2) {
      if (substr(seq_pos, i, i + 2) == "ATG") {
        found <- found + 1
        for (j in (i + 3):nchar(seq_pos) - 2) {
          if (substr(seq_pos, j, j + 2) %in% c("TAG", "TGA", "TAA")) {
            point <- point + 1
            signal <- match(substr(seq_pos, j, j + 2), c("TAG", "TGA", "TAA"))
            break
          } else {
            final_string <- paste0(final_string, substr(seq_pos, j, j))
          }
        }
        if (found > 1) {
          orfs <- append(orfs, paste0("A", final_string))
          final_string <- ''
        }
      }
    }
    
    if (point != 0) {
      for (i in 1:length(orfs)) {
        orf <- orfs[i]
        new_orf <- del_spaces(orf)
        if (nchar(new_orf) > number) {
          signal_str <- c("TAG", "TGA", "TAA")[signal]
          cat(paste("\nThe ORF", i, "is:\n", orf, signal_str, "\n"))
          cat(paste("The amino acid sequence for ORF", i, "is:\n", translation(orf), "\n"))
        }
      }
    } else if (point == 0) {
      cat("No ORFs were found\n")
    } else if (!grepl("ATG", seq_pos)) {
      cat("No reading frames were found\n")
    }
  }
#   return (" ")
}

count <- 0
cat("\n\nMinimal ORF length 'nt' \n 1. 1 \n 2. 30 \n 3. 75 \n 4. 150 \n 5. 300 \n 6. 600 \n")
number <- as.integer(readline(prompt = "Enter the value: "))
cat("\nThe selected minimal length of ORF is:", number)
cat("\n***************************************************************************************************************************************************************************************")

start_time <- proc.time()[3]
print(driver_main(count, result, number))
cat("\n***************************************************************************************************************************************************************************************")

comp <- chartr("ATGC", "TACG", result)
complement <- gsub("\\s+", "", comp)
cat("\nThe complementary sequence is:\n", complement)
count <- 1
print(driver_main(count, complement, number))

# to print execution time by formatting it to 3 decimals
end_time <- proc.time()[3]
execution_time <-sprintf(end_time - start_time, fmt= '%#.3f')
cat(paste("\n\nExecution time:", execution_time, "seconds"))
