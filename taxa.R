library(Rcpp)
library(tidyr)
library(tidyverse)
library(dplyr)
library(dada2)
library(ggplot2)

#ruta

path <- "C:/Users/Admin/Desktop/ptmr/Bioinformatica/secuencias/Proyecto"
fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#filtrado y trim

filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

names(filtFs) <- sample.names


plotQualityProfile(fnFs[1:2])

out <- filterAndTrim(fnFs, filtFs,
                     truncLen=250,  
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=FALSE)
head(out)
#Inferencia de errores y secuencias

errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs[[1]]

#construir tabla de secuencias y eliminar quimeras
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#asignar taxonomía
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Admin/Desktop/ptmr/Bioinformatica/secuencias/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa)

delgados <- c("SRR2912454.fastq.gz", "SRR2912476.fastq.gz", "SRR2912449.fastq.gz", "SRR2912450.fastq.gz", "SRR2912477.fastq.gz")
obesos <- c("SRR2912458.fastq.gz", "SRR2912473.fastq.gz", "SRR2912447.fastq.gz", "SRR2912445.fastq.gz", "SRR2912475.fastq.gz")
perdieron <- c("SRR2912463.fastq.gz", "SRR2912467.fastq.gz", "SRR2912470.fastq.gz", "SRR2912465.fastq.gz", "SRR2912466.fastq.gz")
subieron <- c("SRR2912452.fastq.gz", "SRR2912455.fastq.gz", "SRR2912453.fastq.gz", "SRR2912456.fastq.gz", "SRR2912457.fastq.gz")

muestras <- rownames(seqtab.nochim)
grupo <- ifelse(muestras %in% delgados, "Delgados",
                ifelse(muestras %in% obesos, "Obesos",
                       ifelse(muestras %in% perdieron, "Obesos tratados",
                              ifelse(muestras %in% subieron, "Delgados tratados", "Otro"))))
library(tidyverse)

df <- as.data.frame(seqtab.nochim)
df$Sample <- rownames(df)
df$Grupo <- grupo

df_long <- df %>%
  pivot_longer(cols = -c(Sample, Grupo), names_to = "Secuencia", values_to = "Abundancia")

df_grouped <- df_long %>%
  group_by(Grupo, Secuencia) %>%
  summarise(Total = sum(Abundancia), .groups = "drop")
taxa_df <- as.data.frame(taxa)
taxa_df$Secuencia <- rownames(taxa_df)

df_final <- left_join(df_grouped, taxa_df, by = "Secuencia")
df_final %>%
  group_by(Grupo, Family) %>%
  summarise(Reads = sum(Total), .groups = "drop") %>%
  group_by(Grupo) %>%
  mutate(AbRel = Reads / sum(Reads)) %>%
  ggplot(aes(x = Grupo, y = AbRel, fill = Family)) +
  geom_bar(stat = "identity") +
  ylab("Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


top_familias <- df_final %>%
  group_by(Family) %>%
  summarise(TotalReads = sum(Total), .groups = "drop") %>%
  arrange(desc(TotalReads)) %>%
  slice_head(n = 10) %>%
  pull(Family)

df_top10 <- df_final %>%
  filter(Family %in% top_familias)

df_top10 %>%
  group_by(Grupo, Family) %>%
  summarise(Reads = sum(Total), .groups = "drop") %>%
  group_by(Grupo) %>%
  mutate(AbRel = Reads / sum(Reads)) %>%
  ggplot(aes(x = Grupo, y = AbRel, fill = Family)) +
  geom_bar(stat = "identity") +
  ylab("Abundancia relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_adulto_sano <- data.frame(
  Grupo = rep("Adulto sano", 9),
  Familia = c("Prevotellaceae", "Peptostreptococcaceae", "Pseudohongiellaceae",
              "Bacteroidaceae", "Veillonellaceae", "Rikenellaceae",
              "Lachnospiraceae", "Tannerellaceae", "Ruminococcaceae"),
  Reads = c(509, 39, 37, 33, 20, 15, 15, 10, 10)
)

df_adulto_sano$AbundanciaRelativa <- (df_adulto_sano$Reads / sum(df_adulto_sano$Reads)) * 100

ggplot(df_adulto_sano, aes(x = Grupo, y = AbundanciaRelativa, fill = Familia)) +
  geom_bar(stat = "identity") +
  labs(title = "Composición de familias",
       x = "", y = "Abundancia relativa (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
