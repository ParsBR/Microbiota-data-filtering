setwd("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados")
dir()
getwd ()

#### Human gut species ####
abund_sp <- read.table("250Sample.humanGut_IGC_11M.species.relativeAbun.table")
#Data origin:
#XIE, Hailiang et al. Shotgun Metagenomics of 250 Adult Twins Reveals Genetic and Environmental Impacts on the Gut Microbiome.
#Cell Systems, [s. l.], v. 3, n. 6, p. 572-584.e3, 2016.
#Available at: http://www.cell.com/article/S2405471216303234/fulltext.
#Access em: 8 apr. 2024.

## Deleting rows with value 0

abund_sp [abund_sp == 0] = NA #Transform 0 into NA

abund_no_row_null = abund_sp [rowSums(is.na(abund_sp[,2:251])) != ncol(abund_sp[,2:251]), ] #Delete rows where sum is NA (delete rows where all data is NA)
exist_sp = abund_no_row_null$V1


#######################################################################################################################################
#### Hypertension ####

hypert_sp_original <- read.csv2("Generos_microbiota_humana_hipertensao.csv") #Open csv of genus
#Data origin: Disbiome Database
#JANSSENS, Yorick et al. Disbiome database: Linking the microbiome to disease.
#BMC Microbiology, [s. l.], v. 18, n. 1, p. 1-6, 2018.
#Available at: https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-018-1197-5.
#Access: 8 apr. 2024.

hypert_sp <- hypert_sp_original$genus

increase <- hypert_sp_original[which(hypert_sp_original$status=="+"),]
increase_genus <- increase$genus

decrease <- hypert_sp_original[which(hypert_sp_original$status=="-"),]
decrease_genus <- decrease$genus


#######################################################################################################################################
#### Acetate kinase####

acetate_kinase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Acetato_quinase.csv")

filtered_p1 <- acetate_kinase[acetate_kinase$specie %in% exist_sp, ]
#write.csv2(filtered_p1, "acetate_kinase_filtered.csv", row.names = FALSE)

filtered_hyper_1 <- filtered_p1[filtered_p1$genus %in% hypert_sp, ]

increased_1 <- filtered_p1[filtered_p1$genus %in% increase_genus, ]
decreased_1 <- filtered_p1[filtered_p1$genus %in% decrease_genus, ]


#### Acetyl-CoA sinthetase ####

acetyl_coa_synthetase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Acetil-CoA_sintetase.csv")

filtered_p3 <- acetyl_coa_synthetase[acetyl_coa_synthetase$specie %in% exist_sp, ]
#write.csv2(filtered_p3, "acetyl_coa_synthetase_filtered.csv", row.names = FALSE)

filtered_hyper_3 <- filtered_p3[filtered_p3$genus %in% hypert_sp, ]

increased_3 <- filtered_p3[filtered_p3$genus %in% increase_genus, ]
decreased_3 <- filtered_p3[filtered_p3$genus %in% decrease_genus, ]


#### Butyryl-CoA:acetate_CoA_transferase ####

butyryl_coa_acetate_coa_transferase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Butiril-CoA-acetato_CoA_transferase.csv")
butyryl_coa_acetate_coa_transferase <- butyryl_coa_acetate_coa_transferase[1:65,]

filtered_p5 <- butyryl_coa_acetate_coa_transferase[butyryl_coa_acetate_coa_transferase$specie %in% exist_sp, ]
#write.csv2(filtered_p5, "butyryl_coa_acetate_coa_transferase_filtered.csv", row.names = FALSE)

filtered_hyper_5 <- filtered_p5[filtered_p5$genus %in% hypert_sp, ]

increased_5 <- filtered_p5[filtered_p5$genus %in% increase_genus, ]
decreased_5 <- filtered_p5[filtered_p5$genus %in% decrease_genus, ]


#### Betaine reductase ####

betaine_reductase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Betaina_redutase.csv")
betaine_reductase <- betaine_reductase[1:112,]

filtered_p7 <- betaine_reductase[betaine_reductase$specie %in% exist_sp, ]

filtered_hyper_7 <- filtered_p7[filtered_p7$genus %in% hypert_sp, ]

increased_7 <- filtered_p7[filtered_p7$genus %in% increase_genus, ]
decreased_7 <- filtered_p7[filtered_p7$genus %in% decrease_genus, ]

betaine_beta <-filtered_p7[which(filtered_p7$protein_name=="Betaine_reductase_complex_component_B_sub_beta"),]
betaine_alfa <-filtered_p7[which(filtered_p7$protein_name=="Betaine_reductase_complex_component_B_sub_alpha"),]
#write.csv2(betaine_beta, "betaine_beta_filtered.csv", row.names = FALSE)
#write.csv2(betaine_alfa, "betaine_alfa_filtered.csv", row.names = FALSE)

filtered_hyper_beta <- betaine_beta[betaine_beta$genus %in% hypert_sp, ]
filtered_hyper_alfa <- betaine_alfa[betaine_alfa$genus %in% hypert_sp, ]

increased_beta <- betaine_beta[betaine_beta$genus %in% increase_genus, ]
increased_alfa <- betaine_alfa[betaine_alfa$genus %in% increase_genus, ]
decreased_beta <- betaine_beta[betaine_beta$genus %in% decrease_genus, ]
decreased_alfa <- betaine_alfa[betaine_alfa$genus %in% decrease_genus, ]


#### Butyrate kinase ####

butyrate_kinase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Butirato_quinase.csv")

filtered_p9 <- butyrate_kinase[butyrate_kinase$specie %in% exist_sp, ]
#write.csv2(filtered_p9, "butyrate_kinase_filtered.csv", row.names = FALSE)

filtered_hyper_9 <- filtered_p9[filtered_p9$genus %in% hypert_sp, ]

increased_9 <- filtered_p9[filtered_p9$genus %in% increase_genus, ]
decreased_9 <- filtered_p9[filtered_p9$genus %in% decrease_genus, ]


#### Ergothionase ####

ergothionase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Ergotionase.csv")

filtered_p11 <- ergothionase[ergothionase$specie %in% exist_sp, ]

filtered_hyper_11 <- filtered_p11[filtered_p11$genus %in% hypert_sp, ]

increased_11 <- filtered_p11[filtered_p11$genus %in% increase_genus, ]
decreased_11 <- filtered_p11[filtered_p11$genus %in% decrease_genus, ]


#### Propionate kinase ####

propionate_kinase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Propionato_quinase.csv")

filtered_p13 <- propionate_kinase[propionate_kinase$specie %in% exist_sp, ]
#write.csv2(filtered_p13, "propionate_kinase_filtered.csv", row.names = FALSE)

filtered_hyper_13 <- filtered_p13[filtered_p13$genus %in% hypert_sp, ]

increased_13 <- filtered_p13[filtered_p13$genus %in% increase_genus, ]
decreased_13 <- filtered_p13[filtered_p13$genus %in% decrease_genus, ]


#### Propionyl_CoA_transferase ####

propionyl_coa_transferase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Propionil-CoA_transferase.csv")

filtered_p15 <- propionyl_coa_transferase[propionyl_coa_transferase$specie %in% exist_sp, ]
#write.csv2(filtered_p15, "propionyl_coa_transferase_filtered.csv", row.names = FALSE)

filtered_hyper_15 <- filtered_p15[filtered_p15$genus %in% hypert_sp, ]

increased_15 <- filtered_p15[filtered_p15$genus %in% increase_genus, ]
decreased_15 <- filtered_p15[filtered_p15$genus %in% decrease_genus, ]


#### TMAO reductase ####

TMAO_reductase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/TMAO_redutase.csv")

filtered_p17 <- TMAO_reductase[TMAO_reductase$specie %in% exist_sp, ]

filtered_hyper_17 <- filtered_p17[filtered_p17$genus %in% hypert_sp, ]

TorA <-filtered_p17[which(filtered_p17$protein_name=="Trimethylamine_N-oxide_reductase_TorA "),]
TorZ <-filtered_p17[which(filtered_p17$protein_name=="Trimethylamine_N-oxide_reductase_TorZ"),]
#write.csv2(TorA, "TorA_filtered.csv", row.names = FALSE)
#write.csv2(TorZ, "TorZ_filtered.csv", row.names = FALSE)

filtered_hyper_TorA <- TorA[TorA$genus %in% hypert_sp, ]
filtered_hyper_TorZ <- TorZ[TorZ$genus %in% hypert_sp, ]

increased_17 <- filtered_p17[filtered_p17$genus %in% increase_genus, ]
decreased_17 <- filtered_p17[filtered_p17$genus %in% decrease_genus, ]

increased_TorA <- TorA[TorA$genus %in% increase_genus, ]
increased_TorZ <- TorZ[TorZ$genus %in% increase_genus, ]
decreased_TorA <- TorA[TorA$genus %in% decrease_genus, ]
decreased_TorZ <- TorZ[TorZ$genus %in% decrease_genus, ]


#### Propionyl-CoA:succinate CoA transferase ####

propionyl_coa_succinate_coa_transferase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Propionil-CoA-succinato_CoA_transferase.csv")

filtered_p19 <- propionyl_coa_succinate_coa_transferase[propionyl_coa_succinate_coa_transferase$specie %in% exist_sp, ]
#write.csv2(filtered_p19, "propionyl_coa_succinate_coa_transferase_filtered.csv", row.names = FALSE)

filtered_hyper_19 <- filtered_p19[filtered_p19$genus %in% hypert_sp, ]

increased_19 <- filtered_p19[filtered_p19$genus %in% increase_genus, ]
decreased_19 <- filtered_p19[filtered_p19$genus %in% decrease_genus, ]


#### Choline TMA-lyase ####

choline_tma_lyase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Colina_TMA_liase.csv")

filtered_p21 <- choline_tma_lyase[choline_tma_lyase$specie %in% exist_sp, ]
#write.csv2(filtered_p21, "choline_tma_lyase_filtered.csv", row.names = FALSE)

filtered_hyper_21 <- filtered_p21[filtered_p21$genus %in% hypert_sp, ]

increased_21 <- filtered_p21[filtered_p21$genus %in% increase_genus, ]
decreased_21 <- filtered_p21[filtered_p21$genus %in% decrease_genus, ]


#### Carnitine monooxygenase ####

carnitine_monooxygenase <- read.csv2("D:/Ambiente de trabalho/UFSC Mestrado/Pesquisa/Dados/Proteinas_especies_prontas/Carnitina_monooxigenase.csv")

filtered_p23 <- carnitine_monooxygenase[carnitine_monooxygenase$specie %in% exist_sp, ]

filtered_hyper_23 <- filtered_p23[filtered_p23$genus %in% hypert_sp, ]

increased_23 <- filtered_p23[filtered_p23$genus %in% increase_genus, ]
decreased_23 <- filtered_p23[filtered_p23$genus %in% decrease_genus, ]

CntA <-filtered_p23[which(filtered_p23$protein_name=="Carnitine_monooxygenase_subunit_YeaW_(CntA)_(oxygenase)"),]
CntB <-filtered_p23[which(filtered_p23$protein_name=="Carnitine_monooxygenase_subunit_YeaX_(CntB)_(reductase)"),]
#write.csv2(CntA, "CntA_filtered.csv", row.names = FALSE)
#write.csv2(CntB, "CntB_filtered.csv", row.names = FALSE)

filtered_hyper_CntA <- CntA[CntA$genus %in% hypert_sp, ]
filtered_hyper_CntB <- CntB[CntB$genus %in% hypert_sp, ]

increased_CntA <- CntA[CntA$genus %in% increase_genus, ]
increased_CntB <- CntB[CntB$genus %in% increase_genus, ]
decreased_CntA <- CntA[CntA$genus %in% decrease_genus, ]
decreased_CntB <- CntB[CntB$genus %in% decrease_genus, ]


#####################################################################################################################################

#Merging filtered data
library(tidyverse)
tabela_completa <- bind_rows(filtered_p1,filtered_p3,filtered_p5,filtered_p7,filtered_p9,filtered_p13,filtered_p15,filtered_p17,filtered_p19,filtered_p21,filtered_p23)

#Data Overview
panorama_dados  <- tabela_completa %>%
  filter(!is.na(aa_quantity)) %>% 
  group_by(protein_name) %>%
  summarise(Media_aa = mean(aa_quantity),
            Mediana_aa = median(aa_quantity),
            Min_aa = min(aa_quantity),
            Max_aa = max(aa_quantity),
            Desvio_padrao = sd(aa_quantity),
            Quantidade_de_dados = n())

#Standard deviation indicates how uniform the data set is.

#write.csv2(panorama_dados, "panorama_dados_summarise.csv", row.names = FALSE)

#########################################################################################################################################
#### Table of proteins producing SCFAs and TMA by species ####
library(tidyverse)
species_AGCC <- bind_rows(filtered_p1,filtered_p3,filtered_p5,filtered_p9,filtered_p13,filtered_p15,filtered_p19)
prot_p_sp_AGCC  <- species_AGCC %>%
  group_by(specie) %>%
  summarise(Proteinas_AGCCs = n())


bet <- distinct(filtered_p7,specie, .keep_all = TRUE)
tor <- distinct(filtered_p17,specie, .keep_all = TRUE)
cnt <- distinct(filtered_p23,specie,.keep_all = TRUE)

species_TMA <- bind_rows(bet,tor,filtered_p21,cnt)
prot_p_sp_TMA  <- species_TMA %>%
  group_by(specie) %>%
  summarise(Proteinas_TMA = n())

prot_p_sp_tot <- prot_p_sp_AGCC %>% full_join(prot_p_sp_TMA)
#write.csv2(prot_p_sp_tot, "Proteina_por_especie.csv", row.names = FALSE)

genus_species <- tabela_completa[,c(1,2)] #select specie and genus columns from all filtered data
genus_status_hyper <- hypert_sp_original[,c(1,3)] #selecting columns genus and status

library(dplyr)
genus_species_new <- distinct(genus_species) #removes repeated lines
genus_status_hyper_new <- distinct(genus_status_hyper)

include_genus <- prot_p_sp_tot %>% left_join(genus_species_new)
choose_hyper_genus <- include_genus %>% inner_join(genus_status_hyper_new) #only related to hypertension
choose_hyper_genus_att <- choose_hyper_genus %>%
  arrange(specie)
write.csv2(choose_hyper_genus_att, "Proteina_por_especie_hyper.csv", row.names = FALSE)
#Species of the genus Bacteroides filtered for hypertension, except B. salyersiae and B. thetaiotaomicron,
#were disregarded regarding the change in abundance in hypertension (manual removal). Reason: Only the two
#species mentioned of the genus are registered with reduced abundance in the Disbiome database. The bacteria
#of the other filtered genera were maintained, since the repository reported changes in abundance related to
#the genus, not to determined species.
