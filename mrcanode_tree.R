library(ape)
library(phytools)

tr <- read.nexus("finalbirdsp.tree")


#scaling node age ####
#for our tree

tips_to_keep <- read.table("./tips_to_keep.csv", quote="\"", comment.char="")

tips_to_rm <- setdiff(tr$tip.label, tips_to_keep$V1)

tr <-  drop.tip(tr, tips_to_rm)

scale_factor <- 55.2686 / max(branching.times(tr))
# 
tr$edge.length <- tr$edge.length * scale_factor
# 
max(branching.times(tr))

#write.tree(tr, "full_scaled_ut.tree")

#for the jetz tree

#Plotting fan species-level tree ####
#Getting MRCA nodes for each family for plotting arc labels
tip_source <- list(
  ACANTHIZIDAE=c('Acanthiza_cinerea', 'Pachycare_flavogriseum'),
  PARDALOTIDAE=c('Pardalotus_striatus', 'Pardalotus_rubricatus'),
  DASYORNITHIDAE=c('Dasyornis_broadbenti', 'Dasyornis_brachypterus'),
  MALURIDAE=c('Malurus_alboscapulatus', 'Amytornis_barbatus'),
  POMATOSTOMIDAE=c('Pomatostomus_superciliosus', 'Pomatostomus_ruficeps'),
  PTILONORHYNCHIDAE=c('Ailuroedus_buccoides', 'Ptilonorhynchus_violaceus'),
  CLIMACTERIDAE=c('Climacteris_picumnus','Cormobates_leucophaea'),
  ATRICHORNITHIDAE=c('Atrichornis_rufescens', 'Atrichornis_clamosus'),
  MENURIDAE=c('Menura_novaehollandiae', 'Menura_alberti'),
  PITTIDAE=c('Erythropitta_granatina','Pitta_superba'),#'Pitta_erythrogaster',
  CALYPTOMENIDAE=c('Calyptomena_viridis', 'Calyptomena_hosii'),
  EURYLAIMIDAE=c('Eurylaimus_ochromalus', 'Corydon_sumatranus'),
  CINCLOSOMATIDAE=c('Cinclosoma_punctatum','Ptilorrhoa_leucosticta','Cinclosoma_marginatum'),
  CAMPEPHAGIDAE=c('Edolisoma_tenuirostre','Pericrocotus_divaricatus', 'Pericrocotus_montanus'),
  NEOSITTIDAE=c('Daphoenositta_chrysoptera', 'Daphoenositta_miranda'),
  PSOPHODIDAE=c('Psophodes_cristatus', 'Androphobus_viridis'),
  OREOICIDAE=c('Oreoica_gutturalis', 'Ornorectes_cristatus'),
  ORIOLIDAE=c('Oriolus_chinensis','Pitohui_dichrous','Sphecotheres_hypoleucus'),
  PARAMYTHIIDAE=c('Oreocharis_arfaki', 'Paramythia_olivacea'),
  VIREONIDAE=c('Erpornis_zantholeuca', 'Pteruthius_aeralatus'),
  MACHAERIRHYNCHIDAE=c('Machaerirhynchus_nigripectus','Machaerirhynchus_flaviventer'),
  ARTAMIDAE=c('Artamus_cinereus', 'Melloria_quoyi', 'Peltops_blainvillii', 'Strepera_graculina'),
  AEGITHINIDAE=c('Aegithina_lafresnayei', 'Aegithina_viridissima'),
  VANGIDAE=c('Tephrodornis_virgatus', 'Philentoma_velata'),
  RHIPIDURIDAE=c('Rhipidura_javanica','Chaetorhynchus_papuensis','Rhipidura_tenebrosa'),
  DICRURIDAE=c('Dicrurus_aeneus', 'Dicrurus_sumatranus'),
  MONARCHIDAE=c('Hypothymis_azurea', 'Monarcha_megarhynchus', 'Arses_insularis', 'Myiagra_eichhorni'),
  PARADISAEIDAE=c('Paradisaea_minor', 'Phonygammus_keraudrenii'),
  CORCORACIDAE=c('Struthidea_cinerea', 'Corcorax_melanorhamphos'),
  MELAMPITTIDAE=c('Melampitta_lugubris', 'Megalampitta_gigantea'),
  CNEMOPHILIDAE=c('Cnemophilus_loriae', 'Loboparadisea_sericea'),
  MELANOCHARITIDAE=c('Melanocharis_versteri', 'Oedistoma_iliolophus', 'Toxorhampus_novaeguineae'),
  PETROICIDAE=c('Devioeca_papuana', 'Petroica_multicolor', 'Amalocichla_incerta'),
  STENOSTIRIDAE=c('Chelidorhynx_hypoxanthus', 'Culicicapa_ceylonensis'),
  ORTHONYCHIDAE=c('Orthonyx_temminckii', 'Orthonyx_novaeguineae'),
  ACROCEPHALIDAE=c('Acrocephalus_orientalis', 'Arundinax_aedon'),
  LOCUSTELLIDAE=c('Locustella_lanceolata', 'Megalurulus_grosvenori', 'Poodytes_gramineus'),#'Robsonius_rabori', 
  HIRUNDINIDAE=c('Hirundo_rustica', 'Delichon_lagopodum', 'Cecropis_daurica'),
  PYCNONOTIDAE=c('Pycnonotus_jocosus', 'Alophoixus_pallidus'),
  ZOSTEROPIDAE=c('Zosterops_everetti', 'Yuhina_castaniceps'),
  TIMALIIDAE=c('Timalia_pileata', 'Stachyris_grammiceps'),
  PELLORNEIDAE=c('Trichastoma_pyrrogenys', 'Turdinus_atrigularis'), #'Leonardina_woodi'
  SCOTOCERCIDAE=c('Abroscopus_albogularis','Tesia_cyaniventer', 'Tesia_olivea'),
  DICAEIDAE=c('Dicaeum_aureolimbatum', 'Prionochilus_maculatus'), #'Dicaeum_hypoleucum', 'Dicaeum_bicolor'
  NECTARINIIDAE=c('Leptocoma_sperata', 'Cinnyris_idenburgi'),
  #IRENIDAE=c('Irena_puella'),#'Irena_cyanogastra', 
  CHLOROPSEIDAE=c('Chloropsis_cochinchinensis', 'Chloropsis_sonnerati'),
  PASSERIDAE=c('Passer_domesticus', 'Passer_montanus' ),#'Hypocryptadius_cinnamomeus'
  MOTACILLIDAE=c('Motacilla_alba', 'Anthus_godlewskii'),
  FRINGILLIDAE=c('Fringilla_montifringilla', 'Pyrrhula_nipalensis'),#'Pyrrhula_leucogenis'
  MELIPHAGIDAE=c('Meliphaga_lewinii','Acanthorhynchus_tenuirostris', 'Myzomela_cruentata'),
  PACHYCEPHALIDAE=c('Colluricincla_harmonica','Pachycephala_simplex'),
  CORVIDAE=c('Corvus_coronoides', 'Cissa_chinensis', 'Platysmurus_leucopterus'),
  LANIIDAE=c('Lanius_schach', 'Lanius_cristatus'),
  PARIDAE=c('Parus_cinereus', 'Melanochlora_sultanea'), #'Pardaliparus_elegans'
  ALAUDIDAE=c('Alauda_gulgula', 'Mirafra_javanica'),
  LEIOTHRICHIDAE=c('Leiothrix_laurinae', 'Garrulax_castanotis'),
  PHYLLOSCOPIDAE=c('Phylloscopus_trivirgatus', 'Phylloscopus_trochiloides', 'Phylloscopus_fuscatus'),
  TURDIDAE=c('Cochoa_viridis', 'Zoothera_aurea'), #'Turdus_chrysolaus', 
  MUSCICAPIDAE=c('Muscicapa_sibirica', 'Larvivora_ruficeps', 'Saxicola_maurus'),
  STURNIDAE=c('Sturnia_malabarica', 'Aplonis_minor'), #'Rhabdornis_inornatus'
  SITTIDAE=c('Sitta_frontalis', 'Sitta_neglecta'),
  PLOCEIDAE=c('Ploceus_philippinus', 'Ploceus_manyar'),
  ESTRILDIDAE = c('Poephila_personata', 'Lonchura_leucogastra', 'Erythrura_trichroa'),
  EMBERIZIDAE = c('Emberiza_fucata', 'Emberiza_rutila'),#'Emberiza_sulphurata'
  CISTICOLIDAE = c('Orthotomus_sepium','Cisticola_exilis'),#'Orthotomus_castaneiceps', 
  SYLVIIDAE = c('Cholornis_unicolor', 'Paradoxornis_guttaticollis')
)

mrca_nd <- data.frame(Family = character(), MRCA_Node = integer(), stringsAsFactors = FALSE)

# Iterate through each family
for (family in names(tip_source)) {
  tip_names <- tip_source[[family]]
  
  # Find the MRCA for the current family
  mrca <- getMRCA(tr, tip_names)
  
  # Add the results to the dataframe
  mrca_nd <- rbind(mrca_nd, data.frame(Family = family, MRCA_Node = mrca, stringsAsFactors = FALSE))
}

mrca_nd <- mrca_nd[order(mrca_nd$Family), ]

nodes <- mrca_nd$MRCA_Node
labels <- mrca_nd$Family

tr <- force.ultrametric(tr, method = "extend")
is.ultrametric(tr)
tr <- read.tree("full_scaled_ut.tree")

dev.off()
pdf(file = "bird_tree_new.pdf", paper = "a4r")

plotTree(tr, type = "fan", lwd=0.5, ftype="off", part = 0.94)
axisPhylo(las = 2)

#arc.cladelabels(tree = tr, text = as.character(mrca_nd$Family), node = mrca_nd$MRCA_Node, cex = 0.4)

for (i in 1:length(mrca_nd$Family)) {
  arc.cladelabels(text = mrca_nd$Family[i], node = mrca_nd$MRCA_Node[i], lab.offset = 1.05, cex = 0.4, orientation = "horizontal")
}

dev.off()

#check that grafted and ungrafted branch lengths correlate ########
tree_folder <- "./sub_results"

tree_files <- list.files(path = tree_folder, pattern = "\\.tree$", full.names = TRUE)

phy_trees <- list()

outgroup <- c("Melopsittacus_undulatus", "Casuarius_casuarius", 
              "Dromaius_novaehollandiae", "Apteryx_australis",
              "Chelonia_mydas", "Testudo_graeca", "Alligator_sinensis", 
              "Crocodylus_porosus", "Varanus_komodoensis", "Caretta_caretta")

for (file in tree_files) {
  print(paste("Processing file:", file))
  tree <- read.tree(file)
  tree <- drop.tip(tree, outgroup)
  phy_trees[[basename(file)]] <- tree
}

mrca_ls <- mrca_nd$MRCA_Node

for (i in seq_along(nodes)) {
  tryCatch(
    {
      # Extract the subtree for the current node
      subtree <- extract.clade(tr, nodes[i])
      
      filename <- paste("figures/", tree_name, ".png", sep = "")
      png(filename)
      
      # Create a new plot
      plot_obj <- plot(phy_trees[[i]]$edge.length ~ subtree$edge.length,
                       xlab = "Grafted Edge Length", ylab = "Ungrafted Edge Length")
      
      # Add a title or labels as needed
      tree_name <- gsub("\\.tree$", "", names(phy_trees)[i])
      title(main = tree_name)
      
            print(plot_obj)
      dev.off()
    },
    error = function(e) {
      cat("Error occurred for tree", i, tree_name, ": ", conditionMessage(e), "\n")
    }
  )
}

