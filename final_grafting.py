import ete3
from Bio import Phylo

BASE_TREE = 'fullpass_oliveros100'  # NOTE: you need to name the backbone tree!

# the base 72sp reference tree
# t72sp = Phylo.read(BASE_TREE, 'nexus')
# taxa72sp = set([x.name for x in t72sp.get_terminals()])

# outgroup taxa
outgroup = {'Melopsittacus_undulatus', 'Casuarius_casuarius','Dromaius_novaehollandiae', 'Apteryx_australis','Chelonia_mydas', 'Testudo_graeca', 'Alligator_sinensis','Crocodylus_porosus', 'Varanus_komodoensis', 'Caretta_caretta'}

# replacement rules:
# - each subtree can have one or more rules specifying which clades to copy
# - 'source' specifies taxa whose common ancestor will be copied from subtree
# - 'target' specifies taxa whose common ancestor will be replace in 72sp tree with 'source'
# - 'placeholder' if True (default is False) will treat 'source' as a placeholder tip to insert as sibling of 'target'

direct_on_72sp = {
    'ACANTHIZIDAE': [
        {
            'source': ['Acanthiza_cinerea', 'Pachycare_flavogriseum'],
            'target': ['Acanthiza_cinerea']
        }
    ],
    'PARDALOTIDAE': [
        {
            'source': ['Pardalotus_striatus', 'Pardalotus_rubricatus'],
            'target': ['Pardalotus_striatus']
        }
    ],
    'DASYORNITHIDAE': [
        {
            'source': ['Dasyornis_broadbenti', 'Dasyornis_brachypterus'],
            'target': ['Dasyornis_broadbenti']
        }
    ],
    'MALURIDAE': [
        {
            'source': ['Malurus_alboscapulatus', 'Amytornis_barbatus'],
            'target': ['Malurus_alboscapulatus']
        }
    ],
    'ORTHONYCHIDAE': [
        {
            'source': ['Orthonyx_temminckii', 'Orthonyx_novaeguineae'],
            'target': ['Orthonyx_temminckii']
        }
    ],
    'POMATOSTOMIDAE': [
        {
            'source': ['Pomatostomus_superciliosus', 'Pomatostomus_ruficeps'],
            'target': ['Pomatostomus_superciliosus']
        }
    ],
    'PTILONORHYNCHIDAE': [
        {
            'source': ['Ailuroedus_buccoides', 'Ptilonorhynchus_violaceus'],
            'target': ['Ailuroedus_buccoides', 'Ptilonorhynchus_violaceus']
        }
    ],
    'CLIMACTERIDAE': [
        {
            'source': ['Climacteris_picumnus','Cormobates_leucophaea'],
            'target': ['Climacteris_melanurus','Cormobates_leucophaea']
        }
    ],
    'ATRICHORNITHIDAE': [
        {
            'source': ['Atrichornis_rufescens', 'Atrichornis_clamosus'],
            'target': ['Atrichornis_rufescens']
        }
    ],
    'MENURIDAE': [
        {
            'source': ['Menura_novaehollandiae', 'Menura_alberti'],
            'target': ['Menura_novaehollandiae']
        }
    ],
    'PITTIDAE': [
        {
            'source': ['Pitta_erythrogaster', 'Erythropitta_granatina','Pitta_superba'],
            'target': ['Pitta_erythrogaster','Pitta_cyanea']
        }
    ],
    'CALYPTOMENIDAE': [
        {
            'source': ['Calyptomena_viridis', 'Calyptomena_hosii'],
            'target': ['Calyptomena_viridis']
        }
    ],
    'EURYLAIMIDAE': [
        {
            'source': ['Eurylaimus_ochromalus', 'Corydon_sumatranus'],
            'target': ['Eurylaimus_ochromalus']
        }
    ],
    'CINCLOSOMATIDAE': [
        {
            'source': ['Cinclosoma_punctatum','Ptilorrhoa_leucosticta','Cinclosoma_marginatum'],
            'target': ['Cinclosoma_punctatum','Ptilorrhoa_leucosticta']
        }
    ],
    'CAMPEPHAGIDAE': [
        {
            'source': ['Edolisoma_tenuirostre','Pericrocotus_divaricatus', 'Pericrocotus_montanus'],
            'target': ['Edolisoma_tenuirostre','Pericrocotus_divaricatus']
        }
    ],
    'NEOSITTIDAE': [
        {
            'source': ['Daphoenositta_chrysoptera', 'Daphoenositta_miranda'],
            'target': ['Daphoenositta_chrysoptera']
        }
    ],
    'PSOPHODIDAE': [
        {
            'source': ['Psophodes_cristatus', 'Androphobus_viridis'],
            'target': ['Psophodes_cristatus']
        }
    ],
    'OREOICIDAE': [
        {
            'source': ['Oreoica_gutturalis', 'Ornorectes_cristatus'],
            'target': ['Oreoica_gutturalis']
        }
    ],
    'ORIOLIDAE': [
        {
            'source': ['Oriolus_chinensis','Pitohui_dichrous','Sphecotheres_hypoleucus'],
            'target': ['Oriolus_chinensis','Pitohui_dichrous']
        }
    ],
    'PARAMYTHIIDAE': [
        {
            'source': ['Oreocharis_arfaki', 'Paramythia_olivacea'],
            'target': ['Oreocharis_arfaki']
        }
    ],
    'VIREONIDAE': [
        {
            'source': ['Erpornis_zantholeuca', 'Pteruthius_aeralatus'],
            'target': ['Erpornis_zantholeuca', 'Pteruthius_aeralatus']
        }
    ],
    'MACHAERIRHYNCHIDAE': [
        {
            'source': ['Machaerirhynchus_nigripectus','Machaerirhynchus_flaviventer'],
            'target': ['Machaerirhynchus_nigripectus']
        }
    ],
    'ARTAMIDAE': [
        {
            'source': ['Artamus_cinereus', 'Melloria_quoyi', 'Peltops_blainvillii', 'Strepera_graculina'],
            'target': ['Artamus_cinereus', 'Melloria_quoyi', 'Peltops_blainvillii', 'Strepera_graculina']
        }
    ],
    'AEGITHINIDAE': [
        {
            'source': ['Aegithina_lafresnayei', 'Aegithina_viridissima'],
            'target': ['Aegithina_lafresnayei']
        }
    ],
    'VANGIDAE': [
        {
            'source': ['Tephrodornis_virgatus', 'Philentoma_velata'],
            'target': ['Tephrodornis_virgatus']
        }
    ],
    'RHIPIDURIDAE': [
        {
            'source': ['Rhipidura_javanica','Chaetorhynchus_papuensis','Rhipidura_tenebrosa'],
            'target': ['Rhipidura_javanica']
        }
    ],
    'DICRURIDAE': [
        {
            'source': ['Dicrurus_aeneus', 'Dicrurus_sumatranus'],
            'target': ['Dicrurus_aeneus']
        }
    ],
    'MONARCHIDAE': [
        {
            'source': ['Hypothymis_azurea', 'Monarcha_megarhynchus', 'Arses_insularis', 'Myiagra_eichhorni'],
            'target': ['Hypothymis_azurea']
        }
    ],
    'PARADISAEIDAE': [
        {
            'source': ['Paradisaea_minor', 'Phonygammus_keraudrenii',],
            'target': ['Paradisaea_minor', 'Phonygammus_keraudrenii']
        }
    ],
    'CORCORACIDAE': [
        {
            'source': ['Struthidea_cinerea', 'Corcorax_melanorhamphos'],
            'target': ['Struthidea_cinerea']
        }
    ],
    'MELAMPITTIDAE': [
        {
            'source': ['Melampitta_lugubris', 'Megalampitta_gigantea'],
            'target': ['Melampitta_lugubris']
        }
    ],
    'CNEMOPHILIDAE': [
        {
            'source': ['Cnemophilus_loriae', 'Loboparadisea_sericea'],
            'target': ['Cnemophilus_loriae', 'Loboparadisea_sericea']
        }
    ],
    'MELANOCHARITIDAE': [
        {
            'source': ['Melanocharis_versteri', 'Oedistoma_iliolophus', 'Toxorhampus_novaeguineae'],
            'target': ['Melanocharis_versteri', 'Oedistoma_iliolophus', 'Toxorhamphus_novaeguineae']
        }
    ],
    'PETROICIDAE': [
        {
            'source': ['Devioeca_papuana', 'Petroica_multicolor', 'Amalocichla_incerta'],
            'target': ['Devioeca_papuana', 'Petroica_multicolor']
        }
    ],
    'STENOSTIRIDAE': [
        {
            'source': ['Chelidorhynx_hypoxanthus', 'Culicicapa_ceylonensis'],
            'target': ['Chelidorhynx_hypoxanthus', 'Culicicapa_ceylonensis']
        }
    ],
    'ORTHONYCHIDAE': [
        {
            'source': ['Orthonyx_temminckii', 'Orthonyx_novaeguineae'],
            'target': ['Orthonyx_temminckii']
        }
    ],
    'ACROCEPHALIDAE': [
        {
            'source': ['Acrocephalus_orientalis', 'Arundinax_aedon'],
            'target': ['Acrocephalus_orientalis']
        }
    ],
    'LOCUSTELLIDAE': [
        {
            'source': ['Locustella_lanceolata', 'Megalurulus_grosvenori', 'Robsonius_rabori', 'Poodytes_gramineus'],
            'target': ['Locustella_lanceolata', 'Megalurus_palustris', 'Robsonius_rabori']
        }
    ],
    'HIRUNDINIDAE': [
        {
            'source': ['Hirundo_rustica', 'Delichon_lagopodum', 'Cecropis_daurica'],
            'target': ['Hirundo_rustica']
        }
    ],
    'PYCNONOTIDAE': [
        {
            'source': ['Pycnonotus_jocosus', 'Alophoixus_pallidus'],
            'target': ['Pycnonotus_jocosus']
        }
    ],
    'CHLOROPSEIDAE': [
        {
            'source': ['Chloropsis_cochinchinensis', 'Chloropsis_sonnerati'],
            'target': ['Chloropsis_cochinchinensis', 'Chloropsis_sonnerati']
        }
    ],
    'ZOSTEROPIDAE': [
        {
            'source': ['Zosterops_everetti', 'Yuhina_castaniceps'],
            'target': ['Zosterops_everetti']
        }
    ],
    'TIMALIIDAE': [
        {
            'source': ['Timalia_pileata', 'Stachyris_grammiceps', ],
            'target': ['Timalia_pileata']
        }
    ],
    'PELLORNEIDAE': [
        {
            'source': ['Trichastoma_pyrrogenys', 'Leonardina_woodi'],
            'target': ['Trichastoma_pyrrogenys']
        }
    ],
    'SCOTOCERCIDAE': [
        {
            'source': ['Abroscopus_albogularis','Tesia_cyaniventer', 'Tesia_olivea'],
            'target': ['Abroscopus_albogularis','Tesia_cyaniventer']
        }
    ],
    'DICAEIDAE': [
        {
            'source': ['Dicaeum_hypoleucum', 'Dicaeum_bicolor'],
            'target': ['Dicaeum_hypoleucum']
        }
    ],
    'NECTARINIIDAE': [
        {
            'source': ['Leptocoma_sperata', 'Cinnyris_idenburgi'],
            'target': ['Leptocoma_sperata']
        }
    ],
    'IRENIDAE': [
        {
            'source': ['Irena_cyanogastra', 'Irena_puella'],
            'target': ['Irena_cyanogastra', 'Irena_puella']
        }
    ],
    'CHLOROPSEIDAE': [
        {
            'source': ['Chloropsis_cochinchinensis', 'Chloropsis_sonnerati'],
            'target': ['Chloropsis_cochinchinensis', 'Chloropsis_sonnerati']
        }
    ],
    'PASSERIDAE': [
        {
            'source': ['Passer_domesticus', 'Hypocryptadius_cinnamomeus'],
            'target': ['Passer_domesticus']
        }
    ],
    'MOTACILLIDAE': [
        {
            'source': ['Motacilla_alba', 'Anthus_godlewskii'],
            'target': ['Motacilla_alba']
        }
    ],
    'FRINGILLIDAE': [
        {
            'source': ['Fringilla_montifringilla', 'Pyrrhula_leucogenis'],
            'target': ['Fringilla_montifringilla']
        }
    ],
    'MELIPHAGIDAE': [
        {
            'source': ['Meliphaga_lewinii','Acanthorhynchus_tenuirostris', 'Myzomela_cruentata'],
            'target': ['Meliphaga_montana','Acanthorhynchus_tenuirostris'],
        }
    ],
    'PACHYCEPHALIDAE': [
        {
            'source': ['Colluricincla_harmonica','Pachycephala_simplex'],
            'target': ['Colluricincla_harmonica','Pachycephala_vitiensis']
        }
    ],
    'CORVIDAE': [
        {
            'source': ['Corvus_coronoides', 'Cissa_chinensis', 'Platysmurus_leucopterus'],
            'target': ['Corvus_coronoides']
        }
    ],
    'LANIIDAE': [
        {
            'source': ['Lanius_schach', 'Lanius_cristatus'],
            'target': ['Lanius_schach']
        }
    ],
    'PARIDAE': [
        {
            'source': ['Parus_cinereus', 'Pardaliparus_elegans'],
            'target': ['Parus_cinereus']
        }
    ],
    'ALAUDIDAE': [
        {
            'source': ['Alauda_gulgula', 'Mirafra_javanica'],
            'target': ['Alauda_gulgula']
        }
    ],
    'LEIOTHRICHIDAE': [
        {
            'source': ['Leiothrix_laurinae', 'Garrulax_castanotis'],
            'target': ['Leiothrix_laurinae']
        }
    ],
    'PHYLLOSCOPIDAE': [
        {
            'source': ['Phylloscopus_trivirgatus', 'Phylloscopus_trochiloides', 'Phylloscopus_fuscatus'],
            'target': ['Phylloscopus_trivirgatus']
        }
    ],
    'TURDIDAE': [
        {
            'source': ['Turdus_chrysolaus', 'Cochoa_viridis', 'Zoothera_aurea'],
            'target': ['Turdus_chrysolaus']
        }
    ],
    'MUSCICAPIDAE': [
        {
            'source': ['Muscicapa_sibirica', 'Larvivora_ruficeps', 'Saxicola_maurus'],
            'target': ['Muscicapa_sibirica']
        }
    ],
    'STURNIDAE': [
        {
            'source': ['Sturnia_malabarica','Rhabdornis_inornatus'],
            'target': ['Sturnia_malabarica']
        }
    ],
    'SITTIDAE': [
        {
            'source': ['Sitta_frontalis', 'Sitta_neglecta'],
            'target': ['Sitta_frontalis']
        }
    ],
    'PLOCEIDAE': [
        {
            'source': ['Ploceus_philippinus', 'Ploceus_manyar'],
            'target': ['Ploceus_philippinus']
        }
    ],
    'ESTRILDIDAE': [
        {
            'source': ['Poephila_personata', 'Lonchura_leucogastra','Erythrura_trichroa'],
            'target': ['Poephila_personata']
        }
    ],
    'EMBERIZIDAE': [
        {
            'source': ['Emberiza_fucata', 'Emberiza_sulphurata'],
            'target': ['Emberiza_fucata']
        }
    ],
    'CISTICOLIDAE': [
        {
            'source': ['Orthotomus_castaneiceps', 'Orthotomus_sepium'],
            'target': ['Orthotomus_castaneiceps']
        }
    ],
    'SYLVIIDAE' : [
        {
            'source': ['Cholornis_unicolor', 'Paradoxornis_guttaticollis'],
            'target': ['Cholornis_unicolor']
        }
    ]
}


# # add annotation indicating source of 95% interval (and, therefore, node age)
# for n in t72sp.find_elements():
#     if hasattr(n, 'comment') and n.comment is not None:
#         n.comment = n.comment[:-1] + ',from="72sp"]'

# uncomment if you want separate trees for each step - useful for debugging
# t72sp = Phylo.read(BASE_TREE, 'nexus')
# for n in t72sp.find_elements():
#     if hasattr(n, 'comment') and n.comment is not None:
#         n.comment = n.comment[:-1] + ',from="72sp"]'


def do_one_rule_for_one_tree(backbone_tree: ete3.Tree, subtree: ete3.Tree, graft_rule):
    if graft_rule.get('placeholder', False):
        from Bio.Phylo.BaseTree import Clade
        ca_backbone = backbone_tree.common_ancestor(graft_rule['target'])
        print("ca_backbone is", ca_backbone)
        parent = backbone_tree.get_path(ca_backbone)[-2]
        print("parent is", parent)
        sibling_index = 0 if parent.clades[0] == ca_backbone else 1
        print("sibling index is:", sibling_index)
        ca_subtree = Clade(branch_length=1, name=graft_rule['source'][0])
        ca_subtree.comment = '[&from="placeholder"]'
        new_branch_length = parent.clades[sibling_index].branch_length / 2
        parent.clades[sibling_index].branch_length = parent.clades[sibling_index].branch_length + new_branch_length
        parent.clades[sibling_index].branch_length = new_branch_length
        new_clade = Clade(branch_length=new_branch_length, clades=[parent.clades[sibling_index], ca_subtree])
        new_clade.comment = '[&from="placeholder"]'
        parent.clades[sibling_index] = new_clade
        distance_backbone = backbone_tree.distance(new_clade, graft_rule['target'][0])
        ca_subtree.branch_length = distance_backbone
    else:
        ca_subtree = subtree.common_ancestor(graft_rule['source'])
        ca_backbone = backbone_tree.common_ancestor(graft_rule['target'])
        # check all taxa in 72sp are present in subtree
        taxa_ca_backbone = set([x.name for x in ca_backbone.get_terminals()])
        taxa_ca_subtree = set([x.name for x in ca_subtree.get_terminals()])
        common_taxa = taxa_ca_backbone.intersection(taxa_ca_subtree)
        # if not all the taxa from the 72sp tree are present in the common taxa
        if common_taxa != taxa_ca_backbone:
            print(f'Backbone tree has {taxa_ca_backbone.difference(common_taxa)} but at not in subtree')
        # get the distance from tip to common ancestor in 72sp and subtree
        distance_backbone = backbone_tree.distance(ca_backbone, next(iter(common_taxa)))
        distance_subtree = subtree.distance(ca_subtree, next(iter(common_taxa)))
        print("ca_subtree branch length is", ca_subtree.branch_length)

        # new stem length for 72sp tree
        new_stem_length = ca_backbone.branch_length + (distance_backbone - distance_subtree)
        print({new_stem_length}, "is new stem length for now")

        if new_stem_length < 0 and distance_backbone == 0:
           scaling_factor = ca_backbone.branch_length/(distance_subtree + ca_subtree.branch_length)
           for node in ca_subtree.find_clades(order="postorder"):
               node.branch_length *= scaling_factor
           distance_subtree = subtree.distance(ca_subtree, next(iter(common_taxa)))
           new_stem_length = ca_backbone.branch_length + (distance_backbone - distance_subtree)

#        max_backbone_length = max([clade.branch_length for clade in ca_backbone.find_clades(order="preorder")])

 #       total_subtree_length = sum([clade.branch_length for clade in ca_subtree.find_clades(order="preorder")])

  #      scaling_factor = max_backbone_length / total_subtree_length

   #     print(f"{scaling_factor} is obtained from {max_backbone_length}/{total_subtree_length}")

 #       for node in ca_subtree.find_clades(order="postorder"):
  #          node.branch_length *= scaling_factor
#            if node.branch_length > max_backbone_length:
 #               node.branch_length = max_backbone_length

 #       new_stem_length = ca_backbone.branch_length

        print(f'MRCA age: {distance_backbone} -> {distance_subtree}')
        print(f'Adjusted branch length: {ca_backbone.branch_length} -> {new_stem_length}')
        print(f'Replacing {len(ca_backbone.get_terminals())} taxa with {len(ca_subtree.get_terminals())}')

        ca_backbone.branch_length = new_stem_length
        ca_backbone.comment = ca_subtree.comment
        ca_backbone.clades = ca_subtree.clades
        ca_backbone.name = None

    return backbone_tree

    # prefix names
    # for c in ca_72sp.get_terminals():
    #     c.name = f'{index}_{rule_index}_{c.name}'

    # check distances
    # for rule in rules_so_far:
    #     ca_72sp = t72sp.common_ancestor(rule)
    #     print('Node of MRCA', rule, ':', t72sp.distance(ca_72sp, rule[0]))
    #
    # print(f'Write {index}_{rule_index}')
    # Phylo.write(t72sp, f'5000sp_{index}_{rule_index}.tree', 'nexus')


def do_all_rules_on_tree(backbone_tree_file, prune_graft_rules, gather=False):
    if gather:
        filename = backbone_tree_file
        backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')
        for n in backbone_tree_file.find_elements():
            if hasattr(n, 'comment') and n.comment is not None and "from=" not in n.comment:
                n.comment = n.comment[:-1] + f',from="{filename}"]'

    for tree_file, graft_rules in prune_graft_rules.items():
        backbone_tree_file = do_all_rules_for_one_tree(backbone_tree_file, tree_file, graft_rules, gather)
        print()

    return backbone_tree_file


def do_all_rules_for_one_tree(backbone_tree_file, tree_file, graft_rules, gather):
    subtree = Phylo.read(f'{tree_file}.tree', 'nexus')
    for n in subtree.find_elements():
        if hasattr(n, 'comment') and n.comment is not None and "from=" not in n.comment:
            n.comment = n.comment[:-1] + f',from="{tree_file}"]'

    # fix tilda names in subtree
    for terminal in subtree.get_terminals():
        if terminal.name.endswith('~'):
            terminal.name = terminal.name[:-1]

    total_subtree_taxa = len(subtree.get_terminals()) - len(outgroup)
    print(f'{"-" * 10} Subtree {tree_file}, {total_subtree_taxa} taxa {"-" * 10}')

    if isinstance(backbone_tree_file, str):
        backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')

    rule_index = 0
    for graft_rule in graft_rules:
        if isinstance(backbone_tree_file, str):
            backbone_tree_file = Phylo.read(f'{backbone_tree_file}.tree', 'nexus')
        rule_index += 1
        print(f'- Applying rule {rule_index}')
        grafted_tree = do_one_rule_for_one_tree(backbone_tree_file, subtree, graft_rule)
        out_tree = f'Working_{tree_file}_{rule_index}.tree'
        Phylo.write(grafted_tree, out_tree, 'nexus')
        print()

    return backbone_tree_file

all_mammals = do_all_rules_on_tree('fullpass_oliveros100', direct_on_72sp, gather=True)
Phylo.write(all_mammals, 'finalbirdsp.tree', 'nexus')
#---------------------------------------------------#
# SAC-211110 | Read file and, in case "\[" or "\]"
# are found, then remove them.
# E.g., [\[&95%={0.0694354, 0.280382},from="Marsupialia-time"\]]
# This seems to happen when using the WLS
with open('finalbirdsp.tree', 'r') as file :
  inp_file = file.read()

inp_file = inp_file.replace('\[', '')
inp_file = inp_file.replace('\]', '')

with open('finalbirdsp.tree', 'w') as file:
  file.write(inp_file)
#----------------------------------------------------#

# taxa_in_1 = {t.name for t in Phylo.read('Rodentia_squirrel.tree', 'nexus').get_terminals()}
# taxa_in_2 = {t.name for t in Phylo.read('Rodentia_ctenohystrica_3.tree', 'nexus').get_terminals()}
# taxa_in_3 = {t.name for t in Phylo.read('Rodentia_subtree1.tree', 'nexus').get_terminals()}
# taxa_in_4 = {t.name for t in Phylo.read('Rodentia_subtree2.tree', 'nexus').get_terminals()}

# index = 0
# for tree_file, graft_rules in prune_graft_rules.items():
#     index += 1
#     subtree = Phylo.read(tree_file, 'nexus')
#     for n in subtree.find_elements():
#         if hasattr(n, 'comment') and n.comment is not None:
#             n.comment = n.comment[:-1] + f',from="subtree"]'
#
#     print(f'{index}. {tree_file}')
#
#     # fix tilda names in subtree
#     for t in subtree.get_terminals():
#         if t.name.endswith('~'):
#             t.name = t.name[:-1]
#
#     total_subtree_taxa = len(subtree.get_terminals()) - len(outgroup)
#
#     grafted_taxa_count = 0
#     rule_index = 0
#     rules_so_far = list()
#     for graft_rule in graft_rules:
#         do_graft()
#
#     print(f'{grafted_taxa_count}/{total_subtree_taxa} grafted')
#     print('-' * 80)

# afrotheria - monophyletic replace
# lagomorpha - monophyletic replace
# chiroptera - monophyletic replace
# marsupialia - monophyletic replace
# xenarthra - monophyletic replace
# euarchonta - separate (i) scandentia (ii) primates+dermoptera
# squirrel - single taxon replaced
