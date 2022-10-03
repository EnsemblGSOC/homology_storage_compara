import pytest
from gene_tree import GeneTree, GeneTreeException, HomologyType, IntervalIndex
import pandas as pd

brca_homology_table = pd.read_csv('test/test_data/BRCA_homology.tsv', sep='\t')

def test_gene_tree_load():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    ca = gt.get_common_ancestor(['ENSCING00000002681', 'ENSCING00000016195'])
    assert len(ca) > 0
    assert ca[0].phyloxml_clade.events.duplications != 0


def test_get_lca_type():
    gt = GeneTree()
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
    lca_type, _ = gt._get_lca_type_score(('ENSCING00000002681', 'ENSCSAVG00000006167'))
    assert lca_type == 'duplication'


def test_event_node_annotation():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/HLA_F_gene_tree.xml')
    gt.load_ref_table('test/test_data/HLA_F_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    print('total dup:', len([x for x in gt.trees[0].get_descendants() 
                            if x.phyloxml_clade.events is not None and
                            x.phyloxml_clade.events.duplications == 1]))
    print('total sep:', len([x for x in gt.trees[0].get_descendants() 
                            if x.phyloxml_clade.events is not None and
                            x.phyloxml_clade.events.speciations == 1]))
    print(gt.trees[0])
    assert len([x for x in gt.trees[0].get_descendants()
                if x.phyloxml_clade.events is not None and
                x.phyloxml_clade.events.duplications != 0]) > 0


def test_gene_tree_exception():
    gt = GeneTree()
    with pytest.raises(GeneTreeException):
        gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
        gt.annotate_event_nodes()


def test_export_phylo_xml():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/HOXA10_gene_tree.xml')
    gt.load_ref_table('test/test_data/HOXA10_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    gt.export_phylo_xml('test/test_out/HOXA10_with_events.xml')


def test_export_large_phylo_xml():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/BRCA_gene_tree.xml')
    gt.load_ref_table('test/test_data/BRCA_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    gt.export_phylo_xml('test/test_out/BRCA_gene_tree_with_events.xml')


def test_fast_export_large_phylo_xml():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/HLA_F_gene_tree.xml')
    gt.load_ref_table('test/test_data/HLA_F_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes_fast()
    gt.export_phylo_xml('test/test_out/fasttest_HLA_F_gene_tree_with_events.xml')


def test_homology_inference():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    type1 = gt.get_homology_type('ENSCSAVG00000006167', 'ENSCING00000012110')
    type2 = gt.get_homology_type('ENSCSAVG00000006167', 'ENSCSAVG00000009942')
    assert type1 == HomologyType.ortholog_one2one
    assert type2 == HomologyType.within_species_paralog
    # TODO: test many-to-many ortholog


def test_complex_homology_inference():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_out/BRCA_gene_tree_with_events.xml')
    # gt.load_ref_table('test/test_data/BRCA_gene_tree_ref_table.tsv')
    # gt.annotate_event_nodes()
    assert gt.get_homology_type('ENSSGRG00000017509', 'ENSSGRG00000036407') == \
        _get_actual_homology_type('ENSSGRG00000017509', 'ENSSGRG00000036407')


def test_get_all_orthologs():
    gt = GeneTree()
    gt.load_phylo_xml("pipeline/data_out_without_seq/SNTA1_gene_tree.xml")
    orthologs = gt.get_all_orthologs("MGP_SPRETEiJ_G0018509")
    assert len(orthologs) == 141


def test_get_all_paralogs():
    gt = GeneTree()
    gt.load_phylo_xml("test/test_out/RNU6_92P_with_events.xml")
    gt.trees[0].get_common_ancestor(gt.trees[0].get_leaves_by_name("ENSG00000272393")[0], gt.trees[0].get_leaves_by_name("ENSPEMG00000027410")[0])
    paralogs = gt.get_all_paralogs("ENSG00000272393")
    assert len(paralogs) == 354


def test_node_labeling():
    gt = GeneTree()
    gt.load_phylo_xml("test/test_out/ENSCSAVG00000010931_with_events.xml")
    idx = IntervalIndex(gt)
    print(idx.internal_labels.values())


def test_get_all_orthologs_indexed():
    gt = GeneTree()
    gt.load_phylo_xml("pipeline/data_out_without_seq/SNTA1_gene_tree.xml")
    gt.create_interval_index()
    orthologs = gt.get_all_orthologs_indexed("MGP_SPRETEiJ_G0018509")
    print("label of query gene: {}".format(gt.interval_index.leaf_labels[gt.trees[0].get_leaves_by_name("MGP_SPRETEiJ_G0018509")[0]]))
    print(set(orthologs))
    actual_orthologs = gt.get_all_orthologs("MGP_SPRETEiJ_G0018509")
    actual_ortholog_labels = [gt.interval_index.leaf_labels[x] for x in actual_orthologs]
    print(actual_ortholog_labels)
    assert len(set(orthologs)) == len(gt.get_all_orthologs("MGP_SPRETEiJ_G0018509"))


def test_get_all_paralogs_indexed():
    gt = GeneTree()
    # gt.load_phylo_xml("test/test_out/RNU6_92P_with_events.xml")
    gt.load_phylo_xml("pipeline/data_out/TSHZ2_gene_tree.xml")
    gt.create_interval_index()
    paralogs = gt.get_all_paralogs_indexed("ENSG00000182463")
    print("label of query gene: {}".format(gt.interval_index.leaf_labels[gt.trees[0].get_leaves_by_name("ENSG00000182463")[0]]))
    print(set(paralogs))
    actual_paralogs = gt.get_all_paralogs("ENSG00000182463")
    actual_paralog_labels = [gt.interval_index.leaf_labels[x] for x in actual_paralogs]
    print(actual_paralog_labels)
    assert len(set(paralogs)) == len(gt.get_all_paralogs("ENSG00000182463"))


def _get_actual_homology_type(gene1, gene2):
    homology_type = brca_homology_table[
        (brca_homology_table['stable_id_1'] == gene1) &
        (brca_homology_table['stable_id_2'] == gene2)
    ]['description'].values[0]
    if homology_type == 'ortholog_one2one':
        return HomologyType.ortholog_one2one
    elif homology_type == 'ortholog_one2many':
        return HomologyType.ortholog_one2many
    elif homology_type == 'ortholog_many2many':
        return HomologyType.ortholog_many2many
    elif homology_type == 'within_species_paralog':
        return HomologyType.within_species_paralog
    elif homology_type == 'between_species_paralog':
        return HomologyType.between_species_paralog
    elif homology_type == 'gene_split':
        return HomologyType.gene_split
    else:
        return HomologyType.none
