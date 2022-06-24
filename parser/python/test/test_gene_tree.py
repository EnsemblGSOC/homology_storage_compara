import pytest
from gene_tree import GeneTree, GeneTreeException, HomologyType
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
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
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
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    gt.export_phylo_xml('test/test_out/ENSCSAVG00000010931_with_events.xml')


def test_export_large_phylo_xml():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/BRCA_gene_tree.xml')
    gt.load_ref_table('test/test_data/BRCA_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    gt.export_phylo_xml('test/test_out/BRCA_gene_tree_with_events.xml')


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
