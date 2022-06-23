import pytest
from gene_tree import GeneTree, GeneTreeException, HomologyType


def test_gene_tree_load():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    ca = gt.get_common_ancestor(['ENSCING00000002681', 'ENSCING00000016195'])
    assert len(ca) > 0
    assert ca[0].phyloxml_clade.events.duplications != 0


def test_get_lca_type():
    gt = GeneTree()
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
    lca_type = gt._get_lca_type(('ENSCING00000002681', 'ENSCSAVG00000006167'))
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


def test_homology_inference():
    gt = GeneTree()
    gt.load_phylo_xml('test/test_data/ENSCSAVG00000010931_gene_tree.xml')
    gt.load_ref_table('test/test_data/ENSCSAVG00000010931_gene_tree_ref_table.tsv')
    gt.annotate_event_nodes()
    type1 = gt.get_homology_type('ENSCSAVG00000006167', 'ENSCING00000012110')
    type2 = gt.get_homology_type('ENSCSAVG00000006167', 'ENSCSAVG00000009942')
    assert type1 == HomologyType.ortholog_one2one
    assert type2 == HomologyType.between_species_paralog
    # TODO: test many-to-many ortholog
