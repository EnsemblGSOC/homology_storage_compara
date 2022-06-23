import xml
from ete3 import Tree, Phyloxml, PhyloxmlTree
import ete3
from typing import List, Dict, Tuple, Union, Any
import pandas as pd
from enum import Enum


class HomologyType(Enum):
    """
    Enum for homology type.
    """
    ortholog_one2one = 1
    ortholog_one2many = 2
    ortholog_many2many = 3
    between_species_paralog = 4
    within_species_paralog = 5
    gene_split = 6
    other = 7
    not_homologous = 8


class GeneTree:
    trees: List[PhyloxmlTree] = []
    ref_table: pd.DataFrame = None

    def load_phylo_xml(self, xml_file: str):
        """
        Load the gene tree from a PhyloXML file.
        """
        trees = Phyloxml()
        trees.build_from_file(xml_file)
        for t in trees.get_phylogeny():
            self.trees.append(t)

    def get_common_ancestor(self, taxon_names: list):
        """
        Get the common ancestor of the given taxa.
        """
        if len(self.trees) > 0 and len(taxon_names) > 1:
            return [x.get_common_ancestor(taxon_names) for x in self.trees]
        return []

    def get_homology_type(self, taxon1: str, taxon2: str):
        """
        Get the homology type between two taxa.
        """
        cas = self.get_common_ancestor([taxon1, taxon2])
        if cas is not None and len(cas) > 0:
            ca = cas[0]
            if ca.phyloxml_clade.events is not None:
                if ca.phyloxml_clade.events.speciations is not None:
                    # we have a pair of orthologs
                    children = ca.get_children()
                    leaves = [x for x in children if x.is_leaf()]
                    if len(leaves) == 1:
                        return HomologyType.ortholog_one2many
                    elif len(leaves) == 2:
                        return HomologyType.ortholog_one2one
                    else:
                        return HomologyType.ortholog_many2many
                elif ca.phyloxml_clade.events.duplications is not None:
                    # we have a pair of paralogs
                    descendants = ca.get_descendants()
                    speciation_descendants = [x for x in descendants
                                              if x.phyloxml_clade.events is not None
                                              and x.phyloxml_clade.events.speciations > 0]
                    # TODO: this is not quite correct, we might have to check
                    # the species of the leaves
                    if len(speciation_descendants) == 1:
                        return HomologyType.between_species_paralog
                    elif len(speciation_descendants) == 0:
                        return HomologyType.within_species_paralog
                return HomologyType.other
        else:
            return HomologyType.not_homologous

    def load_ref_table(self, ref_table_file: str, delimiter: str = '\t'):
        """
        Load the reference table for event node annotation.
        """
        self.ref_table = pd.read_csv(ref_table_file, delimiter=delimiter)

    def annotate_event_nodes(self):
        """
        Annotate the event nodes in the gene tree.
        """
        if len(self.trees) > 0:
            if self.ref_table is not None:
                for t in self.trees:
                    leaves = t.get_leaves()
                    # all pairs of leaves
                    pairs = [{x, y} for x in leaves for y in leaves if x != y]
                    for l1, l2 in pairs:
                        lca = t.get_common_ancestor([l1.name, l2.name])
                        lca_type, confidence = self._get_lca_type_score((l1.name, l2.name))
                        if lca is not None and lca_type != 'unknown':
                            if lca_type == 'duplication':
                                if lca.phyloxml_clade.events is None:
                                    lca.phyloxml_clade.events = ete3.phyloxml.Events(duplications=1)
                                else:
                                    lca.phyloxml_clade.events.duplications = 1
                                # set the confidence score for duplication event
                                lca.phyloxml_clade.events.set_confidence(
                                    ete3.phyloxml.Confidence(type_='duplication', valueOf_=confidence))
                            elif lca_type == 'speciation':
                                if lca.phyloxml_clade.events is None:
                                    lca.phyloxml_clade.events = ete3.phyloxml.Events(speciations=1)
                                else:
                                    lca.phyloxml_clade.events.speciations = 1
            else:
                raise GeneTreeException('Reference table is not loaded.')

    def _get_lca_type_score(self, taxon_names: Tuple):
        if self.ref_table is not None:
            id1 = taxon_names[0]
            id2 = taxon_names[1]
            n1 = self.ref_table[self.ref_table['stable_id'] == id1]
            n2 = self.ref_table[self.ref_table['stable_id'] == id2]
            root_id = n1['root_id'].values[0]
            if len(n1) > 0 and len(n2) > 0:
                while id1 != id2 and id1 != root_id and id2 != root_id:
                    n1 = self.ref_table[self.ref_table['node_id'] == n1['parent_id'].values[0]]
                    n2 = self.ref_table[self.ref_table['node_id'] == n2['parent_id'].values[0]]
                    id1 = n1['node_id'].values[0]
                    id2 = n2['node_id'].values[0]
                    # if node_id is equal to root_id, we are at the root
                    if (id1 == n1['root_id'].values[0]):
                        id1 = n1['root_id'].values[0]
                    if (id2 == n2['root_id'].values[0]):
                        id2 = n2['root_id'].values[0]
                return n1['node_type'].values[0], n1['duplication_confidence_score'].values[0]
        else:
            raise GeneTreeException('Reference table is not loaded.')

    def export_phylo_xml(self, xml_file: str):
        """
        Print the gene tree in PhyloXML format.
        """
        if len(self.trees) > 0:
            for t in self.trees:
                phylo_xml = Phyloxml()
                phylo_xml.add_phylogeny(t)
                with open(xml_file, 'w') as f:
                    phylo_xml.export(f)


class GeneTreeException(Exception):
    pass
