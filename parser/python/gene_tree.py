from turtle import update
import xml
from ete3 import Tree, Phyloxml, PhyloxmlTree
import ete3
from typing import List, Dict, Tuple, Union, Any
import pandas as pd
from enum import Enum

from interval_tree import IntervalTree, construct_interval_tree, exclude_intervals


class HomologyType(Enum):
    """
    Enum for homology type.
    """
    ortholog_one2one = 1
    ortholog_one2many = 2
    ortholog_many2many = 3
    between_species_paralog = 4
    within_species_paralog = 5
    other_paralog = 6
    gene_split = 7
    other = 8
    none = 9


class GeneTree:
    trees: List[PhyloxmlTree]
    ref_table: Union[pd.DataFrame, None]

    def __init__(self):
        self.trees = []
        self.ref_table = None

    def load_phylo_xml(self, xml_file: str):
        """
        Load the gene tree from a PhyloXML file.
        """
        tree = Phyloxml()
        tree.build_from_file(xml_file)
        for t in tree.get_phylogeny():
            self.trees.append(t)

    def get_common_ancestor(self, taxon_names: list):
        """
        Get the common ancestor of the given taxa.
        """
        if len(self.trees) > 0 and len(taxon_names) > 1:
            return [x.get_common_ancestor(taxon_names) for x in self.trees]
        return []

    def get_homology_type(self, taxon1: str, taxon2: str,
                          ignore_between_species_paralog: bool = True) -> HomologyType:
        """
        Get the homology type between two taxa. If ignore_between_species_paralog
        is set to True, paralogs will only be reported if they are within the same
        species.
        """
        # TODO: potential improvement: skip common ancestor query
        # if the two nodes are of different species if
        # ignore_between_species_paralog is True and we can
        # determine that the LCA is not a duplication node
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
                    if not ignore_between_species_paralog:
                        species1 = self.trees[0].get_leaves_by_name(taxon1)[0].species
                        species2 = self.trees[0].get_leaves_by_name(taxon2)[0].species
                        if species1 == species2:
                            return HomologyType.within_species_paralog
                        else:
                            return HomologyType.between_species_paralog
                    else:
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
                        if lca is not None:
                            # annotate only if the node is not
                            # already annotated
                            if (lca.phyloxml_clade.events is None
                                or (lca.phyloxml_clade.events.duplications is None
                                    and lca.phyloxml_clade.events.speciations is None)):
                                lca_type, confidence = self._get_lca_type_score((l1.name, l2.name))
                                if lca_type != 'unknown' and lca_type != 'dubious':
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

    def annotate_event_nodes_fast(self):
        """
        Fast annotation of event nodes. This method assumes all unannotated
        nodes are speciation nodes. Only works with XML files downloaded
        from Ensembl website.
        """
        if len(self.trees) > 0:
            for t in self.trees:
                for n in t.get_descendants():
                    if n.phyloxml_clade.events is None and not n.is_leaf():
                        n.phyloxml_clade.events = ete3.phyloxml.Events(speciations=1)
                if t.phyloxml_clade.events is None and not t.is_leaf():
                    t.phyloxml_clade.events = ete3.phyloxml.Events(speciations=1)

    def _get_lca_type_score(self, taxon_names: Tuple):
        if self.ref_table is not None:
            id1 = self.ref_table[self.ref_table['stable_id'] == taxon_names[0]]['node_id'].values[0]
            id2 = self.ref_table[self.ref_table['stable_id'] == taxon_names[1]]['node_id'].values[0]
            n1 = self.ref_table[self.ref_table['node_id'] == id1]
            n2 = self.ref_table[self.ref_table['node_id'] == id2]
            root_id = n1['root_id'].values[0]
            if len(n1) > 0 and len(n2) > 0:
                n1_path = [id1]
                while root_id != n1['node_id'].values[0]:
                    n1 = self.ref_table[self.ref_table['node_id'] == n1['parent_id'].values[0]]
                    n1_path.append(n1['node_id'].values[0])
                n2_path = [id2]
                while root_id != n2['node_id'].values[0]:
                    n2 = self.ref_table[self.ref_table['node_id'] == n2['parent_id'].values[0]]
                    n2_path.append(n2['node_id'].values[0])

                # find first intersection of the two paths
                n1_path = n1_path[::-1]
                n2_path = n2_path[::-1]
                for i in range(min(len(n1_path), len(n2_path)) - 1, -1, -1):
                    if n1_path[i] == n2_path[i]:
                        event_type = self.ref_table[self.ref_table['node_id'] == n1_path[i]]['node_type'].values[0]
                        confidence = self.ref_table[self.ref_table['node_id'] == n1_path[i]]['duplication_confidence_score'].values[0]
                        return event_type, confidence
        else:
            raise GeneTreeException('Reference table is not loaded.')

    def export_phylo_xml(self, xml_file: str):
        """
        Print the gene tree in PhyloXML format.
        """
        if len(self.trees) > 0:
            for t in self.trees:
                t.name = None
                for n in t.get_descendants():
                    if n.name == '':
                        n.name = None
                phylo_xml = Phyloxml()
                phylo_xml.add_phylogeny(t)
                with open(xml_file, 'w') as f:
                    phylo_xml.export(f)

    def get_all_orthologs(self, gene: str):
        """
        Get all orthologs of a given gene
        """
        if len(self.trees) > 0:
            orthologs = []
            genes = []
            for t in self.trees:
                leaves = t.get_leaves_by_name(gene)
                if len(leaves) > 0:
                    genes.extend(leaves)
            for o in genes:
                for t in self.trees:
                    for lf in t.get_leaves():
                        if lf.name != gene:
                            ca = t.get_common_ancestor([lf, o])
                            if ca is not None:
                                if ca.phyloxml_clade.events is not None and \
                                   ca.phyloxml_clade.events.speciations is not None and \
                                   ca.phyloxml_clade.events.speciations > 0:
                                    orthologs.append(lf)
            return orthologs
        else:
            return []

    def get_all_paralogs(self, gene: str):
        """
        Get all paralogs of a given gene
        """
        if len(self.trees) > 0:
            paralogs = []
            genes = []
            for t in self.trees:
                leaves = t.get_leaves_by_name(gene)
                if len(leaves) > 0:
                    genes.extend(leaves)
            for o in genes:
                for t in self.trees:
                    for lf in t.get_leaves():
                        if lf.name != gene:
                            ca = t.get_common_ancestor([lf, o])
                            if ca is not None:
                                if ca.phyloxml_clade.events is not None and \
                                   ca.phyloxml_clade.events.duplications is not None and \
                                   ca.phyloxml_clade.events.duplications > 0:
                                    paralogs.append(lf)
            return paralogs
        else:
            return []

    def height(self, i=0) -> int:
        """
        Get the height of the i-th gene tree
        """
        def _height_rec(t: PhyloxmlTree):
            if t is None:
                return 0
            h = 0
            for n in t.get_children():
                h = max(h, _height_rec(n))
            return h + 1
        if i < len(self.trees):
            return _height_rec(self.trees[i])
        else:
            raise IndexError()

    def create_interval_index(self):
        """
        Create an interval index for the gene trees.
        """
        self.interval_index = IntervalIndex(self)

    def get_all_orthologs_indexed(self, gene: str):
        """
        Get all orthologs of a given gene
        """
        if self.interval_index is not None:
            return self.interval_index.get_all_orthologs(gene)
        else:
            return self.get_all_orthologs(gene)

    def get_all_paralogs_indexed(self, gene: str):
        """
        Get all paralogs of a given gene
        """
        if self.interval_index is not None:
            return self.interval_index.get_all_paralogs(gene)
        else:
            return self.get_all_paralogs(gene)


class GeneTreeException(Exception):
    pass


class IntervalIndex:
    """
    Interval index for GeneTree
    """
    def __init__(self, gt: GeneTree) -> None:
        self.gene_tree = gt
        self.internal_nodes = []
        self.leaf_nodes = []
        self.leaf_labels = {}
        self.internal_labels = {}  # node -> label
        self._label_leaves()
        self._label_internal_nodes()
        self.internal_labels_rev = {v: k for k, v in self.internal_labels.items()} # label -> node
        self.leaf_labels_rev = {v: k for k, v in self.leaf_labels.items()}
        self.interval_tree = construct_interval_tree(self.internal_labels_rev.keys())
        self.speciation_nodes_labels = [x for x in self.internal_labels_rev.keys()
                                        if self.internal_labels_rev[x].phyloxml_clade.events and
                                        self.internal_labels_rev[x].phyloxml_clade.events.speciations and
                                        self.internal_labels_rev[x].phyloxml_clade.events.speciations > 0]
        self.duplication_nodes_labels = [x for x in self.internal_labels_rev.keys()
                                        if self.internal_labels_rev[x].phyloxml_clade.events and
                                        self.internal_labels_rev[x].phyloxml_clade.events.duplications and
                                        self.internal_labels_rev[x].phyloxml_clade.events.duplications > 0]

    def _label_leaves(self):
        self.leaf_nodes = self.gene_tree.trees[0].get_leaves()
        self.leaf_labels = {self.leaf_nodes[x]: x for x in range(len(self.leaf_nodes))}

    def _label_internal_nodes(self):
        self.internal_nodes = [x for x in self.gene_tree.trees[0].get_descendants() if not x.is_leaf()]
        self.internal_nodes.append(self.gene_tree.trees[0])
        self.internal_labels = {self.internal_nodes[x]: (-1, -1) for x in range(len(self.internal_nodes))}
        for i in range(len(self.leaf_nodes)):
            ancestors = self.leaf_nodes[i].get_ancestors()
            for a in ancestors:
                left, right = self.internal_labels[a]
                if left == -1:
                    self.internal_labels[a] = (i, self.internal_labels[a][1])
                if right == -1:
                    self.internal_labels[a] = (self.internal_labels[a][0], 1)
                if i < left:
                    self.internal_labels[a] = (i, self.internal_labels[a][1])
                if i > right:
                    self.internal_labels[a] = (self.internal_labels[a][0], i)

    def get_all_orthologs(self, gene: str):
        label = self.leaf_labels[self.gene_tree.trees[0].get_leaves_by_name(gene)[0]]
        # intervals = self.interval_tree.search(label)
        leaf = self.gene_tree.trees[0].get_leaves_by_name(gene)[0]
        intervals = [self.internal_labels[x] for x in leaf.get_ancestors()]
        intervals = sorted(intervals, key=lambda x: abs(x[1] - x[0]))
        visited = []
        orthologs = []
        for i in intervals:
            # internal node corresponding to the label
            internal_node = self.internal_labels_rev[i]
            # only consider those internal nodes that are speciation nodes
            if i in self.speciation_nodes_labels:
                # exclude those already visited
                for j in range(i[0], i[1]+1):
                    is_excluded = False
                    for x in visited:
                        if j in range(x[0], x[1]+1):
                            is_excluded = True
                    if j != label and not is_excluded:
                        orthologs.append(j)
            visited.append(i)
        return orthologs

    def get_all_paralogs(self, gene: str):
        label = self.leaf_labels[self.gene_tree.trees[0].get_leaves_by_name(gene)[0]]
        intervals = self.interval_tree.search(label)
        intervals = sorted(intervals, key=lambda x: abs(x[1] - x[0]))
        visited = []
        paralogs = []
        for i in intervals:
            # only consider those internal nodes that are duplication nodes
            if i in self.duplication_nodes_labels:
                # exclude those already visited
                for j in range(i[0], i[1]+1):
                    is_excluded = False
                    for x in visited:
                        if j in range(x[0], x[1]+1):
                            is_excluded = True
                    if j != label and not is_excluded:
                        paralogs.append(j)
            visited.append(i)
        return paralogs
