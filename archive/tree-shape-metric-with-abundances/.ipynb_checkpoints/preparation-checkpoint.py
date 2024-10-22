import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
from ete3 import Tree, TreeNode
from gctree import CollapsedTree, CollapsedForest
from typing import List, Union, Optional, Callable


#####---------------------------------------------------------------------------------------#####
##### CONFIGURATIONS
#####---------------------------------------------------------------------------------------#####
feature_names = ['n_nodes',
                 'n_leaves',
                 'p_leaves',
                 'n_internals',
                 'p_internals',
                 'n_ptn',
                 'p_ptn',
                 'n_obs',
                 'p_obs',
                 'n_inf',
                 'p_inf',
                 'n_cells',
                 'min_abund',
                 'max_abund',
                 'avg_abund',
                 'min_depth',
                 'max_depth',
                 'avg_depth',
                 'min_topodepth',
                 'max_topodepth',
                 'avg_topodepth',
                 'trunk',
                 'topotrunk',
                 'min_dlfirst_split_node',
                 'max_dlfirst_split_node',
                 'avg_dlfirst_split_node',
                 'min_od',
                 'max_od',
                 'avg_od',
                 'min_od2',
                 'avg_od2',
                 'min_dlsn',
                 'max_dlsn',
                 'avg_dlsn'
                ]

#####---------------------------------------------------------------------------------------#####
##### HELPER FUNCTIONS
#####---------------------------------------------------------------------------------------#####
def agg_values(vals, base_name: str, agg_funcs: list[str] = None, fun_set: str = 'tree') -> pd.Series:
            """ Aggregate values of a given sequence of integers.
        
            :param vals: A sequence of integers.
            :param base_name: A name for the given sequence for the index of the resulting series.
            :param agg_funcs: Optionally, the names of the functions that should be applied.
                'sum', 'min', 'max', 'avg', 'std' and 'gini' are possible.
            :param fun_set: A predefined set of functions for either tree features ('tree'): 'min', 'max', 'avg'
                or forest features ('forest'): 'sum', 'avg', 'std', 'max', 'gini'.
            :return: A pandas Series of the computed aggregate values.
            """
            vals = pd.Series(vals, dtype=int)
            index = []
            d = []
        
            if agg_funcs is None:
                if fun_set == 'tree':
                    index = [stat + '_' + base_name for stat in ['min', 'max', 'avg']]
                    d = [vals.min(), vals.max(), vals.mean()]
                elif fun_set == 'forest':
                    index = [stat + '_' + base_name for stat in ['sum', 'avg', 'std', 'max', 'gini']]
                    d = [vals.sum(), vals.mean(), vals.std(), vals.max(), gini(vals)]
            else:
                if 'sum' in agg_funcs:
                    index += ['sum_' + base_name]
                    d += [vals.sum()]
                if 'min' in agg_funcs:
                    index += ['min_' + base_name]
                    d += [vals.min()]
                if 'max' in agg_funcs:
                    index += ['max_' + base_name]
                    d += [vals.max()]
                if 'avg' in agg_funcs:
                    index += ['avg_' + base_name]
                    d += [vals.mean()]
                if 'std' in agg_funcs:
                    index += ['std_' + base_name]
                    d += [vals.std()]
                if 'gini' in agg_funcs:
                    index += ['gini_' + base_name]
                    d += [gini(vals)]
        
            return pd.Series(d, index=index)
    
def f_agg_values(vals, base_name: str) -> pd.Series:
    """ A wrapper for the agg_values function, computing aggregate values 'sum', 'avg', 'std', 'max', 'gini',
    used for forest features. """
    return agg_values(vals, base_name, fun_set='forest')

def gini(ll):
    """ Computes the Gini coefficient for a given sequence. """
    ll = sorted(ll)
    n_ = len(ll)
    if n_ == 0 or sum(ll) <= 0:
        return 0.
    coef_ = 2. / n_
    const_ = (n_ + 1.) / n_
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(ll)])
    return coef_*weighted_sum/(sum(ll)) - const_


def get_mid(sample_dir: str) -> int:
    """ Get the unique MID from a sample directory name. """
    return int(re.findall('[0-9]+', sample_dir)[0])


def od(node: TreeNode) -> int:
    """ Computes the outgoing degree of a given TreeNode. """
    return len(node.children)


def trunk(t, topo: bool = False) -> float:
    """ Computes the (topological) trunk length for a given GCtree following the definition in the thesis.
    If the tree does not contain any split nodes, the trunk length is equal to the depth of its only leaf. """
    if t.first_split_node() is not None:
        return t.node_depth(t.first_split_node(), topo=topo)
    else:
        return t.node_depth(t.leaves[0], topo=topo)


def depth(t, topo: bool = False) -> float:
    """ Computes the (topological) depth, aka maximum path length,
    for a given GCtree following the definition in the thesis. """
    return max(t.node_depth(node, topo=topo) for node in t.nodes)
    
#####---------------------------------------------------------------------------------------#####
##### MAIN CLASSES
#####---------------------------------------------------------------------------------------#####

##### GCTREE #####
class GCtree(CollapsedTree):
    """
    A subclass of the gctree CollapsedTree class.
    Provides some handy attributes for the analysis of tree characteristics.
    """
    def __init__(self, tree: TreeNode = None, path: str = None):
        super().__init__(tree=None, allow_repeats=False)
        self.tree = tree
        self.path = path
        self.tree_name = path.split("/")[-1]
        self.sample_name = path.split("/")[-2]
        self.root = tree
        self.nodes = list(self.tree.traverse())
        self.leaves = list(self.tree)
        self.internal_nodes = [node for node in self.nodes if not node.is_leaf() and not node.is_root()]
        self.passthrough_nodes = [node for node in self.internal_nodes if len(node.children) == 1]
        self.split_nodes = [node for node in self.internal_nodes if len(node.children) > 1]
        self.observed_nodes = [node for node in self.nodes if node.abundance > 0]
        self.inferred_nodes = [node for node in self.nodes if node.abundance == 0]

    def node_depth(self, node: TreeNode, topo: bool = False) -> float:
        """ The (topological) path length from the root to a given node.

        :param node: The node whose depth is calculated.
        :param topo: True for topological depth.
        :return: The depth of the given node in the tree.
        """
        return node.get_distance(self.root, topology_only=topo)
        
    def first_split_node(self) -> Optional[TreeNode]:
        """ The first split node in the tree, if existent. """
        if len(self.root.children) > 1:
            return self.root
        elif len(self.split_nodes) > 0:
            return min(self.split_nodes, key=lambda s: self.node_depth(s, topo=True))
        else:
            return None
        
    def t_size_features(self) -> pd.Series:
        """ GCtree size properties as defined in Table 3.1 in the thesis for a given GCtree. """
        n_nodes = len(self.nodes)
        n_leaves = len(self.leaves)
        n_internals = len(self.internal_nodes)
        n_ptn = len(self.passthrough_nodes)
    
        index = ['n_nodes', 'n_leaves', 'p_leaves', 'n_internals', 'p_internals', 'n_ptn', 'p_ptn']
        d = [n_nodes, n_leaves, n_leaves / n_nodes, n_internals, n_internals / n_nodes, n_ptn, n_ptn / n_nodes]
        return pd.Series(data=d, index=index)

    
    def t_abundance_features(self) -> pd.Series:
        """ GCtree abundance properties as defined in Table 3.1 in the thesis for a given GCtree. """
        obs_abundances = [node.abundance for node in self.observed_nodes]
    
        n_nodes = len(self.nodes)
        n_obs = len(obs_abundances)
        n_inf = n_nodes - n_obs
        n_cells = sum(obs_abundances)
    
        index = ['n_obs', 'p_obs', 'n_inf', 'p_inf', 'n_cells']
        d = [n_obs, n_obs / n_nodes, n_inf, n_inf / n_nodes, n_cells]
    
        ser = pd.Series(data=d, index=index)
        agg_abundance = agg_values(obs_abundances, base_name='abund')
    
        return pd.concat([ser, agg_abundance], axis=0)
        
    def t_len_features(self) -> pd.Series:
        """ GCtree length properties as defined in Table 3.1 in the thesis for a given GCtree. """
        depths = [self.node_depth(leaf) for leaf in self.leaves]
        topo_depths = [self.node_depth(leaf, topo=True) for leaf in self.leaves]
        dlfirst_split_node = []
        if self.first_split_node() is not None:
            dlfirst_split_node = [leaf.get_distance(target=self.first_split_node(), topology_only=False) for leaf in self.leaves]
        
        trunk_features = pd.Series(index=['trunk', 'topotrunk'], data=[trunk(self, topo=False), trunk(self, topo=True)])
    
        return pd.concat([agg_values(depths, 'depth'), agg_values(topo_depths, 'topodepth'), trunk_features,
                          agg_values(dlfirst_split_node, 'dlfirst_split_node')], axis=0)
    
    
    def t_bushiness_features(self) -> pd.Series:
        """ GCtree bushiness properties as defined in Table 3.1 in the thesis for a given GCtree. """
        ods = [len(node.children) for node in self.nodes if not node.is_leaf()]
        ods2 = [len(split.children) for split in self.split_nodes]
        dlsn = []
        if len(ods2) > 0:
            dlsn = [min(leaf.get_distance(target=split, topology_only=False) for split in self.split_nodes) for leaf in self.leaves]
    
        return pd.concat([agg_values(ods, 'od'), agg_values(ods2, 'od2', ['min', 'avg']),
                          agg_values(dlsn, 'dlsn')], axis=0)
        
    def summarize_tree_features(self) -> pd.Series:
        """ GCtree features as defined in Table 3.1 in the thesis for a given GCtree.
        Additionally, sample MID and name of the tree. """
        d = [self.sample_name, self.tree_name]
        index = ["Sample", "Tree"]
        info_features = pd.Series(d, index=index)
    
        return pd.concat([info_features, self.t_size_features(), self.t_abundance_features(), self.t_len_features(),
                          self.t_bushiness_features()], axis=0)

def read_tree_from_path(treedir):
    treedir = str(treedir)
    nk_path = os.path.join(treedir, "gctree.out.inference.1.nk")
    if os.path.exists(nk_path) == True:
        ab_dict_path = os.path.join(treedir, "abund.csv")
        abund_df = pd.read_csv(ab_dict_path, index_col=0, names=['val'])
        ab_dict = abund_df.to_dict().get('val')
        tree_path = treedir
        tree = Tree(newick=nk_path, format=1)
        if ab_dict is not None:
            for node in tree.traverse():
                node.add_feature('abundance', ab_dict.get(node.name, 0))
        treeobj = GCtree(tree = tree, path = tree_path)
        return treeobj
#####---------------------------------------------------------------------------------------#####
##### GCFOREST #####

class GCforest(CollapsedForest):
    """
    A subclass of the gctree CollapsedForest class.
    """

    def __init__(
        self,
        forest: List[GCtree] = None
    ):
        CollapsedForest.__init__(self, forest=None)
        if forest is None:
            forest = list()
        self.tree_list = forest

    @property
    def tree_list(self) -> list[GCtree]:
        """ A list of all GCtrees in the forest. """
        return self._ctrees

    @tree_list.setter
    def tree_list(self, value: list[GCtree]):
        self._ctrees = value
        self.n_trees = len(value)

    def sample_tree(self, cond: Callable[[GCtree], bool] = None) -> GCtree:
        """ Sample a random GCtree from the forest, optionally matching the given condition. """
        if cond is None:
            return random.choice(self.tree_list)
        else:
            return random.choice(self.choose_trees(cond))

    def choose_tree(self, cond: Callable[[GCtree], bool]) -> GCtree:
        """ Return the first GCtree matching the given condition. """
        return self.choose_trees(cond)[0]

    def choose_trees(self, cond: Callable[[GCtree], bool]) -> List[GCtree]:
        """ Return a list of every GCtree in the forest matching the given condition. """
        return list(filter(cond, self.tree_list))

class LabForest(GCforest):
    """
    A subclass of the above defined GCforest class.
    Represents a forest of imported GCtrees that all belong to the same lab sample.
    """

    def __init__(
        self,
        forest: list[GCtree] = None,
        path: str = None
    ):
        GCforest.__init__(self, forest)
        self._path = path
        self._missing_trees = list()

    @property
    def path(self) -> str:
        """ The absolute path to the forest directory in which the tree directories are located. """
        return self._path

    @property
    def mid(self) -> Optional[int]:
        """ The unique MID that identifies the lab sample on which phylogenetic inference has been performed. """
        if self.path is not None:
            sample_dir = self.path.split('/')[-1]
            return get_mid(sample_dir)
        else:
            return None

    @property
    def missing_trees(self) -> list[int]:
        """ A list of cell counts for each genotype with missing tree files. """
        return self._missing_trees

    @missing_trees.setter
    def missing_trees(self, value):
        self._missing_trees = value

    @property
    def n_trees_plus_missing(self):
        """ Number of trees in the forest, including genotypes without tree files. """
        return self.n_trees + len(self.missing_trees)


def lab_forest_features(f: LabForest, MID: str) -> pd.Series:
    """ Computes a feature vector for a given forest of imported GCtrees.
    The choice of features is based on insights gathered during the work on the thesis. """
    node_counts = [len(t.nodes) for t in f.tree_list]
    leaf_counts = [len(t.leaves) for t in f.tree_list]
    # ods for non-leaf nodes:
    ods = [od_ for t in f.tree_list for od_ in [od(v) for v in [t.root] + t.internal_nodes]]
    # ods for split nodes:
    ods2 = [od_ for od_ in ods if od_ > 1]
    node_abunds = [ab_ for t in f.tree_list for ab_ in [v.abundance for v in t.nodes]]
    # total cellular abundance for each tree, including missing tree's abundances:
    tree_abunds = [sum(v.abundance for v in t.nodes) for t in f.tree_list] + f.missing_trees
    # cellular abundance for observed genotypes:
    obs_abunds = [ab_ for ab_ in node_abunds if ab_ > 0]
    topo_depths = [depth(t, topo=True) for t in f.tree_list]
    depths = [depth(t, topo=False) for t in f.tree_list]
    trunks = [trunk(t, topo=False) for t in f.tree_list]
    # topo_trunks not included as they are trivial for regarded data (1 for almost every tree).

    tot_node_count = sum(node_counts)

    misc_index = ['MID', 
                  'n_trees_plus_missing', 
                  'n_trees',
                  'p_singletons',
                  'p_leaves',
                  'p_observed',
                  'avg_node_abund']
    misc_d = [MID, 
              f.n_trees_plus_missing, 
              f.n_trees,
              tree_abunds.count(1) / f.n_trees_plus_missing,
              sum(leaf_counts) / tot_node_count,
              len(obs_abunds) / len(node_abunds),
              np.mean(obs_abunds)]

    misc_features = pd.Series(misc_d, index=misc_index)
    node_counts = f_agg_values(node_counts, 'n_nodes')
    tree_abunds = f_agg_values(tree_abunds, 'tree_abund')
    ods = f_agg_values(ods, 'od')
    ods2 = f_agg_values(ods2, 'od2')
    topo_depths = f_agg_values(topo_depths, 'topodepth')
    depths = f_agg_values(depths, 'depth')
    trunks = f_agg_values(trunks, 'trunk')
    outputdf = pd.concat([misc_features, 
                          pd.concat([node_counts, tree_abunds, ods, ods2, topo_depths, depths, trunks], axis = 0).round(decimals=4)])
    return pd.DataFrame(outputdf)