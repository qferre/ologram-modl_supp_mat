"""
CL-MAX approximate itemset miner.

Source: https://github.com/MohsenFatemii/CL-MAX


Example:

>>> name_of_dataset = 'chess'
>>> min_support = 0.9
>>> num_cluster = 10
>>> rounding_threshold = 0.9

>>> cl_max = CL_MAX(min_support, num_cluster, rounding_threshold)
>>> path = cl_max.find_path()
>>> cl_max.read_dataset_from_file(name_of_dataset, ',')
>>> cl_max.load_one_hot_dataset(name_of_dataset)
>>> cl_max.remove_non_frequent_single_items()
>>> MFIs = cl_max._CL_MAX()
>>> print(MFIs)

"""

from functools import lru_cache
import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
import json
import os

from functools import lru_cache
import numpy as np
import json


# ---------------------------------------------------------------------------- #


class TransVerticalBitmaps:
    """ Stores transactions as vertical bitmap of individual items and helps count itemset supports efficiently 
    Parameters
    ----------
    transactions : list of sets
            The list of transactions
    Attributes
    ----------
    transactions : list of sets
            The list of transactions
    n_transactions : int
            The number of transactions
    items : list
            The list of all items in sorted order
    items_vertical_bitmaps : dict
            The dictionary of vertical bitmap representation of transactions indexed by item
    """

    def __init__(self, transactions):

        self.transactions = transactions
        self.n_transactions = len(self.transactions)

        # Extract the list of items in the transactions
        items = set()
        for transaction in self.transactions:
            items.update(transaction)
        self.items = sorted(items)

        self.items_vertical_bitmaps = {item: np.zeros(
            shape=(self.n_transactions,), dtype=np.bool) for item in self.items}

        for i_transaction, transaction in enumerate(self.transactions):
            for item in transaction:
                self.items_vertical_bitmaps[item][i_transaction] = True

    @lru_cache(maxsize=32)
    def compVerticalBitmap(self, itemset):
        """ Compute the vertical bitmap of the given itemset

        Parameters
        ----------
        itemset : tuple
                The tuple of items (or itemset) for which support is to be counted
        Returns
        -------
        np.array(bool)
                Vertical bitmap of transactions in the which itemset is present
        """

        if len(itemset) == 1:
            item = itemset[0]
            return self.items_vertical_bitmaps[item]

        else:
            last_item = itemset[-1]
            return self.compVerticalBitmap(itemset[:-1]) & self.items_vertical_bitmaps[last_item]

    def countSupport(self, itemset):
        """ Count the support of the itemset in the transactions 
        Parameters
        ----------
        itemset : tuple
                The tuple of items for which support is to be counted
        Returns
        -------
        int
                The support count of the itemset in the transactions
        """

        itemset_vertical_bitmap = self.compVerticalBitmap(itemset)
        itemset_support_count = np.count_nonzero(itemset_vertical_bitmap)

        return itemset_support_count


class MafiaNode:
    """ A node in the MAFIA candidate itemset tree 

    Parameters
    ----------
    head : set
            The set of all items in head of the node
    tail : tuple
            The tuple of all items in tail of the node
    support_count : int optional(default=None)
            The support count of head of the node
    Attributes
    ----------
    head : set
            The set of all items in head of the node
    tail : tuple
            The tuple of all items in tail of the node
    support_count : int
            The support count of head of the node
    Notes
    -----
            The :attr:`head` is of type tuple as opposed to set or list since it is passed to :meth:`compVerticalBitmap`
            of :class:`transactions` that caches results and therefore requires the inputs to **hashable**.
    """

    def __init__(self, head, tail, support_count=None):
        self.head = head
        self.tail = tail.copy()
        self.support_count = support_count


def _mafiaAlgorithm(current_node, MFIs, transactions, min_support_count):

    # HUTMFI Pruning - Prune the subtree, if the HUT of the node has any superset in MFI
    head_union_tail = current_node.head + tuple(current_node.tail)
    if any(all(item in mfi for item in head_union_tail) for mfi in MFIs):
        return

    # Count the support of all children of the node
    node_children_support_cnts = [(item, transactions.countSupport(
        current_node.head + (item,))) for item in current_node.tail]
    # Extract the frequent children of the node and their support counts
    node_freq_children_sup_cnts = [
        (item, support_count) for item, support_count in node_children_support_cnts if support_count >= min_support_count]

    # The items in tail with same support as parent
    node_children_items_parent_eq = []
    # The items in node's tail (except parent equivalence items) sorted by desc support
    node_tail_items_sup_cnts = []

    for item, support_count in node_freq_children_sup_cnts:
        if support_count == current_node.support_count:
            node_children_items_parent_eq.append(item)
        else:
            node_tail_items_sup_cnts.append((item, support_count))

    # Sort the items in the trimmed tail by increasing support
    node_tail_items_sup_cnts.sort(key=lambda x: x[1])
    node_tail_items = [item for item, support in node_tail_items_sup_cnts]

    current_node.head += tuple(node_children_items_parent_eq)
    current_node.tail = node_tail_items

    is_leaf = not bool(current_node.tail)

    for i, item in enumerate(current_node.tail):
        new_node_head = current_node.head + (item,)
        new_node_tail = current_node.tail[i+1:]
        new_node_support_cnt = node_tail_items_sup_cnts[i][1]
        new_node = MafiaNode(new_node_head, new_node_tail,
                             new_node_support_cnt)

        is_hut = (i == 0)  # if i is the first element in the tail
        _mafiaAlgorithm(new_node, MFIs, transactions, min_support_count)

    # if current node is a leaf and no superset of current node head in MFIs
    if is_leaf and current_node.head and not any(all(item in mfi for item in current_node.head) for mfi in MFIs):
        MFIs.append(set(current_node.head))


def mafiaAlgorithm(transactions, min_support_count, possible_candidates, seen_dict):
    """ Extract the MFIs (Maximal Frequent Itemsets) from transactions with min support count using MAFIA Algorithm
    Parameters
    ----------
    transactions : list of sets
            The list of transactions
    min_support_count : int
            The minimum support count threshold
    possible_candidates : list of lists
            The Possible Candidates for being Maximal Frequent Itemsets
    seen_dict : dictionary
            A dictionary for prevent sending repetitive patterns to MAFIA
    Returns
    -------
    list of sets
            The list of all maximal frequent itemsets
    """

    transactions_vertical_bitmaps = TransVerticalBitmaps(transactions)
    MFIs = []
    i = 1
    for item in possible_candidates:
        # random.shuffle(item)
        # Create the root node of MAFIA candidate itemset tree
        str_item = json.dumps(item)
        #print(f"Extacting Frequent Itemsets from cluster number {i}")
        if seen_dict[str_item] == 0:
            seen_dict[str_item] = 1
            mafia_cand_itemset_root = MafiaNode(
                tuple(), item, transactions_vertical_bitmaps.n_transactions)
            # Perform the MAFIA algorithm
            _mafiaAlgorithm(mafia_cand_itemset_root, MFIs,
                            transactions_vertical_bitmaps, min_support_count)
        i += 1
    return MFIs








# ---------------------------------------------------------------------------- #




class CL_MAX():
    def __init__(self, min_support, num_of_clusters, rounding_threshold, batch_size=100):
        self.min_support = min_support
        self.num_of_clusters = num_of_clusters
        self.threshold = rounding_threshold
        self.seen_dict = {}
        self.path = ''
        self.batch_size = batch_size

    def find_path(self):
        self.path = os.getcwd()+'/Datasets'
        return self.path

    def read_dataset_from_file(self, name_of_dataset, delimiter):
        df = []
        maximum = 0
        with open(self.path+'/Actual Datasets/'+name_of_dataset+'.csv') as file:
            for line in file:
                transaction = [int(item)
                               for item in line.replace('\n', '').split(delimiter)]
                df.append(transaction)
                maximum = np.max([maximum, np.max(transaction)])
        self.maximum = maximum
        self.transactions = df

    def load_one_hot_dataset(self, name_of_dataset):
        df = pd.read_csv(self.path+'/One-hot/'+name_of_dataset+'01.csv')
        self.dataset = df

    def remove_non_frequent_single_items(self):
        # This step replaces sorted_items with an object like [ [count of item, item] for item in items]
        count = np.zeros(self.maximum+1)
        for transaction in self.transactions:
            for item in transaction:
                count[item] += 1

        for i in range(len(count)):
            if count[i] < self.min_support:
                self.dataset[i] = 0
        temp = [[count[i], i] for i in range(self.maximum)]
        self.sorted_items = sorted(temp)

    def cluster_transactions(self):
        kmeans = MiniBatchKMeans(n_clusters=self.num_of_clusters,
                                 batch_size=self.batch_size, max_iter=20).fit(self.dataset)
        labels = kmeans.labels_
        cnt_labels = np.zeros(self.num_of_clusters)
        clusters = np.zeros((self.num_of_clusters, self.dataset.shape[1]))
        for i in range(len(labels)):
            clusters[labels[i]] += np.array(self.dataset.iloc[i])
            cnt_labels[labels[i]] += 1
        for i in range(len(clusters)):
            clusters[i] /= cnt_labels[i]
        for i in range(len(clusters)):
            clusters[i][clusters[i] >= self.threshold] = 1
            clusters[i][clusters[i] < self.threshold] = 0
        clusters = np.array(clusters, dtype=np.int)
        self.clusters = clusters
        return clusters

    def convert_clusters_to_itemset(self):
        items = []
        for cluster in self.clusters:
            s = []
            for i in range(len(cluster)):
                if cluster[i] == 1:
                    s.append(i)
            if len(s) > 0:
                items.append(s)
        return items

    def order_correction(self, possible_candidates):
        # Return all candidates in correct order
        temp_candidates = []
        for candidate in possible_candidates:
            temp_can = []
            for i in range(self.maximum):
                if self.sorted_items[i][1] in candidate:
                    temp_can.append(self.sorted_items[i][1])
            temp_candidates.append(temp_can)
        return temp_candidates



    def _CL_MAX(self):
        clusters = self.cluster_transactions()
        possible_candidates = self.convert_clusters_to_itemset()
        print("Candidates from clustering:", possible_candidates)

        possible_candidates = self.order_correction(possible_candidates)
        for item in possible_candidates:
            str_item = json.dumps(item)
            self.seen_dict[str_item] = 0
        min_support_count = self.min_support*len(self.transactions)
        MFIs = mafiaAlgorithm(self.transactions, min_support_count,
                              possible_candidates, self.seen_dict)

        print("After MAFIA:", MFIs)
        return MFIs

