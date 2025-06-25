from typing import Set, Dict, List


class Node:
    """Node
    Store basic information about an entry in GFF file
    and parents and/or children of the entry
    """

    def __init__(self, idx: str) -> None:
        """
        Initializes a new `Node` object.

        Parameters:
            idx (str): Identifier, from "ID" attribute of GFF file.

        Attributes:
            _parents (set[str]): A set of parent indices.
            _children (set[str]): A set of child indices.

        """
        self.idx = idx
        self._parents: Set[str] = set()
        self._children: Set[str] = set()

    def add_parent(self, parent: str) -> None:
        self._parents.add(parent)

    def add_child(self, child: str) -> None:
        self._children.add(child)

    @property
    def parents(self) -> Set[str]:
        return self._parents

    @property
    def children(self) -> Set[str]:
        return self._children

    @property
    def is_root(self) -> bool:
        return len(self._parents) == 0

    @property
    def is_leaf(self) -> bool:
        return len(self._children) == 0

    @property
    def is_singleton(self) -> bool:
        return len(self._parents) == 0 and len(self._children) == 0

    def __repr__(self) -> str:
        return f"Node(id: {self.idx}. parents: ({', '.join(self._parents)}) children: ({','.join(self._children)}))"


def get_genes(features: Dict[str, Node], idx: str) -> Set[str]:
    """get_genes get genes
    Given a dictionary of Nodes and an id, return the set of genes that are parents of the node
    Args:
        features: dict[str, Node], dictionary of Nodes (features)
        idx: str, unique id of the node

    Returns:
        Set[str], set of genes
    """
    genes: Set[str] = set()
    stack: List[Node] = [features[idx]]
    while stack:
        node: Node = stack.pop()
        if node.is_root:
            genes.add(node.idx)
        else:
            for parent in node.parents:
                stack.append(features[parent])
    return genes
