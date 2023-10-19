"""treedraw.py
===========.

Adapted from:
https://github.com/cvzi/py_treedraw

This module implements a layout algorithm for a tree.
Each node in the tree has a layout property that contains an x and y
coordinate.
The layout is calculated by calling the method "walker(distance)".
The distance property indicated the distance between nodes on the same level.
If the tree is build up using the addChild/removeChild methods, the layout
will be calculated in linear time.
The algoithm is a python implemenation of this publication
<a href="http://citeseer.ist.psu.edu/buchheim02improving.html"> "Improving
Walker's Algorithm to Run in Linear Time"</a> by Christoph Buchheim, Michael
Junger, Sebastian Leipert"
"""
from __future__ import annotations

__all__ = ["Tree", "Node"]

__version__ = "1.4"

import math

try:
    xrange(5)
    myrange = xrange
except NameError:
    myrange = range


class Tree:
    def __init__(self, data) -> None:
        # Start a tree with a root node
        self.root = Node(data)
        self.root.tree = self
        self.root.leftSibling = None
        self.nodes = [self.root]

    def addNode(self, node):
        # add a child to the root
        return self.root.addNode(node)

    def addChild(self, data):
        # add a child to the root
        return self.root.addChild(data)

    def removeChild(self, node):
        # remove a child from the root
        return self.root.removeChild(node)

    def walker(self, distance=1.0):
        # Init layout algorithm
        self._firstWalk(self.root, distance)
        self._secondWalk(self.root, -self.root.layout.prelim)

    def _add(self, node):
        if node not in self.nodes:
            self.nodes.append(node)
        return node

    def _firstWalk(self, v, distance):
        if v.children:
            defaultAncestor = v.children[0]
            for w in v.children:
                self._firstWalk(w, distance)
                self._apportion(w, defaultAncestor, distance)
            self._executeShifts(v)
            midpoint = 0.5 * (v.children[0].layout.prelim + v.children[-1].layout.prelim)
            w = self._leftSibling(v)
            if w is not None:
                v.layout.prelim = w.layout.prelim + distance
                v.layout.mod = v.layout.prelim - midpoint
            else:
                v.layout.prelim = midpoint
        else:
            ls = self._leftSibling(v)
            if ls:
                v.layout.prelim = ls.layout.prelim + distance

    def _secondWalk(self, v, m):
        v.layout.x(v.layout.prelim + m)
        for w in v.children:
            self._secondWalk(w, m + v.layout.mod)

    def _apportion(self, v, defaultAncestor, distance):
        w = self._leftSibling(v)
        if w is not None:
            v_p_o = v
            v_p_i = v
            v_m_i = w
            v_m_o = v_p_i.parent.children[0]
            s_p_i = v_p_i.layout.mod
            s_p_o = v_p_o.layout.mod
            s_m_i = v_m_i.layout.mod
            s_m_o = v_m_o.layout.mod
            while v_m_i.nextRight() and v_p_i.nextLeft():
                v_m_i = v_m_i.nextRight()
                v_p_i = v_p_i.nextLeft()
                v_m_o = v_m_o.nextLeft()
                v_p_o = v_p_o.nextRight()
                v_p_o.layout.ancestor = v
                shift = v_m_i.layout.prelim + s_m_i - (v_p_i.layout.prelim + s_p_i) + distance
                if shift > 0:
                    self._moveSubtree(self._ancestor(v_m_i, v, defaultAncestor), v, shift)
                    s_p_i = s_p_i + shift
                    s_p_o = s_p_o + shift
                s_m_i = s_m_i + v_m_i.layout.mod
                s_p_i = s_p_i + v_p_i.layout.mod
                s_m_o = s_m_o + v_m_o.layout.mod
                s_p_o = s_p_o + v_p_o.layout.mod
            if v_m_i.nextRight() and v_p_o.nextRight() is None:
                v_p_o.layout.thread = v_m_i.nextRight()
                v_p_o.layout.mod = v_p_o.layout.mod + s_m_i - s_p_o
            if v_p_i.nextLeft() and v_m_o.nextLeft() is None:
                v_m_o.layout.thread = v_p_i.nextLeft()
                v_m_o.layout.mod = v_m_o.layout.mod + s_p_i - s_m_o
                defaultAncestor = v
        return defaultAncestor

    @staticmethod
    def _leftSibling(v):
        if v.leftSibling != -1:
            return v.leftSibling
        else:
            if v.parent is None or not v.parent.children:
                return None
            last = None
            for w in v.parent.children:
                if w == v:
                    if last is not None:
                        return last
                    return None
                last = w
            return None

    @staticmethod
    def _moveSubtree(w_m, w_p, shift):
        subtrees = w_p.number() - w_m.number()
        if subtrees == 0:
            subtrees = 0.0000001
        w_p.layout.change = w_p.layout.change - shift / subtrees
        w_p.layout.shift = w_p.layout.shift + shift
        w_m.layout.change = w_m.layout.change + shift / subtrees
        w_p.layout.prelim = w_p.layout.prelim + shift
        w_p.layout.mod = w_p.layout.mod + shift

    @staticmethod
    def _executeShifts(v):
        shift = 0
        change = 0
        i = len(v.children)
        for i in myrange(len(v.children) - 1, -1, -1):
            w = v.children[i]
            w.layout.prelim = w.layout.prelim + shift
            w.layout.mod = w.layout.mod + shift
            change = change + w.layout.change
            shift = shift + w.layout.shift + change

    @staticmethod
    def _ancestor(v_i, v, defaultAncestor):
        if v_i.layout.ancestor.parent == v.parent:
            return v_i.layout.ancestor
        return defaultAncestor


class Node:
    class Layout:
        def __init__(self, v) -> None:
            self.node = v
            self.mod = 0
            self.thread = None
            self.ancestor = v
            self.prelim = 0
            self.shift = 0
            self.change = 0
            self.pos = [None, None]
            self.number = -1  # undefined

        def x(self, value=None):
            if value is not None:
                self.pos[0] = value
                return None
            else:
                return self.pos[0]

        def y(self, value=None):
            if value is not None:
                self.pos[1] = value
                return None
            else:
                if self.pos[1] is None:
                    self.pos[1] = self.node.level()
                return self.pos[1]

    def __init__(self, data) -> None:
        self.tree = None
        self.data = data
        self.leftSibling = -1  # undefined, outdated
        self.children = []
        self.parent = None
        self.layout = self.Layout(self)

    def addNode(self, node):
        # Add an existing tree/node as a child

        # Set left sibling
        if self.children:
            node.leftSibling = self.children[-1]
            node.layout.number = node.leftSibling.layout.number + 1
        else:
            node.leftSibling = None
            node.layout.number = 0

        # Append to node
        node.parent = self
        self.children.append(node)

        # Add to tree
        root = self
        i = 0
        while root.parent is not None:
            root = root.parent
            i += 1
        root.tree._add(node)
        node.tree = root.tree
        self.layout.pos[1] = i  # Level
        return node

    def addChild(self, data):
        # Create a new node and add it as a child
        return self.addNode(Node(data))

    def removeChild(self, v):
        j = -1
        for i in myrange(len(self.children)):
            if self.children[i] == v:
                del self.children[i]
                j = i
                break

        for i in myrange(len(self.tree.nodes)):
            if self.tree.nodes[i] == v:
                del self.tree.nodes[i]
                break

        # Update left sibling
        if j == 0:
            self.children[0].leftSibling = None
        elif j > 0:
            self.children[j].leftSibling = self.children[j - 1]
        else:  # j == -1
            return

        # Update numbers
        for i in myrange(j, len(self.children)):
            self.children[i].layout.number = i

        # Remove children of the deleted node
        i = 0
        while i < len(self.tree.nodes):
            if self.tree.nodes[i] in v.children:
                del self.tree.nodes[i]
            else:
                i += 1

        v.children = []

    def nextLeft(self):
        if self.children:
            return self.children[0]
        else:
            return self.layout.thread

    def nextRight(self):
        if self.children:
            return self.children[-1]
        else:
            return self.layout.thread

    def level(self):
        if self.layout.pos[1] is not None:
            return self.layout.pos[1]
        n = self.parent
        i = 0
        while n is not None:
            n = n.parent
            i += 1
        return i

    def position(self, origin=(0, 0), scalex=1.0, scaley=None):
        """Return position as integer
        Examples:
        position(origin)
        position(origin, 10)
        position(origin, 10, 15)
        position(origin, (10, 15)).
        """
        if scaley is None:
            if hasattr(scalex, "__getitem__"):
                scaley = scalex[1]
                scalex = scalex[0]
            else:
                scaley = scalex
        return (
            origin[0] + int(math.ceil(self.layout.x() * scalex)),
            origin[1] + int(math.ceil(self.layout.y() * scaley)),
        )

    def positionf(self, origin=(0.0, 0.0), scalex=1.0, scaley=None):
        """Return position as floating point
        Examples:
        position(origin)
        position(origin, 10)
        position(origin, 10, 15)
        position(origin, (10, 15)).
        """
        if scaley is None:
            if hasattr(scalex, "__getitem__"):
                scaley = scalex[1]
                scalex = scalex[0]
            else:
                scaley = scalex
        return (origin[0] + (self.layout.x() * scalex), origin[1] + (self.layout.y() * scaley))

    def number(self):
        if self.layout.number != -1:
            return self.layout.number
        else:
            if self.parent is None:
                return 0
            else:
                i = 0
                for node in self.parent.children:
                    if node == self:
                        return i
                    i += 1
            msg = "Error in number(self)!"
            raise Exception(msg)
