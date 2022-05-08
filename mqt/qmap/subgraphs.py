import os.path

import networkx as nx
import matplotlib.pyplot as plt
import itertools
import json
from typing import List, Set

path_prefix = '../../test/'


def compute_unique_subgraphs(graph, filename: str):
    sgs = {}

    with open(path_prefix + "subgraphs/" + filename + ".txt", encoding='utf-8', mode='w+') as f:
        for k in range(1, graph.number_of_nodes() + 1):
            sgs[k] = []
            print("---", k, "---")
            for sg in (graph.subgraph(selected_nodes) for selected_nodes in itertools.combinations(graph, k)):
                if nx.is_connected(sg):
                    new_class = True
                    for g in sgs[k]:
                        if nx.is_isomorphic(sg, g):
                            new_class = False
                            break
                    if new_class:
                        f.write(str(sg.nodes) + '\n')
                        sgs[k].append(sg)
            f.write("---\n")
            print(len(sgs[k]))
        ncols = max([len(v) for v in sgs.values()])
        print("Max number of unique subgraphs:", ncols)


def load_subgraphs(graph, filename: str):
    sgs = {}
    idx = 1
    with open(path_prefix + 'subgraphs/' + filename + ".txt", encoding='utf-8') as f:
        sgs[idx] = []
        for line in f:
            line = line.rstrip('\n')
            if line == "---":
                idx += 1
                sgs[idx] = []
            else:
                sgs[idx].append(graph.subgraph(json.loads(line)))
    return sgs


def load_subgraphs_from_file(filename: str, nqubits: int) -> List[Set[int]]:
    sgs = []
    with open(filename, encoding='utf-8') as f:
        q = 1
        for line in f:
            line = line.rstrip('\n')
            if q == nqubits:
                if line == "---":
                    return sgs
                else:
                    sgs.append(set(eval(line)))
                    print(sgs[-1])
            else:
                if line == "---":
                    q += 1
    return sgs


def draw_subgraphs(graph, filename: str):
    sgs = load_subgraphs(graph, filename)

    ncols = max([len(v) for v in sgs.values()])
    nrows = graph.number_of_nodes()

    aspect_ratio = float(ncols) / float(nrows)
    fig_height = 5. * nrows

    width = fig_height * aspect_ratio
    height = fig_height
    if width > 655.36:
        ratio = 655.35 / width
        width *= ratio
        height *= ratio

    plt.subplots(figsize=(width, height))
    for row in range(nrows):
        graphs = sgs[row + 1]
        for col in range(len(graphs)):
            graph = graphs[col]
            plt.subplot2grid((nrows, ncols), (row, col))
            nx.draw_networkx(graph, pos=nx.kamada_kawai_layout(graph), with_labels=False)

    plt.tight_layout()
    plt.savefig(fname=arch + '.pdf', format='pdf')


def create_architecture_files(graph, filename: str, max_nodes: int = 10):
    sgs = load_subgraphs(graph, filename)

    if not os.path.exists(path_prefix + "subgraphs/" + filename):
        os.makedirs(path_prefix + "subgraphs/" + filename)

    for k in range(1, min(max_nodes, graph.number_of_nodes()) + 1):
        if not os.path.exists(path_prefix + "subgraphs/" + filename + "/" + str(k)):
            os.makedirs(path_prefix + "subgraphs/" + filename + "/" + str(k))

        graphs = sgs[k]
        idx = 0
        for subgraph in graphs:
            idx += 1
            with open(path_prefix + "subgraphs/" + filename + "/" + str(k) + "/" + filename + '_' + str(k) + '_' + str(idx) + ".arch", 'w+') as f:
                f.write(str(subgraph.number_of_nodes()) + '\n')
                for edge in subgraph.edges():
                    f.write("{} {}\n".format(edge[0], edge[1]))
                    f.write("{} {}\n".format(edge[1], edge[0]))


if __name__ == '__main__':
    archs = {
        # 'ibm_lagos': ["0 1", "1 2", "1 3", "3 5", "4 5", "5 6"],
        # 'ibmq_quadalupe': ["0 1", "1 2", "2 3", "3 5", "1 4", "5 8",
        #                    "4 7", "6 7", "8 9", "7 10", "8 11", "10 12",
        #                    "12 15", "12 13", "13 14", "11 14"],
        # 'rigetti_16q': ["0 1", "1 2", "2 3", "3 4", "4 5", "5 6", "6 7", "7 0",
        #                 "10 11", "11 12", "12 13", "13 14", "14 15", "15 16", "16 17", "17 10",
        #                 "2 15", "1 16"],
        # 'sycamore_12q': ["0 4", "1 4", "1 5", "2 5", "2 6", "3 6", "3 7",
        #                  "4 8", "4 9", "5 9", "5 10", "6 10", "6 11", "7 11"],
        # 'sycamore_16q': ["0 4", "1 4", "1 5", "2 5", "2 6", "3 6", "3 7",
        #                  "4 8", "4 9", "5 9", "5 10", "6 10", "6 11", "7 11",
        #                  "8 12", "9 12", "9 13", "10 13", "10 14", "11 14", "11 15"],
        'ibmq_ehningen': ["0 1", "1 2", "2 3", "3 5", "1 4", "5 8",
                          "4 7", "6 7", "8 9", "7 10", "8 11", "10 12",
                          "12 15", "12 13", "13 14", "11 14",
                          "12 15", "14 16", "15 18", "16 19", "17 18",
                          "19 20", "18 21", "19 22", "21 23", "22 25",
                          "23 24", "24 25", "25 26"]
    }

    for arch, edgelist in archs.items():
        G = nx.parse_edgelist(edgelist, nodetype=int)

        # compute_unique_subgraphs(G, arch)

        # draw_subgraphs(G, arch)

        create_architecture_files(G, arch)
