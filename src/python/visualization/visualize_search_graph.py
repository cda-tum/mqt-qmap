from distinctipy import distinctipy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import json
from typing import List, Tuple, Dict, Optional, Union, Callable, Any, Literal
import os
from mqt.qmap.visualization.treedraw import Tree
from ipywidgets import interactive, VBox, HBox, IntSlider, Play, jslink, Layout, Widget
from copy import deepcopy
from _plotly_utils.basevalidators import ColorValidator, ColorscaleValidator
import re
from dataclasses import dataclass
from random import shuffle

@dataclass
class TwoQbitMultiplicity:
    q0: int
    q1: int
    forward: int
    backward: int

@dataclass
class SearchNode:
    id: int
    parent: Optional[int]
    fixed_cost: float
    heuristic_cost: float
    lookahead_penalty: float
    is_valid_mapping: bool
    final: bool
    depth: int
    layout: Tuple[int, ...]
    swaps: Tuple[Tuple[int, int], ...]
    
    def total_cost(self):
        return self.fixed_cost + self.heuristic_cost + self.lookahead_penalty
    
    def total_fixed_cost(self):
        return self.fixed_cost + self.lookahead_penalty

Position = Tuple[float, float]
Colorscale = Union[str, List[str], List[Tuple[float, str]]]

def _is_len_iterable(obj) -> bool:
    if isinstance(obj, str):
        return False
    # this is actually the safest way to do this:
    # https://stackoverflow.com/questions/1952464/in-python-how-do-i-determine-if-an-object-is-iterable
    try:
        iter(obj)
        len(obj)
        return True
    except TypeError:
        return False

def _is_number(x):
    return type(x) is float or type(x) is int

def _remove_first_lines(string: str, n: int) -> str:
    return string.split("\n", n)[n]

def _get_avg_min_distance(seq):
    arr = sorted(seq)
    sum_dist = 0
    for i in range(1, len(arr)):
        sum_dist += arr[i] - arr[i-1]
    return sum_dist / (len(arr) - 1)

def _reverse_layout(seq):
    r = [-1]*len(seq)
    for i, v in enumerate(seq):
        r[v] = i
    return r

def _copy_to_dict(target: Dict[Any, Any], source: Dict[Any, Any], recursive: bool = True):
    for key, value in source.items():
        if recursive and isinstance(value, dict):
            if key not in target or not isinstance(target[key], dict):
                target[key] = value
            else:
                _copy_to_dict(target[key], value)
        else:
            target[key] = value
    return target

def _parse_search_graph(file_path: str, final_node_id: int) -> Tuple[nx.Graph, int]:
    graph = nx.Graph()
    root = None
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split(";")
            node = int(line[0])
            parent = int(line[1])
            data = SearchNode(
                node, 
                parent if parent != node else None,
                float(line[2]),
                float(line[3]),
                float(line[4]),
                line[5].strip() == "1",
                node == final_node_id,
                int(line[6]),
                tuple(int(q) for q in line[7].strip().split(",")),
                tuple() if len(line[8].strip()) == 0 else tuple(tuple(int(q) for q in swap.split(" ")) for swap in line[8].strip().split(","))
            )
            graph.add_node(node, data=data)
            if parent == node:
                root = node
            else:
                graph.add_edge(node, parent)
    return graph, root

def _layout_search_graph(
        search_graph: nx.Graph, 
        root: int,
        method: Literal['walker', 'dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'], 
        walker_factor: float, # 0-1
        tapered_layer_heights: bool) -> Dict[int, Position]:
    
    if method == 'walker':
        tree = Tree(root)
        nodes = {root: tree.root}
        for node in search_graph.nodes:
            if node != root:
                nodes[node] = nodes[search_graph.nodes[node]['data'].parent].addChild(node)
        tree.walker(walker_factor)
        pos = {}
        for node in tree.nodes:
            pos[node.data] = node.position(origin=(0, 0), scalex=100, scaley=-1)
    else:
        pos = graphviz_layout(search_graph, prog=method, root=search_graph.nodes[0])
        
    if not tapered_layer_heights:
        return pos
    
    layers_x = {}
    layer_spacings = {}
    for n in pos:
        x, y = pos[n]
        if y not in layers_x:
            layers_x[y] = [x]
        layers_x[y].append(x)
    if len(layers_x) <= len(search_graph.nodes):
        return pos
    
    avg_spacing = 0
    for y in layers_x:
        layer_spacings[y] = _get_avg_min_distance(layers_x[y])
        avg_spacing += layer_spacings[y]
    avg_spacing /= len(layers_x)
    layer_spacings_items = sorted(layer_spacings.items(), key=lambda x: x[0])
    for i in range(len(layer_spacings_items)):
        before = layer_spacings_items[i-1][1] if i != 0 else 0
        after = layer_spacings_items[i+1][1] if i != len(layer_spacings_items)-1 else before
        if before == 0:
            before = after
        y, spacing = layer_spacings_items[i]
        if spacing == 0 or (spacing < (before+after)/3):
            layer_spacings[y] = (before+after)/2
            if layer_spacings[y] == 0:
                layer_spacings[y] = avg_spacing
    current_y = 0
    layers_new_y = {}
    for old_y in sorted(layers_x.keys()):
        layers_new_y[old_y] = current_y
        current_y += layer_spacings[old_y]
    for n in pos:
        x, y = pos[n]
        pos[n] = (x, layers_new_y[y])
    
    return pos

def _prepare_search_graph_scatters(
        number_of_scatters: int,
        color_scale: List[Colorscale],
        invert_color_scale: List[bool],
        search_node_colorbar_title: List[Optional[str]],
        search_node_colorbar_spacing: float,
        use3d: bool,
        draw_stems: bool,
        draw_edges: bool,
        plotly_settings: Dict[str, Dict[str, Any]]) -> Tuple[List[Union[go.Scatter, go.Scatter3d]], List[Union[go.Scatter, go.Scatter3d]], Optional[go.Scatter3d]]: # nodes, edges, stems
    if not use3d:
        _copy_to_dict(plotly_settings['search_nodes'], dict(
            marker=dict(
                colorscale=color_scale[0],
                reversescale=invert_color_scale[0],
                colorbar=dict(
                    title=search_node_colorbar_title[0],
                )
            )
        ))
        ps = plotly_settings['search_nodes']
        if search_node_colorbar_title[0] is None:
            ps = deepcopy(ps)
            del ps['marker']['colorbar']
        node_scatter = go.Scatter(**ps)
        node_scatter.x = []
        node_scatter.y = []
        if draw_edges:
            edge_scatter = go.Scatter(**plotly_settings['search_edges'])
            edge_scatter.x = []
            edge_scatter.y = []
            edge_scatters = [edge_scatter]
        else:
            edge_scatters = []
        return [node_scatter], edge_scatters, None
    else:
        node_scatters = []
        n_colorbars = 0
        for i in range(number_of_scatters):
            _copy_to_dict(plotly_settings['search_nodes'], dict(
                marker=dict(
                    colorscale=color_scale[i],
                    reversescale=invert_color_scale[i],
                    colorbar=dict(
                        y=1-search_node_colorbar_spacing*n_colorbars,
                        title=search_node_colorbar_title[i],
                    )
                )
            ))
            ps = plotly_settings['search_nodes']
            if search_node_colorbar_title[i] is None:
                ps = deepcopy(ps)
                del ps['marker']['colorbar']
            else:
                n_colorbars += 1
            scatter = go.Scatter3d(**ps)
            scatter.x = []
            scatter.y = []
            scatter.z = []
            node_scatters.append(scatter)
        if draw_edges:
            edge_scatters = []
            for _ in range(number_of_scatters):
                scatter = go.Scatter3d(**plotly_settings['search_edges'])
                scatter.x = []
                scatter.y = []
                scatter.z = []
                edge_scatters.append(scatter)
        else:
            edge_scatters = []
        stem_scatter = None
        if draw_stems:
            stem_scatter = go.Scatter3d(**plotly_settings['search_node_stems'])
            stem_scatter.x = []
            stem_scatter.y = []
            stem_scatter.z = []
        return node_scatters, edge_scatters, stem_scatter

def _prepare_search_graph_scatter_data(
        number_of_scatters: int,
        search_graph: nx.Graph, 
        search_pos: Dict[int, Position], 
        use3d: bool,
        search_node_color: List[Union[str, Callable[[SearchNode], float]]], 
        # 'total_cost' | 'fixed_cost' | 'heuristic_cost' | 'lookahead_penalty' | static HTML color (e.g. 'blue' or '#0000FF') | 
        # function that takes a node parameter dict and returns a float (e.g. lambda node: node['total_cost'] + node['fixed_cost'])
        # node parameter dict: {"fixed_cost": float, "heuristic_cost": float, "lookahead_penalty": float, "is_valid_mapping": bool, 
        #                       "final": bool, "depth": int, "layout": Tuple[int, ...], "swaps": Tuple[Tuple[int, int], ...]}
        # or list of the above if 3d graph is used and multiple points per node are defined in search_node_height (lengths need to match)
        # if in that case no is list is provided all points per node will be the same color
        prioritize_search_node_color: List[bool],
        search_node_height: List[Callable[[SearchNode], float]], 
        color_valid_mapping: Optional[str], # static HTML color (e.g. 'blue' or '#0000FF')
        color_final_node: Optional[str], # static HTML color (e.g. 'blue' or '#0000FF')
        draw_stems: bool, # only applicable for 3D plots
        draw_edges: bool) -> Tuple[
                List[Optional[float]], 
                List[Optional[float]], 
                Tuple[List[Optional[float]], ...],
                List[float],
                List[float],
                Tuple[List[float], ...],
                List[Union[float, str]],
                List[Optional[float]], 
                List[Optional[float]], 
                List[Optional[float]], 
                float,
                float,
                float,
                float,
                float,
                float
            ]: #edge_x, edge_y, edge_z, node_x, node_y, node_z, node_color, stem_x, stem_y, stem_z, min_x, max_x, min_y, max_y, min_z, max_z
    edge_x = []
    edge_y = []
    edge_z = tuple([] for _ in range(number_of_scatters))
    node_x = []
    node_y = []
    node_z = tuple([] for _ in range(number_of_scatters))
    node_color = tuple([] for _ in range(number_of_scatters))
    stem_x = []
    stem_y = []
    stem_z = []
    min_x = None
    max_x = None
    min_y = None
    max_y = None
    min_z = None
    max_z = None
    
    for node in search_graph.nodes():
        node_params = search_graph.nodes[node]['data']
        nx, ny = search_pos[node]
        min_nz = 0
        max_nz = 0
        
        node_x.append(nx)
        node_y.append(ny)
        if use3d:
            for i in range(len(search_node_height)):
                nz = search_node_height[i](node_params)
                node_z[i].append(nz)
                min_nz = min(min_nz, nz)
                max_nz = max(max_nz, nz)
            if draw_stems:
                stem_x.append(nx)
                stem_y.append(ny)
                stem_z.append(min_nz)
                stem_x.append(nx)
                stem_y.append(ny)
                stem_z.append(max_nz)
                stem_x.append(None)
                stem_y.append(None)
                stem_z.append(None)
        if min_x is None:
            min_x = nx
            max_x = nx
            min_y = ny
            max_y = ny
            min_z = min_nz
            max_z = max_nz
        else:
            min_x = min(min_x, nx)
            max_x = max(max_x, nx)
            min_y = min(min_y, ny)
            max_y = max(max_y, ny)
            if use3d:
                min_z = min(min_z, min_nz)
                max_z = max(max_z, max_nz)
        
        ncolor = None
        if len(search_node_color) == 1:
            if callable(search_node_color[0]):
                ncolor = search_node_color[0](node_params)
            else:
                ncolor = search_node_color[0]
        for i in range(number_of_scatters):
            prio_color = prioritize_search_node_color[i] if len(prioritize_search_node_color) > 1 else prioritize_search_node_color[0]
            if not prio_color and color_final_node is not None and node_params.final:
                node_color[i].append(color_final_node)
            elif not prio_color and color_valid_mapping is not None and node_params.is_valid_mapping:
                node_color[i].append(color_valid_mapping)
            elif ncolor is not None:
                node_color[i].append(ncolor)
            elif callable(search_node_color[i]):
                node_color[i].append(search_node_color[i](node_params))
            else:
                node_color[i].append(search_node_color[i])
    
    if draw_edges:
        nodes_indices = dict()
        for i, n in enumerate(list(search_graph.nodes())):
            nodes_indices[n] = i
        for n0, n1 in search_graph.edges():
            n0_i = nodes_indices[n0]
            n1_i = nodes_indices[n1]
            edge_x.append(node_x[n0_i])
            edge_x.append(node_x[n1_i])
            edge_x.append(None)
            edge_y.append(node_y[n0_i])
            edge_y.append(node_y[n1_i])
            edge_y.append(None)
            if use3d:
                for i in range(len(node_z)):
                    edge_z[i].append(node_z[i][n0_i])
                    edge_z[i].append(node_z[i][n1_i])
                    edge_z[i].append(None)
    
    return edge_x, edge_y, edge_z, node_x, node_y, node_z, node_color, stem_x, stem_y, stem_z, min_x, max_x, min_y, max_y, min_z, max_z

def _draw_search_graph_nodes(
        scatters: List[Union[go.Scatter, go.Scatter3d]],
        x: List[float],
        y: List[float],
        z: Tuple[List[float], ...],
        color: List[Union[float, str]],
        use3d: bool):
    for i in range(len(scatters)):
        scatters[i].x = x
        scatters[i].y = y
        if use3d:
            scatters[i].z = z[i]
        scatters[i].marker.color = color[i] if len(color) > 1 else color[0]

def _draw_search_graph_stems(
        scatter: go.Scatter3d,
        x: List[float],
        y: List[float],
        z: List[float]):
    scatter.x = x
    scatter.y = y
    scatter.z = z

def _draw_search_graph_edges(
        scatters: List[Union[go.Scatter, go.Scatter3d]],
        x: List[float],
        y: List[float],
        z: Tuple[List[float], ...],
        use3d: bool):
    for i in range(len(scatters)):
        scatters[i].x = x
        scatters[i].y = y
        if use3d:
            scatters[i].z = z[i]

def _parse_arch_graph(file_path: str) -> nx.Graph:
    arch = None
    with open(file_path, "r") as file:
        arch = json.load(file)
    fidelity = None
    if 'fidelity' in arch:
        fidelity = arch['fidelity']

    edges = set()
    nqbits = 0
    for q0, q1 in arch['coupling_map']:
        edge = (q0, q1) if q0 < q1 else (q1, q0)
        if fidelity is not None:
            edge += (fidelity['swap_fidelity_costs'][q0][q1],)
        else:
            edge += (float(30),) # TODO: remove static value, replace by something more sensible
        edges.add(edge)
        nqbits = max(nqbits, q0, q1)
    nqbits += 1

    graph = nx.Graph()
    graph.add_nodes_from(list(range(nqbits)))
    graph.add_weighted_edges_from(edges)
    return graph

def _draw_architecture_edges(
        arch_graph: nx.Graph, 
        arch_pos: Dict[int, Position], 
        plotly_settings: Dict[str, Dict[str, Any]]) -> Tuple[go.Scatter, go.Scatter]:
    edge_x = []
    edge_y = []
    edge_label_x = []
    edge_label_y = []
    edge_label_text = []
    for n1, n2 in arch_graph.edges():
        x0, y0 = arch_pos[n1]
        x1, y1 = arch_pos[n2]
        mid_point = ((x0+x1)/2, (y0+y1)/2)
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
        edge_label_x.append(mid_point[0])
        edge_label_y.append(mid_point[1])
        edge_label_text.append("{:.3f}".format(arch_graph[n1][n2]['weight']))
    
    _copy_to_dict(plotly_settings['architecture_edges'], dict(
        x=edge_x, y=edge_y
    ))
    _copy_to_dict(plotly_settings['architecture_edge_labels'], dict(
        x=edge_label_x, y=edge_label_y, text=edge_label_text
    ))
    
    return (
        go.Scatter(**plotly_settings['architecture_edges']), 
        go.Scatter(**plotly_settings['architecture_edge_labels'])
    )

def _draw_architecture_nodes(
        scatter: go.Scatter,
        arch_pos: Dict[int, Position],
        considered_qubit_colors: Dict[int, str],
        initial_qubit_position: List[int],
        single_qubit_multiplicity: List[int],
        two_qubit_individual_multiplicity: List[int]):
    x = []
    y = []
    color = []
    text = []
    for log_qubit, phys_qubit in enumerate(initial_qubit_position):
        if phys_qubit == -1:
            continue
        if log_qubit not in considered_qubit_colors:
            continue
        nx, ny = arch_pos[phys_qubit]
        x.append(nx)
        y.append(ny)
        color.append(considered_qubit_colors[log_qubit])
        text.append(f'q{log_qubit}<br>1-gates: {single_qubit_multiplicity[log_qubit]}x<br>2-gates: {two_qubit_individual_multiplicity[log_qubit]}x')
    
    scatter.x = x
    scatter.y = y
    scatter.text = text
    scatter.marker.color = color

def _draw_swap_arrows(
        arch_pos: Dict[int, Position],
        initial_layout: List[int],
        swaps: Tuple[Tuple[int, int], ...],
        considered_qubit_colors: Dict[int, str],
        arrow_offset: float,
        arrow_spacing_x: float,
        arrow_spacing_y: float,
        shared_swaps: bool,
        plotly_settings: Dict[str, Dict[str, Any]]) -> List[go.layout.Annotation]:
    layout = initial_layout.copy()

    swap_arrow_props = {} # (q0, q1) -> [{color: str, straight: bool, color2: Optional[str]}, ...]    \\ q0 < q1, color2 is color at shaft side of arrow for a shared swap
    for sw in swaps:
        edge = (sw[0], sw[1]) if sw[0] < sw[1] else (sw[1], sw[0])
        if edge not in swap_arrow_props:
            swap_arrow_props[edge] = []
        props_list = swap_arrow_props[edge]
        
        sw0_considered = layout[sw[0]] in considered_qubit_colors
        sw1_considered = layout[sw[1]] in considered_qubit_colors
        if shared_swaps and sw0_considered and sw1_considered:
            props_list.append({
                'color': considered_qubit_colors[layout[sw[0]]],
                'straight': sw[0] < sw[1],
                'color2': considered_qubit_colors[layout[sw[1]]]
            })
        else:
            if sw0_considered:
                props_list.append({
                    'color': considered_qubit_colors[layout[sw[0]]],
                    'straight': sw[0] < sw[1],
                    'color2': None
                })
            if sw1_considered:
                props_list.append({
                    'color': considered_qubit_colors[layout[sw[1]]],
                    'straight': sw[0] >= sw[1],
                    'color2': None
                })
        layout[sw[0]], layout[sw[1]] = layout[sw[1]], layout[sw[0]]

    list_of_arch_arrows = []
    for edge, props in swap_arrow_props.items():
        x0, y0 = arch_pos[edge[0]]
        x1, y1 = arch_pos[edge[1]]
        v = (x1-x0, y1-y0)
        n = -v[1], v[0]
        norm = (v[0]**2 + v[1]**2)**0.5
        
        x0, y0 = x0+arrow_offset*v[0], y0+arrow_offset*v[1]
        x1, y1 = x1-arrow_offset*v[0], y1-arrow_offset*v[1]
        n = arrow_spacing_x*n[0]/norm, arrow_spacing_y*n[1]/norm
        
        n_offsets = (len(props)-1)/2
        for p in props:
            a_x0, a_y0, a_x1, a_y1 = (x0, y0, x1, y1) if p['straight'] else (x1, y1, x0, y0)
            a_x0, a_y0, a_x1, a_y1 = a_x0+n_offsets*n[0], a_y0+n_offsets*n[1], a_x1+n_offsets*n[0], a_y1+n_offsets*n[1]
            n_offsets -= 1
            if p['color2'] is not None:
                mx, my = (a_x0+a_x1)/2, (a_y0+a_y1)/2
                _copy_to_dict(plotly_settings['arrows'], dict(
                    x=a_x1, #arrow end
                    y=a_y1,
                    ax=mx, #arrow start
                    ay=my,
                    arrowcolor=p['color'],
                ))
                arrow = go.layout.Annotation(plotly_settings['arrows'])
                list_of_arch_arrows.append(arrow)
                _copy_to_dict(plotly_settings['arrows'], dict(
                    x=a_x0, #arrow end
                    y=a_y0,
                    ax=mx, #arrow start
                    ay=my,
                    arrowcolor=p['color2'],
                ))
                arrow = go.layout.Annotation(plotly_settings['arrows'])
                list_of_arch_arrows.append(arrow)
            else:
                _copy_to_dict(plotly_settings['arrows'], dict(
                    x=a_x1, #arrow end
                    y=a_y1,
                    ax=a_x0, #arrow start
                    ay=a_y0,
                    arrowcolor=p['color'],
                ))
                arrow = go.layout.Annotation(plotly_settings['arrows'])
                list_of_arch_arrows.append(arrow)
    return list_of_arch_arrows

def _visualize_layout(
        fig: go.Figure,
        search_node: SearchNode,
        arch_node_trace: go.Scatter,
        arch_node_positions: Dict[int, Position],
        initial_layout: List[int],
        considered_qubit_colors: Dict[int, str],
        swap_arrow_offset: float,
        arch_x_arrow_spacing: float,
        arch_y_arrow_spacing: float,
        show_shared_swaps: bool,
        layout_node_trace_index: int, # current_node_layout_visualized
        search_node_trace: Optional[go.Scatter],
        plotly_settings: Dict[str, Dict[str, Any]]):
    
    layout = search_node.layout
    swaps = search_node.swaps
    _copy_to_dict(plotly_settings['stats_legend'], dict(
        text=f'Node:   <b>{search_node.id}</b><br>'+
             f'Cost:     <b>{search_node.total_cost():.3f}</b> = {search_node.fixed_cost:.3f} + '+
                f'{search_node.heuristic_cost:.3f} + {search_node.lookahead_penalty:.3f} (fixed + heuristic + lookahead)<br>'+
             f'Depth:  {search_node.depth}<br>'+
             f'Valid / Final:  {"yes" if search_node.is_valid_mapping else "no"} / {"yes" if search_node.final else "no"}',
    ))
    stats = go.layout.Annotation(**plotly_settings['stats_legend'])
    annotations = []
    if len(swaps) > 0:
        annotations = _draw_swap_arrows(
            arch_node_positions, 
            initial_layout, 
            swaps, 
            considered_qubit_colors, 
            swap_arrow_offset, 
            arch_x_arrow_spacing, 
            arch_y_arrow_spacing,
            show_shared_swaps,
            plotly_settings
        )
    annotations.append(stats)
    
    arch_node_x = []
    arch_node_y = []
    for log_qubit, phys_qubit in enumerate(_reverse_layout(layout)):
        if phys_qubit == -1:
            continue
        if log_qubit not in considered_qubit_colors:
            continue
        x, y = arch_node_positions[phys_qubit]
        arch_node_x.append(x)
        arch_node_y.append(y)
    
    with fig.batch_update():
        arch_node_trace.x = arch_node_x
        arch_node_trace.y = arch_node_y
        fig.layout.annotations = annotations
        if search_node_trace is not None:
            marker_line_widths = [0]*len(search_node_trace.x)
            marker_line_widths[layout_node_trace_index] = 2
            search_node_trace.marker.line.width = marker_line_widths

def _load_layer_data(
        data_logging_path: str,
        layer: int,
        layout_method: Literal['walker', 'dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'], 
        layout_walker_factor: float, # 0-1
        tapered_layer_heights: bool,
        number_of_node_traces: int,
        use3d: bool, 
        node_color: List[Union[str, Callable[[SearchNode], float]]],
        prioritize_node_color: List[bool],
        node_height: List[Callable[[SearchNode], float]], 
        color_valid_mapping: Optional[str], 
        color_final_node: Optional[str],
        draw_stems: bool,
        draw_edges: bool) -> Tuple[
            nx.Graph, # search_graph
            List[int], # initial_layout
            List[int], # initial_qbit_positions
            Dict[int, str], # considered_qubit_colors
            List[int], # single_qbit_multiplicities
            List[TwoQbitMultiplicity], # two_qbit_multiplicities
            List[int], # individual_two_qbit_multiplicities
            int, # final_node_id
            List[float], # search_node_scatter_data_x
            List[float], # search_node_scatter_data_y
            List[List[float]], # search_node_scatter_data_z
            List[List[Union[float, str]]], # search_node_scatter_data_color
            List[float], # search_node_stem_scatter_data_x
            List[float], # search_node_stem_scatter_data_y
            List[float], # search_node_stem_scatter_data_z
            List[float], # search_edge_scatter_data_x
            List[float], # search_edge_scatter_data_y
            List[List[float]], # search_edge_scatter_data_z
            float, # search_min_x
            float, # search_max_x
            float, # search_min_y
            float, # search_max_y
            float, # search_min_z
            float, # search_max_z
        ]:
    
    if not os.path.exists(f"{data_logging_path}layer_{layer}.json"):
        raise FileNotFoundError(f"No data at {data_logging_path}layer_{layer}.json")
    if not os.path.exists(f"{data_logging_path}nodes_layer_{layer}.csv"):
        raise FileNotFoundError(f"No data at {data_logging_path}nodes_layer_{layer}.csv")
    
    circuit_layer = None
    with open(f"{data_logging_path}layer_{layer}.json", "r") as circuit_layer_file:
        circuit_layer = json.load(circuit_layer_file)

    single_q_mult = circuit_layer['single_qubit_multiplicity']
    two_q_mult_raw = circuit_layer['two_qubit_multiplicity']
    two_q_mult: List[TwoQbitMultiplicity] = []
    for mult in two_q_mult_raw:
        two_q_mult.append(TwoQbitMultiplicity(mult['q1'], mult['q2'], mult['backward'], mult['forward']))
    initial_layout = circuit_layer['initial_layout']
    initial_positions = _reverse_layout(initial_layout)
    final_node_id = circuit_layer['final_node_id']

    graph, graph_root = _parse_search_graph(f"{data_logging_path}nodes_layer_{layer}.csv", final_node_id)
    
    pos = _layout_search_graph(
        graph, 
        graph_root, 
        layout_method, 
        layout_walker_factor,
        tapered_layer_heights
    )
    
    (
        edge_x, edge_y, edge_z, 
        node_x, node_y, node_z, node_color_data, 
        stem_x, stem_y, stem_z, 
        min_x, max_x, min_y, max_y, min_z, max_z
    ) = _prepare_search_graph_scatter_data(
        number_of_node_traces, 
        graph,
        pos, 
        use3d, 
        node_color,
        prioritize_node_color,
        node_height, 
        color_valid_mapping, 
        color_final_node,
        draw_stems,
        draw_edges
    )
    
    considered_qubit_colors = {}
    considered_qubit_ngroups = 0
    two_q_mult_individual = [0]*len(single_q_mult)
    for mult in two_q_mult:
        two_q_mult_individual[mult.q0] += mult.backward + mult.forward
        two_q_mult_individual[mult.q1] += mult.backward + mult.forward
        considered_qubit_colors[mult.q0] = considered_qubit_ngroups
        considered_qubit_colors[mult.q1] = considered_qubit_ngroups
        considered_qubit_ngroups += 1
    for i, q in enumerate(single_q_mult):
        if q != 0 and i not in considered_qubit_colors:
            considered_qubit_colors[i] = considered_qubit_ngroups
            considered_qubit_ngroups += 1
    considered_qubits_color_codes = [distinctipy.get_hex(c) for c in distinctipy.get_colors(max(10, considered_qubit_ngroups))]
    shuffle(considered_qubits_color_codes)
    for q in considered_qubit_colors:
        considered_qubit_colors[q] = considered_qubits_color_codes[considered_qubit_colors[q]]
    
    return (
        graph,
        initial_layout,
        initial_positions,
        considered_qubit_colors,
        single_q_mult,
        two_q_mult,
        two_q_mult_individual,
        final_node_id,
        node_x,
        node_y,
        node_z,
        node_color_data,
        stem_x,
        stem_y,
        stem_z,
        edge_x,
        edge_y,
        edge_z,
        min_x,
        max_x,
        min_y,
        max_y,
        min_z,
        max_z
    )

_total_cost_lambda: Callable[[SearchNode], float] = lambda n: n.total_cost()
_total_fixed_cost_lambda: Callable[[SearchNode], float] = lambda n: n.total_fixed_cost()
_fixed_cost_lambda: Callable[[SearchNode], float] = lambda n: n.fixed_cost
_heuristic_cost_lambda: Callable[[SearchNode], float] = lambda n: n.heuristic_cost
_lookahead_penalty_lambda: Callable[[SearchNode], float] = lambda n: n.lookahead_penalty
def _cost_string_to_lambda(cost_string: str) -> Optional[Callable[[SearchNode], float]]:
    if cost_string == 'total_cost':
        return _total_cost_lambda
    elif cost_string == 'total_fixed_cost':
        return _total_fixed_cost_lambda
    elif cost_string == 'fixed_cost':
        return _fixed_cost_lambda
    elif cost_string == 'heuristic_cost':
        return _heuristic_cost_lambda
    elif cost_string == 'lookahead_penalty':
        return _lookahead_penalty_lambda
    else:
        return None

default_plotly_settings = dict(
    layout=dict(
        autosize=False,
        showlegend=False,
        hovermode='closest',
        coloraxis_colorbar_x=-0.15,
    ),
    
    arrows=dict(
        text="",
        showarrow=True,
        arrowhead=5,
        arrowwidth=2,
    ),
    
    stats_legend=dict(
        align='left',
        showarrow=False,
        xref='paper',
        yref='paper',
        x=1,
        y=1.175,
        bordercolor='black',
        borderwidth=1,
    ),
    
    search_nodes=dict(
        mode='markers',
        hoverinfo='none',
        marker=dict(
            showscale=True,
            size=10,
            colorbar=dict(
                thickness=15,
                orientation='h',
                lenmode='fraction',
                len=0.5,
                xref='paper',
                x=0.21,
                yref='container',
            )
        )
    ),
    
    search_node_stems=dict(
        hoverinfo='none',
        mode='lines'
    ),
    
    search_edges=dict(
        hoverinfo='none',
        mode='lines'
    ),
    
    architecture_nodes=dict(
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
        )
    ),
    
    architecture_edges=dict(
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'
    ),
    
    architecture_edge_labels=dict(
        hoverinfo='none',
        mode='text'
    ),
    
    search_xaxis=dict(
        showgrid=False, 
        zeroline=False, 
        showticklabels=False,
        title_text=''
    ),
    
    search_yaxis=dict(
        showgrid=False, 
        zeroline=False, 
        showticklabels=False,
        title_text=''
    ),
    search_zaxis=dict(
        title_text=''
    ),
    architecture_xaxis=dict(
        showgrid=False, 
        zeroline=False, 
        showticklabels=False,
        title_text=''
    ),
    architecture_yaxis=dict(
        showgrid=False, 
        zeroline=False, 
        showticklabels=False,
        title_text=''
    )
)

def _visualize_search_graph_check_parameters(
        data_logging_path: str, 
        layer: Union[int, Literal['interactive']],
        architecture_node_positions: Optional[Dict[int, Position]],
        architecture_layout_method: Literal['dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'],
        search_node_layout_method: Literal['walker', 'dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'],
        search_node_layout_walker_factor: float,
        search_graph_border: float,
        architecture_border: float,
        swap_arrow_spacing: float,
        swap_arrow_offset: float,
        use3d: bool,
        projection: Literal['orthographic', 'perspective'],
        width: int,
        height: int,
        draw_search_edges: bool,
        search_edges_width: float,
        search_edges_color: str,
        search_edges_dash: str,
        tapered_search_layer_heights: bool,
        show_layout: Optional[Literal['hover', 'click']],
        show_swaps: bool,
        show_shared_swaps: bool,
        color_valid_mapping: Optional[str],
        color_final_node: Optional[str],
        search_node_color: Union[str, Callable[[SearchNode], float], List[Union[str, Callable[[SearchNode], float]]]],
        prioritize_search_node_color: Union[bool, List[bool]],
        search_node_color_scale: Union[Colorscale, List[Colorscale]],
        search_node_invert_color_scale: Union[bool, List[bool]],
        search_node_colorbar_title: Optional[Union[str, List[Optional[str]]]],
        search_node_colorbar_spacing: float,
        search_node_height: Union[str, Callable[[SearchNode], float], List[Union[str, Callable[[SearchNode], float]]]], 
        draw_stems: bool,
        stems_width: float,
        stems_color: str,
        stems_dash: str,
        show_search_progression: bool,
        search_progression_step: int,
        search_progression_speed: float,
        plotly_settings: Dict[str, Dict[str, Any]]
            ) -> Tuple[
                str, # data_logging_path
                bool, # hide_layout
                bool, # draw_stems
                int, # number_of_node_traces
                List[Callable[[SearchNode], float]], # search_node_height
                List[Union[str, Callable[[SearchNode], float]]], # search_node_color
                List[Colorscale], # search_node_color_scale
                List[bool], # search_node_invert_color_scale
                List[bool], # prioritize_search_node_color
                List[str], # search_node_colorbar_title
                Dict[str, Dict[str, Any]] # plotly_settings
            ]:
    
    if type(data_logging_path) is not str:
        raise ValueError('data_logging_path must be a string')
    if data_logging_path[-1] != "/":
        data_logging_path += "/"
    if not os.path.exists(data_logging_path):
        raise FileNotFoundError(f"Path {data_logging_path} does not exist.")
    
    if type(layer) is not int and layer != 'interactive':
        raise ValueError('layer must be an integer or string literal "interactive"')
    
    if architecture_node_positions is not None:
        if type(architecture_node_positions) is not dict:
            raise ValueError('architecture_node_positions must be a dict of the form {qubit_index: (x: float, y: float)}')
        for i in architecture_node_positions:
            if type(i) is not int:
                raise ValueError('architecture_node_positions must be a dict of the form {qubit_index: (x: float, y: float)}')
            if type(architecture_node_positions[i]) is not tuple:
                raise ValueError('architecture_node_positions must be a dict of the form {qubit_index: (x: float, y: float)}')
            if len(architecture_node_positions[i]) != 2:
                raise ValueError('architecture_node_positions must be a dict of the form {qubit_index: (x: float, y: float)}')
            if not _is_number(architecture_node_positions[i][0]) or not _is_number(architecture_node_positions[i][1]):
                raise ValueError('architecture_node_positions must be a dict of the form {qubit_index: (x: float, y: float)}')
    
    if architecture_layout_method not in ['dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork']:
        raise ValueError('architecture_layout_method must be one of "dot", "neato", "fdp", "sfdp", "circo", "twopi", "osage", "patchwork"')
    
    if search_node_layout_method not in ['walker', 'dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork']:
        raise ValueError('search_node_layout_method must be one of "walker", "dot", "neato", "fdp", "sfdp", "circo", "twopi", "osage", "patchwork"')
    
    if search_node_layout_method == 'walker' and (not _is_number(search_node_layout_walker_factor) or search_node_layout_walker_factor < 0):
        raise ValueError('search_node_layout_walker_factor must be a non-negative float')
    
    if not _is_number(search_graph_border) or search_graph_border < 0:
        raise ValueError('search_graph_border must be a non-negative float')
    
    if not _is_number(architecture_border) or architecture_border < 0:
        raise ValueError('architecture_border must be a non-negative float')
    
    if not _is_number(swap_arrow_spacing) or swap_arrow_spacing < 0:
        raise ValueError('swap_arrow_spacing must be a non-negative float')
    
    if not _is_number(swap_arrow_offset) or swap_arrow_offset < 0 or swap_arrow_offset >= 0.5:
        raise ValueError('swap_arrow_offset must be a float between 0 and 0.5')
    
    if type(use3d) is not bool:
        raise ValueError('use3d must be a boolean')
    
    if projection != 'orthographic' and projection != 'perspective':
        raise ValueError('projection must be either "orthographic" or "perspective"')
    
    if type(width) is not int or width < 1:
        raise ValueError('width must be a positive integer')
    
    if type(height) is not int or height < 1:
        raise ValueError('height must be a positive integer')
    
    if type(draw_search_edges) is not bool:
        raise ValueError('draw_search_edges must be a boolean')
    
    if type(search_edges_width) is not float or search_edges_width <= 0:
        raise ValueError('search_edges_width must be a positive float')
    
    if ColorValidator.perform_validate_coerce(search_edges_color, allow_number=False) is None:
        raise ValueError(ColorValidator('search_edges_color', 'visualize_search_graph').description())
    
    if search_edges_dash not in ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'] and not re.match(r'^(\d+(px|%)?(\s*\d+(px|%)?)*)$', search_edges_dash):
        raise ValueError('search_edges_dash must be one of "solid", "dot", "dash", "longdash", "dashdot", "longdashdot" or a string containing a dash length list in '+
                         'pixels or percentages (e.g. "5px 10px 2px 2px", "5, 10, 2, 2", "10\% 20\% 40\%")')
    
    if type(tapered_search_layer_heights) is not bool:
        raise ValueError('tapered_search_layer_heights must be a boolean')
    
    if show_layout not in ['hover', 'click'] and show_layout is not None:
        raise ValueError('show_layout must be one of "hover", "click" or None')
    hide_layout = (show_layout is None)
    
    if type(show_swaps) is not bool:
        raise ValueError('show_swaps must be a boolean')
    
    if type(show_shared_swaps) is not bool:
        raise ValueError('show_shared_swaps must be a boolean')
    
    if color_valid_mapping is not None and ColorValidator.perform_validate_coerce(color_valid_mapping, allow_number=False) is None:
        raise ValueError('color_valid_mapping must be None or a color specified as:\n' + 
                         _remove_first_lines(ColorValidator('color_valid_mapping', 'visualize_search_graph').description(), 1))
    
    if color_final_node is not None and ColorValidator.perform_validate_coerce(color_final_node, allow_number=False) is None:
        raise ValueError('color_final_node must be None or a color specified as:\n' + 
                         _remove_first_lines(ColorValidator('color_final_node', 'visualize_search_graph').description(), 1))
    
    if type(draw_stems) is not bool:
        raise ValueError('draw_stems must be a boolean')
    
    if type(stems_width) is not float or stems_width <= 0:
        raise ValueError('stems_width must be a positive float')
    
    if ColorValidator.perform_validate_coerce(stems_color, allow_number=False) is None:
        raise ValueError(ColorValidator('stems_color', 'visualize_search_graph').description())
    
    if stems_dash not in ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot'] and not re.match(r'^(\d+(px|%)?(\s*\d+(px|%)?)*)$', stems_dash):
        raise ValueError('stems_dash must be one of "solid", "dot", "dash", "longdash", "dashdot", "longdashdot" or a string containing a dash length list in '+
                         'pixels or percentages (e.g. "5px 10px 2px 2px", "5, 10, 2, 2", "10\% 20\% 40\%")')
    
    if type(show_search_progression) is not bool:
        raise ValueError('show_search_progression must be a boolean')
    
    if type(search_progression_step) is not int or search_progression_step < 1:
        raise ValueError('search_porgression_step must be a positive integer')
    
    if not _is_number(search_progression_speed) or search_progression_speed <= 0:
        raise ValueError('search_progression_speed must be a positive float')
    
    if type(plotly_settings) is not dict and any([(
            type(plotly_settings[key]) is not dict or 
            key not in [
                'layout', 'arrows', 'stats_legend', 'search_nodes', 'search_edges', 'architecture_nodes', 'architecture_edges',
                'architecture_edge_labels', 'search_xaxis', 'search_yaxis', 'search_zaxis', 'architecture_xaxis', 'architecture_yaxis'
            ]) for key in plotly_settings.keys()]):
        raise ValueError("plotly_settings must be a dict with any of these entries:" + 
            "{\n" + 
            "    'layout': settings for plotly.graph_objects.Layout (of subplots figure)\n" + 
            "    'arrows': settings for plotly.graph_objects.layout.Annotation\n" + 
            "    'stats_legend': settings for plotly.graph_objects.layout.Annotation\n" + 
            "    'search_nodes': settings for plotly.graph_objects.Scatter resp. ...Scatter3d\n" + 
            "    'search_edges': settings for plotly.graph_objects.Scatter resp. ...Scatter3d\n" + 
            "    'architecture_nodes': settings for plotly.graph_objects.Scatter\n" + 
            "    'architecture_edges': settings for plotly.graph_objects.Scatter\n" + 
            "    'architecture_edge_labels': settings for plotly.graph_objects.Scatter\n" + 
            "    'search_xaxis': settings for plotly.graph_objects.layout.XAxis resp. ...layout.scene.XAxis\n" + 
            "    'search_yaxis': settings for plotly.graph_objects.layout.YAxis resp. ...layout.scene.YAxis\n" + 
            "    'search_zaxis': settings for plotly.graph_objects.layout.scene.ZAxis\n" + 
            "    'architecture_xaxis': settings for plotly.graph_objects.layout.XAxis\n" + 
            "    'architecture_yaxis': settings for plotly.graph_objects.layout.YAxis\n" + 
            "}")
    
    plotly_set = deepcopy(default_plotly_settings)
    plotly_set['layout']['width'] = width
    plotly_set['layout']['height'] = height
    plotly_set['search_node_stems']['line'] = dict(width=stems_width, color=stems_color, dash=stems_dash)
    plotly_set['search_edges']['line'] = dict(width=search_edges_width, color=search_edges_color, dash=search_edges_dash)
    if use3d:
        plotly_set['layout']['scene'] = dict(
            camera=dict(
                projection=dict(
                    type=projection
                )
            )
        )
        plotly_set['arrows']['xref'] = 'x'
        plotly_set['arrows']['yref'] = 'y'
        plotly_set['arrows']['axref'] = 'x'
        plotly_set['arrows']['ayref'] = 'y'
    else:
        draw_stems = False
        plotly_set['arrows']['xref'] = 'x2'
        plotly_set['arrows']['yref'] = 'y2'
        plotly_set['arrows']['axref'] = 'x2'
        plotly_set['arrows']['ayref'] = 'y2'
    _copy_to_dict(plotly_set, plotly_settings)
    
    number_of_node_traces = len(search_node_height) if use3d and _is_len_iterable(search_node_height) else 1
    
    if not _is_number(search_node_colorbar_spacing) or search_node_colorbar_spacing <= 0 or search_node_colorbar_spacing >= 1:
        raise ValueError('search_node_colorbar_spacing must be a float between 0 and 1')
    
    if search_node_colorbar_title is None:
        search_node_colorbar_title = [None] * number_of_node_traces
    if _is_len_iterable(search_node_colorbar_title):
        if len(search_node_colorbar_title) > 1 and not use3d:
            raise ValueError("search_node_colorbar_title can only be a list in a 3D plot.")
        if len(search_node_colorbar_title) != number_of_node_traces:
            raise ValueError(f"Length of search_node_colorbar_title ({len(search_node_colorbar_title)}) does not match length of search_node_height ({number_of_node_traces}).")
        untitled = 1
        for i, title in enumerate(search_node_colorbar_title):
            if title is None:
                color = None
                if _is_len_iterable(search_node_color):
                    if len(search_node_color) == len(search_node_colorbar_title):
                        color = search_node_color[i]
                    else:
                        color = search_node_color[0]
                else:
                    color = search_node_color
                if type(color) is not str:
                    search_node_colorbar_title[i] = f'Untitled{untitled}'
                    untitled += 1
                elif color == 'total_cost':
                    search_node_colorbar_title[i] = 'Total cost'
                elif color == 'total_fixed_cost':
                    search_node_colorbar_title[i] = 'Total fixed cost'
                elif color == 'fixed_cost':
                    search_node_colorbar_title[i] = 'Fixed cost'
                elif color == 'heuristic_cost':
                    search_node_colorbar_title[i] = 'Heuristic cost'
                elif color == 'lookahead_penalty':
                    search_node_colorbar_title[i] = 'Lookahead penalty'
                else:
                    search_node_colorbar_title[i] = color
            elif type(title) is not str:
                raise ValueError(f"search_node_colorbar_title must be None, a string, or list of strings and None.")
    elif type(search_node_colorbar_title) is str:
        search_node_colorbar_title = [search_node_colorbar_title] * number_of_node_traces
    else:
        raise ValueError(f"search_node_colorbar_title must be None, a string, or list of strings and None.")
    
    if _is_len_iterable(search_node_color_scale):
        if len(search_node_color_scale) > 1 and not use3d:
            raise ValueError("search_node_color_scale can only be a list in a 3D plot.")
        if len(search_node_color_scale) != number_of_node_traces:
            raise ValueError(f"Length of search_node_color_scale ({len(search_node_color_scale)}) does not match length of search_node_height ({number_of_node_traces}).")
        cs_validator = ColorscaleValidator('search_node_color_scale', 'visualize_search_graph')
        for cs in search_node_color_scale:
            try:
                cs_validator.validate_coerce(cs)
            except:
                raise ValueError('search_node_color_scale must be a list of colorscales or a colorscale, specified as:\n' + 
                                _remove_first_lines(cs_validator.description(), 1))
    else:
        cs_validator = ColorscaleValidator('search_node_color_scale', 'visualize_search_graph')
        try:
            cs_validator.validate_coerce(search_node_color_scale)
        except:
            raise ValueError('search_node_color_scale must be a list of colorscales or a colorscale, specified as:\n' + 
                             _remove_first_lines(cs_validator.description(), 1))
        search_node_color_scale = [search_node_color_scale] * number_of_node_traces
    
    if _is_len_iterable(search_node_invert_color_scale):
        if len(search_node_invert_color_scale) > 1 and not use3d:
            raise ValueError("search_node_invert_color_scale can only be a list in a 3D plot.")
        if len(search_node_invert_color_scale) != number_of_node_traces:
            raise ValueError(f"Length of search_node_invert_color_scale ({len(search_node_invert_color_scale)}) does not match length of search_node_height ({number_of_node_traces}).")
        for i, invert in enumerate(search_node_invert_color_scale):
            if type(invert) is not bool:
                raise ValueError(f"search_node_invert_color_scale must be a boolean or list of booleans.")
    elif type(search_node_invert_color_scale) is not bool:
        raise ValueError(f"search_node_invert_color_scale must be a boolean or list of booleans.")
    else:
        search_node_invert_color_scale = [search_node_invert_color_scale] * number_of_node_traces
    
    if _is_len_iterable(prioritize_search_node_color):
        if len(prioritize_search_node_color) > 1 and not use3d:
            raise ValueError("prioritize_search_node_color can only be a list in a 3D plot.")
        if len(prioritize_search_node_color) != number_of_node_traces:
            raise ValueError(f"Length of prioritize_search_node_color ({len(prioritize_search_node_color)}) does not match length of search_node_height ({number_of_node_traces}).")
        for i, invert in enumerate(prioritize_search_node_color):
            if type(invert) is not bool:
                raise ValueError(f"prioritize_search_node_color must be a boolean or list of booleans.")
    elif type(prioritize_search_node_color) is not bool:
        raise ValueError(f"prioritize_search_node_color must be a boolean or list of booleans.")
    else:
        prioritize_search_node_color = [prioritize_search_node_color] * number_of_node_traces
    
    if _is_len_iterable(search_node_color):
        if len(search_node_color) > 1 and not use3d:
            raise ValueError("search_node_color can only be a list in a 3D plot.")
        if len(search_node_color) != number_of_node_traces:
            raise ValueError(f"Length of search_node_color ({len(search_node_color)}) does not match length of search_node_height ({number_of_node_traces}).")
        for i, c in enumerate(search_node_color):
            if isinstance(c, str):
                l = _cost_string_to_lambda(c)
                if l is not None:
                    search_node_color[i] = l
                elif ColorValidator.perform_validate_coerce(c, allow_number=False) is None:
                    raise ValueError(f'search_node_color[{i}] is neither a valid cost function preset ("total_cost", "total_fixed_cost", "fixed_cost", '+
                                     '"heuristic_cost", "lookahead_penalty") nor a respective callable, nor a valid color string specified as:\n' + 
                                     _remove_first_lines(ColorValidator('search_node_color', 'visualize_search_graph').description(), 1))
                else: # static color
                    search_node_colorbar_title[i] = None
    elif isinstance(search_node_color, str):
        l = _cost_string_to_lambda(search_node_color)
        if l is not None:
            search_node_color = [l]
        elif ColorValidator.perform_validate_coerce(search_node_color, allow_number=False) is None:
            raise ValueError(f'search_node_color is neither a list nor a valid cost function preset ("total_cost", "total_fixed_cost", "fixed_cost", '+
                              '"heuristic_cost", "lookahead_penalty") nor a respective callable, nor a valid color string specified as:\n' + 
                              _remove_first_lines(ColorValidator('search_node_color', 'visualize_search_graph').description(), 1))
        else: # static color
            search_node_color = [search_node_color]
            search_node_colorbar_title = [None]*number_of_node_traces
    else:
        raise ValueError('search_node_color must be a cost function preset ("total_cost", "total_fixed_cost", "fixed_cost", "heuristic_cost", '+
                         '"lookahead_penalty") or a respective callable, or a valid color string, or a list of the above.')
    
    if use3d:
        if _is_len_iterable(search_node_height):
            lambdas = set()
            for i, c in enumerate(search_node_height):
                if isinstance(c, str):
                    l = _cost_string_to_lambda(c)
                    if l is None:
                        raise ValueError(f"Unkown cost function preset search_node_height[{i}]: {c}")
                    elif l in lambdas:
                        raise ValueError(f"search_node_height must not contain the same cost function multiple times: {c}")
                    else:
                        search_node_height[i] = l
                        lambdas.add(l)
                elif not callable(c):
                    raise ValueError(f"search_node_height must be a cost function preset ('total_cost', 'total_fixed_cost', 'fixed_cost', 'heuristic_cost', "+
                                     "'lookahead_penalty') or a respective callable, or a list of the above.")
        elif isinstance(search_node_height, str):
            l = _cost_string_to_lambda(search_node_height)
            if l is not None:
                search_node_height = [l]
            else:
                raise ValueError(f"Unkown cost function preset search_node_height: {search_node_height}")
        else:
            raise ValueError("search_node_height must be a list of cost functions or a single cost function.")
    else:
        search_node_height = []
        
    return (data_logging_path, hide_layout, draw_stems, number_of_node_traces, search_node_height, search_node_color, search_node_color_scale, 
        search_node_invert_color_scale, prioritize_search_node_color, search_node_colorbar_title, plotly_set)




def visualize_search_graph(
        data_logging_path: str, 
        layer: Union[int, Literal['interactive']] = 'interactive', # 'interactive' (slider menu) | index
        architecture_node_positions: Optional[Dict[int, Position]] = None,
        architecture_layout_method: Literal['dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'] = 'sfdp',
        search_node_layout_method: Literal['walker', 'dot', 'neato', 'fdp', 'sfdp', 'circo', 'twopi', 'osage', 'patchwork'] = 'walker',
        search_node_layout_walker_factor: float = 0.6,
        search_graph_border: float = 0.05,
        architecture_border: float = 0.05,
        swap_arrow_spacing: float = 0.05,
        swap_arrow_offset: float = 0.05,
        use3d: bool = True,
        projection: Literal['orthographic', 'perspective'] = 'perspective',
        width: int = 1400,
        height: int = 700,
        draw_search_edges: bool = True,
        search_edges_width: float = 0.5,
        search_edges_color: str = '#888',
        search_edges_dash: str = 'solid', # 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot', string containing a dash length list in pixels or percentages (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%')
        tapered_search_layer_heights: bool = True,
        show_layout: Optional[Literal['hover', 'click']] = 'hover',
        show_swaps: bool = True,
        show_shared_swaps: bool = True, # if one swap moves 2 considered qubits -> combine into one two-colored and two-headed arrow
        color_valid_mapping: Optional[str] = 'green', # static HTML color (e.g. 'blue' or '#0000FF')
        color_final_node: Optional[str] = 'red', # static HTML color (e.g. 'blue' or '#0000FF')
        search_node_color: Union[str, Callable[[SearchNode], float], List[Union[str, Callable[[SearchNode], float]]]] = 'total_cost',
        # 'total_cost' | 'fixed_cost' | 'heuristic_cost' | 'lookahead_penalty' | static HTML color (e.g. 'blue' or '#0000FF') | 
        # function that takes a SearchNode and returns a float (e.g. lambda n: n.fixed_cost + n.heuristic_cost)
        # node fields: {"fixed_cost": float, "heuristic_cost": float, "lookahead_penalty": float, "is_valid_mapping": bool, 
        #                       "final": bool, "depth": int, "layout": Tuple[int, ...], "swaps": Tuple[Tuple[int, int], ...]}
        # or list of the above if 3d graph is used and multiple points per node are defined in search_node_height (lengths need to match)
        # if in that case no is list is provided all points per node will be the same color
        prioritize_search_node_color: Union[bool, List[bool]] = False, # if True, search_node_color will be prioritized over color_valid_mapping and color_final_node
        search_node_color_scale: Union[str, List[str]] = 'YlGnBu', # https://plotly.com/python/builtin-colorscales/
        # aggrnyl     agsunset    blackbody   bluered     blues       blugrn      bluyl       brwnyl
        # bugn        bupu        burg        burgyl      cividis     darkmint    electric    emrld
        # gnbu        greens      greys       hot         inferno     jet         magenta     magma
        # mint        orrd        oranges     oryel       peach       pinkyl      plasma      plotly3
        # pubu        pubugn      purd        purp        purples     purpor      rainbow     rdbu
        # rdpu        redor       reds        sunset      sunsetdark  teal        tealgrn     turbo
        # viridis     ylgn        ylgnbu      ylorbr      ylorrd      algae       amp         deep
        # dense       gray        haline      ice         matter      solar       speed       tempo
        # thermal     turbid      armyrose    brbg        earth       fall        geyser      prgn
        # piyg        picnic      portland    puor        rdgy        rdylbu      rdylgn      spectral
        # tealrose    temps       tropic      balance     curl        delta       oxy         edge
        # hsv         icefire     phase       twilight    mrybm       mygbm       armylg      falllg
        search_node_invert_color_scale: Union[bool, List[bool]] = True,
        search_node_colorbar_title: Optional[Union[str, List[Optional[str]]]] = None,
        search_node_colorbar_spacing: float = 0.06,
        search_node_height: Union[str, Callable[[SearchNode], float], List[Union[str, Callable[[SearchNode], float]]]] = 'total_cost', 
        # just as with search_node_color (without color strings), but possible to specify a list to draw multiple point per node on different heights
        # only applicable if use3d is True
        draw_stems: bool = False,
        stems_width: float = 0.7,
        stems_color: str = '#444',
        stems_dash: str = 'solid', # 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot', string containing a dash length list in pixels or percentages (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%')
        show_search_progression: bool = True,
        search_progression_step: int = 10,
        search_progression_speed: float = 2, # steps per second
        plotly_settings: Dict[str, Dict[str, Any]] = dict(),
        # {
        #   'layout': settings for plotly.graph_objects.Layout (of subplots figure)
        #   'arrows': settings for plotly.graph_objects.layout.Annotation
        #   'stats_legend': settings for plotly.graph_objects.layout.Annotation
        #   'search_nodes': settings for plotly.graph_objects.Scatter resp. ...Scatter3d
        #   'search_edges': settings for plotly.graph_objects.Scatter resp. ...Scatter3d
        #   'architecture_nodes': settings for plotly.graph_objects.Scatter
        #   'architecture_edges': settings for plotly.graph_objects.Scatter
        #   'architecture_edge_labels': settings for plotly.graph_objects.Scatter
        #   'search_xaxis': settings for plotly.graph_objects.layout.XAxis resp. ...layout.scene.XAxis
        #   'search_yaxis': settings for plotly.graph_objects.layout.YAxis resp. ...layout.scene.YAxis
        #   'search_zaxis': settings for plotly.graph_objects.layout.scene.ZAxis
        #   'architecture_xaxis': settings for plotly.graph_objects.layout.XAxis
        #   'architecture_yaxis': settings for plotly.graph_objects.layout.YAxis
        # }
        
        # TODO: show archticture edge labels (and control text?)
        # TODO: control hover text of search (especially for multiple points per node!) and architecture nodes?
    ) -> Widget:
    
    # check and process all parameters
    (
        data_logging_path, hide_layout, draw_stems, number_of_node_traces, 
        search_node_height, search_node_color, search_node_color_scale, 
        search_node_invert_color_scale, prioritize_search_node_color, 
        search_node_colorbar_title, plotly_settings
    ) = _visualize_search_graph_check_parameters(
        data_logging_path, 
        layer,
        architecture_node_positions,
        architecture_layout_method,
        search_node_layout_method,
        search_node_layout_walker_factor,
        search_graph_border,
        architecture_border,
        swap_arrow_spacing,
        swap_arrow_offset,
        use3d,
        projection,
        width,
        height,
        draw_search_edges,
        search_edges_width,
        search_edges_color,
        search_edges_dash,
        tapered_search_layer_heights,
        show_layout,
        show_swaps,
        show_shared_swaps,
        color_valid_mapping,
        color_final_node,
        search_node_color,
        prioritize_search_node_color,
        search_node_color_scale,
        search_node_invert_color_scale,
        search_node_colorbar_title,
        search_node_colorbar_spacing,
        search_node_height, 
        draw_stems,
        stems_width,
        stems_color,
        stems_dash,
        show_search_progression,
        search_progression_step,
        search_progression_speed,
        plotly_settings
    )
    
    # function-wide variables
    search_graph: Optional[nx.Graph] = None
    arch_graph: Optional[nx.Graph] = None
    initial_layout: List[int] = []
    considered_qubit_colors: Dict[int, str] = dict()
    current_node_layout_visualized: Optional[int] = None
    current_layer: Optional[int] = None
    
    arch_x_spacing: Optional[float] = None
    arch_y_spacing: Optional[float] = None
    arch_x_arrow_spacing: Optional[float] = None
    arch_y_arrow_spacing: Optional[float] = None
    arch_x_min: Optional[float] = None
    arch_x_max: Optional[float] = None
    arch_y_min: Optional[float] = None
    arch_y_max: Optional[float] = None
    
    search_node_traces: List[Union[go.Scatter, go.Scatter3d]] = []
    search_edge_traces: List[Union[go.Scatter, go.Scatter3d]] = []
    search_node_stem_trace: Optional[go.Scatter3d] = None
    arch_node_trace: Optional[go.Scatter] = None
    arch_edge_trace: Optional[go.Scatter] = None
    arch_edge_label_trace: Optional[go.Scatter] = None
    
    # one possible entry per layer (if not present, layer was not loaded yet)
    search_graphs: Dict[int, nx.Graph] = dict()
    initial_layouts:  Dict[int, List[int]] = dict()
    initial_qbit_positions:  Dict[int, List[int]] = dict() # reverses of initial_layouts
    layers_considered_qubit_colors:  Dict[int, Dict[int, str]] = dict()
    single_qbit_multiplicities: Dict[int, List[int]] = dict()
    individual_two_qbit_multiplicities: Dict[int, List[int]] = dict()
    two_qbit_multiplicities: Dict[int, List[TwoQbitMultiplicity]] = dict()
    final_node_ids: Dict[int, int] = dict()
    search_node_scatter_data_x: Dict[int, List[float]] = dict()
    search_node_scatter_data_y: Dict[int, List[float]] = dict()
    search_node_scatter_data_z: Dict[int, List[List[float]]] = dict()
    search_node_scatter_data_color: Dict[int, List[List[Union[float, str]]]] = dict()
    search_node_stem_scatter_data_x: Dict[int, List[float]] = dict()
    search_node_stem_scatter_data_y: Dict[int, List[float]] = dict()
    search_node_stem_scatter_data_z: Dict[int, List[float]] = dict()
    search_edge_scatter_data_x: Dict[int, List[float]] = dict()
    search_edge_scatter_data_y: Dict[int, List[float]] = dict()
    search_edge_scatter_data_z: Dict[int, List[List[float]]] = dict()
    search_min_x: Dict[int, float] = dict()
    search_max_x: Dict[int, float] = dict()
    search_min_y: Dict[int, float] = dict()
    search_max_y: Dict[int, float] = dict()
    search_min_z: Dict[int, float] = dict()
    search_max_z: Dict[int, float] = dict()
    
    sub_plots: Optional[go.Figure] = None
    
    number_of_layers = 0
    
    # parse general mapping info
    with open(f"{data_logging_path}mapping_result.json", "r") as result_file:
        number_of_layers = json.load(result_file)['statistics']['layers']
    
    if type(layer) is int:
        if layer >= number_of_layers:
            raise ValueError(f"Invalid layer {layer}. There are only {number_of_layers} layers in the data log.")
    
    # prepare search graph traces
    search_node_traces, search_edge_traces, search_node_stem_trace = _prepare_search_graph_scatters(
        number_of_node_traces, 
        search_node_color_scale,
        search_node_invert_color_scale,
        search_node_colorbar_title,
        search_node_colorbar_spacing,
        use3d,
        draw_stems,
        draw_search_edges,
        plotly_settings
    )
    
    # parse architecture info and prepare respective traces
    if not hide_layout:
        arch_graph = _parse_arch_graph(f"{data_logging_path}architecture.json")
        
        if architecture_node_positions is None:
            for node in arch_graph.nodes:
                print(node, arch_graph.nodes[node])
            for edge in arch_graph.edges:
                print(edge, arch_graph.edges[edge])
            architecture_node_positions = graphviz_layout(arch_graph, prog=architecture_layout_method)
        elif len(architecture_node_positions) != len(arch_graph.nodes):
            raise ValueError(f"architecture_node_positions must contain positions for all {len(arch_graph.nodes)} architecture nodes.")
        
        arch_x_spacing = _get_avg_min_distance([p[0] for p in architecture_node_positions.values()])
        arch_y_spacing = _get_avg_min_distance([p[1] for p in architecture_node_positions.values()])
        arch_x_arrow_spacing = (arch_x_spacing if arch_x_spacing > 0 else arch_y_spacing) * swap_arrow_spacing
        arch_y_arrow_spacing = (arch_y_spacing if arch_y_spacing > 0 else arch_x_spacing) * swap_arrow_spacing
        arch_x_min = min(architecture_node_positions.values(), key=lambda p: p[0])[0]
        arch_x_max = max(architecture_node_positions.values(), key=lambda p: p[0])[0]
        arch_y_min = min(architecture_node_positions.values(), key=lambda p: p[1])[1]
        arch_y_max = max(architecture_node_positions.values(), key=lambda p: p[1])[1]

        arch_edge_trace, arch_edge_label_trace = _draw_architecture_edges(
            arch_graph, architecture_node_positions, plotly_settings)
    
        arch_node_trace = go.Scatter(**plotly_settings['architecture_nodes'])
        arch_node_trace.x = []
        arch_node_trace.y = []
    
    # create figure and add traces
    if not hide_layout:
        sub_plots = make_subplots(rows=1, cols=2,specs=[[{'is_3d': use3d}, {'is_3d': False}]])
    else:
        sub_plots = go.Figure()
    sub_plots.update_layout(plotly_settings['layout'])
    
    active_trace_indices = {}
    
    active_trace_indices['search_edges'] = []
    for trace in search_edge_traces:
        active_trace_indices['search_edges'].append(len(sub_plots.data))
        sub_plots.add_trace(trace, row=1 if not hide_layout else None, col=1 if not hide_layout else None)
    
    if draw_stems:
        active_trace_indices['search_node_stems'] = len(sub_plots.data)
        sub_plots.add_trace(search_node_stem_trace, row=1 if not hide_layout else None, col=1 if not hide_layout else None)
    
    active_trace_indices['search_nodes'] = []
    for trace in search_node_traces:
        active_trace_indices['search_nodes'].append(len(sub_plots.data))
        sub_plots.add_trace(trace, row=1 if not hide_layout else None, col=1 if not hide_layout else None)
    
    xaxis1 = sub_plots.layout.scene.xaxis if use3d else (sub_plots.layout.xaxis if hide_layout else sub_plots.layout.xaxis1)
    yaxis1 = sub_plots.layout.scene.yaxis if use3d else (sub_plots.layout.yaxis if hide_layout else sub_plots.layout.yaxis1)
    xaxis1.update(**plotly_settings['search_xaxis'])
    yaxis1.update(**plotly_settings['search_yaxis'])
    if use3d:
        sub_plots.layout.scene.zaxis.update(**plotly_settings['search_zaxis'])
    
    if not hide_layout:
        active_trace_indices['arch_edges'] = len(sub_plots.data)
        sub_plots.add_trace(arch_edge_trace, row=1, col=2)
        active_trace_indices['arch_edge_labels'] = len(sub_plots.data)
        sub_plots.add_trace(arch_edge_label_trace, row=1, col=2)
        active_trace_indices['arch_nodes'] = len(sub_plots.data)
        sub_plots.add_trace(arch_node_trace, row=1, col=2)
        arch_node_trace = sub_plots.data[-1]
        
        xaxis2 = sub_plots['layout']['xaxis'] if use3d else sub_plots['layout']['xaxis2']
        yaxis2 = sub_plots['layout']['yaxis'] if use3d else sub_plots['layout']['yaxis2']
        
        plotly_settings['architecture_xaxis']['range'] = [
            arch_x_min - abs(arch_x_max-arch_x_min)*architecture_border, 
            arch_x_max + abs(arch_x_max-arch_x_min)*architecture_border
        ]
        plotly_settings['architecture_yaxis']['range'] = [
            arch_y_min - abs(arch_y_max-arch_y_min)*architecture_border,
            arch_y_max + abs(arch_y_max-arch_y_min)*architecture_border
        ]
        x_diff = plotly_settings['architecture_xaxis']['range'][1] - plotly_settings['architecture_xaxis']['range'][0]
        y_diff = plotly_settings['architecture_yaxis']['range'][1] - plotly_settings['architecture_yaxis']['range'][0]
        if x_diff == 0:
            mid = plotly_settings['architecture_xaxis']['range'][0]
            plotly_settings['architecture_xaxis']['range'] = [mid-y_diff/2, mid+y_diff/2]
        if y_diff == 0:
            mid = plotly_settings['architecture_yaxis']['range'][0]
            plotly_settings['architecture_yaxis']['range'] = [mid-x_diff/2, mid+x_diff/2]
        
        xaxis2.update(**plotly_settings['architecture_xaxis'])
        yaxis2.update(**plotly_settings['architecture_yaxis'])

    fig = go.FigureWidget(sub_plots)
    
    # update trace variables to active traces (traces are internally copied when creating fig)
    search_node_traces = []
    for i in active_trace_indices['search_nodes']:
        search_node_traces.append(fig.data[i])
    search_edge_traces = []
    for i in active_trace_indices['search_edges']:
        search_edge_traces.append(fig.data[i])
    if draw_stems:
        search_node_stem_trace = fig.data[active_trace_indices['search_node_stems']]

    if not hide_layout:
        
        arch_node_trace = fig.data[active_trace_indices['arch_nodes']]
        arch_edge_trace = fig.data[active_trace_indices['arch_edges']]
        arch_edge_label_trace = fig.data[active_trace_indices['arch_edge_labels']]
        
    # define interactive callbacks
    def visualize_search_node_layout(trace, points, selector):
        nonlocal current_node_layout_visualized
        
        if current_layer is None:
            return
        point_inds = []
        try:
            point_inds = points.point_inds
        except:
            point_inds = points
        if len(point_inds) == 0:
            return
        if current_node_layout_visualized == point_inds[0]:
            return
        current_node_layout_visualized = point_inds[0]
        node_index = list(search_graph.nodes())[current_node_layout_visualized]
        
        _visualize_layout(
            fig, 
            search_graph.nodes[node_index]['data'],
            arch_node_trace,
            architecture_node_positions,
            initial_layout,
            considered_qubit_colors,
            swap_arrow_offset,
            arch_x_arrow_spacing,
            arch_y_arrow_spacing,
            show_shared_swaps,
            current_node_layout_visualized,
            search_node_traces[0] if not use3d else None,
            plotly_settings
        )
    
    if not hide_layout:
        for trace in search_node_traces:
            if show_layout == 'hover':
                trace.on_hover(visualize_search_node_layout)
            elif show_layout == 'click':
                trace.on_click(visualize_search_node_layout)
    
    def update_timestep(timestep):
        timestep = timestep['new']
        if current_layer is None:
            return
        with fig.batch_update():
            data_x = search_node_scatter_data_x[current_layer][:timestep]
            data_y = search_node_scatter_data_y[current_layer][:timestep]
            if draw_stems:
                # line trace data: begin, end, None, begin, end, None, ...
                search_node_stem_trace.x = search_node_stem_scatter_data_x[current_layer][:timestep*3]
                search_node_stem_trace.y = search_node_stem_scatter_data_y[current_layer][:timestep*3]
                search_node_stem_trace.z = search_node_stem_scatter_data_z[current_layer][:timestep*3]
            for i, trace in enumerate(search_node_traces):
                trace.x = data_x
                trace.y = data_y
                if use3d:
                    trace.z = search_node_scatter_data_z[current_layer][i][:timestep]
    
    timestep_play = Play(
        value=1,
        min=1,
        max=2,
        step=search_progression_step,
        interval=int(1000/search_progression_speed),
        disabled=False
    )
    timestep_slider = IntSlider(min=timestep_play.min, max=timestep_play.max, step=timestep_play.step, layout=Layout(width=f'{width-232}px'))
    jslink((timestep_play, 'value'), (timestep_slider, 'value'))
    timestep_play.observe(update_timestep, names='value')
    timestep = HBox([timestep_play, timestep_slider])
    
    def update_layer(new_layer: int):
        nonlocal current_layer, search_graph, initial_layout, considered_qubit_colors, current_node_layout_visualized
        
        if current_layer == new_layer:
            return
        current_layer = new_layer
        
        if current_layer not in search_graphs:
            (
                search_graphs[current_layer],
                initial_layouts[current_layer],
                initial_qbit_positions[current_layer],
                layers_considered_qubit_colors[current_layer],
                single_qbit_multiplicities[current_layer],
                two_qbit_multiplicities[current_layer],
                individual_two_qbit_multiplicities[current_layer],
                final_node_ids[current_layer],
                search_node_scatter_data_x[current_layer],
                search_node_scatter_data_y[current_layer],
                search_node_scatter_data_z[current_layer],
                search_node_scatter_data_color[current_layer],
                search_node_stem_scatter_data_x[current_layer],
                search_node_stem_scatter_data_y[current_layer],
                search_node_stem_scatter_data_z[current_layer],
                search_edge_scatter_data_x[current_layer],
                search_edge_scatter_data_y[current_layer],
                search_edge_scatter_data_z[current_layer],
                search_min_x[current_layer],
                search_max_x[current_layer],
                search_min_y[current_layer],
                search_max_y[current_layer],
                search_min_z[current_layer],
                search_max_z[current_layer]
            ) = _load_layer_data(
                data_logging_path,
                current_layer,
                search_node_layout_method,
                search_node_layout_walker_factor,
                tapered_search_layer_heights,
                number_of_node_traces,
                use3d,
                search_node_color,
                prioritize_search_node_color,
                search_node_height,
                color_valid_mapping,
                color_final_node,
                draw_stems,
                draw_search_edges
            )
        search_graph = search_graphs[current_layer]
        initial_layout = initial_layouts[current_layer]
        considered_qubit_colors = layers_considered_qubit_colors[current_layer]
        with fig.batch_update():
            xaxis1 = fig.layout.scene.xaxis if use3d else (fig.layout.xaxis if hide_layout else fig.layout.xaxis1)
            yaxis1 = fig.layout.scene.yaxis if use3d else (fig.layout.yaxis if hide_layout else fig.layout.yaxis1)
            xaxis1.range = [
                search_min_x[current_layer] - abs(search_max_x[current_layer]-search_min_x[current_layer])*search_graph_border, 
                search_max_x[current_layer] + abs(search_max_x[current_layer]-search_min_x[current_layer])*search_graph_border
            ]
            yaxis1.range = [
                search_min_y[current_layer] - abs(search_max_y[current_layer]-search_min_y[current_layer])*search_graph_border,
                search_max_y[current_layer] + abs(search_max_y[current_layer]-search_min_y[current_layer])*search_graph_border
            ]
            if use3d:
                fig.layout.scene.zaxis.range = [
                    search_min_z[current_layer] - abs(search_max_z[current_layer]-search_min_z[current_layer])*search_graph_border,
                    search_max_z[current_layer] + abs(search_max_z[current_layer]-search_min_z[current_layer])*search_graph_border
                ]
            
            _draw_search_graph_nodes(
                search_node_traces, 
                search_node_scatter_data_x[current_layer], 
                search_node_scatter_data_y[current_layer], 
                search_node_scatter_data_z[current_layer], 
                search_node_scatter_data_color[current_layer], 
                use3d
            )
            if draw_stems:
                _draw_search_graph_stems(
                    search_node_stem_trace, 
                    search_node_stem_scatter_data_x[current_layer], 
                    search_node_stem_scatter_data_y[current_layer], 
                    search_node_stem_scatter_data_z[current_layer]
                )
            if draw_search_edges:
                _draw_search_graph_edges(
                    search_edge_traces, 
                    search_edge_scatter_data_x[current_layer], 
                    search_edge_scatter_data_y[current_layer], 
                    search_edge_scatter_data_z[current_layer], 
                    use3d
                )
            
            if not hide_layout:
                _draw_architecture_nodes(
                    arch_node_trace,
                    architecture_node_positions, 
                    considered_qubit_colors, 
                    initial_qbit_positions[current_layer], 
                    single_qbit_multiplicities[current_layer], 
                    individual_two_qbit_multiplicities[current_layer]
                )
        
        current_node_layout_visualized = None
        visualize_search_node_layout(None, [0], None)
        cl = current_layer
        current_layer = None # to prevent triggering redraw in update_timestep
        timestep_play.max = len(search_graph.nodes)
        timestep_slider.max = len(search_graph.nodes)
        timestep_play.value = len(search_graph.nodes)
        timestep_slider.value = len(search_graph.nodes)
        current_layer = cl
    
    if type(layer) is int:
        update_layer(layer)
    else:
        update_layer(0)

    layer_slider = interactive(update_layer, new_layer=IntSlider(min=0, max=number_of_layers, step=1, value=0, description='Layer:', layout=Layout(width=f'{width-80}px')))
    
    vbox_elems = [fig]
    if show_search_progression:
        vbox_elems.append(timestep)
    if layer == 'interactive':
        vbox_elems.append(layer_slider)
    return VBox(vbox_elems)