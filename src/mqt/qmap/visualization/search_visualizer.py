from __future__ import annotations

from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any, Callable, Literal

from mqt.qmap.visualization.visualize_search_graph import Position, SearchNode, visualize_search_graph

if TYPE_CHECKING:
    from ipywidgets import Widget


class SearchVisualizer:
    def __init__(self, data_logging_path: str | None = None) -> None:
        if data_logging_path is not None:
            self.data_logging_path = data_logging_path
            self.data_logging_tmp_dir = None
        else:
            self.data_logging_tmp_dir = TemporaryDirectory()
            self.data_logging_path = self.data_logging_tmp_dir.name

    def close(self):
        if self.data_logging_tmp_dir is not None:
            self.data_logging_tmp_dir.cleanup()
            self.data_logging_path = None

    def visualize_search_graph(
        self,
        layer: int | Literal["interactive"] = "interactive",  # 'interactive' (slider menu) | index
        architecture_node_positions: dict[int, Position] | None = None,
        architecture_layout_method: Literal[
            "dot", "neato", "fdp", "sfdp", "circo", "twopi", "osage", "patchwork"
        ] = "sfdp",
        search_node_layout_method: Literal[
            "walker", "dot", "neato", "fdp", "sfdp", "circo", "twopi", "osage", "patchwork"
        ] = "walker",
        search_node_layout_walker_factor: float = 0.6,
        search_graph_border: float = 0.05,
        architecture_border: float = 0.05,
        swap_arrow_spacing: float = 0.05,
        swap_arrow_offset: float = 0.05,
        use3d: bool = True,
        projection: Literal["orthographic", "perspective"] = "perspective",
        width: int = 1400,
        height: int = 700,
        draw_search_edges: bool = True,
        search_edges_width: float = 0.5,
        search_edges_color: str = "#888",
        search_edges_dash: str = "solid",  # 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot', string containing a dash length list in pixels or percentages (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%')
        tapered_search_layer_heights: bool = True,
        show_layout: Literal["hover", "click"] | None = "hover",
        show_swaps: bool = True,
        show_shared_swaps: bool = True,  # if one swap moves 2 considered qubits -> combine into one two-colored and two-headed arrow
        color_valid_mapping: str | None = "green",  # static HTML color (e.g. 'blue' or '#0000FF')
        color_final_node: str | None = "red",  # static HTML color (e.g. 'blue' or '#0000FF')
        search_node_color: str
        | (Callable[[SearchNode], float] | list[str | Callable[[SearchNode], float]]) = "total_cost",
        # 'total_cost' | 'fixed_cost' | 'heuristic_cost' | 'lookahead_penalty' | static HTML color (e.g. 'blue' or '#0000FF') |
        # function that takes a SearchNode and returns a float (e.g. lambda n: n.fixed_cost + n.heuristic_cost)
        # node fields: {"fixed_cost": float, "heuristic_cost": float, "lookahead_penalty": float, "is_valid_mapping": bool,
        #                       "final": bool, "depth": int, "layout": Tuple[int, ...], "swaps": Tuple[Tuple[int, int], ...]}
        # or list of the above if 3d graph is used and multiple points per node are defined in search_node_height (lengths need to match)
        # if in that case no is list is provided all points per node will be the same color
        prioritize_search_node_color: bool
        | list[
            bool
        ] = False,  # if True, search_node_color will be prioritized over color_valid_mapping and color_final_node
        search_node_color_scale: str | list[str] = "YlGnBu",  # https://plotly.com/python/builtin-colorscales/
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
        search_node_invert_color_scale: bool | list[bool] = True,
        search_node_colorbar_title: str | list[str | None] | None = None,
        search_node_colorbar_spacing: float = 0.06,
        search_node_height: str
        | (Callable[[SearchNode], float] | list[str | Callable[[SearchNode], float]]) = "total_cost",
        # just as with search_node_color (without color strings), but possible to specify a list to draw multiple point per node on different heights
        # only applicable if use3d is True
        draw_stems: bool = False,
        stems_width: float = 0.7,
        stems_color: str = "#444",
        stems_dash: str = "solid",  # 'solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot', string containing a dash length list in pixels or percentages (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%')
        show_search_progression: bool = True,
        search_progression_step: int = 10,
        search_progression_speed: float = 2,  # steps per second
        plotly_settings: dict[str, dict[str, Any]] | None = None,
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
    ) -> Widget:
        if plotly_settings is None:
            plotly_settings = {}
        if self.data_logging_path is None:
            msg = "SearchVisualizer has already been closed and data logs have been discarded."
            raise ValueError(msg)
        return visualize_search_graph(
            data_logging_path=self.data_logging_path,
            layer=layer,
            architecture_node_positions=architecture_node_positions,
            architecture_layout_method=architecture_layout_method,
            search_node_layout_method=search_node_layout_method,
            search_node_layout_walker_factor=search_node_layout_walker_factor,
            search_graph_border=search_graph_border,
            architecture_border=architecture_border,
            swap_arrow_spacing=swap_arrow_spacing,
            swap_arrow_offset=swap_arrow_offset,
            use3d=use3d,
            projection=projection,
            width=width,
            height=height,
            draw_search_edges=draw_search_edges,
            search_edges_width=search_edges_width,
            search_edges_color=search_edges_color,
            search_edges_dash=search_edges_dash,
            tapered_search_layer_heights=tapered_search_layer_heights,
            show_layout=show_layout,
            show_swaps=show_swaps,
            show_shared_swaps=show_shared_swaps,
            color_valid_mapping=color_valid_mapping,
            color_final_node=color_final_node,
            search_node_color=search_node_color,
            prioritize_search_node_color=prioritize_search_node_color,
            search_node_color_scale=search_node_color_scale,
            search_node_invert_color_scale=search_node_invert_color_scale,
            search_node_colorbar_title=search_node_colorbar_title,
            search_node_colorbar_spacing=search_node_colorbar_spacing,
            search_node_height=search_node_height,
            draw_stems=draw_stems,
            stems_width=stems_width,
            stems_color=stems_color,
            stems_dash=stems_dash,
            show_search_progression=show_search_progression,
            search_progression_step=search_progression_step,
            search_progression_speed=search_progression_speed,
            plotly_settings=plotly_settings,
        )
