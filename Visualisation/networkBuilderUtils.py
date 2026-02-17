import networkx as nx
import pandas as pd
import dash
from dash import html
import dash_cytoscape as cyto
from networkx.classes import Graph


def export_network_to_html(
    graph: nx.Graph,
    filename: str = "network.html",
    truncate_len: int = 20,
    directed: bool | None = None
):
    """
    Export a NetworkX graph as a standalone HTML file with draggable nodes using vis.js.
    Negative edges are red, positive edges are gray.
    Node and edge labels are truncated and rounded to 2 decimal places.
    Includes info whether graph is directed or not.

    Parameters
    ----------
    graph : nx.Graph
        The NetworkX graph to export.
    filename : str, optional
        Output HTML filename.
    truncate_len : int, optional
        Maximum length for string attributes in labels/tooltips.
    directed : bool or None, optional
        Whether to treat the graph as directed. If None, will use graph.is_directed().
    """
    nodes = []
    edges = []

    # Decide directedness
    if directed is None:
        directed = graph.is_directed()

    # Nodes
    for node, attrs in graph.nodes(data=True):
        title_parts = []
        for k, v in attrs.items():
            if isinstance(v, (int, float)):
                value = round(v, 2)
            else:
                value = str(v)[:truncate_len]
            title_parts.append(f"{k}: {value}")
        title = "<br>".join(title_parts)
        nodes.append({
            "id": str(node),
            "label": str(node),
            "title": title
        })

    # Edges
    for u, v, attrs in graph.edges(data=True):
        weight = attrs.get("weight", 1)
        if isinstance(weight, (int, float)):
            weight_display = round(weight, 2)
        else:
            weight_display = str(weight)[:truncate_len]
        color = attrs.get("color", "#888")
        edge_data = {
            "from": str(u),
            "to": str(v),
            "label": str(weight_display),
            "color": color
        }
        if directed:
            edge_data["arrows"] = "to"
        edges.append(edge_data)

    graph_type = "Directed Graph" if directed else "Undirected Graph"

    html_content = f"""
<!DOCTYPE html>
<html>
<head>
  <title>Network Graph</title>
  <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
  <style type="text/css">
    body {{
      font-family: Arial, sans-serif;
    }}
    #title {{
      font-size: 18px;
      margin: 10px;
      font-weight: bold;
    }}
    #network {{
      width: 100%;
      height: 800px;
      border: 1px solid lightgray;
    }}
  </style>
</head>
<body>
  <div id="title">{filename}</div>
  <div id="network"></div>
  <script type="text/javascript">
    var nodes = new vis.DataSet({nodes});
    var edges = new vis.DataSet({edges});

    var container = document.getElementById('network');
    var data = {{
      nodes: nodes,
      edges: edges
    }};
    var options = {{
      nodes: {{
        shape: 'dot',
        size: 20
      }},
      edges: {{
        font: {{ align: 'top' }},
        smooth: {{ type: 'dynamic' }}
      }},
      physics: {{
        enabled: false
      }},
      interaction: {{
        hover: true,
        tooltipDelay: 100,
        dragNodes: true
      }}
    }};
    var network = new vis.Network(container, data, options);
  </script>
</body>
</html>
"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Network saved to {filename} ({graph_type})")


