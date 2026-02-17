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
      description = attrs.get("description", "")

      # Format hover tooltip
      if description:
          tooltip = f"{description.replace(chr(10), '\n')}"
      else:
          tooltip = f"{node}"

      nodes.append({
          "id": str(node),
          "label": str(node),   # stays visible on graph
          "title": tooltip      # shown only on hover
      })


    # Edges
    for u, v, attrs in graph.edges(data=True):
      color = attrs.get("color", "#888")

      description = attrs.get("description", "")

      if description:
          tooltip = f"{description.replace(chr(10), '\n')}"
      else:
          tooltip = f"{u} ↔ {v}"

      edge_data = {
          "from": str(u),
          "to": str(v),
          "color": color,
          "title": tooltip
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
    #legend {{
      position: absolute;
      top: 20px;
      right: 20px;
      background: white;
      padding: 12px;
      border: 1px solid #ccc;
      border-radius: 8px;
      box-shadow: 2px 2px 8px rgba(0,0,0,0.1);
      font-size: 14px;
    }}

    .legend-item {{
      display: flex;
      align-items: center;
      margin-bottom: 6px;
    }}

    .legend-color {{
      width: 16px;
      height: 16px;
      margin-right: 8px;
      border-radius: 3px;
    }}

  </style>
</head>
<body>
  <div id="title">{filename}</div>
  <div id="network"></div>
  <div id="legend">
    <strong>Edge Legend</strong>
    <div class="legend-item">
      <div class="legend-color" style="background: green;"></div>
      Present in both networks
    </div>
    <div class="legend-item">
      <div class="legend-color" style="background: red;"></div>
      Missing in provided network
    </div>
    <div class="legend-item">
      <div class="legend-color" style="background: orange;"></div>
      Only in reference network
    </div>
  </div>

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


