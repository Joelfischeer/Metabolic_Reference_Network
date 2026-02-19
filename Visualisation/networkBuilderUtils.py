import networkx as nx
from networkx.classes import Graph

def export_network_to_html(
    graph: nx.Graph,
    filename: str = "network.html",
    truncate_len: int = 20,
    directed: bool | None = None
):
    nodes = []
    edges = []

    if directed is None:
        directed = graph.is_directed()

    # --- Nodes ---
    for node, attrs in graph.nodes(data=True):
        description = attrs.get("description", "")
        tooltip = description.replace('\n', '<br>') if description else ""
        nodes.append({
            "id": str(node),
            "label": str(node),
            "customTooltip": tooltip  # only use our custom floating tooltip
        })

    # --- Edges ---
    for u, v, attrs in graph.edges(data=True):
        color = attrs.get("color", "#888")
        description = attrs.get("description", "")
        tooltip = description.replace('\n', '<br>') if description else ""
        edge_data = {
            "from": str(u),
            "to": str(v),
            "color": color,
            "customTooltip": tooltip
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
  <style>
    body {{ font-family: Arial, sans-serif; }}
    #title {{ font-size: 18px; margin: 10px; font-weight: bold; }}
    #network {{ width: 100%; height: 800px; border: 1px solid lightgray; }}
    #legend {{
      position: absolute; top: 20px; right: 20px; background: white;
      padding: 12px; border: 1px solid #ccc; border-radius: 8px;
      box-shadow: 2px 2px 8px rgba(0,0,0,0.1); font-size: 14px;
    }}
    .legend-item {{ display: flex; align-items: center; margin-bottom: 6px; }}
    .legend-color {{ width: 16px; height: 16px; margin-right: 8px; border-radius: 3px; }}

    #tooltip {{
        position: absolute;
        display: none;
        background: white;
        border: 1px solid #333;
        padding: 6px;
        border-radius: 5px;
        box-shadow: 2px 2px 6px rgba(0,0,0,0.2);
        pointer-events: none;
        max-width: 300px;
        font-size: 13px;
        line-height: 1.3;
    }}
  </style>
</head>
<body>
  <div id="title">{filename}</div>
  <div id="network"></div>
  <div id="legend">
    <strong>Edge Legend</strong>
    <div class="legend-item"><div class="legend-color" style="background: green;"></div>Present in both networks</div>
    <div class="legend-item"><div class="legend-color" style="background: red;"></div>Missing in provided network</div>
    <div class="legend-item"><div class="legend-color" style="background: orange;"></div>Only in reference network</div>
  </div>

  <div id="tooltip"></div>

  <script type="text/javascript">
    var nodes = new vis.DataSet({nodes});
    var edges = new vis.DataSet({edges});
    var container = document.getElementById('network');
    var data = {{ nodes: nodes, edges: edges }};
    var options = {{
      nodes: {{
        shape: 'dot',
        size: 20
      }},
      edges: {{
        smooth: {{ type: 'dynamic' }}
      }},
      physics: {{ enabled: false }},
      interaction: {{
        hover: true,
        dragNodes: true
      }}
    }};
    var network = new vis.Network(container, data, options);

    // Floating tooltip div
    var tooltipDiv = document.getElementById('tooltip');

    // Node hover
    network.on("hoverNode", function(params) {{
        var node = nodes.get(params.node);
        tooltipDiv.innerHTML = node.customTooltip;
        tooltipDiv.style.left = params.event.pageX + 10 + "px";
        tooltipDiv.style.top = params.event.pageY + 10 + "px";
        tooltipDiv.style.display = "block";
    }});
    network.on("blurNode", function(params) {{
        tooltipDiv.style.display = "none";
    }});

    // Edge hover
    network.on("hoverEdge", function(params) {{
        var edge = edges.get(params.edge);
        tooltipDiv.innerHTML = edge.customTooltip;
        tooltipDiv.style.left = params.event.pageX + 10 + "px";
        tooltipDiv.style.top = params.event.pageY + 10 + "px";
        tooltipDiv.style.display = "block";
    }});
    network.on("blurEdge", function(params) {{
        tooltipDiv.style.display = "none";
    }});
  </script>
</body>
</html>
"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"Network saved to {filename} ({graph_type})")
