import json
import networkx as nx


def export_network_to_cytoscape_dashboard(
    graph: nx.Graph,
    filename: str = "network_dashboard.html",
    directed: bool | None = None
):
    """
    Export a NetworkX graph to an interactive Cytoscape.js dashboard HTML file.

    Features:
    - Draggable nodes
    - Click node → highlight connected edges & neighbors (colors preserved)
    - Background click resets highlight
    - Force and Circle layouts
    - Edge color legend
    - Hover tooltips
    """

    if directed is None:
        directed = graph.is_directed()

    elements = []

    # --- Nodes ---
    for node, attrs in graph.nodes(data=True):
        elements.append({
            "data": {
                "id": str(node),
                "label": str(node),
                "description": attrs.get("description", ""),
                "color": attrs.get("color", "#2563eb")
            }
        })

    # --- Edges ---
    for u, v, attrs in graph.edges(data=True):
        elements.append({
            "data": {
                "id": f"{u}_{v}",
                "source": str(u),
                "target": str(v),
                "color": attrs.get("color", "#9ca3af"),
                "description": attrs.get("description", "")
            }
        })

    elements_json = json.dumps(elements)
    arrow_shape = "triangle" if directed else "none"

    html_content = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{filename}</title>
  <script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>

  <style>
    body {{
      margin: 0;
      font-family: Inter, system-ui, sans-serif;
      background: #f3f4f6;
    }}

    #header {{
      padding: 16px 24px;
      font-size: 18px;
      font-weight: 600;
      background: white;
      border-bottom: 1px solid #e5e7eb;
    }}

    #controls {{
      padding: 12px 24px;
      background: white;
      border-bottom: 1px solid #e5e7eb;
    }}

    #cy {{
      width: 100%;
      height: 80vh;
      background: white;
    }}

    #tooltip {{
      position: absolute;
      display: none;
      background: white;
      border: 1px solid #d1d5db;
      padding: 8px 10px;
      border-radius: 8px;
      font-size: 13px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.12);
      max-width: 260px;
      pointer-events: none;
      z-index: 999;
    }}

    #legend {{
      position: absolute;
      top: 20px;
      right: 20px;
      background: white;
      border: 1px solid #d1d5db;
      border-radius: 8px;
      padding: 10px;
      font-size: 13px;
      box-shadow: 2px 2px 8px rgba(0,0,0,0.1);
    }}

    .legend-item {{
      display: flex;
      align-items: center;
      margin-bottom: 4px;
    }}

    .legend-color {{
      width: 16px;
      height: 16px;
      margin-right: 6px;
      border-radius: 3px;
    }}
  </style>
</head>
<body>

<div id="header">{filename}</div>

<div id="controls">
  <button onclick="runLayout('cose')">Force Layout</button>
  <button onclick="runLayout('circle')">Circle Layout</button>
  <button onclick="fitGraph()">Fit Graph</button>
</div>

<div id="cy"></div>
<div id="tooltip"></div>

<div id="legend">
  <strong>Edge Color Legend</strong>
  <div class="legend-item"><div class="legend-color" style="background: green;"></div>Present in both networks</div>
  <div class="legend-item"><div class="legend-color" style="background: red;"></div>Missing in provided network</div>
  <div class="legend-item"><div class="legend-color" style="background: orange;"></div>Only in reference network</div>
</div>

<script>

  const elements = {elements_json};

  const cy = cytoscape({{
    container: document.getElementById('cy'),
    elements: elements,

    style: [
      {{
        selector: 'node',
        style: {{
          'background-color': 'data(color)',
          'label': 'data(label)',
          'color': '#111827',
          'text-valign': 'center',
          'text-halign': 'center',
          'font-size': '12px',
          'width': 45,
          'height': 45,
          'opacity': 1
        }}
      }},
      {{
        selector: 'edge',
        style: {{
          'width': 2,
          'line-color': 'data(color)',
          'target-arrow-color': 'data(color)',
          'target-arrow-shape': '{arrow_shape}',
          'curve-style': 'bezier',
          'opacity': 1
        }}
      }},
      {{
        selector: '.faded',
        style: {{
          'opacity': 0.1
        }}
      }}
    ],

    layout: {{
      name: 'cose',
      animate: true,
      padding: 100,
      nodeRepulsion: 12000,
      idealEdgeLength: 180,
      edgeElasticity: 200,
      gravity: 0.2
    }}
  }});

  function runLayout(name) {{
    cy.layout({{
      name: name,
      animate: true,
      padding: 50
    }}).run();
  }}

  function fitGraph() {{
    cy.fit();
  }}

  // --- TOOLTIP ---
  const tooltip = document.getElementById('tooltip');

  cy.on('mouseover', 'node, edge', function(evt) {{
    const desc = evt.target.data('description');
    if (desc) {{
      tooltip.innerHTML = desc.replace(/\\n/g, "<br>");
      tooltip.style.display = 'block';
    }}
  }});

  cy.on('mousemove', function(evt) {{
    tooltip.style.left = evt.originalEvent.pageX + 10 + 'px';
    tooltip.style.top = evt.originalEvent.pageY + 10 + 'px';
  }});

  cy.on('mouseout', 'node, edge', function() {{
    tooltip.style.display = 'none';
  }});

  // --- CLICK HIGHLIGHT LOGIC (COLORS PRESERVED) ---
  cy.on('tap', 'node', function(evt) {{
    const node = evt.target;

    // Fade all elements
    cy.elements().addClass('faded');

    // Keep selected node visible
    node.removeClass('faded');

    // Connected edges
    const connectedEdges = node.connectedEdges();
    connectedEdges.removeClass('faded');

    // Neighbor nodes
    const neighbors = node.connectedNodes();
    neighbors.removeClass('faded');
  }});

  // Background click → reset
  cy.on('tap', function(evt) {{
    if (evt.target === cy) {{
      cy.elements().removeClass('faded');
    }}
  }});

</script>

</body>
</html>
"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html_content)

    graph_type = "Directed Graph" if directed else "Undirected Graph"
    print(f"Dashboard saved to {filename} ({graph_type})")
