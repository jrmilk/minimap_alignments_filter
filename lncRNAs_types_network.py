# =========================================================
# Interactive + SVG lncRNA similarity network
# Using Jaccard similarity matrix
#
# Outputs:
#   - lncrna_jaccard_network.html
#   - lncrna_jaccard_network.svg
#
# Google Colab / Jupyter compatible
# =========================================================

!pip install pyvis networkx pandas matplotlib -q

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pyvis.network import Network
from IPython.display import IFrame


# =========================================================
# Jaccard similarity matrix
# =========================================================

data = {
    "lincRNA": {
        "lincRNA": 1.0,
        "exonic_sense": 0.3715,
        "exonic_antisense": 0.0997,
        "intronic_sense": 0.2637,
        "intronic_antisense": 0.1292,
    },

    "exonic_sense": {
        "lincRNA": 0.3715,
        "exonic_sense": 1.0,
        "exonic_antisense": 0.1122,
        "intronic_sense": 0.4469,
        "intronic_antisense": 0.1521,
    },

    "exonic_antisense": {
        "lincRNA": 0.0997,
        "exonic_sense": 0.1122,
        "exonic_antisense": 1.0,
        "intronic_sense": 0.0854,
        "intronic_antisense": 0.2539,
    },

    "intronic_sense": {
        "lincRNA": 0.2637,
        "exonic_sense": 0.4469,
        "exonic_antisense": 0.0854,
        "intronic_sense": 1.0,
        "intronic_antisense": 0.1648,
    },

    "intronic_antisense": {
        "lincRNA": 0.1292,
        "exonic_sense": 0.1521,
        "exonic_antisense": 0.2539,
        "intronic_sense": 0.1648,
        "intronic_antisense": 1.0,
    }
}


# =========================================================
# Convert dataframe
# =========================================================

df = pd.DataFrame(data)

print("\nJaccard similarity matrix:\n")
print(df)


# =========================================================
# Create graph
# =========================================================

G = nx.Graph()

classes = list(df.index)


# =========================================================
# Add nodes
# =========================================================

for cls in classes:

    G.add_node(cls)


# =========================================================
# Add edges
# =========================================================

for i in range(len(classes)):

    for j in range(i + 1, len(classes)):

        a = classes[i]
        b = classes[j]

        similarity = float(df.loc[a, b])

        if similarity > 0:

            G.add_edge(
                a,
                b,
                value=similarity
            )


# =========================================================
# Interactive HTML network
# =========================================================

net = Network(
    height="900px",
    width="100%",
    bgcolor="white",
    font_color="black",
    notebook=True,
    cdn_resources="in_line"
)

net.from_nx(G)


# =========================================================
# Customize nodes
# =========================================================

for node in net.nodes:

    node["size"] = 40
    node["label"] = node["id"]


# =========================================================
# Customize edges
# =========================================================

for edge in net.edges:

    similarity = edge.get("value", 0)

    edge["width"] = similarity * 30

    edge["label"] = f"{similarity:.2f}"

    edge["title"] = f"Jaccard similarity: {similarity:.4f}"


# =========================================================
# Physics/layout
# =========================================================

net.set_options("""
var options = {

  "physics": {

    "enabled": true,

    "barnesHut": {

      "gravitationalConstant": -4000,
      "centralGravity": 0.15,
      "springLength": 200,
      "springConstant": 0.03,
      "damping": 0.09

    }
  },

  "interaction": {

    "hover": true,
    "dragNodes": true,
    "dragView": true,
    "zoomView": true
  }
}
""")


# =========================================================
# Save interactive HTML
# =========================================================

html_output = "lncrna_jaccard_network.html"

net.save_graph(html_output)

print(f"\nSaved HTML:")
print(html_output)


# =========================================================
# Static SVG version
# =========================================================

plt.figure(figsize=(10, 10))

# layout
pos = nx.spring_layout(
    G,
    seed=42,
    k=1.2
)

# edge widths
edge_widths = [
    G[u][v]["value"] * 15
    for u, v in G.edges()
]

# draw network
nx.draw_networkx_nodes(
    G,
    pos,
    node_size=4000
)

nx.draw_networkx_labels(
    G,
    pos,
    font_size=12
)

nx.draw_networkx_edges(
    G,
    pos,
    width=edge_widths
)

# edge labels
edge_labels = {
    (u, v): f"{d['value']:.2f}"
    for u, v, d in G.edges(data=True)
}

nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=edge_labels,
    font_size=10
)

plt.axis("off")

svg_output = "lncrna_jaccard_network.svg"

plt.savefig(
    svg_output,
    format="svg",
    bbox_inches="tight"
)

plt.close()

print(f"\nSaved SVG:")
print(svg_output)


# =========================================================
# Show interactive network
# =========================================================

IFrame(
    src=html_output,
    width="100%",
    height=900
)
