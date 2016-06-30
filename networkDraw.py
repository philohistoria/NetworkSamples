import networkx as nx
import codecs
import matplotlib.pyplot as plt
from networkx.readwrite import json_graph
f = codecs.open("romantic.csv",'r',encoding="utf-8")
G=nx.read_adjlist(f,delimiter=",")
#out = json_graph.node_link_data(G)
nx.write_pajek(G,"pajek.txt")
#print out
#nx.draw_spring(G)
#plt.show()