pacman::p_load(
  DiagrammeR, # for flow diagrams
  networkD3, # For alluvial/Sankey diagrams
  tidyverse
) # data management and visualization




DiagrammeR::grViz("               # All instructions are within a large character string
digraph immune_genes_data_cleaning_pipeline {    # 'digraph' means 'directional graph', then the graph name 
  
  # graph statement
  #################
  graph [layout = dot,
         #splines=line,
         rankdir = TB,            # layout top-to-bottom
         fontsize = 12,
         ranksep = 0.7,  # Increase the rank separation value
         fontname = 'Times New Roman']
  

  # nodes
  #################
  node [shape = box,           # shape = circle
       fixedsize = true,
       width = 1.3,
      style='filled',
      fillcolor='WhiteSmoke',
      bgcolor='WhiteSmoke',
      color='DimGray'
]                      
  
A   [label = 'Housekeeping\ngene'] #, fillcolor ='Gray70'
B1   [label = 'Gene by gene\nCt'] #, fillcolor ='Gray70'
B2  [label = 'Discard\nall genes',
fontcolor = red] 
D4   [label = 'Diff.Ct'] #, fillcolor ='Gray70'
E5   [label = 'Valid:\nuse mean Ct', fontcolor = red]
E6   [label = 'Invalid\n(discard)',
fontcolor = red]
IMP1   [label = 'Impute\nrelative conc.',fontcolor = red]
J3 [label = '',fillcolor=white,color=white] #empty node used for spacing

#Notes [label = 'Ct: Cycle Threshold; LOD: Limit of Detection; T.E.: Technical Error ', color = white]
  
  # edges
  #######
edge [color = 'DimGray', arrowhead = 'vee']
A   -> B1 [label = 'valid']
A   -> B2 [label = 'invalid']
B1 -> IMP1 [xlabel = 'mean≥LOD']
B1 -> D4 [xlabel = 'mean<LOD']
B1 -> IMP2 [label = '2 NA']
B1 -> D2 [label = '1 NA']
D2 -> IMP2 [label = '≥(LOD-T.E.)']
D2 -> F1 [label = '<(LOD-T.E.)']
D4 -> E5 [xlabel = '≤pipetting error\nthreshold']
D4 -> E6 [label = '>pipetting error\nthreshold']
IMP1 -> J3 [color = white]

  # grouped edge
  #{E3} -> IMP [label = '']


 subgraph cluster_group {
    label = 'Group'
    #style = rounded;
    color = black;
    D2 [label = 'Non NA Ct'];
    F1 [label = 'Discard NA Ct', fontcolor = red];
    IMP2   [label = 'Impute\nrelative conc.',fontcolor = red]

  }

#ranking of parts of graph
{rank = same; D2; D4;}
{rank = same; IMP1; IMP2; J3; F1; E5; E6} # 
#{rank = sink; Notes; }
}
")


# Export the graph as a PNG image
export_graph(graph, file_name = "graph.png")


# # valid:\n mean-2sd<Ct<mean+2sd\nΔCt<0.5\nneither NA


## #F3   [label = 'rerun RTqPCR\ndiscard\nuse mean', fontcolor = darkgreen, fontsize=9]

## #E4   [label = 'discard higher\nCt value', fontcolor = red]


#E3 -> F2
#D4 -> E4 [label = '>T.E.']
#D4 -> E6 [label = 'Poisson threshold < ΔCt ≤ T.E.']


#############################
