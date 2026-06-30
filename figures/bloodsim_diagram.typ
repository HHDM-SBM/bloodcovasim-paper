#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#set page(width: auto, height: auto, margin: 5mm, fill: white)

#diagram(
  debug: true,
  node-stroke: black,
  node-corner-radius: 5pt,

  node((0, 0), [Conducting\ a Covasim Simulation], name: <aa>),
  node((0, 1), [Initializing\ BloodSim agents], name: <ab>),
  node((0, 2), [Selecting\ a Testing Strategy], name: <ac>),
  node((0, 3), [Aggregation\ of laboratory tests], name: <ad>),

  node((1, 1), [Initialization], name: <bb>),
  node((1, 3), [Conducting\ Statistical tests], name: <bd>),

  node((2, 0), [Determining\ the presence of\ a comorbid condition], name: <ca>),
  node((2, 1), [Determination of\ the baseline CRP level], name: <cb>),
  node((2, 2), [Determining when\ infectious state change], name: <cc>),
  node((2, 3), [Metrics Computation], name: <cd>),

  node((3, 1), [Modeling\ CRP Dynamics], name: <db>),
  node((3, 3), [Visualization], name: <dd>),

  


  edge(<aa>, <ab>, "-|>"),
  edge(<ab>, <ac>, "-|>"),
  edge(<ac>, <ad>, "-|>"),
  edge(<ad>, <bd>, "-|>"),
  edge(<bd>, <cd>, "-|>"),
  edge(<cd>, <dd>, "-|>"),
  edge(<ab>, <bb>, [For Every agent], "-|>", bend: 40deg),
  edge(<bb>, <ab>, "-|>", bend: 40deg),

  node((1.35, 1), [], name: <splitL>, stroke: none),
  node((1.35, 0), [], name: <splitLU>, stroke: none),
  node((1.35, 2), [], name: <splitLD>, stroke: none),
  node((2.6, 1), [], name: <splitR>, stroke: none),
  node((2.6, 0), [], name: <splitRU>, stroke: none),
  node((2.6, 2), [], name: <splitRD>, stroke: none),

  edge(<bb>, <splitL>, "-"),
  edge(<splitL>, <splitLU>, "-"),
  edge(<splitL>, <splitLD>, "-"),
  edge(<splitLU>, <ca>, "-|>"),
  edge(<splitL>, <cb>, "-|>"),
  edge(<splitLD>, <cc>, "-|>"),

  edge(<ca>, <splitRU>, "-"),
  edge(<cb>, <splitR>, "-"),
  edge(<cc>, <splitRD>, "-"),
  edge(<splitRU>, <splitR>, "-"),
  edge(<splitRD>, <splitR>, "-"),
  edge(<splitR>, <db>, "-|>"),
)