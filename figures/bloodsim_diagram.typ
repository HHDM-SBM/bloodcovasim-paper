#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#import fletcher.shapes: diamond
#set page(width: auto, height: auto, margin: 5mm, fill: white)
#set text(font: "Inter")

#diagram(
  debug: false,
  node-stroke: black,
  node-corner-radius: 5pt,

  node((1, 0), [Run\ Covasim simulation], name: <sim>),
  node((1, 1), [Initialize\ BloodSim agents], name: <agents>),
  node((1, 2), [More\ agents?], name: <loop>, shape: diamond),
  node((1, 3), [Select\ testing strategy], name: <strategy>),
  node((1, 4), [Aggregate\ lab tests], name: <labs>),

  node((2, 2), [Run\ detectors], name: <stats>),
  node((2, 3), [Compute\ metrics], name: <metrics>),
  node((2, 4), [Visualize\ results], name: <viz>),
  
  node((0, 1), [Assign\ comorbidity], name: <comorb>),
  node((0, 2), [Set baseline\ CRP], name: <baseline>),
  node((0, 3), [Detect infection-\ state transitions], name: <infection>),
  node((0, 4), [Model\ CRP dynamics], name: <crp>),

  node((0.5, 2), name: <from_loop1>, stroke: none),
  node((0.5, 1), name: <from_loop2>, stroke: none),
  node((0.5, 4), name: <to_loop1>, stroke: none),
  node((0.5, 2.5), name: <to_loop2>, stroke: none),
  node((1.5, 4), name: <to_stats1>, stroke: none),
  node((1.5, 2), name: <to_stats2>, stroke: none),


  edge(<sim>, <agents>, "-|>"),
  edge(<agents>, <loop>, "-|>"),
  edge(<strategy>, <labs>, "-|>"),
  edge(<stats>, <metrics>, "-|>"),
  edge(<metrics>, <viz>, "-|>"),
  edge(<comorb>, <baseline>, "-|>"),
  edge(<baseline>, <infection>, "-|>"),
  edge(<infection>, <crp>, "-|>"),

  edge(<loop>, <strategy>, [no], "-|>"),
  edge(<loop>, <from_loop1>, [yes], "-"),
  edge(<from_loop1>, <from_loop2>, "-"),
  edge(<from_loop2>, <comorb>, "-|>"),

  edge(<crp>, <to_loop1>, "-"),
  edge(<to_loop1>, <to_loop2>, "-"),
  edge(<to_loop2>, <loop>, "-|>"),

  edge(<labs>, <to_stats1>, "-"),
  edge(<to_stats1>, <to_stats2>, "-"),
  edge(<to_stats2>, <stats>, "-|>"),
)