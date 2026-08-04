#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#set page(width: auto, height: auto, margin: 5mm, fill: white)
#set text(font: "Inter")

#diagram(
  debug: false,
  node-stroke: black,
  node-corner-radius: 5pt,

  node((0, -1), [Run Covasim\ simulation], name: <covasim>),
  node((0, 0), [Initialize\ BloodSim agents], name: <agents>),
  node((1, 0.5), [Select random\ alive agents], name: <strategy>),
  node((1, 1.5), [Aggregate\ lab tests], name: <labs>),

  node((1, 2.5), [Run\ detectors], name: <stats>),
  node((1, 3.5), [Compute\ metrics], name: <metrics>),
  node((1, 4.5), [Visualize\ results], name: <viz>),
  
  node((0, 1.5), [Assign\ comorbidity], name: <comorb>),
  node((0, 2.5), [Set baseline\ CRP], name: <baseline>),
  node((0, 3.5), [Detect infection-\ state transitions], name: <infection>),
  node((0, 4.5), [Model\ CRP dynamics], name: <crp>),
  
  node((0,1), align(top+center)[for each agent], stroke: none),
  node(enclose: ((0,1), (0,4.5)), inset: 10pt, name: <agentcycle>),
  node((1,0), align(top+center)[for each day], stroke: none),
  node(enclose: ((1,0), (1,1.5)), inset: 10pt, name: <bloodsimcycle>),
  
  node((0, 5.5), name: <tmp1>),
  node((0.5, 5.5), name: <tmp2>),
  node((0.5, -1), name: <tmp3>),
  node((1, -1), name: <tmp4>),

  edge(<covasim>, <agents>, "-|>"),
  edge(<agents>, <agentcycle>, "-|>"),
  edge(<strategy>, <labs>, "-|>"),
  edge(<stats>, <metrics>, "-|>"),
  edge(<metrics>, <viz>, "-|>"),
  edge(<comorb>, <baseline>, "-|>"),
  edge(<baseline>, <infection>, "-|>"),
  edge(<infection>, <crp>, "-|>"),
  edge(<bloodsimcycle>, <stats>, "-|>"),
  
  edge(<agentcycle>, <tmp1>, "-"),
  edge(<tmp1>, <tmp2>, "-"),
  edge(<tmp2>, <tmp3>, "-"),
  edge(<tmp3>, <tmp4>, "-"),
  edge(<tmp4>, <bloodsimcycle>, "-|>"),
)