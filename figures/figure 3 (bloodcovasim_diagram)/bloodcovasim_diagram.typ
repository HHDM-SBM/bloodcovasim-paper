#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#set page(width: auto, height: auto, margin: 5mm, fill: white)
#set text(font: "Inter")

#diagram(
  debug: false,
  node-stroke: black,
  node-corner-radius: 5pt,

  node((0, 0), [Run Covasim\ simulation], name: <covasim>),
  node((0, 1), [Initialize\ BloodCovasim agents], name: <agents>),
  node((1, 0), [Select random\ alive agents], name: <strategy>),
  node((1, 1), [Aggregate\ lab tests], name: <labs>),

  node((1, 2), [Run\ detectors], name: <stats>),
  node((1, 3), [Compute\ metrics], name: <metrics>),
  node((1, 4), [Visualize\ results], name: <viz>),
  
  node((0, 2), [Set baseline\ CRP], name: <baseline>),
  node((0, 3), [Detect infection-\ state transitions], name: <infection>),
  node((0, 4), [Model\ CRP dynamics], name: <crp>),
  
  node((0,1.625), align(top+center)[for each agent], stroke: none),
  node(enclose: ((0,2), (0,4)), inset: 15pt, name: <agentcycle>),
  node((1,-0.475), align(top+center)[for each day], stroke: none),
  node(enclose: ((1,0), (1,1)), inset: 15pt, name: <bloodcovasimcycle>),
  
  node((0, 5), name: <tmp1>),
  node((0.5, 5), name: <tmp2>),
  node((0.5, -1), name: <tmp3>),
  node((1, -1), name: <tmp4>),

  edge(<covasim>, <agents>, "-|>"),
  edge(<agents>, <agentcycle>, "-|>"),
  edge(<strategy>, <labs>, "-|>"),
  edge(<stats>, <metrics>, "-|>"),
  edge(<metrics>, <viz>, "-|>"),
  edge(<baseline>, <infection>, "-|>"),
  edge(<infection>, <crp>, "-|>"),
  edge(<bloodcovasimcycle>, <stats>, "-|>"),
  
  edge(<agentcycle>, <tmp1>, "-"),
  edge(<tmp1>, <tmp2>, "-"),
  edge(<tmp2>, <tmp3>, "-"),
  edge(<tmp3>, <tmp4>, "-"),
  edge(<tmp4>, <bloodcovasimcycle>, "-|>"),
)