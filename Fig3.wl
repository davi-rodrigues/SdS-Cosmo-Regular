(* ::Package:: *)

(* ::Title:: *)
(*Fig. 3*)


(* ::Text:: *)
(*From "Schwarzschild-de Sitter spacetime in regular coordinates with cosmological time" (2025) by L. Lima and D.C. Rodrigues.*)


(* ::Author:: *)
(*Davi C. Rodrigues*)


SetDirectory[NotebookDirectory[]];

$PlotTheme = {
  "Scientific", 
  "FrameGrid", 
  "BoldColor", 
  "SerifLabels", 
  "SizeScale", 
  {"BackgroundColor", White}
};

optionsPlot = {FrameStyle-> Directive[{Black, FontFamily -> "Times", FontSize -> 15}]};

Clear[\[CapitalLambda], m]
f[x_]= 1 - 2 m/x - \[CapitalLambda]/3 x^2;

a[\[Eta]_] = -(Sqrt[3]/(Sqrt[\[CapitalLambda]] \[Eta])); (* a'>0 *)



(* ::Section:: *)
(*Negative f regions*)


\[CapitalLambda] = 10^-2;
m = 1;

\[ScriptL][\[Eta]t_, rt_] = a[Tan[\[Eta]t \[Pi]/2]] Tan[rt \[Pi]/2];

frameticks = {
  {
    {{-1,-1,{0,0.02}},{-0.75,-0.75,{0,0.02}},{-0.5,-0.5,{0,0.02}}, {-0.25,-0.25,{0,0.02}}, {0,0,{0,0.02}}} (*left*), 
    None  (*right*)
  },
  {
    {{0,0,{0,0.02}},{0.25,0.25,{0,0.02}},{0.5,0.5,{0,0.02}}, {0.75,0.75,{0,0.02}}, {1.0,1,{0,0.02}}}, (*Bottom*)
    None (*top*) 
  }  
};

region = RegionPlot[
  {f[\[ScriptL][\[Eta]t, rt]]< 0}, 
  {rt, 0, 1}, 
  {\[Eta]t, -1, 0}, 
  FrameStyle->Directive[12, Black],
  BoundaryStyle->None, 
  Evaluate@optionsPlot, 
  PlotPoints->60, 
  PlotRangePadding->None,
  FrameTicks -> frameticks
];

contour = ContourPlot[
  Evaluate@Table[a[Tan[\[Eta]t \[Pi]/2]] Tan[rt \[Pi]/2] == i, 
  {i, {10^-1, 10^-0.5, 1, 10^0.5, 10, 10^1.5, 10^2, 10^2.5, 10^3}}], 
  {rt, 0, 1}, 
  {\[Eta]t, -1, 0}, 
  Evaluate@optionsPlot, 
  ContourStyle-> Directive[Darker@Gray, Thickness[0.0015], Dashed],
  PlotPoints->50,
  MaxRecursion->4
];

show = Show[region, contour, FrameLabel-> {"\!\(\*OverscriptBox[\(r\), \(-\)]\)", "\!\(\*OverscriptBox[\(\[Eta]\), \(-\)]\)"}]


(* ::Section:: *)
(*SdS: Outgoing geodesics*)


(* ::Text:: *)
(*Outgoing light geodesics must satisfy*)


(* ::DisplayFormula:: *)
(*\[Eta]'[r] = SubPlus[A]^-1 = -(1/\[ScriptL]^2) 3/\[CapitalLambda] SubMinus[A] = 1/\[ScriptL]^2 3/\[CapitalLambda] 1/2 (Sqrt[f[\[ScriptL]]^2 + 4 \[ScriptL]^2 \[CapitalLambda]/3] - f[\[ScriptL]])*)


Clear[\[Eta]sol, \[Eta], \[ScriptL], \[ScriptL]sol, AmAbs, AmAbssol, \[Eta]prime, \[Eta]primesol, rsolDomainMax, rsolMin, rsolMax, plotGeodesics];

\[ScriptL][r_] = r a[\[Eta][r]];
\[ScriptL]sol[\[Eta]1_][r_] = r a[\[Eta]sol[\[Eta]1][r]];

AmAbs[r_] = 1/2 (Sqrt[f[\[ScriptL][r]]^2 + 4 \[ScriptL][r]^2 \[CapitalLambda]/3] - f[\[ScriptL][r]]); (*This is the absolute value of the negative A, Subscript[A, -]*)
AmAbssol[\[Eta]1_][r_] = 1/2 (Sqrt[f[\[ScriptL]sol[\[Eta]1][r]]^2 + 4 \[ScriptL]sol[\[Eta]1][r]^2 \[CapitalLambda]/3] - f[\[ScriptL]sol[\[Eta]1][r]]);

\[Eta]prime[r_] = 1/\[ScriptL][r]^2 3/\[CapitalLambda] AmAbs[r];
\[Eta]primesol[\[Eta]1_][r_] = 1/\[ScriptL]sol[\[Eta]1][r]^2 3/\[CapitalLambda] AmAbssol[\[Eta]1][r];

\[Eta]sol[\[Eta]1_?NumberQ] := \[Eta]sol[\[Eta]1] =  Block[{},
  Off[NDSolveValue::precw];
  Off[NDSolveValue::ndsz];
  N @ NDSolveValue[
    {
      \[Eta]'[r] == \[Eta]prime[r](*1/\[ScriptL][r]^23/\[CapitalLambda] AmAbs[r]*), 
      \[Eta][1] == \[Eta]1
    }, \[Eta], {r, 0.1, 100}, 
    MaxStepSize -> 0.005, 
    MaxSteps -> Infinity, 
    WorkingPrecision->30, 
    PrecisionGoal -> 4, 
    AccuracyGoal -> \[Infinity]
  ] // Return;
  On[NDSolveValue::precw];
  On[NDSolveValue::ndsz];
];

rsolMin[\[Eta]1_] := \[Eta]sol[\[Eta]1]["Domain"][[1,1]]
rsolDomainMax[\[Eta]1_] := rsolDomainMax[\[Eta]1] = \[Eta]sol[\[Eta]1]["Domain"][[1,2]]; 
rsolMax[\[Eta]1_] := rsolMax[\[Eta]1] = Block[{r}, 
  Min[
    rsolDomainMax[\[Eta]1],
    r /. FindRoot[\[Eta]sol[\[Eta]1][r] == 0, {r, Evaluate[rsolMin[\[Eta]1] + 0.001], Evaluate[rsolDomainMax[\[Eta]1]]}]
  ]
];

plotGeodesics[\[Eta]1_] := plotGeodesics[\[Eta]1] = Block[
  {rtMin, rtMax, tableArrowPositions, ff},
  rtMin = rt /. FindRoot[Tan[rt \[Pi]/2] == rsolMin[\[Eta]1], {rt, 0, 1}];
  rtMax = rt /. FindRoot[Tan[rt \[Pi]/2] == rsolMax[\[Eta]1], {rt, 0, 1}];
  ff[rt_] = 2/\[Pi] ArcTan @ \[Eta]sol[\[Eta]1][Tan[rt \[Pi]/2]];
  Off[InterpolatingFunction::dprec];
  plotCurves = Plot[
    ff[rt], 
    {rt, rtMin, rtMax}, 
    PlotRange->All,
    PlotStyle->Thickness[0.006]
  ]; 
 tableArrowPositions = Table[{rt, ff[rt]}, {rt, rtMin, rtMax, 0.005}];
 listplotArrows = ListPlot[
   tableArrowPositions, 
   Joined -> True, 
   PlotStyle->Thickness[0.006]
 ] /. Line[x_] :> {Arrowheads[Table[.04, {3} (*Number of arrow heads*)]], Arrow[x]};
 Show[plotCurves, listplotArrows]
];

plotOut = Show[show, plotGeodesics[-10^-10], plotGeodesics[-2], plotGeodesics[-10^10]]


(* ::Section:: *)
(*SdS: Ingoing geodesics*)


(* ::Text:: *)
(*Ingoing light geodesics must satisfy*)


(* ::DisplayFormula:: *)
(*\[Eta]'[r] = - SubPlus[A] =  1/\[ScriptL]^2 3/\[CapitalLambda] SubMinus[A] = - (1/\[ScriptL]^2) 3/\[CapitalLambda] 1/2 (Sqrt[f[\[ScriptL]]^2 + 4 \[ScriptL]^2 \[CapitalLambda]/3] - f[\[ScriptL]])*)


Clear[\[Eta]sol, \[Eta], \[ScriptL], \[ScriptL]sol, AmAbs, AmAbssol, \[Eta]prime, \[Eta]primesol, rsolDomainMax, rsolMin, rsolMax, plotGeodesicsIn];

\[CapitalLambda] = 10^-2;
m = 1;

\[ScriptL][r_] = r a[\[Eta][r]];
\[ScriptL]sol[\[Eta]1_][r_] = r a[\[Eta]sol[\[Eta]1][r]];

AmAbs[r_] = 1/2 (Sqrt[f[\[ScriptL][r]]^2 + 4 \[ScriptL][r]^2 \[CapitalLambda]/3] - f[\[ScriptL][r]]);
AmAbssol[\[Eta]1_][r_] = 1/2 (Sqrt[f[\[ScriptL]sol[\[Eta]1][r]]^2 + 4 \[ScriptL]sol[\[Eta]1][r]^2 \[CapitalLambda]/3] - f[\[ScriptL]sol[\[Eta]1][r]]);

\[Eta]prime[r_] = - (1/\[ScriptL][r]^2) 3/\[CapitalLambda] AmAbs[r];
\[Eta]primesol[\[Eta]1_][r_] = - (1/\[ScriptL]sol[\[Eta]1][r]^2) 3/\[CapitalLambda] AmAbssol[\[Eta]1][r];

\[Eta]sol[\[Eta]1_?NumberQ] := \[Eta]sol[\[Eta]1] =  Block[{},
  Off[NDSolveValue::precw];
  Off[NDSolveValue::ndsz];
  N @ NDSolveValue[
    {
      \[Eta]'[r] == \[Eta]prime[r](*1/\[ScriptL][r]^23/\[CapitalLambda] AmAbs[r]*), 
      \[Eta][1] == \[Eta]1
    }, \[Eta], {r, 1. 10^-10, 100}, 
    MaxStepSize -> 0.001, 
    MaxSteps -> Infinity, 
    WorkingPrecision->30, 
    PrecisionGoal -> 4, 
    AccuracyGoal -> \[Infinity]
  ] // Return;
  On[NDSolveValue::precw];
  On[NDSolveValue::ndsz];
];

rsolDomainMin[\[Eta]1_] := \[Eta]sol[\[Eta]1]["Domain"][[1,1]];
rsolDomainMax[\[Eta]1_] := \[Eta]sol[\[Eta]1]["Domain"][[1,2]]; 

rsolMax[\[Eta]1_] := rsolDomainMax[\[Eta]1];
rsolMin[\[Eta]1_] := rsolMin[\[Eta]1] = Block[{r}, 
  Max[
    rsolDomainMin[\[Eta]1],
    Quiet[r /. FindRoot[\[Eta]sol[\[Eta]1][r] == 0, {r, Evaluate[rsolDomainMin[\[Eta]1]], Evaluate[rsolMax[\[Eta]1] -  0.001]}]]
  ]
];

plotGeodesicsIn[\[Eta]1_] := plotGeodesicsIn[\[Eta]1] = Block[
  {rtMin, rtMax, tableArrowPositions, ff},
  rtMin = rt /. FindRoot[Tan[rt \[Pi]/2] == rsolMin[\[Eta]1], {rt, 0, 1}];
  rtMax = rt /. FindRoot[Tan[rt \[Pi]/2] == rsolMax[\[Eta]1], {rt, 0, 1}];
  ff[rt_] = 2/\[Pi] ArcTan @ \[Eta]sol[\[Eta]1][Tan[rt \[Pi]/2]];
  Off[InterpolatingFunction::dprec];
  plotCurves = Plot[
    ff[rt], 
    {rt, rtMin, rtMax}, 
    PlotRange->All,
    PlotStyle->Directive[Thickness[0.006], Purple]
  ]; 
 tableArrowPositions = Reverse @ Table[{rt, ff[rt]}, {rt, rtMin, rtMax, 0.005}];
 listplotArrows = ListPlot[
   tableArrowPositions, 
   Joined -> True, 
   PlotStyle->Directive[Thickness[0.006], Purple]
 ] /. Line[x_] :> {Arrowheads[Table[.04, {3} (*Number of arrow heads*)]], Arrow[x]};
 Show[plotCurves, listplotArrows]
];

plotEtaRgeodesics = Show[plotOut, plotGeodesicsIn[-0.00001], plotGeodesicsIn[-2],plotGeodesicsIn[-10^5], FrameLabel-> {}]


Export["plotEtaRgeodesics.pdf", plotEtaRgeodesics]
