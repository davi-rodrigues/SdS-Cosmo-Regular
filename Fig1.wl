
(* ::Title:: *)
(*Fig. 1*)


(* ::Text:: *)
(*From "Schwarzschild-de Sitter spacetime in regular coordinates with cosmological time" (2025) by L. Lima and D.C. Rodrigues.*)


(* ::Author:: *)
(*Davi C. Rodrigues*)


SetDirectory[NotebookDirectory[]];

Clear[m, CC]
Am = A /. Solve[A == 1 - (2 m)/r - r^2 CC/3 (1 - A^-1), A][[1]];
Ap = A /. Solve[A == 1 - (2 m)/r - r^2 CC/3 (1 - A^-1), A][[2]];
f = 1 - 2 m/r -  CC/3 r^2; 

Echo[Expand[Am],  "Am = \!\(\*FormBox[SubscriptBox[\(A\), \(-\)],
TraditionalForm]\) = "];
Echo[Expand[Ap],  "Ap = \!\(\*FormBox[SubscriptBox[\(A\), \(+\)],
TraditionalForm]\) = "];
Echo[Expand[f],  "f = "];


m = 1;
CC = 1/100;

legend = LineLegend[
  {"\!\(\*SubscriptBox[\(A\), \(+\)]\)",  "\!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\)", "\!\(\*SubscriptBox[\(A\), \(-\)]\)"},
  LegendFunction-> (Framed[#, RoundingRadius -> 5, FrameStyle-> Thin]&),
  LegendMarkerSize->{{30, 15}}
];
    
plot = Plot[
  {Ap, f,  Am, 1, - (CC/3) r^2}, 
  {r, 0, 21},
  PlotRange-> {Automatic, {-2.0, 1.1}}, 
  PlotTheme->"Scientific",
  Background->White,
  FrameStyle->Directive[Black, 15],
  PlotLegends-> Placed[legend, Scaled[{0.25, 0.27}]],
  LabelStyle->{FontFamily -> "Times", FontSize-> 15},
  PlotStyle -> {
   {RGBColor[0.1, 0.2, 0.6], Thick}, 
   {RGBColor[0.9, 0.6, 0.1], Thick}, 
   {RGBColor[0.1, 0.6, 0.5], Thick},
   {RGBColor[0.1, 0.2, 0.6], Dashed, Thickness[0.004]} ,
   {RGBColor[0.1, 0.6, 0.5], DotDashed, Thickness[0.004]}
  },
  GridLinesStyle->Directive[Darker@Gray, Dashed],
  PlotRangePadding-> 0
]

Export["plotAlphas.pdf", plot]



