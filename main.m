(* ::Package:: *)

BeginPackage["Mirion`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*Usage*)


(* ::Section:: *)
(*Function*)


IBP::usage = "\
IBP[u,x] \:7ed9\:51fa \[Integral]u \[DifferentialD]x \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,v,x] \:7ed9\:51fa \[Integral]u \[DifferentialD]v \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,{x,a,b}] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)u \[DifferentialD]x \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,v,{x,a,b}] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)u \[DifferentialD]v \:7684\:5206\:90e8\:79ef\:5206.";  


IBS::usage = "\
IBS[f,ex\[Rule]et] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex\[Rule]et \:540e\:7684\:7ed3\:679c.
IBS[f,x\[Rule]et,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 x\[Rule]et \:540e\:7684\:7ed3\:679c.
IBS[f,ex\[Rule]t,x] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex\[Rule]t \:540e\:7684\:7ed3\:679c.
IBS[f,ex\[Rule]et,x,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex\[Rule]et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex\[Rule]et] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:6362\:5143 ex\[Rule]et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},x\[Rule]et,t] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:6362\:5143 x\[Rule]et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex\[Rule]t,x] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:6362\:5143 ex\[Rule]t \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex\[Rule]et,x,t] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:6362\:5143 ex\[Rule]et \:540e\:7684\:7ed3\:679c.";


ApartArcTan::usage = "\
ApartArcTan[expr] \:7ed9\:51fa arctan(expr) \:7684\:88c2\:9879.
ApartArcTan[expr,x] \:7ed9\:51fa arctan(expr) \:5173\:4e8e x \:7684\:88c2\:9879."; 


ComplexFactor::usage = "\
ComplexFactor[poly] \:7ed9\:51fa poly \:5728\:590d\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3.
ComplexFactor[poly,x] \:7ed9\:51fa poly \:5173\:4e8e x \:5728\:590d\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3."; 


ComplexApart::usage = "\
ComplexApart[expr] \:7ed9\:51fa expr \:5728\:590d\:6570\:57df\:4e0a\:7684\:88c2\:9879.
ComplexApart[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:590d\:6570\:57df\:4e0a\:7684\:88c2\:9879."; 


RealFactor::usage = "\
RealFactor[poly] \:7ed9\:51fa poly \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3.
RealFactor[poly,x] \:7ed9\:51fa poly \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3."; 


RealApart::usage = "\
RealApart[expr] \:7ed9\:51fa expr \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879.
RealApart[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879."; 


FindIdentities::usage = "\
FindIdentities[expr1,expr2,x] \:7ed9\:51fa\:5173\:4e8e expr1,expr2 \:7684\:6052\:7b49\:5f0f.
";


CorrectionTest::usage = "\
CorrectionTest[org,res,{x,\!\(\*SubscriptBox[\(x\), \(min\)]\),\!\(\*SubscriptBox[\(x\), \(max\)]\)}] \:ff08\:5b9e\:9a8c\:6027\:ff09\:68c0\:9a8c res \:5173\:4e8e org \:7684\:4fee\:6b63\:662f\:5426\:6b63\:786e."; 


BivariablePlot::usage = "\
BivariablePlot[list,x] \:ff08\:5b9e\:9a8c\:6027\:ff09\:7ed8\:5236\:591a\:5143 list \:7684\:5173\:7cfb\:56fe.
";


(* ::Section:: *)
(*Option*)


GenerateConstant::usage = "GenerateConstant \:5173\:4e8e\:662f\:5426\:751f\:6210\:5e38\:6570\:7684\:9009\:9879"; 
SimplifyFunction::usage = "SimplifyFunction \:5173\:4e8e\:5316\:7b80\:51fd\:6570\:7684\:9009\:9879"; 
VerifyDiscontinuities::usage = "VerifyDiscontinuities \:5173\:4e8e\:662f\:5426\:9a8c\:8bc1\:95f4\:65ad\:70b9\:7684\:9009\:9879"; 


Begin["`Private`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*Definition*)


(* ::Section:: *)
(*Private*)


(* ::Subsection:: *)
(*TEST*)


Attributes[TEST] = {HoldAll}; 
TEST[code_, n_:1] := Module[{}, ClearSystemCache[]; AbsoluteTiming[Do[code, n - 1]; code]]


(* ::Subsection:: *)
(*ExpressionPivot*)


Attributes[ExpressionPivot] = {Listable}; 
ExpressionPivot[expr_] := FirstCase[expr, _Symbol?(Not @* NumericQ), Symbol, {-1}]; 


(* ::Subsection:: *)
(*CoefficientSeparation*)


Attributes[CoefficientSeparation] = {Listable}; 
CoefficientSeparation[expr_, x_] /; FreeQ[expr, x] := {expr, 1}; 
CoefficientSeparation[expr_, x_] := Replace[expr, Longest[c_.]*(r_) /; FreeQ[c, x] -> {c, r}]; 


(* ::Subsection:: *)
(*FirstCoefficient*)


Attributes[FirstCoefficient] = {Listable}; 
FirstCoefficient[poly_, x_] := Coefficient[poly, x, Exponent[poly, x]]; 


(* ::Subsection:: *)
(*PolynomialRoots*)


Attributes[PolynomialRoots] = {Listable}; 
Options[PolynomialRoots] = {ToRadicals -> False}; 
PolynomialRoots[a_, x_, OptionsPattern[]] /; FreeQ[a, x] := {}; 
PolynomialRoots[(a_.)*(x_) + (b_.), x_, OptionsPattern[]] /; FreeQ[{a, b}, x] := {-(b/a)}; 
PolynomialRoots[poly_, x_, OptionsPattern[]] := List @@ (Roots[poly == 0, x, Cubics -> OptionValue[ToRadicals], 
      Quartics -> OptionValue[ToRadicals]] /. Equal -> (#2 & )); 


(* ::Subsection:: *)
(*ArcTanLimitAtPositiveInfinity*)


Attributes[ArcTanLimitAtPositiveInfinity] = {Listable}; 
ArcTanLimitAtPositiveInfinity[expr_, x_] := Module[{nd, e, r}, nd = NumeratorDenominator[Together[expr]]; e = Exponent[nd, x]; 
     r = Divide @@ Coefficient[nd, x, e]; Piecewise[{{Sign[r]*(Pi/2), Subtract @@ e > 0}, {0, Subtract @@ e < 0}, {ArcTan[r], True}}]]; 


(* ::Section:: *)
(*Public*)


(* ::Subsection:: *)
(*IBP*)


Options[IBP] = {Assumptions -> $Assumptions}; 
IBP[u_, v_, x_Symbol, OptionsPattern] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; 
     u*v - coef*Inactive[Integrate][rem, x]]; 
IBP[f_, x_Symbol, OptionsPattern[]] := IBP[f, x, x]; 
IBP[u_, v_, {x_Symbol, a_, b_}, OptionsPattern[]] := Module[{assum = OptionValue[Assumptions], c, r}, 
    {c, r} = CoefficientSeparation[v*D[u, x], x]; Limit[u*v, x -> b, Assumptions -> assum] - Limit[u*v, x -> a, Assumptions -> assum] - 
      c*Inactive[Integrate][r, {x, a, b}]]; 
IBP[f_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := IBP[f, x, {x, a, b}, opts]; 


(* ::Subsection:: *)
(*IBS*)


Options[IBS] = {Assumptions -> $Assumptions, Inactive -> False}; 
IBS[f_, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{assum = OptionValue[Assumptions], inactive = OptionValue[Inactive], u, result, c, r}, 
    result = If[ex === x, f /. x -> et, IntegrateChangeVariables[Inactive[Integrate][f, x], u, u == ex, Assumptions -> assum][[1]] /. 
          C[_] -> 0 /. u -> et]*D[et, t]; If[inactive, {c, r} = CoefficientSeparation[Together[result], t]; c*Inactive[Integrate][r, t], 
      result]]; 
IBS[f_, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, x -> et, x, t, opts]; 
IBS[f_, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, ex -> t, x, t, opts]; 
IBS[f_, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{assum = OptionValue[Assumptions], u, temp, result, aa, bb, c, r}, 
    temp = If[ex === x, Inactive[Integrate][f, {x, a, b}] /. x -> u, IntegrateChangeVariables[Inactive[Integrate][f, {x, a, b}], u, u == ex, 
        Assumptions -> assum]]; {result, {t, aa, bb}} = List @@ If[et === t, temp /. u -> t, IntegrateChangeVariables[temp, t, u == et, 
         Assumptions -> assum]]; {c, r} = CoefficientSeparation[Together[result], t]; c*Inactive[Integrate][r, {t, a, b}]]; 
IBS[f_, {a_, b_}, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, x -> et, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> t, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 


(* ::Subsection:: *)
(*ApartArcTan*)


Attributes[ApartArcTan] = {Listable}; 
Options[ApartArcTan] = {GenerateConstant -> True, SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := Module[{result, count = 0, poles, sample, c = {}}, 
    result = Sum[ArcTan[OptionValue[SimplifyFunction][(x - Re[r])/((count += Sign[#1]; #1) & )[Im[r]]]], {r, SolveValues[expr - I == 0, x]}]; 
     If[OptionValue[GenerateConstant], result += ArcTanLimitAtPositiveInfinity[expr, x] - count*(Pi/2); 
       poles = DeleteDuplicates[SolveValues[Denominator[Together[expr]] == 0, x, Reals]]; If[Length[poles] =!= 0, 
        sample = (((1/2)*Table[#1[[i]] + #1[[i + 1]], {i, Length[#1] - 1}] & )[Prepend[#1, First[#1] - 1]] & )[N[poles]]; 
         c = Pi*Round[(1/Pi)*Table[ArcTan[expr] - result /. x -> i, {i, sample}]]; 
         result += Piecewise[Prepend[Table[{c[[i + 1]], poles[[i]] < x < poles[[i + 1]]}, {i, Length[poles] - 1}], 
            {First[c], x < First[poles]}]]]; ]; result]; 
ApartArcTan[expr_, opts:OptionsPattern[]] := ApartArcTan[expr, ExpressionPivot[expr], opts]; 


(* ::Subsection:: *)
(*ComplexFactor*)


Attributes[ComplexFactor] = {Listable}; 
Options[ComplexFactor] = Options[PolynomialRoots]; 
ComplexFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := FirstCoefficient[poly, x]*
    Times @@ (x - PolynomialRoots[poly, x, ToRadicals -> OptionValue[ToRadicals]]); 
ComplexFactor[poly_, opts:OptionsPattern[]] := ComplexFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection:: *)
(*ComplexApart*)


Attributes[ComplexApart] = {Listable}; 
Options[ComplexApart] = Options[ComplexFactor]; 
ComplexApart[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/ComplexFactor[#2, x, ToRadicals -> OptionValue[ToRadicals]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
ComplexApart[poly_, opts:OptionsPattern[]] := ComplexApart[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection:: *)
(*RealFactor*)


Attributes[RealFactor] = {Listable}; 
Options[RealFactor] = {SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
RealFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := 
  FirstCoefficient[poly, x]*Times @@ Table[OptionValue[SimplifyFunction][
      Piecewise[{{x - r, Element[r, Reals]}, {(x - Re[r])^2 + Im[r]^2, True}}]], 
     {r, DeleteDuplicates[PolynomialRoots[poly, x], N[Conjugate[#1] == #2] & ]}]
RealFactor[poly_, opts:OptionsPattern[]] := RealFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection:: *)
(*RealApart*)


Attributes[RealApart] = {Listable}; 
Options[RealApart] = Options[RealFactor]; 
RealApart[expr_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/RealFactor[#2, x, SimplifyFunction -> OptionValue[SimplifyFunction]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
RealApart[expr_, opts:OptionsPattern[]] := RealApart[expr, ExpressionPivot[expr], opts]; 


(* ::Subsection:: *)
(*FindIdentities*)


FindIdentities[expr1_, expr2_, x_Symbol] /; RationalExpressionQ[expr1, x] && RationalExpressionQ[expr2, x] := 
   Module[{p1, p2, roots, limit}, p1 = Numerator[Simplify[expr1]]*Denominator[Simplify[expr2]]; 
     p2 = Numerator[Simplify[expr2]]*Denominator[Simplify[expr1]]; roots = DeleteDuplicates[SolveValues[D[p1, x]*p2 == D[p2, x]*p1, x]]; 
     DeleteDuplicates[Flatten[Reap[Do[limit = Simplify[Limit[p1/p2, x -> i]]; If[limit =!= 0 && Element[limit, Rationals], 
            Sow[Defer[Evaluate[p1]] - limit*p2 == Factor[p1 - limit*p2]]]; , {i, roots}]][[2]]]]]; 


(* ::Section:: *)
(*Experimental*)


(* ::Subsection:: *)
(*CorrectionTest*)


CorrectionTest::overtime = "`` \:8ba1\:7b97\:8d85\:65f6"; 
Options[CorrectionTest] = {Skip -> {}, TimeConstraint -> Infinity, VerifyDiscontinuities -> {}, WorkingPrecision -> $MinPrecision}; 
CorrectionTest[org_, res_, {x_, xmin_, xmax_}, OptionsPattern[]] := 
   Module[{table = {{Style["\:5b9a\:4e49\:57df", Bold], Style["\:95f4\:65ad\:70b9", Bold]}, {Style[True, Bold, RGBColor["#4081ff"]], 
        Style[False, Bold, RGBColor["#4081ff"]]}, {Style[True, Bold, RGBColor["#eb7100"]], Style[False, Bold, RGBColor["#eb7100"]]}, 
       {Style["\:56fe\:50cf", Bold], SpanFromLeft}, {Null, SpanFromLeft}}, skip = OptionValue[Skip], cons = OptionValue[TimeConstraint], 
     veri = OptionValue[VerifyDiscontinuities], prec = OptionValue[WorkingPrecision], dom = True}, 
    If[ !ListQ[skip], skip = {skip}]; If[ !ListQ[veri], veri = {veri}]; If[ !AssociationQ[cons], 
      cons = Association["\:5b9a\:4e49\:57df" -> cons, "\:95f4\:65ad\:70b9" -> cons, "\:56fe\:50cf" -> cons]]; If[ !MemberQ[skip, "\:5b9a\:4e49\:57df"], 
      TimeConstrained[table[[2,1]] = Item[Style[dom = Reduce[FunctionDomain[org, x], Element[x, Reals]], Bold, RGBColor["#4081ff"]], 
           ItemSize -> Automatic]; table[[3,1]] = Item[Style[Reduce[FunctionDomain[res, x] && dom, Element[x, Reals]], Bold, 
            RGBColor["#eb7100"]], ItemSize -> Automatic]; , cons["\:5b9a\:4e49\:57df"], Message[CorrectionTest::overtime, "\:5b9a\:4e49\:57df"]; ]; ]; 
     If[ !MemberQ[skip, "\:95f4\:65ad\:70b9"], 
      TimeConstrained[table[[2,2]] = Item[Style[Reduce[FunctionDiscontinuities[org, x] && dom, Element[x, Reals]], Bold, RGBColor["#4081ff"]], 
           ItemSize -> Automatic]; table[[3,2]] = Item[Style[Reduce[FunctionDiscontinuities[res, x] && dom, Element[x, Reals]], Bold, 
            RGBColor["#eb7100"]], ItemSize -> Automatic]; , cons["\:95f4\:65ad\:70b9"], Message[CorrectionTest::overtime, "\:95f4\:65ad\:70b9"]; ]; ]; 
     If[ !MemberQ[skip, "\:56fe\:50cf"], 
      TimeConstrained[table[[5,1]] = Plot[{org, res, Evaluate[D[res /. {Abs[p_] -> Piecewise[{{p, p > 0}, {-p, p < 0}}], Sign[p_] -> 
                Piecewise[{{1, p > 0}, {-1, p < 0}}]}, x]]}, {x, xmin, xmax}, PlotStyle -> 96, ImageSize -> Medium, 
          Epilog -> {Black, PointSize[Medium], Point[Block[{$MinPrecision = prec}, N[Table[{i, res /. x -> i}, {i, veri}]]]]}], cons["\:56fe\:50cf"], 
        Message[CorrectionTest::overtime, "\:56fe\:50cf"]; ]; ]; Grid[table, Alignment -> {Center, Center}, Spacings -> {1, 1}, 
      ItemSize -> {Scaled[0.5], Automatic}, Frame -> All]]; 


(* ::Subsection:: *)
(*BivariablePlot*)


Options[BivariablePlot] = {PlotLabels -> None}; 
BivariablePlot[list_List, x_Symbol, OptionsPattern[]] := Module[{labels, isValid, const, constQ, edges, vertexes, usedVertexes, edgeStylize, 
     vertexStylize}, labels = OptionValue[PlotLabels]; isValid = labels =!= None && Length[list] === Length[labels]; 
     constQ[value_] := FreeQ[value, x] &&  !PossibleZeroQ[value]; edgeStylize[edge_, value_] := 
      Labeled[edge, Placed[Style[value, 14, FontFamily -> "CMU Serif"], 0.5]]; vertexStylize[value_] := 
      Framed[Style[value, Black, 14, FontFamily -> "CMU Serif"], FrameStyle -> None]; 
     edges = Table[Piecewise[{{{UndirectedEdge[bivars[[1]], bivars[[2]]], const}, 
          constQ[const = Simplify[list[[bivars[[1]]]]^2 + list[[bivars[[2]]]]^2]]}, 
         {Piecewise[{{{DirectedEdge[bivars[[2]], bivars[[1]]], -const}, TrueQ[const < 0]}, {{DirectedEdge[bivars[[1]], bivars[[2]]], const}, 
             True}}], constQ[const = Simplify[list[[bivars[[1]]]]^2 - list[[bivars[[2]]]]^2]]}, {Nothing, True}}], 
       {bivars, Subsets[Range[Length[list]], {2}]}]; usedVertexes = DeleteDuplicates[Cases[edges[[All,1]], _Integer, {2}]]; 
     vertexes = Table[i -> (Evaluate[Inset[vertexStylize[Piecewise[{{labels[[i]] == list[[i]], isValid}, {list[[i]], True}}]], #1]] & ), 
       {i, usedVertexes}]; edges = edgeStylize @@@ edges; Graph[edges, PlotTheme -> "DiagramBlack", VertexLabels -> None, 
      VertexShapeFunction -> vertexes, PlotRangePadding -> Scaled[0.1]]]; 


End[];


EndPackage[];
