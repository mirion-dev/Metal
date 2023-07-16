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
IBS[f,x\[Rule]et,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:5173\:4e8e x\[Rule]et \:7684\:79ef\:5206\:6362\:5143.
IBS[f,ex\[Rule]t,x] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:5173\:4e8e ex\[Rule]t \:7684\:79ef\:5206\:6362\:5143.
IBS[f,ex\[Rule]et,x,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:5173\:4e8e ex\[Rule]et \:7684\:79ef\:5206\:6362\:5143.
IBS[f,{a,b},x\[Rule]et,t] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:5173\:4e8e x\[Rule]et \:7684\:79ef\:5206\:6362\:5143.
IBS[f,{a,b},ex\[Rule]t,x] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:5173\:4e8e ex\[Rule]t \:7684\:79ef\:5206\:6362\:5143.
IBS[f,{a,b},ex\[Rule]et,x,t] \:7ed9\:51fa \!\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]\)f \[DifferentialD]x \:5173\:4e8e ex\[Rule]et \:7684\:79ef\:5206\:6362\:5143.";


ApartArcTan::usage = "\
ApartArcTan[expr] \:7ed9\:51fa arctan(expr) \:7684\:88c2\:9879.
ApartArcTan[expr,x] \:7ed9\:51fa arctan(expr) \:5173\:4e8e x \:7684\:88c2\:9879."; 


RealFactor::usage = "\
RealFactor[expr] \:7ed9\:51fa expr \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3.
RealFactor[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3."; 


RealApart::usage = "\
RealApart[expr] \:7ed9\:51fa expr \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879.
RealApart[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879."; 


CorrectionTest::usage = "\
CorrectionTest[org,res,{x,\!\(\*SubscriptBox[\(x\), \(min\)]\),\!\(\*SubscriptBox[\(x\), \(max\)]\)}] \:68c0\:9a8c res \:5173\:4e8e org \:7684\:4fee\:6b63\:662f\:5426\:6b63\:786e."; 


FindIdentities::usage = "\
FindIdentities[expr1,expr2,x] \:7ed9\:51fa\:5173\:4e8e expr1,expr2 \:7684\:6052\:7b49\:5f0f.
";


BivariablePlot::usage = "\
BivariablePlot[list,x] \:7ed8\:5236\:591a\:5143 list \:7684\:5173\:7cfb\:56fe.
";


(* ::Subsection:: *)
(*Beta*)


Li2Transform::usage = "Beta";
Ti2Transform::usage = "Beta";


(* ::Section:: *)
(*Option*)


GenerateConstant::usage = "GenerateConstant \:5173\:4e8e\:751f\:6210\:5e38\:6570\:7684\:9009\:9879"; 
SimplifyFunction::usage = "SimplifyFunction \:5173\:4e8e\:5316\:7b80\:51fd\:6570\:7684\:9009\:9879"; 
VerifyDiscontinuities::usage = "VerifyDiscontinuities \:5173\:4e8e\:9a8c\:8bc1\:95f4\:65ad\:70b9\:7684\:9009\:9879"; 


Begin["`Private`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*Definition*)


(* ::Section:: *)
(*Private*)


(* ::Subsection:: *)
(*ArcTanLimit*)


Attributes[ArcTanLimit] = {Listable}; 
ArcTanLimit[expr_, x_] := Module[{nd, en, ed, r}, nd = NumeratorDenominator[Together[expr]]; {en, ed} = Exponent[nd, x]; 
     r = Divide @@ Coefficient[nd, x, {en, ed}]; Piecewise[{{Sign[r]*(Pi/2), en > ed}, {0, en < ed}, {ArcTan[r], True}}]]; 


(* ::Section:: *)
(*Public*)


(* ::Subsection:: *)
(*ExpressionPivot*)


Attributes[ExpressionPivot] = {Listable}; 
ExpressionPivot[expr_, default_:Missing[]] := FirstCase[expr, _Symbol?(Not @* NumericQ), default, {-1}]; 


(* ::Subsection:: *)
(*CoefficientSeparation*)


Attributes[CoefficientSeparation] = {Listable}; 
CoefficientSeparation[expr_, x_Symbol] := If[FreeQ[expr, x], {expr, 1}, Replace[expr, Longest[c_.]*(r_) /; FreeQ[c, x] -> {c, r}]]; 


(* ::Subsection:: *)
(*IBP*)


IBP[u_, v_, x_Symbol] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; u*v - coef*Inactive[Integrate][rem, x]]; 
IBP[f_, x_Symbol] := IBP[f, x, x]; 
IBP[u_, v_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; 
     Limit[u*v, x -> b, opts] - Limit[u*v, x -> a, opts] - coef*Inactive[Integrate][rem, {x, a, b}]]; 
IBP[f_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := IBP[f, x, {x, a, b}, opts]; 


(* ::Subsection:: *)
(*IBS*)


IBS[f_, ex_ -> et_, x_Symbol, t_Symbol, opts:OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{u}, If[ex === x, f /. x -> et, IntegrateChangeVariables[Inactive[Integrate][f, x], u, u == ex, opts][[1]] /. u -> et]*D[et, t] /. 
     C[_] -> 0]; 
IBS[f_, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, x -> et, x, t, opts]; 
IBS[f_, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, ex -> t, x, t, opts]; 
IBS[f_, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, x_Symbol, t_Symbol, opts:OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{u, temp}, temp = If[ex === x, Inactive[Integrate][f, {x, a, b}] /. x -> u, IntegrateChangeVariables[
        Inactive[Integrate][f, {x, a, b}], u, u == ex, opts]]; If[et === t, temp /. u -> t, IntegrateChangeVariables[temp, t, u == et, opts]]]; 
IBS[f_, {a_, b_}, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, x -> et, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> t, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 


(* ::Subsection:: *)
(*ApartArcTan*)


Attributes[ApartArcTan] = {Listable}; 
Options[ApartArcTan] = {GenerateConstant -> True, SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := Module[{result, count = 0, poles, sample, c = {}}, 
    result = Sum[ArcTan[OptionValue[SimplifyFunction][(x - Re[r])/((count += Sign[#1]; #1) & )[Im[r]]]], {r, SolveValues[expr == I, x]}]; 
     If[OptionValue[GenerateConstant], result += ArcTanLimit[expr, x] - count*(Pi/2); 
       poles = DeleteDuplicates[SolveValues[Denominator[Together[expr]] == 0, x, Reals]]; If[Length[poles] =!= 0, 
        sample = (((1/2)*Table[#1[[i]] + #1[[i + 1]], {i, Length[#1] - 1}] & )[Prepend[#1, First[#1] - 1]] & )[N[poles]]; 
         c = Pi*Round[(1/Pi)*Table[ArcTan[expr] - result /. x -> i, {i, sample}]]; 
         result += Piecewise[Prepend[Table[{c[[i + 1]], poles[[i]] < x < poles[[i + 1]]}, {i, Length[poles] - 1}], 
            {First[c], x < First[poles]}]]]; ]; result]; 
ApartArcTan[expr_, opts:OptionsPattern[]] := ApartArcTan[expr, ExpressionPivot[expr, Symbol], opts]; 


(* ::Subsection:: *)
(*RealFactor*)


Attributes[RealFactor] = {Listable}; 
Options[RealFactor] = {SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
RealFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := 
   Module[{roots, fcoef}, If[FreeQ[poly, x], Return[poly]]; fcoef = Coefficient[poly, x, Exponent[poly, x]]; 
     roots = DeleteDuplicates[SolveValues[poly == 0, x], N[Conjugate[#1] == #2] & ]; 
     fcoef*Times @@ Table[OptionValue[SimplifyFunction][Piecewise[{{x - r, Element[r, Reals]}, {(x - Re[r])^2 + Im[r]^2, True}}]], 
        {r, roots}]]; 
RealFactor[poly_, opts:OptionsPattern[]] := RealFactor[poly, ExpressionPivot[poly, Symbol], opts]; 


Attributes[RealApart] := {Listable}; 
Options[RealApart] := Options[RealFactor]; 
RealApart[poly_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[poly, x] := 
   Apart[Numerator[poly]/RealFactor[Denominator[poly], opts]]; 
RealApart[poly_, opts:OptionsPattern[]] := RealApart[poly, ExpressionPivot[poly, Symbol], opts]; 


(* ::Subsection:: *)
(*CorrectionTest*)


CorrectionTest::overtime = "`` \:8ba1\:7b97\:8d85\:65f6"; 
Options[CorrectionTest] := {Skip -> {}, TimeConstraint -> Infinity, VerifyDiscontinuities -> {}, WorkingPrecision -> $MinPrecision}; 
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
        Message[CorrectionTest::overtime, "\:56fe\:50cf"]; ]; ]; Print[Grid[table, Alignment -> {Center, Center}, Spacings -> {1, 1}, 
       ItemSize -> {Scaled[0.5], Automatic}, Frame -> All]]; ]; 


(* ::Subsection:: *)
(*FindIdentities*)


FindIdentities[expr1_, expr2_, x_Symbol] /; RationalExpressionQ[expr1, x] && RationalExpressionQ[expr2, x] := 
   Module[{p1, p2, roots, limit}, p1 = Numerator[Simplify[expr1]]*Denominator[Simplify[expr2]]; 
     p2 = Numerator[Simplify[expr2]]*Denominator[Simplify[expr1]]; roots = SolveValues[D[p1/p2, x] == 0, x]; If[ !ListQ[roots], Return[{}]]; 
     DeleteDuplicates[Flatten[Reap[Do[limit = Simplify[Limit[p1/p2, x -> i]]; If[limit =!= 0 && Element[limit, Rationals], 
            Sow[Defer[Evaluate[p1]] - limit*p2 == Factor[p1 - limit*p2]]]; , {i, roots}]][[2]]]]]; 


(* ::Subsection:: *)
(*BivariablePlot*)


Options[BivariablePlot] := {PlotLabels -> None}; 
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


(* ::Section:: *)
(*Experimental*)


(* ::Subsection:: *)
(*Li2Transform*)


Li2Transform[z_, f_Function] /; Module[{t, ft}, ft = Together[f[t]]; PolynomialQ[ft, t] && CoefficientList[ft, t] === {1, -1}] := 
   -PolyLog[2, 1 - z] - Log[z]*Log[1 - z] + Pi^2/6; 
Li2Transform[z_, f_Function] /; Module[{t, ft}, ft = Together[f[t]]; RationalExpressionQ[ft, t] && CoefficientList[Denominator[ft], t] === 
        {0, 1} && CoefficientList[Numerator[ft], t] === {1}] := -PolyLog[2, 1/z] - Log[-z]*Log[z] + Log[z]^2/2 + Pi^2/3; 
Li2Transform[z_, f_Function] /; Module[{t, ft}, ft = Together[f[t]]; RationalExpressionQ[ft, t] && CoefficientList[Denominator[ft], t] === 
        {0, 1} && CoefficientList[Numerator[ft], t] === {-1, 1}] := PolyLog[2, 1 - 1/z] + Log[1 - 1/z]*Log[1/z] - Log[-z]*Log[z] + Log[z]^2/2 + 
    Pi^2/6; 
Li2Transform[z_, f_Function] /; Module[{t, ft}, ft = Together[f[t]]; RationalExpressionQ[ft, t] && CoefficientList[Denominator[ft], t] === 
        {1, -1} && CoefficientList[Numerator[ft], t] === {1}] := PolyLog[2, 1/(1 - z)] + (1/2)*Log[1 - z]^2 - Log[-z]*Log[1 - z] - Pi^2/6; 
Li2Transform[z_, f_Function] /; Module[{t, ft}, ft = Together[f[t]]; RationalExpressionQ[ft, t] && CoefficientList[Denominator[ft], t] === 
        {-1, 1} && CoefficientList[Numerator[ft], t] === {0, 1}] := -PolyLog[2, z/(-1 + z)] + (1/2)*Log[1 - z]^2 - Log[1 - z]*Log[-z] - 
    Log[1/(1 - z)]*Log[z/(-1 + z)]; 
Li2Transform[z1_, z2_] /; PossibleZeroQ[z1 + z2] := (1/2)*PolyLog[2, z1^2]; 


(* ::Subsection:: *)
(*Ti2Transform*)


Options[Ti2Transform] := {Defer -> True}; 
Ti2Transform[z_, "\:5012\:6570", OptionsPattern[]] := Module[{ti2 = If[OptionValue[Defer], Defer, Identity][Ti2]}, 
    ti2[1/z] + Sign[z]*(Pi/2)*Log[Abs[z]]]; 
Ti2Transform[z_, a_, "\:5012\:6570", OptionsPattern[]] := Module[{ti2 = If[OptionValue[Defer], Defer, Identity][Ti2]}, 
    -ti2[1/z, 1/a] + ti2[z] - ti2[a] + ArcTan[a]*Log[Abs[a]] - Sign[a]*(Pi/2)*Log[(z*Sqrt[a^2 + 1])/(z + a)]]; 
Ti2Transform[z_, a_, "\:4ea4\:6362", OptionsPattern[]] := Module[{ti2 = If[OptionValue[Defer], Defer, Identity][Ti2]}, 
    ti2[a, z] + ti2[z] - ti2[a] + ArcTan[a]*Log[(a*Sqrt[z^2 + 1])/(z + a)] - ArcTan[z]*Log[(z*Sqrt[a^2 + 1])/(z + a)]]; 


End[];


EndPackage[];
