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
(*CheckVersion*)


CheckVersion::vtl = "Mathematica \:7248\:672c\:8fc7\:4f4e. \:5f53\:524d\:7684\:7248\:672c\:4e3a ``, \:8981\:6c42\:7684\:7248\:672c\:4e3a ``."; 
CheckVersion[(v_)?NumericQ] := If[TrueQ[$VersionNumber < v], Message[CheckVersion::vtl, $VersionNumber, v]]; 


(* ::Subsection:: *)
(*ExpressionPivot*)


Attributes[ExpressionPivot] := {Listable}; 
ExpressionPivot[expr_] := FirstCase[expr, _Symbol, False, {0, Infinity}]; 


(* ::Subsection:: *)
(*CoefficientSeparation*)


Attributes[CoefficientSeparation] := {Listable}; 
CoefficientSeparation[expr_, x_Symbol] := Module[{coef}, If[FreeQ[expr, x], Return[{expr, 1}]]; 
     Simplify[({#1, expr/#1} & )[Replace[(c_.)*_ /; FreeQ[c, x] -> c][Simplify[expr]]]]]; 


(* ::Section:: *)
(*Public*)


(* ::Subsection:: *)
(*IBP*)


IBP[u_, v_, x_Symbol] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; 
     Simplify[u*v - coef*Inactive[Integrate][rem, x]]]; 
IBP[f_, x_Symbol] := IBP[f, x, x]; 
IBP[u_, v_, {x_Symbol, a_, b_}] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; 
     Simplify[Limit[u*v, x -> b] - Limit[u*v, x -> a] - coef*Inactive[Integrate][rem, {x, a, b}]]]; 
IBP[f_, {x_Symbol, a_, b_}] := IBP[f, x, {x, a, b}]; 


(* ::Subsection:: *)
(*IBS*)


Options[IBS] := {Assumptions -> $Assumptions}; 
IBS[f_, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; x =!= t && FreeQ[ex, t] && FreeQ[et, x] := 
   Module[{assum = OptionValue[Assumptions], u}, CheckVersion[13.1]; 
     Simplify[If[ex === x, f /. x -> et, IntegrateChangeVariables[Inactive[Integrate][f, x], u, u == ex, Assumptions -> assum][[1]] /. u -> et]*
        D[et, t] /. C[_] -> 0, assum]]; 
IBS[f_, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] /; x =!= t && FreeQ[et, x] := IBS[f, x -> et, x, t, opts]; 
IBS[f_, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] /; x =!= t && FreeQ[ex, t] := IBS[f, ex -> t, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; x =!= t && FreeQ[ex, t] && FreeQ[et, x] := 
   Module[{assum = OptionValue[Assumptions], u, temp}, CheckVersion[13.1]; temp = If[ex === x, Inactive[Integrate][f, {x, a, b}] /. x -> u, 
       IntegrateChangeVariables[Inactive[Integrate][f, {x, a, b}], u, u == ex, Assumptions -> assum]]; 
     Simplify[If[et === t, temp /. u -> t, IntegrateChangeVariables[temp, t, u == et, Assumptions -> assum]], assum]]; 
IBS[f_, {a_, b_}, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] /; x =!= t && FreeQ[et, x] := IBS[f, {a, b}, x -> et, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] /; x =!= t && FreeQ[ex, t] := IBS[f, {a, b}, ex -> t, x, t, opts]; 


(* ::Subsection:: *)
(*ApartArcTan*)


Attributes[ApartArcTan] := {Listable}; 
Options[ApartArcTan] := {GenerateConstant -> False, SimplifyFunction -> Simplify}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   Module[{simp = OptionValue[SimplifyFunction], gene = OptionValue[GenerateConstant], roots, result, zeros, offset}, 
    roots = x /. Solve[expr == I, x]; result = Sum[ArcTan[simp[(x - Re[r])/Im[r]]], {r, roots}]; If[ !gene, Return[result]]; 
     zeros = x /. DeleteDuplicates[Solve[Denominator[Together[expr]] == 0, x]]; 
     result += If[Length[zeros] === 0, ArcTan[expr] - result /. x -> 0, 
       offset = Append[Table[Limit[ArcTan[expr], x -> i, Direction -> "FromBelow"] - (result /. x -> i), {i, zeros}], 
          Limit[ArcTan[expr], x -> Last[zeros], Direction -> "FromAbove"] - (result /. x -> Last[zeros])]; 
        Piecewise[Append[Prepend[Table[{offset[[i + 1]], zeros[[i]] < x < zeros[[i + 1]]}, {i, Length[zeros] - 1}], 
           {First[offset], x < First[zeros]}], {Last[offset], x > Last[zeros]}]]]; Simplify[result]]; 
ApartArcTan[expr_, opts:OptionsPattern[]] := ApartArcTan[expr, ExpressionPivot[expr], opts]; 


(* ::Subsection:: *)
(*RealFactor*)


Attributes[RealFactor] := {Listable}; 
Options[RealFactor] := {SimplifyFunction -> RootReduce /* ToRadicals /* Simplify}; 
RealFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := Module[{simp = OptionValue[SimplifyFunction], roots, result = 1}, 
    If[FreeQ[poly, x], Return[poly]]; roots = DeleteDuplicates[x /. Solve[poly == 0, x], N[Conjugate[#1] == #2] & ]; 
     Coefficient[poly, x, Exponent[poly, x]]*
      Times @@ Table[simp[Piecewise[{{x - root, Element[root, Reals]}, {(x - Re[root])^2 + Im[root]^2, True}}]], {root, roots}]]; 
RealFactor[poly_, opts:OptionsPattern[]] := RealFactor[poly, ExpressionPivot[poly], opts]; 


Attributes[RealApart] := {Listable}; 
Options[RealApart] := Options[RealFactor]; 
RealApart[poly_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[poly, x] := Apart[Numerator[poly]/RealFactor[Denominator[poly], opts]]; 
RealApart[poly_, opts:OptionsPattern[]] := RealApart[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection:: *)
(*CorrectionTest*)


CorrectionTest::overtime = "`` \:8ba1\:7b97\:8d85\:65f6"; 
Options[CorrectionTest] := {Skip -> {}, TimeConstraint -> Infinity, VerifyDiscontinuities -> {}, WorkingPrecision -> $MinPrecision}; 
CorrectionTest[org_, res_, {x_, xmin_, xmax_}, OptionsPattern[]] := 
   Module[{table = {{Style["\:5b9a\:4e49\:57df", Bold], Style["\:95f4\:65ad\:70b9", Bold]}, {Style[True, Bold, RGBColor["#4081ff"]], 
        Style[False, Bold, RGBColor["#4081ff"]]}, {Style[True, Bold, RGBColor["#eb7100"]], Style[False, Bold, RGBColor["#eb7100"]]}, 
       {Style["\:56fe\:50cf", Bold], SpanFromLeft}, {Null, SpanFromLeft}}, skip = OptionValue[Skip], cons = OptionValue[TimeConstraint], 
     veri = OptionValue[VerifyDiscontinuities], prec = OptionValue[WorkingPrecision], dom = True}, 
    CheckVersion[12.2]; If[ !ListQ[skip], skip = {skip}]; If[ !ListQ[veri], veri = {veri}]; 
     If[ !AssociationQ[cons], cons = Association["\:5b9a\:4e49\:57df" -> cons, "\:95f4\:65ad\:70b9" -> cons, "\:56fe\:50cf" -> cons]]; 
     If[ !MemberQ[skip, "\:5b9a\:4e49\:57df"], TimeConstrained[table[[2,1]] = Item[Style[dom = Reduce[FunctionDomain[org, x], Element[x, Reals]], Bold, 
            RGBColor["#4081ff"]], ItemSize -> Automatic]; table[[3,1]] = Item[Style[Reduce[FunctionDomain[res, x] && dom, Element[x, Reals]], 
            Bold, RGBColor["#eb7100"]], ItemSize -> Automatic]; , cons["\:5b9a\:4e49\:57df"], Message[CorrectionTest::overtime, "\:5b9a\:4e49\:57df"]; ]; ]; 
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


ClearAll[FindIdentities]; 
FindIdentities[expr1_, expr2_, x_Symbol] /; RationalExpressionQ[expr1, x] && RationalExpressionQ[expr2, x] := 
   Module[{p1, p2, roots, limit}, p1 = Numerator[Simplify[expr1]]*Denominator[Simplify[expr2]]; 
     p2 = Numerator[Simplify[expr2]]*Denominator[Simplify[expr1]]; roots = x /. Solve[D[p1/p2, x] == 0, x]; If[ !ListQ[roots], Return[{}]]; 
     DeleteDuplicates[Flatten[Reap[Do[limit = Simplify[Limit[p1/p2, x -> i]]; If[limit =!= 0 && Element[limit, Rationals], 
            Sow[Defer[Evaluate[p1]] - limit*p2 == Factor[p1 - limit*p2]]]; , {i, roots}]][[2]]]]]; 


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
