(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Mirion`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*Usage*)


(* ::Section:: *)
(*Public*)


TIMING::usage = "\
TIMING[code] \:8fd4\:56de\:6267\:884c code \:6240\:82b1\:8d39\:7684\:65f6\:95f4\:548c\:7ed3\:679c.
TIMING[code,n] \:8fd4\:56de\:6267\:884c code n \:6b21\:6240\:82b1\:8d39\:7684\:65f6\:95f4\:548c\:7ed3\:679c.";


ExpressionPivot::usage = "\
ExpressionPivot[expr] \:8fd4\:56de expr \:4e2d\:9996\:4e2a\:975e\:6570\:503c\:7b26\:53f7.";


CoefficientSeparation::usage = "\
CoefficientSeparation[expr,x] \:8fd4\:56de expr \:7cfb\:6570\:548c\:4f59\:4e0b\:90e8\:5206.";


FirstCoefficient::usage = "\
FirstCoefficient[poly,x] \:8fd4\:56de poly \:6700\:9ad8\:6b21\:9879\:7684\:7cfb\:6570.";


PolynomialRoots::usage = "\
PolynomialRoots[poly,x] \:8fd4\:56de poly \:7684\:6839.";


PolynomialRootApproximant::usage = "\
PolynomialRootApproximant[poly,x] \:8fd4\:56de poly \:7684\:6839\:8fd1\:4f3c.";


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
ApartArcTan[expr,x] \:7ed9\:51fa arctan(expr) \:7684\:88c2\:9879."; 


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


ContinuedFractionExpand::usage = "\
ContinuedFractionExpand[f,{x,n}] \:7ed9\:51fa\:51fd\:6570 f \:7684\:524d n \:9879\:8fde\:5206\:6570\:5c55\:5f00.
";


FromContinuedFractionExpand::usage = "\
FromContinuedFractionExpand[list] \:6839\:636e list \:6784\:9020\:8fde\:5206\:6570.
";


ContinuedFractionExpandPeriod::usage = "\
ContinuedFractionExpandPeriod[\!\(\*SqrtBox[\(r\)]\),x] \:7ed9\:51fa\:51fd\:6570 \!\(\*SqrtBox[\(r\)]\) \:7684\:8fde\:5206\:6570\:5c55\:5f00\:6700\:5c0f\:5468\:671f.
";


PolynomialFit::usage = "\
PolynomialFit[expr,x,n] \:7ed9\:51fa\:8868\:8fbe\:5f0f expr \:7684 n \:6b21\:591a\:9879\:5f0f\:8fd1\:4f3c.
"


PolynomialReverse::usage = "\
PolynomialReverse[poly,x] \:7ed9\:51fa\:591a\:9879\:5f0f poly \:7684\:53cd\:5411\:591a\:9879\:5f0f.
"


IntegrateCF::usage = "\
IntegrateCF[expr,x] \:4f7f\:7528\:8fde\:5206\:6570\:5c55\:5f00\:6cd5\:6c42 expr \:5173\:4e8e x \:7684\:79ef\:5206.
"


(* ::Section:: *)
(*Experimental*)


FindIdentities::usage = "\
FindIdentities[expr1,expr2,x] \:ff08\:5b9e\:9a8c\:6027\:ff09\:7ed9\:51fa\:5173\:4e8e expr1,expr2 \:7684\:6052\:7b49\:5f0f.
";


BivariablePlot::usage = "\
BivariablePlot[list,x] \:ff08\:5b9e\:9a8c\:6027\:ff09\:7ed8\:5236\:591a\:5143 list \:7684\:5173\:7cfb\:56fe.
";


(* ::Section:: *)
(*Option*)


GenerateConstant::usage = "GenerateConstant \:5173\:4e8e\:662f\:5426\:751f\:6210\:5e38\:6570\:7684\:9009\:9879"; 
SimplifyFunction::usage = "SimplifyFunction \:5173\:4e8e\:5316\:7b80\:51fd\:6570\:7684\:9009\:9879"; 


(* ::Subsection::Closed:: *)
(*BeginPrivate*)


Begin["`Private`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*Definition*)


(* ::Section:: *)
(*Public*)


(* ::Subsection::Closed:: *)
(*TIMING*)


Attributes[TIMING] = {HoldAll}; 
TIMING[code_, n_Integer:1] := Module[{}, ClearSystemCache[]; AbsoluteTiming[Do[code, n - 1]; code]]


(* ::Subsection::Closed:: *)
(*ExpressionPivot*)


Attributes[ExpressionPivot] = {Listable}; 
ExpressionPivot[expr_] := FirstCase[expr, _Symbol?(Not @* NumericQ), Symbol, {-1}]; 


(* ::Subsection::Closed:: *)
(*CoefficientSeparation*)


Attributes[CoefficientSeparation] = {Listable}; 
CoefficientSeparation[expr_, x_Symbol] /; FreeQ[expr, x] := {expr, 1}; 
CoefficientSeparation[expr_, x_Symbol] := Replace[expr, Longest[c_.]*(r_) /; FreeQ[c, x] -> {c, r}]; 


(* ::Subsection::Closed:: *)
(*FirstCoefficient*)


Attributes[FirstCoefficient] = {Listable}; 
FirstCoefficient[poly_, x_Symbol] := Coefficient[poly, x, Exponent[poly, x]]; 


(* ::Subsection::Closed:: *)
(*PolynomialRoots*)


Attributes[PolynomialRoots] = {Listable}; 
PolynomialRoots[a_, x_Symbol, OptionsPattern[]] /; FreeQ[a, x] := {}; 
PolynomialRoots[(a_.)*(x_) + (b_.), x_Symbol, OptionsPattern[]] /; FreeQ[{a, b}, x] := {-(b/a)}; 
PolynomialRoots[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := List @@ (Roots[poly == 0, x, Cubics -> False, Quartics -> False] /. 
     Equal -> (#2 & )); 


(* ::Subsection::Closed:: *)
(*PolynomialRootApproximant*)


Attributes[PolynomialRootApproximant] = {Listable}; 
PolynomialRootApproximant[poly_, x_Symbol] /; PolynomialQ[poly, x] := FromDigits[RootApproximant /@ Reverse[CoefficientList[poly, x]], x]; 


(* ::Subsection::Closed:: *)
(*IBP*)


IBP[u_, v_, x_Symbol, OptionsPattern[]] := Module[{coef, rem}, {coef, rem} = CoefficientSeparation[v*D[u, x], x]; u*v - coef*Inactive[Integrate][rem, x]]; 
IBP[f_, x_Symbol, OptionsPattern[]] := IBP[f, x, x]; 
IBP[u_, v_, {x_Symbol, a_, b_}, OptionsPattern[]] := Module[{assum = OptionValue[Assumptions], c, r}, {c, r} = CoefficientSeparation[v*D[u, x], x]; 
     Limit[u*v, x -> b, Assumptions -> assum] - Limit[u*v, x -> a, Assumptions -> assum] - c*Inactive[Integrate][r, {x, a, b}]]; 
IBP[f_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := IBP[f, x, {x, a, b}, opts]; 


(* ::Subsection::Closed:: *)
(*IBS*)


Options[IBS] = {Assumptions -> $Assumptions, Inactive -> False}; 
IBS[f_, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{assum = OptionValue[Assumptions], inactive = OptionValue[Inactive], u, result, c, r}, 
    result = If[ex === x, f /. x -> et, IntegrateChangeVariables[Inactive[Integrate][f, x], u, u == ex, Assumptions -> assum][[1]] /. C[_] -> 0 /. u -> et]*
       D[et, t]; If[inactive, {c, r} = CoefficientSeparation[Together[result], t]; c*Inactive[Integrate][r, t], result]]; 
IBS[f_, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, x -> et, x, t, opts]; 
IBS[f_, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, ex -> t, x, t, opts]; 
IBS[f_, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, x_Symbol, t_Symbol, OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{assum = OptionValue[Assumptions], u, temp, result, aa, bb, c, r}, 
    temp = If[ex === x, Inactive[Integrate][f, {x, a, b}] /. x -> u, IntegrateChangeVariables[Inactive[Integrate][f, {x, a, b}], u, u == ex, 
        Assumptions -> assum]]; {result, {t, aa, bb}} = List @@ If[et === t, temp /. u -> t, IntegrateChangeVariables[temp, t, u == et, Assumptions -> assum]]; 
     {c, r} = CoefficientSeparation[Together[result], t]; c*Inactive[Integrate][r, {t, a, b}]]; 
IBS[f_, {a_, b_}, x_Symbol -> et_, t_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, x -> et, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> t_Symbol, x_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> t, x, t, opts]; 
IBS[f_, {a_, b_}, ex_ -> et_, opts:OptionsPattern[]] := IBS[f, {a, b}, ex -> et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 


(* ::Subsection::Closed:: *)
(*ApartArcTan*)


Attributes[ApartArcTan] = {Listable}; 
Options[ApartArcTan] = {GenerateConstant -> True, SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := Module[{result, count = 0, poles, sample, c = {}}, 
    result = Sum[ArcTan[OptionValue[SimplifyFunction][(x - Re[r])/((count += Sign[#1]; #1) & )[Im[r]]]], {r, SolveValues[expr - I == 0, x]}]; 
     If[OptionValue[GenerateConstant], result += Block[{nd, deg, ratio}, nd = NumeratorDenominator[Together[expr]]; deg = Exponent[nd, x]; 
           If[Subtract @@ deg < 0, 0, ratio = Divide @@ Coefficient[nd, x, deg]; If[Subtract @@ deg > 0, Sign[ratio]*(Pi/2), ArcTan[ratio]]]] - count*(Pi/2); 
       poles = DeleteDuplicates[SolveValues[Denominator[Together[expr]] == 0, x, Reals]]; If[Length[poles] =!= 0, 
        sample = (((1/2)*Table[#1[[i]] + #1[[i + 1]], {i, Length[#1] - 1}] & )[Prepend[#1, First[#1] - 1]] & )[N[poles]]; 
         c = Pi*Round[(1/Pi)*Table[ArcTan[expr] - result /. x -> i, {i, sample}]]; 
         result += Piecewise[Prepend[Table[{c[[i + 1]], poles[[i]] < x < poles[[i + 1]]}, {i, Length[poles] - 1}], {First[c], x < First[poles]}]]]; ]; result]; 


(* ::Subsection::Closed:: *)
(*ComplexFactor*)


Attributes[ComplexFactor] = {Listable}; 
Options[ComplexFactor] = {ToRadicals -> False}; 
ComplexFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := FirstCoefficient[poly, x]*
    Times @@ (x - If[OptionValue[ToRadicals], ToRadicals, Identity][PolynomialRoots[poly, x]]); 
ComplexFactor[poly_, opts:OptionsPattern[]] := ComplexFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection::Closed:: *)
(*ComplexApart*)


Attributes[ComplexApart] = {Listable}; 
Options[ComplexApart] = Options[ComplexFactor]; 
ComplexApart[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/ComplexFactor[#2, x, ToRadicals -> OptionValue[ToRadicals]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
ComplexApart[poly_, opts:OptionsPattern[]] := ComplexApart[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection::Closed:: *)
(*RealFactor*)


Attributes[RealFactor] = {Listable}; 
Options[RealFactor] = {SimplifyFunction -> ToRadicals /* ComplexExpand /* FunctionExpand /* Simplify}; 
RealFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := FirstCoefficient[poly, x]*
   Times @@ Table[OptionValue[SimplifyFunction][Piecewise[{{x - r, Element[r, Reals]}, {(x - Re[r])^2 + Im[r]^2, True}}]], 
     {r, DeleteDuplicates[PolynomialRoots[poly, x], N[Conjugate[#1] == #2] & ]}]
RealFactor[poly_, opts:OptionsPattern[]] := RealFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection::Closed:: *)
(*RealApart*)


Attributes[RealApart] = {Listable}; 
Options[RealApart] = Options[RealFactor]; 
RealApart[expr_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/RealFactor[#2, x, SimplifyFunction -> OptionValue[SimplifyFunction]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
RealApart[expr_, opts:OptionsPattern[]] := RealApart[expr, ExpressionPivot[expr], opts]; 


(* ::Subsection::Closed:: *)
(*ContinuedFractionExpand*)


ContinuedFractionExpand[f_, {x_Symbol, n_Integer}] := Module[{a, b = f}, Quiet[Table[a = Normal[Series[b, {x, Infinity, 0}]]; b = FullSimplify[1/(b - a)]; a, 
      {i, n}], Series::sbyc]]; 


(* ::Subsection::Closed:: *)
(*FromContinuedFractionExpand*)


FromContinuedFractionExpand[list_List] := Fold[#2 + 1/#1 & , Reverse[list]]; 


(* ::Subsection::Closed:: *)
(*ContinuedFractionExpandPeriod*)


Options[ContinuedFractionExpandPeriod] = {MaxIterations -> 100}; 
ContinuedFractionExpandPeriod::lim = "Iteration limit of `1` exceeded."; 
ContinuedFractionExpandPeriod[Sqrt[rad_], x_Symbol, OptionsPattern[]] /; PolynomialQ[rad, x] && SquareFreeQ[rad, x] := 
   Module[{lim = OptionValue[MaxIterations], r, n, a, b, c}, n = Exponent[rad, x]; r = PolynomialReverse[rad, x, n]; {a, b} = {0, 1}; 
     Reap[Do[If[i == lim, Message[ContinuedFractionExpandPeriod::lim, lim]; Abort[]]; c = Normal[Series[a + b*Sqrt[r], {x, 0, n/2}]]; 
         If[Head[c] === Piecewise, c = SelectFirst[c[[1,All,1]], Not @* PossibleZeroQ, c[[2]]]; ]; If[i > 1 &&  !PossibleZeroQ[Coefficient[c, x, 0]], Break[], 
          Sow[PolynomialReverse[c, x, n/2]]]; {a, b} = Simplify[(x^n*{a - c, -b})/((a - c)^2 - b^2*r)]; , {i, Infinity}]][[2,1]]]; 


(* ::Subsection::Closed:: *)
(*PolynomialFit*)


Options[PolynomialFit] = {InterpolatingFunction -> (RandomReal[#1, WorkingPrecision -> 100] & )}; 
PolynomialFit[expr_, x_Symbol, n_Integer, OptionsPattern[]] := Module[{intf = OptionValue[InterpolatingFunction]}, 
    InterpolatingPolynomial[Table[({#1, expr /. x -> #1} & )[intf[i]], {i, n + 1}], x]]; 


(* ::Subsection::Closed:: *)
(*PolynomialReverse*)


PolynomialReverse[poly_, x_Symbol] /; PolynomialQ[poly, x] := FromDigits[CoefficientList[poly, x], x]; 
PolynomialReverse[poly_, x_Symbol, n_Integer] /; PolynomialQ[poly, x] := FromDigits[PadRight[CoefficientList[poly, x], n + 1], x]; 


(* ::Subsection:: *)
(*IntegrateCF*)


Options[IntegrateCFImpl] = Options[IntegrateCF] = {WorkingPrecision -> Automatic, Evaluate -> False, Sequence @@ Options[ContinuedFractionExpandPeriod], 
     Sequence @@ Options[PolynomialFit]}; 
IntegrateCFReplace /: MakeBoxes[IntegrateCFReplace[expr_, repl_], form_] := 
    (InterpretationBox[FrameBox[RowBox[{#1, "/.", #2}], RoundingRadius -> 5], expr /. repl] & ) @@ {MakeBoxes[expr, form], MakeBoxes[repl, form]}; 
IntegrateCFImpl[poly_, rad_, x_, opts:OptionsPattern[]] := Module[{prec = OptionValue[WorkingPrecision], 
     opts1 = FilterRules[{opts}, Options[ContinuedFractionExpandPeriod]], opts2 = FilterRules[{opts}, Options[PolynomialFit]], r, g, p, q, alg, coef}, 
    If[prec === Automatic, prec = If[FreeQ[{poly, rad}, Root], Infinity, 100]]; r = N[rad, prec]; g = Exponent[r, x]/2 - 1; 
     {p, q} = NumeratorDenominator[Together[FromContinuedFractionExpand[ContinuedFractionExpandPeriod[Sqrt[r], x, opts1]]]]; 
     alg = If[prec === Infinity, Cancel, PolynomialFit[#1, x, g, opts2] & ][D[p, x]/q]; coef = Coefficient[poly, x, g]/Coefficient[alg, x, g]; 
     If[prec === Infinity, Simplify, PolynomialRootApproximant[#1, x] & ][{coef, p, q, poly - coef*alg}]]; 
IntegrateCF[0, x_Symbol] := 0; 
IntegrateCF[(poly_)/Sqrt[rad_], x_Symbol, opts:OptionsPattern[]] /; PolynomialQ[poly, x] && PolynomialQ[rad, x] && SquareFreeQ[rad, x] := 
   Module[{n, g}, n = Exponent[rad, x]; g = n/2 - 1; (#1*Log[#2 + #3*Sqrt[rad]] + IntegrateCF[#4/Sqrt[rad], x] & ) @@ IntegrateCFImpl[poly, rad, x, opts] /; 
       !PossibleZeroQ[Coefficient[poly, x, g]]]; 
IntegrateCF[(rat_)/Sqrt[rad_], x_Symbol, opts:OptionsPattern[]] /;  !PolynomialQ[rat, x] && RationalExpressionQ[rat, x] && PolynomialQ[rad, x] && 
     SquareFreeQ[rad, x] := Module[{eval = OptionValue[Evaluate], nume, deno, n, g, alg, trans, coef, r}, 
    {nume, deno} = NumeratorDenominator[Together[rat]]; n = Exponent[rad, x]; g = Ceiling[n/2] - 1; {alg, trans} = {0, PolynomialQuotient[nume, deno, x]}; 
     Do[coef = Residue[rat, {x, y}]; r = Collect[Sum[x^(2*(g + 1) - k)*Coefficient[rad, x, k]*(1 + y*x)^k, {k, 0, n}], x]; 
       {alg, trans} += ({(-RootReduce[coef*#1])*If[eval, ReplaceAll, IntegrateCFReplace][Log[#2 + #3*Sqrt[r]], x -> 1/(x - y)], 
           coef*Sum[Coefficient[#4, x, k]*(x - y)^(g - 1 - k), {k, 0, g - 1}]} & ) @@ IntegrateCFImpl[x^g, r, x, opts]; , {y, PolynomialRoots[deno, x]}]; 
     alg + IntegrateCF[RootReduce[Collect[trans, x]]/Sqrt[rad], x]]; 


(* ::Section:: *)
(*Experimental*)


(* ::Subsection::Closed:: *)
(*FindIdentities*)


FindIdentities[expr1_, expr2_, x_Symbol] /; RationalExpressionQ[expr1, x] && RationalExpressionQ[expr2, x] := 
   Module[{p1, p2, roots, limit}, p1 = Numerator[Simplify[expr1]]*Denominator[Simplify[expr2]]; p2 = Numerator[Simplify[expr2]]*Denominator[Simplify[expr1]]; 
     roots = DeleteDuplicates[SolveValues[D[p1, x]*p2 == D[p2, x]*p1, x]]; 
     DeleteDuplicates[Reap[Do[limit = Simplify[Limit[p1/p2, x -> i]]; If[limit =!= 0 && Element[limit, Rationals], 
           Sow[Defer[Evaluate[p1]] - limit*p2 == Factor[p1 - limit*p2]]]; , {i, roots}]][[2,1]]]]; 


(* ::Subsection::Closed:: *)
(*BivariablePlot*)


Options[BivariablePlot] = {PlotLabels -> None}; 
BivariablePlot[list_List, x_Symbol, OptionsPattern[]] := Module[{labels, isValid, const, constQ, edges, vertexes, usedVertexes, edgeStylize, vertexStylize}, 
    labels = OptionValue[PlotLabels]; isValid = labels =!= None && Length[list] === Length[labels]; constQ[value_] := FreeQ[value, x] &&  !PossibleZeroQ[value]; 
     edgeStylize[edge_, value_] := Labeled[edge, Placed[Style[value, 14, FontFamily -> "Times"], 0.5]]; 
     vertexStylize[value_] := Framed[Style[value, Black, 14, FontFamily -> "Times"], FrameStyle -> None]; 
     edges = Table[Piecewise[{{{UndirectedEdge[bivars[[1]], bivars[[2]]], const}, constQ[const = Simplify[list[[bivars[[1]]]]^2 + list[[bivars[[2]]]]^2]]}, 
         {Piecewise[{{{DirectedEdge[bivars[[2]], bivars[[1]]], -const}, TrueQ[const < 0]}, {{DirectedEdge[bivars[[1]], bivars[[2]]], const}, True}}], 
          constQ[const = Simplify[list[[bivars[[1]]]]^2 - list[[bivars[[2]]]]^2]]}, {Nothing, True}}], {bivars, Subsets[Range[Length[list]], {2}]}]; 
     usedVertexes = DeleteDuplicates[Cases[edges[[All,1]], _Integer, {2}]]; 
     vertexes = Table[i -> (Evaluate[Inset[vertexStylize[Piecewise[{{labels[[i]] == list[[i]], isValid}, {list[[i]], True}}]], #1]] & ), {i, usedVertexes}]; 
     edges = edgeStylize @@@ edges; Graph[edges, ImageSize -> Medium, PlotTheme -> "DiagramBlack", VertexLabels -> None, VertexShapeFunction -> vertexes, 
      PlotRangePadding -> Scaled[0.1]]]; 


(* ::Subsection::Closed:: *)
(*EndPrivate*)


End[];


(* ::Subsection::Closed:: *)
(*EndPackage*)


EndPackage[];
