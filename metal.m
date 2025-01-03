(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*[BeginPackage]*)


BeginPackage["Metal`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*\:8bf4\:660e*)


(* ::Section::Closed:: *)
(*\:5b9e\:7528\:51fd\:6570*)


TimingTest::usage = "\
TimingTest[code] \:8fd4\:56de\:6267\:884c code \:6240\:82b1\:8d39\:7684\:65f6\:95f4\:548c\:7ed3\:679c.
TimingTest[code,n] \:8fd4\:56de\:6267\:884c code n \:6b21\:6240\:82b1\:8d39\:7684\:65f6\:95f4\:548c\:7ed3\:679c.";


PassOptions::usage = "\
PassOptions[from,to,opts] \:5728 from \:548c to \:95f4\:4f20\:9012 opts.";


(* ::Section::Closed:: *)
(*\:8868\:8fbe\:5f0f*)


AlgebraicExpressionQ::usage = "\
AlgebraicExpressionQ[expr,x] \:8fd4\:56de expr \:662f\:5426\:4e3a\:4ee3\:6570\:8868\:8fbe\:5f0f.";


ExpressionPivot::usage = "\
ExpressionPivot[expr] \:8fd4\:56de expr \:4e2d\:9996\:4e2a\:975e\:6570\:503c\:7b26\:53f7.";


CoefficientSeparation::usage = "\
CoefficientSeparation[expr,x] \:8fd4\:56de expr \:7cfb\:6570\:548c\:4f59\:4e0b\:90e8\:5206.";


(* ::Section::Closed:: *)
(*\:591a\:9879\:5f0f*)


FirstCoefficient::usage = "\
FirstCoefficient[poly,x] \:8fd4\:56de poly \:6700\:9ad8\:6b21\:9879\:7684\:7cfb\:6570.";


PolynomialReverse::usage = "\
PolynomialReverse[poly,x] \:7ed9\:51fa poly \:7684\:53cd\:5411\:591a\:9879\:5f0f.";


PolynomialRoots::usage = "\
PolynomialRoots[poly,x] \:8fd4\:56de poly \:7684\:6839.";


PolynomialRootApproximant::usage = "\
PolynomialRootApproximant[poly,x] \:8fd4\:56de poly \:7684\:6839\:8fd1\:4f3c.";


ComplexFactor::usage = "\
ComplexFactor[poly] \:7ed9\:51fa poly \:5728\:590d\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3.
ComplexFactor[poly,x] \:7ed9\:51fa poly \:5173\:4e8e x \:5728\:590d\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3."; 


RealFactor::usage = "\
RealFactor[poly] \:7ed9\:51fa poly \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3.
RealFactor[poly,x] \:7ed9\:51fa poly \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:56e0\:5f0f\:5206\:89e3."; 


(* ::Section::Closed:: *)
(*\:521d\:7b49\:51fd\:6570*)


PolynomialFit::usage = "\
PolynomialFit[expr,x,n] \:7ed9\:51fa expr \:7684 n \:6b21\:591a\:9879\:5f0f\:8fd1\:4f3c.";


ApartArcTan::usage = "\
ApartArcTan[expr] \:7ed9\:51fa arctan(expr) \:7684\:88c2\:9879.
ApartArcTan[expr,x] \:7ed9\:51fa arctan(expr) \:5173\:4e8e x \:7684\:88c2\:9879."; 


ComplexApart::usage = "\
ComplexApart[expr] \:7ed9\:51fa expr \:5728\:590d\:6570\:57df\:4e0a\:7684\:88c2\:9879.
ComplexApart[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:590d\:6570\:57df\:4e0a\:7684\:88c2\:9879."; 


RealApart::usage = "\
RealApart[expr] \:7ed9\:51fa expr \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879.
RealApart[expr,x] \:7ed9\:51fa expr \:5173\:4e8e x \:5728\:5b9e\:6570\:57df\:4e0a\:7684\:88c2\:9879."; 


(* ::Section::Closed:: *)
(*\:8fde\:5206\:6570*)


ContinuedFractionExpand::usage = "\
ContinuedFractionExpand[f,{x,n}] \:7ed9\:51fa\:51fd\:6570 f \:7684\:524d n \:9879\:8fde\:5206\:6570\:5c55\:5f00.";


ContinuedFractionExpandPeriod::usage = "\
ContinuedFractionExpandPeriod[\!\(\*FormBox[SqrtBox[\(r\)],
TraditionalForm]\),x] \:7ed9\:51fa\:51fd\:6570 \!\(\*FormBox[SqrtBox[\(r\)],
TraditionalForm]\) \:7684\:8fde\:5206\:6570\:5c55\:5f00\:6700\:5c0f\:5468\:671f.";


FromContinuedFractionExpand::usage = "\
FromContinuedFractionExpand[list] \:6839\:636e list \:6784\:9020\:8fde\:5206\:6570.";


(* ::Section::Closed:: *)
(*\:4e0d\:5b9a\:79ef\:5206*)


IBP::usage = "\
IBP[u,x] \:7ed9\:51fa \[Integral]u\[DifferentialD]v \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,v,x] \:7ed9\:51fa \[Integral]u\[DifferentialD]x \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,{x,a,b}] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]u \[DifferentialD]x\),
TraditionalForm]\) \:7684\:5206\:90e8\:79ef\:5206.
IBP[u,v,{x,a,b}] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]u \[DifferentialD]v\),
TraditionalForm]\) \:7684\:5206\:90e8\:79ef\:5206.";  


IBS::usage = "\
IBS[f,ex==et] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex==et \:540e\:7684\:7ed3\:679c.
IBS[f,x==et,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 x==et \:540e\:7684\:7ed3\:679c.
IBS[f,ex==t,x] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex==t \:540e\:7684\:7ed3\:679c.
IBS[f,ex==et,x,t] \:7ed9\:51fa \[Integral]f \[DifferentialD]x \:6362\:5143 ex==et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex==et] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]f \[DifferentialD]x\),
TraditionalForm]\) \:6362\:5143 ex==et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},x==et,t] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]f \[DifferentialD]x\),
TraditionalForm]\) \:6362\:5143 x==et \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex==t,x] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]f \[DifferentialD]x\),
TraditionalForm]\) \:6362\:5143 ex==t \:540e\:7684\:7ed3\:679c.
IBS[f,{a,b},ex==et,x,t] \:7ed9\:51fa \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(a\), \(b\)]f \[DifferentialD]x\),
TraditionalForm]\) \:6362\:5143 ex==et \:540e\:7684\:7ed3\:679c.";


IntegrateCF::usage = "\
IntegrateCF[expr,x] \:4f7f\:7528\:8fde\:5206\:6570\:5c55\:5f00\:6cd5\:6c42 expr \:5173\:4e8e x \:7684\:4e0d\:5b9a\:79ef\:5206.";


(* ::Section::Closed:: *)
(*\:4ee3\:6570\:66f2\:7ebf*)


IgusaClebschInvariants::usage = "\
IgusaClebschInvariants[poly,x] \:8fd4\:56de y^2=poly \:7684 Igusa-Clebsch \:4e0d\:53d8\:91cf."


(* ::Section::Closed:: *)
(*\:7ed8\:56fe*)


FastComplexPlot3D::usage = "\
FastComplexPlot3D[f,{x,xmin,xmax},n] \:751f\:6210\:51fd\:6570 f \:5173\:4e8e x \:7684\:4e09\:7ef4\:590d\:5e73\:9762\:7ed8\:56fe.";


(* ::Section::Closed:: *)
(*\:9009\:9879*)


GenerateConstant::usage = "GenerateConstant \:5173\:4e8e\:662f\:5426\:751f\:6210\:5e38\:6570\:7684\:9009\:9879"; 
SimplifyFunction::usage = "SimplifyFunction \:5173\:4e8e\:5316\:7b80\:51fd\:6570\:7684\:9009\:9879"; 


(* ::Subsection::Closed:: *)
(*[BeginPrivate]*)


Begin["`Private`"]; 
ClearAll["`*"]; 


(* ::Chapter:: *)
(*\:5b9a\:4e49*)


(* ::Section:: *)
(*\:5b9e\:7528\:51fd\:6570*)


(* ::Subsection::Closed:: *)
(*TimingTest*)


Attributes[TimingTest] = {HoldAll}; 
TimingTest[code_, n_Integer:1] := Module[{}, ClearSystemCache[]; AbsoluteTiming[Do[code, n - 1]; code]]; 


(* ::Subsection::Closed:: *)
(*PassOptions*)


PassOptions[from_, to_, opts:OptionsPattern[]] := Sequence @@ FilterRules[GatherBy[{opts, Sequence @@ Options[from]}, First][[All,1]], Options[to]]; 


(* ::Section:: *)
(*\:8868\:8fbe\:5f0f*)


(* ::Subsection::Closed:: *)
(*AlgebraicExpressionQ*)


Attributes[AlgebraicExpressionQ] = {Listable}; 
AlgebraicExpressionQ[expr_, x_Symbol] := Catch[RationalExpressionQ[expr /. Sqrt[a_] :> a /. (a_)^(b_) :> (If[ !FreeQ[b, x], Throw[False]]; a), x]]; 


(* ::Subsection::Closed:: *)
(*ExpressionPivot*)


Attributes[ExpressionPivot] = {Listable}; 
ExpressionPivot[expr_] := Catch[Do[If[ !FreeQ[expr, i], Throw[i]], {i, ToExpression["{x,y,z,u,v,w,t}"]}]; FirstCase[expr, _Symbol?( !NumericQ[#1] & ), Symbol, 
      {-1}]]; 


(* ::Subsection::Closed:: *)
(*CoefficientSeparation*)


Attributes[CoefficientSeparation] = {Listable}; 
CoefficientSeparation[expr_, x_Symbol] /; FreeQ[expr, x] := {expr, 1}; 
CoefficientSeparation[expr_, x_Symbol] := Replace[expr, Longest[c_.]*(r_) /; FreeQ[c, x] -> {c, r}]; 


(* ::Section:: *)
(*\:591a\:9879\:5f0f*)


(* ::Subsection::Closed:: *)
(*FirstCoefficient*)


Attributes[FirstCoefficient] = {Listable}; 
FirstCoefficient[poly_, x_Symbol] := Coefficient[poly, x, Exponent[poly, x]]; 


(* ::Subsection::Closed:: *)
(*PolynomialReverse*)


PolynomialReverse[poly_, x_Symbol] /; PolynomialQ[poly, x] := FromDigits[CoefficientList[poly, x], x]; 
PolynomialReverse[poly_, x_Symbol, n_Integer] /; PolynomialQ[poly, x] := FromDigits[PadRight[CoefficientList[poly, x], n + 1], x]; 


(* ::Subsection::Closed:: *)
(*PolynomialRoots*)


Attributes[PolynomialRoots] = {Listable}; 
Options[PolynomialRoots] = {Cubics -> False, Quartics -> False}; 
PolynomialRoots[poly_, x_Symbol, opts:OptionsPattern[]] /; PolynomialQ[poly, x] := 
   Last @@@ Join[ToRules[Roots[poly == 0, x, PassOptions[PolynomialRoots, Roots, opts]]]]; 


(* ::Subsection::Closed:: *)
(*PolynomialRootApproximant*)


Attributes[PolynomialRootApproximant] = {Listable}; 
PolynomialRootApproximant[poly_, x_Symbol] /; PolynomialQ[poly, x] := FromDigits[RootApproximant /@ Reverse[CoefficientList[poly, x]], x]; 


(* ::Subsection::Closed:: *)
(*ComplexFactor*)


Attributes[ComplexFactor] = {Listable}; 
Options[ComplexFactor] = Options[PolynomialRoots]; 
ComplexFactor[poly_, x_Symbol, opts:OptionsPattern[]] /; PolynomialQ[poly, x] := 
   FirstCoefficient[poly, x]*Times @@ (x - PolynomialRoots[poly, x, PassOptions[ComplexFactor, PolynomialRoots, opts]]); 
ComplexFactor[poly_, opts:OptionsPattern[]] := ComplexFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection::Closed:: *)
(*RealFactor*)


Attributes[RealFactor] = {Listable}; 
Options[RealFactor] = {WorkingPrecision -> 100}; 
RealFactor[poly_, x_Symbol, OptionsPattern[]] /; PolynomialQ[poly, x] := FirstCoefficient[poly, x]*
    Product[ToRadicals[If[Element[r, Reals], x - r, PolynomialRootApproximant[N[(x - Re[r])^2 + Im[r]^2, OptionValue[WorkingPrecision]], x]]], 
     {r, DeleteDuplicates[PolynomialRoots[poly, x], N[Conjugate[#1] == #2] & ]}]; 
RealFactor[poly_, opts:OptionsPattern[]] := RealFactor[poly, ExpressionPivot[poly], opts]; 


(* ::Section:: *)
(*\:521d\:7b49\:51fd\:6570*)


(* ::Subsection::Closed:: *)
(*PolynomialFit*)


Options[PolynomialFit] = {InterpolatingFunction -> (RandomReal[#1, WorkingPrecision -> 100] & )}; 
PolynomialFit[expr_, x_Symbol, n_Integer, OptionsPattern[]] := Module[{func = OptionValue[InterpolatingFunction]}, 
    InterpolatingPolynomial[Table[({#1, expr /. x -> #1} & )[func[i]], {i, n + 1}], x]]; 


(* ::Subsection::Closed:: *)
(*ApartArcTan*)


Attributes[ApartArcTan] = {Listable}; 
Options[ApartArcTan] = {GenerateConstant -> True, SimplifyFunction -> (Simplify[FunctionExpand[ComplexExpand[ToRadicals[#1]]]] & )}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] /; RationalExpressionQ[expr, x] := Module[{result, re, im, count = 0, poles, sample, c = {}}, 
    result = Sum[{re, im} = ReIm[r]; count += Sign[im]; ArcTan[OptionValue[SimplifyFunction][(x - re)/im]], 
       {r, PolynomialRoots[(#1 - #2*I & ) @@ NumeratorDenominator[Together[expr]], x]}]; If[OptionValue[GenerateConstant], 
      result += Block[{nd, deg, ratio}, nd = NumeratorDenominator[Together[expr]]; deg = Exponent[nd, x]; If[Subtract @@ deg < 0, 0, 
            ratio = Divide @@ Coefficient[nd, x, deg]; If[Subtract @@ deg > 0, Sign[ratio]*(Pi/2), ArcTan[ratio]]]] - count*(Pi/2); 
       poles = DeleteDuplicates[SolveValues[Denominator[Together[expr]] == 0, x, Reals]]; If[Length[poles] =!= 0, 
        sample = (((1/2)*Table[#1[[i]] + #1[[i + 1]], {i, Length[#1] - 1}] & )[Prepend[#1, First[#1] - 1]] & )[N[poles]]; 
         c = Pi*Round[(1/Pi)*Table[ArcTan[expr] - result /. x -> i, {i, sample}]]; 
         result += Piecewise[Prepend[Table[{c[[i + 1]], poles[[i]] < x < poles[[i + 1]]}, {i, Length[poles] - 1}], {First[c], x < First[poles]}]]]; ]; result]; 
ApartArcTan[expr_, opts:OptionsPattern[]] := ApartArcTan[expr, ExpressionPivot[expr], opts]; 


(* ::Subsection::Closed:: *)
(*ComplexApart*)


Attributes[ComplexApart] = {Listable}; 
Options[ComplexApart] = Options[ComplexFactor]; 
ComplexApart[expr_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/ComplexFactor[#2, x, PassOptions[ComplexFactor, ComplexApart, opts]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
ComplexApart[poly_, opts:OptionsPattern[]] := ComplexApart[poly, ExpressionPivot[poly], opts]; 


(* ::Subsection::Closed:: *)
(*RealApart*)


Attributes[RealApart] = {Listable}; 
Options[RealApart] = Options[RealFactor]; 
RealApart[expr_, x_Symbol, opts:OptionsPattern[]] /; RationalExpressionQ[expr, x] := 
   (Apart[#1/RealFactor[#2, x, PassOptions[RealApart, RealFactor, opts]], x] & ) @@ NumeratorDenominator[Together[expr]]; 
RealApart[expr_, opts:OptionsPattern[]] := RealApart[expr, ExpressionPivot[expr], opts]; 


(* ::Section:: *)
(*\:8fde\:5206\:6570*)


(* ::Subsection::Closed:: *)
(*ContinuedFractionExpand*)


ContinuedFractionExpand[f_, {x_Symbol, n_Integer}] := Module[{a, b = f}, Quiet[Table[a = Normal[Series[b, {x, Infinity, 0}]]; b = FullSimplify[1/(b - a)]; a, 
      {i, n}], Series::sbyc]]; 


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
(*FromContinuedFractionExpand*)


FromContinuedFractionExpand[list_List] := Fold[#2 + 1/#1 & , Reverse[list]]; 


(* ::Section:: *)
(*\:4e0d\:5b9a\:79ef\:5206*)


(* ::Subsection::Closed:: *)
(*IBP*)


Options[IBP] = {Assumptions -> $Assumptions}; 
IBP[u_, v_, x_Symbol] := Module[{c, r}, {c, r} = CoefficientSeparation[v*D[u, x], x]; u*v - c*Inactive[Integrate][r, x]]; 
IBP[f_, x_Symbol] := IBP[f, x, x]; 
IBP[u_, v_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := Module[{c, r, optLimit = PassOptions[IBP, Limit, opts]}, 
    {c, r} = CoefficientSeparation[v*D[u, x], x]; Limit[u*v, x -> b, optLimit] - Limit[u*v, x -> a, optLimit] - c*Inactive[Integrate][r, {x, a, b}]]; 
IBP[f_, {x_Symbol, a_, b_}, opts:OptionsPattern[]] := IBP[f, x, {x, a, b}, opts]; 


(* ::Subsection::Closed:: *)
(*IBS*)


Options[IBS] = {Assumptions -> $Assumptions, Method -> Automatic}; 
IBS[f_, (ex_) == (et_), x_Symbol, t_Symbol, opts:OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{u, result, g, optMethod = OptionValue[Method]}, 
    result = If[ex === x, f /. x -> et, If[(optMethod === Automatic && AlgebraicExpressionQ[ex, x]) || optMethod === GroebnerBasis, 
          First[SolveValues[Factor[GroebnerBasis[{ex == u, g*D[ex, x] == f}, {x, u}, {x}]] == 0, g]], 
          First[IntegrateChangeVariables[Inactive[Integrate][f, x], u, u == ex, PassOptions[IBS, IntegrateChangeVariables, opts]]] /. C[_] -> 0] /. u -> et]*
       D[et, t]; result]; 
IBS[f_, (x_Symbol) == (et_), t_Symbol, opts:OptionsPattern[]] := IBS[f, x == et, x, t, opts]; 
IBS[f_, (ex_) == (t_Symbol), x_Symbol, opts:OptionsPattern[]] := IBS[f, ex == t, x, t, opts]; 
IBS[f_, (ex_) == (et_), opts:OptionsPattern[]] := IBS[f, ex == et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 
IBS[f_, {a_, b_}, (ex_) == (et_), x_Symbol, t_Symbol, opts:OptionsPattern[]] /; FreeQ[ex, t] &&  !FreeQ[ex, x] && FreeQ[et, x] &&  !FreeQ[et, t] := 
   Module[{u, temp, result, newA, newB, c, r, optICV = PassOptions[IBS, IntegrateChangeVariables, opts]}, 
    temp = If[ex === x, Inactive[Integrate][f, {x, a, b}] /. x -> u, IntegrateChangeVariables[Inactive[Integrate][f, {x, a, b}], u, u == ex, optICV]]; 
     {result, {t, newA, newB}} = List @@ If[et === t, temp /. u -> t, IntegrateChangeVariables[temp, t, u == et, optICV]]; 
     {c, r} = CoefficientSeparation[result, t]; c*Inactive[Integrate][r, {t, newA, newB}]]; 
IBS[f_, {a_, b_}, (x_Symbol) == (et_), t_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, x == et, x, t, opts]; 
IBS[f_, {a_, b_}, (ex_) == (t_Symbol), x_Symbol, opts:OptionsPattern[]] := IBS[f, {a, b}, ex == t, x, t, opts]; 
IBS[f_, {a_, b_}, (ex_) == (et_), opts:OptionsPattern[]] := IBS[f, {a, b}, ex == et, ExpressionPivot[ex], ExpressionPivot[et], opts]; 


(* ::Subsection::Closed:: *)
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
(*\:4ee3\:6570\:66f2\:7ebf*)


(* ::Subsection::Closed:: *)
(*IgusaClebschInvariants*)


{IgusaClebschInvariantsI2, IgusaClebschInvariantsI4, IgusaClebschInvariantsI6, IgusaClebschInvariantsI10} = 
   Uncompress["1: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"]; 
Attributes[IgusaClebschInvariants] = {Listable}; 
IgusaClebschInvariants[poly_, x_Symbol] /; PolynomialQ[poly, x] := Module[{v}, v = PadRight[CoefficientList[poly, x], 7]; 
     {IgusaClebschInvariantsI2 . v . v, IgusaClebschInvariantsI4 . v . v . v . v, IgusaClebschInvariantsI6 . v . v . v . v . v . v, 
      IgusaClebschInvariantsI10 . v . v . v . v . v . v . v . v . v . v}]; 


(* ::Section:: *)
(*\:7ed8\:56fe*)


(* ::Subsection::Closed:: *)
(*FastComplexPlot3D*)


Options[FastComplexPlot3D] = Join[{Mesh -> None, AxesLabel -> {Re, Im}, InterpolationOrder -> 0, WorkingPrecision -> MachinePrecision, Chop -> 10^(-10), 
     DiscretePlot -> True}, Options[ListPlot3D]]; 
FastComplexPlot3D[fn_, {z_Symbol, zmin_, zmax_}, n_Integer, opts:OptionsPattern[]] /; n > 0 := 
   Module[{delta, prec, rmin, rmax, imin, imax, f, data, colorf, cnt, cur, colorfd, ticks}, delta = OptionValue[Chop]; prec = OptionValue[WorkingPrecision]; 
     {rmin, imin} = ReIm[N[zmin, prec]]; {rmax, imax} = ReIm[N[zmax, prec]]; f = Function[z, (If[NumericQ[#1], Chop[#1, delta], Indeterminate] & )[N[fn, prec]]]; 
     data = Quiet[Table[f[a + b*I], {a, rmin, rmax, (rmax - rmin)/n}, {b, imin, imax, (imax - imin)/n}]]; 
     colorf = Function[{x, y}, Hue[Arg[data[[Round[n*x + 1],Round[n*y + 1]]]]/(2*Pi), 0.75]]; cnt = 0; 
     colorfd = Function[{x, y}, If[Mod[cnt++, 4] == 0, cur = data[[Round[n*x + 1],Round[n*y + 1]]]]; Hue[Arg[cur]/(2*Pi), 0.75]]; 
     ticks = Table[{(n*i)/4 + 1, (If[Abs[#1 - Round[#1]] < 0.01, Round[#1], #1] & )[Round[rmin + ((rmax - rmin)/4)*i, 0.01]]}, {i, 0, 4}]; 
     ListPlot3D[Map[Abs, Transpose[data], {2}], ColorFunction -> If[OptionValue[DiscretePlot], colorfd, colorf], Ticks -> {ticks, ticks, Automatic}, 
      PassOptions[FastComplexPlot3D, ListPlot3D, opts]]]; 


(* ::Subsection::Closed:: *)
(*[EndPrivate]*)


End[];


(* ::Subsection::Closed:: *)
(*[EndPackage]*)


EndPackage[];
