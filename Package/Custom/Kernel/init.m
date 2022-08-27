(* ::Package:: *)

(* ::Subsection:: *)
(*\:5206\:90e8\:79ef\:5206*)


IntegrateByParts::usage = "\:5206\:90e8\:79ef\:5206"; 
IntegrateByParts[u_, v_, x_Symbol] := Simplify[u*v - Inactive[Integrate][v*D[u, x], x]]; 
IntegrateByParts[f_, x_Symbol] := IntegrateByParts[f, x, x]; 
IntegrateByParts[u_, v_, {x_Symbol, a_, b_}] := Simplify[Limit[u*v, x -> b] - Limit[u*v, x -> a] - Inactive[Integrate][v*D[u, x], {x, a, b}], 
    a < x < b]; 
IntegrateByParts[f_, {x_Symbol, a_, b_}] := IntegrateByParts[f, x, {x, a, b}]; 
IBP = IntegrateByParts; 


(* ::Subsection:: *)
(*\:79ef\:5206\:6362\:5143(bug)*)


IntegrateBySubstitution::usage = "\:79ef\:5206\:6362\:5143"; 
IntegrateBySubstitution[f_, x_Symbol -> {et_, t_Symbol}, cond_:Null] := If[cond =!= Null, FullSimplify[#1, cond] & , Simplify][f /. x -> et]*
    D[et, t]; 
IntegrateBySubstitution[f_, {ex_, x_Symbol} -> t_Symbol, cond_:Null] := 
   Simplify[IntegrateBySubstitution[f, x -> {First[x /. Solve[ex == t, x]], t}, cond]]; 
IntegrateBySubstitution[f_, {ex_, x_Symbol} -> {et_, t_Symbol}, cond_:Null] := 
   Simplify[IntegrateBySubstitution[f, x -> {First[x /. Solve[ex == et, x]], t}, cond]]; 
IBS = IntegrateBySubstitution; 


(* ::Subsection:: *)
(*arctan \:88c2\:9879*)


ApartArcTan::usage = "arctan \:88c2\:9879\n\:4f7f\:7528 Constant->True \:8ba1\:7b97\:5e38\:6570\n\:4f7f\:7528 SimplifyFunction->_ \:8bbe\:7f6e\:5316\:7b80\:51fd\:6570"; 
SimplifyFunction = SimplifyFunction; 
Options[ApartArcTan] = {Constant -> False, SimplifyFunction -> Simplify}; 
ApartArcTan[expr_, x_Symbol, OptionsPattern[]] := Module[{root, simp, res, pole, diff}, root = x /. Solve[expr == I, x]; 
     simp = OptionValue[SimplifyFunction]; res = Sum[ArcTan[simp[(x - Re[root[[k]]])/Im[root[[k]]]]], {k, 1, Length[root]}]; 
     If[ !OptionValue[Constant], Return[res]]; pole = (#1[[1]] & ) /@ FunctionPoles[{expr, Element[x, Reals]}, x]; 
     res += If[Length[pole] === 0, ArcTan[expr] - res /. x -> 0, 
       diff = Append[Table[Limit[ArcTan[expr], x -> i, Direction -> "FromBelow"] - (res /. x -> i), {i, pole}], 
          Limit[ArcTan[expr], x -> Last[pole], Direction -> "FromAbove"] - (res /. x -> Last[pole])]; 
        Piecewise[Append[Prepend[Table[{diff[[i + 1]], pole[[i]] < x < pole[[i + 1]]}, {i, Length[pole] - 1}], {First[diff], x < First[pole]}], 
          {Last[diff], x > Last[pole]}]]]; Simplify[res]]; 


(* ::Subsection:: *)
(*\:5b9e\:56e0\:5f0f\:5206\:89e3*)


RealFactor::usage = "\:5b9e\:56e0\:5f0f\:5206\:89e3\n\:4f7f\:7528 SimplifyFunction->_ \:8bbe\:7f6e\:5316\:7b80\:51fd\:6570"; 
SimplifyFunction = SimplifyFunction; 
Options[RealFactor] = {SimplifyFunction -> (Simplify[ToRadicals[RootReduce[#1]]] & )}; 
RealFactor[(poly_)?PolynomialQ, x_Symbol, OptionsPattern[]] := Module[{coef, simp, sol}, coef = CoefficientList[poly, x]; 
     simp = OptionValue[SimplifyFunction]; If[ !Element[coef, Reals], Return[HoldForm[RealFactor[poly, x]]]]; 
     sol = Sort[x /. Solve[poly == 0, x], Less @@ Re[N[{#1, #2}]] & ]; 
     Do[If[i > Length[sol], Break[]]; If[Im[sol[[i]]] === 0, sol[[i]] = x - simp[sol[[i]]], 
        sol[[i]] = x^2 - simp[sol[[i]] + sol[[i + 1]]]*x + simp[sol[[i]]*sol[[i + 1]]]; sol = Delete[sol, i + 1]; ]; , {i, Infinity}]; 
     Last[coef]*Times @@ sol]; 
