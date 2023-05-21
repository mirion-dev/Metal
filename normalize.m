(* ::Package:: *)

BoxForm`MakeTraditionalBoxes[Log] := "ln"; 
BoxForm`MakeTraditionalBoxes[ArcSin] := "arcsin"; 
BoxForm`MakeTraditionalBoxes[ArcCos] := "arccos"; 
BoxForm`MakeTraditionalBoxes[ArcTan] := "arctan"; 
BoxForm`MakeTraditionalBoxes[ArcCot] := "arccot"; 
BoxForm`MakeTraditionalBoxes[ArcSec] := "arcsec"; 
BoxForm`MakeTraditionalBoxes[ArcCsc] := "arccsc"; 
BoxForm`MakeTraditionalBoxes[ArcSinh] := "arsinh"; 
BoxForm`MakeTraditionalBoxes[ArcCosh] := "arcosh"; 
BoxForm`MakeTraditionalBoxes[ArcTanh] := "artanh"; 
BoxForm`MakeTraditionalBoxes[ArcCoth] := "arcoth"; 
BoxForm`MakeTraditionalBoxes[ArcSech] := "arsech"; 
BoxForm`MakeTraditionalBoxes[ArcCsch] := "arcsch"; 
BoxForm`MakeTraditionalExpression["arsinh"] := ArcSinh; 
BoxForm`MakeTraditionalExpression["arcosh"] := ArcCosh; 
BoxForm`MakeTraditionalExpression["artanh"] := ArcTanh; 
BoxForm`MakeTraditionalExpression["arcoth"] := ArcCoth; 
BoxForm`MakeTraditionalExpression["arsech"] := ArcSech; 
BoxForm`MakeTraditionalExpression["arcsch"] := ArcCsch; 


Unprotect[StirlingS1, StirlingS2]; 
Clear[StirlingS1, StirlingS2]; 
Protect[StirlingS1, StirlingS2]; 
BoxForm`MakeTraditionalBoxes[StirlingS1] := InterpretationBox["s", StirlingS1]; 
MakeBoxes[StirlingS2[a_, b_], TraditionalForm] := TemplateBox[{MakeBoxes[a, TraditionalForm], MakeBoxes[b, TraditionalForm]}, "StirlingS2"]; 
MakeBoxes[StirlingS2[a_, b_], TraditionalForm] := 
   TagBox[RowBox[{"{", GridBox[{{TagBox[MakeBoxes[a, TraditionalForm], Identity, Editable -> True, Selectable -> True]}, 
        {TagBox[MakeBoxes[b, TraditionalForm], Identity, Editable -> True, Selectable -> True]}}], "}"}], InterpretTemplate[StirlingS2[#1, #2] & ], 
    Editable -> False, Selectable -> False]; 
