/* 
===========================================================================
Standalone Microbenchmark: Native Intrinsic (k*P) vs Interpreted Jacobi wNAF
Curve Target: BLS12-479+
===========================================================================
*/

p_val := 1040582850105330743815570414239668597581975900197275493095847470017583224665441180782578049215009464174278151947833029365685675660572997417552727;
Fp := FiniteField(p_val);

// Base Weierstrass curve for the native intrinsic test
E := EllipticCurve([Fp | 0, 1]);
P_native := Random(E);
while P_native eq E!0 do P_native := Random(E); end while;

// Jacobi Quartic constants and dummy variables
a_2 := Fp!12345; 
DummyPrj := [Fp!0, Fp!0, Fp!0]; 
P_jacobi := [Fp!123, Fp!456, Fp!1]; // Dummy point to load the arithmetic pipeline

DedicatedDblJq_Prj := function(Dummy, a_2_val, PrjPoint)
    X_1:=PrjPoint[1]; Y_1:=PrjPoint[2]; Z_1:=PrjPoint[3];
    A:=2*Z_1^2; B:=X_1^2; C:=Y_1^2; E:=2*X_1*Y_1; 
    F:=a_2_val*B; G:=C-F; H:=A-G;
    return [E*H, 2*C*G - H^2, H^2];
end function;

DedicatedAddJq_Prj := function(Dummy, a_2_val, PrjPoint1, PrjPoint2)
    X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
    X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; // Z_2 is strictly 1 here
    A:=Z_1; B:=Y_1*Y_2; C:=X_1*X_2;
    D:=X_1*Y_2 - X_2*Y_1; F:=X_1^2; G:=(Z_1*X_2)^2;
    return [D*(F-G), (B-a_2_val*C)*(F+G)-2*C*(A^2+C^2), A*D^2];
end function;

NAF_w := procedure(~ArrayNAF, Scalar, w)
    Rise2_w := 2^w; 
    Rise2_wP1 := 2^(w+1);
    Tmp := Scalar;
    while Tmp gt 0 do
       if IsOdd(Tmp) then
          Val := Tmp mod Rise2_wP1;
          if Val ge Rise2_w then Val -:= Rise2_wP1; end if;
          Tmp -:= Val;
       else
          Val := 0;
       end if;
       Append(~ArrayNAF, Val);
       Tmp div:= 2;
    end while;
end procedure;

CreateLookUpTable_wNAF := function(DummyPrj, a_2_val, w, PrjPoint)
    Rise2_w := 2^w;
    Table := [];
    for i in [1..Rise2_w] do Append(~Table, DummyPrj); end for;
    
    Table[1] := PrjPoint; 
    
    // Compute P2 and normalize it on the fly (Z=1)
    P2 := DedicatedDblJq_Prj(DummyPrj, a_2_val, PrjPoint);
    Zinv2 := P2[3]^-1;
    P2 := [P2[1]*Zinv2, P2[2]*Zinv2, Fp!1];
    
    // Populate the remainder of the table using fast dedicated addition
    for i in [3 .. Rise2_w-1 by 2] do
        P_next := DedicatedAddJq_Prj(DummyPrj, a_2_val, Table[i-2], P2);
        Zinv := P_next[3]^-1;
        Table[i] := [P_next[1]*Zinv, P_next[2]*Zinv, Fp!1];
    end for;
    return Table;
end function;

ScalarMulJqwNAF_Prj := function(DummyPrj, a_2_val, w, Scalar, PrjPoint)
    if Scalar eq 0 then return DummyPrj; end if;
    
    ArrayNAF := [];
    NAF_w(~ArrayNAF, Scalar, w);
    Table := CreateLookUpTable_wNAF(DummyPrj, a_2_val, w, PrjPoint);
    
    ResultPoint := DummyPrj;
    started := false;
    
    for i := #ArrayNAF to 1 by -1 do
        if started then
            ResultPoint := DedicatedDblJq_Prj(DummyPrj, a_2_val, ResultPoint);
        end if;
        
        Val := ArrayNAF[i];
        if Val ne 0 then
            P_add := Val gt 0 select Table[Val] else [-Table[-Val][1], Table[-Val][2], Table[-Val][3]];
            if not started then
                ResultPoint := P_add;
                started := true;
            else
                ResultPoint := DedicatedAddJq_Prj(DummyPrj, a_2_val, ResultPoint, P_add);
            end if;
        end if;
    end for;
    return ResultPoint;
end function;

iterations := 100000;
BlindScalar := Random(2^478, 2^479 - 1);
w_obf := 5;

printf "Running Standalone Microbenchmark (%o iterations)...\n\n", iterations;

// --- TEST 1: Native Magma Intrinsic (Compiled C) ---
t := ClockCycles();
for i in [1..iterations] do
    res := BlindScalar * P_native;
end for;
native_cycles := (ClockCycles() - t) div iterations;
printf "Native Magma Intrinsic (k * P_weierstrass): %o cycles\n", native_cycles;

// --- TEST 2: Interpreted Jacobi wNAF Pipeline ---
t := ClockCycles();
for i in [1..iterations] do
    res := ScalarMulJqwNAF_Prj(DummyPrj, a_2, w_obf, BlindScalar, P_jacobi);
end for;
jacobi_cycles := (ClockCycles() - t) div iterations;
printf "Interpreted Jacobi wNAF (ScalarMulJqwNAF):  %o cycles\n", jacobi_cycles;

printf "\nRatio (Jacobi / Native): %o%%\n", RealField(5)!((jacobi_cycles / native_cycles) * 100);
quit;
