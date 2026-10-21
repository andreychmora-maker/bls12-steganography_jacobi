/* 
===========================================================================
Standalone Benchmark Script for BLS12-479+ and BLS12-539+ Jacobi Pipelines 
Features: Pure Algebraic Cycle Counter (wNAF + On-the-fly Normalization)
===========================================================================
*/

function RunBenchmark(p_val, blinding_bits, clearing_bits, iterations)
    Fp := FiniteField(p_val);
    
    // We use d=1 for strict compliance with your 9M+5S formula
    a_2 := Fp!12345; 
    DummyPrj := [Fp!0, Fp!0, Fp!0]; 
    
    // Addition formulas for P(1,1,1) from your code
    DedicatedDblJq_Prj := function(Dummy, a_2_val, PrjPoint)
        X_1:=PrjPoint[1]; Y_1:=PrjPoint[2]; Z_1:=PrjPoint[3];
        A:=2*Z_1^2; B:=X_1^2; C:=Y_1^2; E:=2*X_1*Y_1; 
        F:=a_2_val*B; G:=C-F; H:=A-G;
        return [E*H, 2*C*G - H^2, H^2];
    end function;

    UnifiedAddJq_Prj := function(Dummy, a_2_val, PrjPoint1, PrjPoint2)
        X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
        X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; Z_2:=PrjPoint2[3];
        A:=Z_1*Z_2; B:=Y_1*Y_2; C:=X_1*X_2; 
        D:=(X_1+Y_1)*(X_2+Y_2)-B-C;
        F:=(X_1*Z_2)^2; G:=(Z_1*X_2)^2; H:=A^2; J:=C^2; 
        L:=H-J; M:=L^2; N:=(A+L)^2-H-M; P:=(A+C)^2-H-J;
        return [D*N, 2*A*(B+a_2_val*C)*(H+J)+2*P*(F+G), 2*M];
    end function;

    // Ultra-fast addition (9M+5S), requiring Z2=1
    DedicatedAddJq_Prj := function(Dummy, a_2_val, PrjPoint1, PrjPoint2)
        X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
        X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; 
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
        
        // Doubling and normalization
        P2 := DedicatedDblJq_Prj(DummyPrj, a_2_val, PrjPoint);
        Zinv2 := P2[3]^-1;
        P2 := [P2[1]*Zinv2, P2[2]*Zinv2, Fp!1];
        
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

    w_obf := 5;
    w_deobf := 2;

    ObfuscateCls := 0; 
    DeObfuscateCls := 0; 

    // Dummy base points to load the Jacobi arithmetic
    P_base := [Fp!123, Fp!456, Fp!1]; 
    P_u := [Fp!789, Fp!101, Fp!1];
    P_v := [Fp!112, Fp!131, Fp!1];

    BlindScalar := Random(2^(blinding_bits-1), 2^blinding_bits - 1);
    FastCofactor := Random(2^(clearing_bits-1), 2^clearing_bits - 1);

    for i in [1..iterations] do
        // === OBFUSCATION ===
        t := ClockCycles(); 
        R_cof_J := ScalarMulJqwNAF_Prj(DummyPrj, a_2, w_obf, BlindScalar, P_base);
        // Emulating result normalization
        Zinv := R_cof_J[3]^-1;
        R_cof_norm := [R_cof_J[1]*Zinv, R_cof_J[2]*Zinv, Fp!1];
        // Emulating the final addition
        P_isog_J := DedicatedAddJq_Prj(DummyPrj, a_2, P_base, R_cof_norm);
        ObfuscateCls +:= (ClockCycles() - t);

        // === DEOBFUSCATION ===
        t := ClockCycles(); 
        P_sum_J := DedicatedAddJq_Prj(DummyPrj, a_2, P_u, P_v);
        ZinvSum := P_sum_J[3]^-1;
        P_sum_J := [P_sum_J[1]*ZinvSum, P_sum_J[2]*ZinvSum, Fp!1];
        P_cleared_J := ScalarMulJqwNAF_Prj(DummyPrj, a_2, w_deobf, FastCofactor, P_sum_J);
        DeObfuscateCls +:= (ClockCycles() - t);
    end for;

    return ObfuscateCls div iterations, DeObfuscateCls div iterations;
end function;

iterations := 100000;

p_479 := 1040582850105330743815570414239668597581975900197275493095847470017583224665441180782578049215009464174278151947833029365685675660572997417552727;
printf "Running raw algebraic benchmarks for BLS12-479+ (%o iterations)...\n", iterations;
obf_479, deobf_479 := RunBenchmark(p_479, 479, 80, iterations);
printf "BLS12-479+ -> Obfuscation: %o cycles | Deobfuscation: %o cycles\n\n", obf_479, deobf_479;

p_539 := 1199710345211519035416219429741786876319311421933497483626991418658478319088544596915796254702839442454587420641007826572438113715342206999017869775578752368205297; 
printf "Running raw algebraic benchmarks for BLS12-539+ (%o iterations)...\n", iterations;
obf_539, deobf_539 := RunBenchmark(p_539, 539, 80, iterations);
printf "BLS12-539+ -> Obfuscation: %o cycles | Deobfuscation: %o cycles\n", obf_539, deobf_539;

quit;
