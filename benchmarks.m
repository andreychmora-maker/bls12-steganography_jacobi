/* 
===========================================================================
Standalone Benchmark Script for BLS12-479+ and BLS12-539+ Jacobi Pipelines 
Features: Asymmetric wNAF, Explicit 2-Isogenies, Billet-Joye Fusion
===========================================================================
*/

function RunBenchmark(p_val, iterations)
    Fp := FiniteField(p_val);

    // Dynamically calculate precise BLS12 parameters from 'p'
    z_approx := Iroot(3 * p_val, 6);
    z_val := 0;
    for cand in [-z_approx-2 .. -z_approx+2] cat [z_approx-2 .. z_approx+2] do
        if (cand-1)^2 * (cand^4 - cand^2 + 1) div 3 + cand eq p_val then
            z_val := cand;
            break;
        end if;
    end for;

    q := z_val^4 - z_val^2 + 1;
    Cofactor := (z_val - 1)^2 div 3;

    B_target := Fp!1;
    A_isog   := Fp!-15;
    B_isog   := Fp!22;

    E_target := EllipticCurve([Fp | 0, B_target]);
    E_isog   := EllipticCurve([Fp | A_isog, B_isog]);
    Identity := E_target!0;

    // Jacobi Quartic Constants (Strictly Aligned with Hisil Arithmetic)
    theta := Fp!2; // 2-torsion point on E_isog
    a_bj  := Fp!6; // Shifted curve parameter 'a'
    b_bj  := Fp!-3; // Shifted curve parameter 'b'

    a_1  := -a_bj / 4; 
    a_2  := 2 * a_1;   
    d_jq := (a_bj^2 - 4 * b_bj) / 16; 

    DummyPrj := [Fp!0, Fp!0, Fp!0];

    // Birational Mappings: Exact Homogeneous Hisil Coordinates
    W2J := function(P)
        if P eq E_isog!0 then return [Fp!0, Fp!1, Fp!1]; end if;
        x := P[1]; y := P[2];
        if y eq 0 and x eq theta then return [Fp!0, Fp!-1, Fp!1]; end if;
        x_shift := x - theta;
        return [2 * x_shift * y, x_shift * (x_shift^2 - b_bj), y^2];
    end function;

    J2W := function(PJ)
        XJ := PJ[1]; YJ := PJ[2]; ZJ := PJ[3];
        if ZJ eq 0 then return E_isog!0; end if;
        if XJ eq 0 and YJ eq ZJ then return E_isog!0; end if;
        if XJ eq 0 and YJ eq -ZJ then return E_isog![theta, 0]; end if;

        num_u := 2*ZJ^2 + 2*YJ*ZJ - XJ^2;
        den_u := XJ^2;
        num_v := 2*ZJ*(2*ZJ^2 + 2*YJ*ZJ - 3*XJ^2);
        den_v := XJ^3;
        
        if den_u eq 0 or den_v eq 0 then return E_isog!0; end if;
        return E_isog![num_u / den_u, num_v / den_v];
    end function;

    // Jacobi Quartic Dedicated Doubling (3M+4S+1C+3add)
    DedicatedDblJq_Prj := function(Result, a_2_val, PrjPoint)
        X_1:=PrjPoint[1]; Y_1:=PrjPoint[2]; Z_1:=PrjPoint[3];
        A:=2*Z_1^2; B:=X_1^2; C:=Y_1^2; E:=2*X_1*Y_1; 
        F:=a_2_val*B; G:=C-F; H:=A-G;
        Result[3]:=H^2;
        Result[1]:=E*H;
        Result[2]:=2*C*G-Result[3];
        return Result;
    end function;

    // Jacobi Quartic Unified Addition (9M+7S+1C+15add)
    UnifiedAddJq_Prj := function(Result, d, a_2_val, PrjPoint1, PrjPoint2)
        X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
        X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; Z_2:=PrjPoint2[3];
        A:=Z_1*Z_2; B:=Y_1*Y_2; C:=X_1*X_2; 
        D:=(X_1+Y_1)*(X_2+Y_2)-B-C;
        F:=(X_1*Z_2)^2; G:=(Z_1*X_2)^2; H:=A^2; J:=C^2; 
        K:=d*J; 
        L:=H-K; M:=L^2; N:=(A+L)^2-H-M; P:=(A+C)^2-H-J;
        Result[1]:=D*N;
        Result[2]:=2*A*(B+a_2_val*C)*(H+K)+2*d*P*(F+G); 
        Result[3]:=2*M;
        return Result;
    end function;

    // Generate wNAF array
    NAF_w := procedure(~ArrayNAF, Scalar, w)
        Rise2_w := 2^w; 
        Rise2_wP1 := 2^(w+1);
        Tmp := Scalar;
        while Tmp gt 0 do
           if IsOdd(Tmp) then
              Val := Tmp mod Rise2_wP1;
              if Val ge Rise2_w then
                Val -:= Rise2_wP1;
              end if;
              Tmp -:= Val;
           else
              Val := 0;
           end if;
           Append(~ArrayNAF, Val);
           Tmp div:= 2;
        end while;
    end procedure;

    // Create lookup table for odd multiples
    CreateLookUpTable_wNAF := function(DummyPrj, d, a_2_val, w, PrjPoint)
        Rise2_w := 2^w;
        Table := [];
        for i in [1..Rise2_w] do Append(~Table, DummyPrj); end for;
        
        Table[1] := PrjPoint;
        P2 := DedicatedDblJq_Prj(DummyPrj, a_2_val, PrjPoint);
        
        for i in [3 .. Rise2_w-1 by 2] do
            Table[i] := UnifiedAddJq_Prj(DummyPrj, d, a_2_val, Table[i-2], P2);
        end for;
        return Table;
    end function;

    // Main wNAF scalar multiplication function
    ScalarMulJqwNAF_Prj := function(DummyPrj, d, a_2_val, w, Scalar, PrjPoint)
        if Scalar eq 0 then return [Fp!0, Fp!1, Fp!1]; end if;
        
        ArrayNAF := [];
        NAF_w(~ArrayNAF, Scalar, w);
        Table := CreateLookUpTable_wNAF(DummyPrj, d, a_2_val, w, PrjPoint);
        
        ResultPoint := DummyPrj;
        started := false;
        
        for i := #ArrayNAF to 1 by -1 do
            if started then
                // Use ultra-fast dedicated doubling (3M+4S)
                ResultPoint := DedicatedDblJq_Prj(DummyPrj, a_2_val, ResultPoint);
            end if;
            
            Val := ArrayNAF[i];
            if Val ne 0 then
                P_add := Val gt 0 select Table[Val] else [-Table[-Val][1], Table[-Val][2], Table[-Val][3]];
                if not started then
                    ResultPoint := P_add;
                    started := true;
                else
                    ResultPoint := UnifiedAddJq_Prj(DummyPrj, d, a_2_val, ResultPoint, P_add);
                end if;
            end if;
        end for;
        return ResultPoint;
    end function;

    // Forward Map & Backward Map (2-Isogenies)
    Isogeny2Target := function(P_isog)
        if P_isog eq E_isog!0 then return Identity; end if;
        x := P_isog[1]; y := P_isog[2];
        x_den := x - 2; 
        if x_den eq 0 then return Identity; end if;
        X_isog := x - 3 / x_den;
        Y_isog := y + 3 * y / (x_den^2);
        return E_target![X_isog / 4, Y_isog / 8];
    end function;

    EvaluateBackwardMap := function(P_target)
        if P_target eq Identity then return E_isog!0; end if;
        x := P_target[1]; y := P_target[2];
        if x eq -1 then return E_isog!0; end if; 
        X_hat := x + 3 / (x + 1);
        Y_hat := y * (x^2 + 2*x - 2) / ((x + 1)^2);
        return E_isog![X_hat, Y_hat];
    end function;

    // Setup Scalars
    FastCofactor := Abs(z_val - 1);
    L_val := 2;
    fc_mod_q := FastCofactor mod q;
    ScalarFactor := (L_val * fc_mod_q) mod q;
    InvScalar := Integers()!( (Integers(q)!ScalarFactor)^-1 );
    BlindScalar := q * (FastCofactor div 3); 
    
    // wNAF window sizes (Asymmetric approach)
    w_obf := 5;   // Wide window for heavy full-scalar blinding
    w_deobf := 2; // Minimal window for lightweight Fast Cofactor Clearing

    // SENDER PIPELINE (OBFUSCATION) - JACOBI OPTIMIZED
    ObfuscatePublicKey_Jacobi := function(PublicKey_target)
        P_scaled := InvScalar * PublicKey_target;
        P_isog_base := EvaluateBackwardMap(P_scaled); 
        P_isog_base_J := W2J(P_isog_base);
        
        R_base := Random(E_isog);
        R_base_J := W2J(R_base);
        R_cof_J := ScalarMulJqwNAF_Prj(DummyPrj, d_jq, a_2, w_obf, BlindScalar, R_base_J);
        
        P_isog_J := UnifiedAddJq_Prj(DummyPrj, d_jq, a_2, P_isog_base_J, R_cof_J);
        P_isog := J2W(P_isog_J);
        
        P_u_isog := Random(E_isog); 
        P_v_isog := P_isog - P_u_isog;
        return P_u_isog, P_v_isog;
    end function;

    // RECEIVER PIPELINE (DEOBFUSCATION) - JACOBI OPTIMIZED
    DeobfuscatePublicKey_Jacobi := function(P_u_isog, P_v_isog)
        P_u_J := W2J(P_u_isog);
        P_v_J := W2J(P_v_isog);
        
        P_sum_J := UnifiedAddJq_Prj(DummyPrj, d_jq, a_2, P_u_J, P_v_J);
        P_cleared_J := ScalarMulJqwNAF_Prj(DummyPrj, d_jq, a_2, w_deobf, FastCofactor, P_sum_J);
        P_cleared_W := J2W(P_cleared_J);
        
        return Isogeny2Target(P_cleared_W);
    end function;

    ObfuscateCls := 0; 
    DeObfuscateCls := 0; 

    // Generate the public key once outside the loop to drastically reduce absolute test execution time
    SK := Random(1, q-1);
    G_target := Cofactor * Random(E_target);
    while G_target eq Identity do G_target := Cofactor * Random(E_target); end while;
    PK := SK * G_target;

    for i in [1..iterations] do
        t := ClockCycles(); 
        pu_b, pv_b := ObfuscatePublicKey_Jacobi(PK);
        ObfuscateCls +:= (ClockCycles() - t);

        t := ClockCycles(); 
        _ := DeobfuscatePublicKey_Jacobi(pu_b, pv_b);
        DeObfuscateCls +:= (ClockCycles() - t);
    end for;

    return ObfuscateCls div iterations, DeObfuscateCls div iterations;
end function;

iterations := 1000;

p_479 := 1040582850105330743815570414239668597581975900197275493095847470017583224665441180782578049215009464174278151947833029365685675660572997417552727;
printf "Running benchmarks for BLS12-479+ (%o iterations)...\n", iterations;
obf_479, deobf_479 := RunBenchmark(p_479, iterations);
printf "BLS12-479+ -> Obfuscation: %o cycles | Deobfuscation: %o cycles\n\n", obf_479, deobf_479;

p_539 := 1199710345211519035416219429741786876319311421933497483626991418658478319088544596915796254702839442454587420641007826572438113715342206999017869775578752368205297;
printf "Running benchmarks for BLS12-539+ (%o iterations)...\n", iterations;
obf_539, deobf_539 := RunBenchmark(p_539, iterations);
printf "BLS12-539+ -> Obfuscation: %o cycles | Deobfuscation: %o cycles\n", obf_539, deobf_539;

quit;
