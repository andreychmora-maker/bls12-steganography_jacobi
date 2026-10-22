/* 
===========================================================================
Standalone Benchmark Script for BLS12-479+ and BLS12-539+ Jacobi Pipelines 
Features: Pure Algebraic Cycle Counter (Montgomery Ladder + Normalization)
===========================================================================
*/

function RunBenchmark(p_val, blinding_bits, clearing_bits, iterations)
    Fp := FiniteField(p_val);
    
    // We use d=1 for strict compliance with your 9M+5S formula
    a_2 := Fp!12345; 
    DummyPrj := [Fp!0, Fp!0, Fp!0]; 
    
    // Doubling formulas for P(X,Y,Z)
    DedicatedDblJq_Prj := function(Dummy, a_2_val, PrjPoint)
        X_1:=PrjPoint[1]; Y_1:=PrjPoint[2]; Z_1:=PrjPoint[3];
        A:=2*Z_1^2; B:=X_1^2; C:=Y_1^2; E:=2*X_1*Y_1; 
        F:=a_2_val*B; G:=C-F; H:=A-G;
        return [E*H, 2*C*G - H^2, H^2];
    end function;

    // Unified Addition (Handles arbitrary Z1, Z2) - Used inside Montgomery Ladder
    UnifiedAddJq_Prj := function(Dummy, a_2_val, PrjPoint1, PrjPoint2)
        X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
        X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; Z_2:=PrjPoint2[3];
        A:=Z_1*Z_2; B:=Y_1*Y_2; C:=X_1*X_2; 
        D:=(X_1+Y_1)*(X_2+Y_2)-B-C;
        F:=(X_1*Z_2)^2; G:=(Z_1*X_2)^2; H:=A^2; J:=C^2; 
        L:=H-J; M:=L^2; N:=(A+L)^2-H-M; P:=(A+C)^2-H-J;
        return [D*N, 2*A*(B+a_2_val*C)*(H+J)+2*P*(F+G), 2*M];
    end function;

    // Ultra-fast addition (9M+5S), requiring Z2=1 - Used for outside operations
    DedicatedAddJq_Prj := function(Dummy, a_2_val, PrjPoint1, PrjPoint2)
        X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
        X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; 
        A:=Z_1; B:=Y_1*Y_2; C:=X_1*X_2;
        D:=X_1*Y_2 - X_2*Y_1; F:=X_1^2; G:=(Z_1*X_2)^2;
        return [D*(F-G), (B-a_2_val*C)*(F+G)-2*C*(A^2+C^2), A*D^2];
    end function;

    // Strict Constant-Time Montgomery Ladder
    ScalarMulJqMontgomery_Prj := function(DummyPrj, a_2_val, Scalar, PrjPoint, num_bits)
        if Scalar eq 0 then return DummyPrj; end if;
        
        bits := IntegerToSequence(Scalar, 2);
        Reverse(~bits); // Ensure bits are processed MSB to LSB
        
        // Since Scalar is generated as Random(2^(num_bits-1), ...), bits[1] is always 1.
        // We initialize R0 = P, R1 = 2P and start loop from the second MSB.
        R0 := PrjPoint;
        R1 := DedicatedDblJq_Prj(DummyPrj, a_2_val, PrjPoint);
        
        for i in [2..num_bits] do
            if bits[i] eq 0 then
                R1 := UnifiedAddJq_Prj(DummyPrj, a_2_val, R0, R1);
                R0 := DedicatedDblJq_Prj(DummyPrj, a_2_val, R0);
            else
                R0 := UnifiedAddJq_Prj(DummyPrj, a_2_val, R0, R1);
                R1 := DedicatedDblJq_Prj(DummyPrj, a_2_val, R1);
            end if;
        end for;
        
        return R0;
    end function;

    ObfuscateCls := 0; 
    DeObfuscateCls := 0; 

    // Dummy base points to load the Jacobi arithmetic
    P_base := [Fp!123, Fp!456, Fp!1]; 
    P_u := [Fp!789, Fp!101, Fp!1];
    P_v := [Fp!112, Fp!131, Fp!1];

    // Generate scalars. MSB is strictly 1 to enforce exact bit-length execution.
    BlindScalar := Random(2^(blinding_bits-1), 2^blinding_bits - 1);
    FastCofactor := Random(2^(clearing_bits-1), 2^clearing_bits - 1);

    for i in [1..iterations] do
        // === OBFUSCATION ===
        t := ClockCycles(); 
        R_cof_J := ScalarMulJqMontgomery_Prj(DummyPrj, a_2, BlindScalar, P_base, blinding_bits);
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
        P_cleared_J := ScalarMulJqMontgomery_Prj(DummyPrj, a_2, FastCofactor, P_sum_J, clearing_bits);
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
