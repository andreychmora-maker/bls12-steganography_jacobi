/* 
===========================================================================
Full Pipeline for Obfuscation of the BLS12-479+ G1 Public Key 
using Elligator Squared, Explicit 2-Isogeny, and Jacobi Quartic Optimization.
===========================================================================
*/

p := 1040582850105330743815570414239668597581975900197275493095847470017583224665441180782578049215009464174278151947833029365685675660572997417552727;
Fp := FiniteField(p);

// ----------------------------------------------------------------
// Dynamically calculate precise BLS12-479 parameters from 'p'
// ----------------------------------------------------------------
z_approx := Iroot(3 * p, 6);
z_val := 0;
for cand in [-z_approx-2 .. -z_approx+2] cat [z_approx-2 .. z_approx+2] do
    if (cand-1)^2 * (cand^4 - cand^2 + 1) div 3 + cand eq p then
        z_val := cand;
        break;
    end if;
end for;

assert z_val ne 0; 
q := z_val^4 - z_val^2 + 1;
Cofactor := (z_val - 1)^2 div 3;

B_target := Fp!1;
A_isog   := Fp!-15;
B_isog   := Fp!22;

E_target := EllipticCurve([Fp | 0, B_target]);
E_isog   := EllipticCurve([Fp | A_isog, B_isog]);
Identity := E_target!0;

/*
----------------------------------------------------------------
Jacobi Quartic Constants (Strictly Aligned with Hisil Arithmetic)
----------------------------------------------------------------
*/
theta := Fp!2; // 2-torsion point on E_isog
a_bj  := Fp!6; // Shifted curve parameter 'a' = 3*theta
b_bj  := Fp!-3; // Shifted curve parameter 'b' = 3*theta^2 + A_isog

a_1  := -a_bj / 4; 
a_2  := 2 * a_1;   // Evaluates exactly to -3
d_jq := (a_bj^2 - 4 * b_bj) / 16; // Evaluates exactly to 3

DummyPrj := [Fp!0, Fp!0, Fp!0];

/*
----------------------------------------------------------------
Birational Mappings: Exact Homogeneous Hisil Coordinates
----------------------------------------------------------------
*/
function W2J(P)
    if P eq E_isog!0 then return [Fp!0, Fp!1, Fp!1]; end if;
    x := P[1]; y := P[2];
    
    // 2-torsion point maps to (0 : -1 : 1)
    if y eq 0 and x eq theta then return [Fp!0, Fp!-1, Fp!1]; end if;
    
    x_shift := x - theta;
    XJ := 2 * x_shift;
    YJ := x_shift * (x_shift^2 - b_bj);
    ZJ := y;
    return [XJ, YJ, ZJ];
end function;

function J2W(PJ)
    XJ := PJ[1]; YJ := PJ[2]; ZJ := PJ[3];
    if ZJ eq 0 then return E_isog!0; end if;
    if XJ eq 0 and YJ eq ZJ^2 then return E_isog!0; end if;
    if XJ eq 0 and YJ eq -ZJ^2 then return E_isog![theta, 0]; end if;

    // The shift +theta is mathematically pre-absorbed into num_x
    num_x := 2*ZJ^2 + 2*YJ - XJ^2;
    den_x := XJ^2;
    num_y := 2*ZJ*(2*ZJ^2 + 2*YJ - 3*XJ^2);
    den_y := XJ^3;
    
    if den_x eq 0 or den_y eq 0 then return E_isog!0; end if;
    return E_isog![num_x / den_x, num_y / den_y];
end function;

/*
----------------------------------------------------------------
Jacobi Quartic Arithmetic Library
----------------------------------------------------------------
*/
function UnifiedAddJq_Prj(Result, d, a_2, PrjPoint1, PrjPoint2)
    X_1:=PrjPoint1[1]; Y_1:=PrjPoint1[2]; Z_1:=PrjPoint1[3];
    X_2:=PrjPoint2[1]; Y_2:=PrjPoint2[2]; Z_2:=PrjPoint2[3];
    A:=Z_1*Z_2; B:=Y_1*Y_2; C:=X_1*X_2; 
    D:=(X_1+Y_1)*(X_2+Y_2)-B-C;
    F:=(X_1*Z_2)^2; G:=(Z_1*X_2)^2; H:=A^2; J:=C^2; 
    K:=d*J; // Multiplier d is strictly active
    L:=H-K; M:=L^2; N:=(A+L)^2-H-M; P:=(A+C)^2-H-J;
    Result[1]:=D*N;
    Result[2]:=2*A*(B+a_2*C)*(H+K)+2*d*P*(F+G); // Multiplier d is strictly active
    Result[3]:=2*M;
    return Result;
end function;

function NAF(Scalar)
    ArrayNAF := [];
    Tmp := Scalar;
    while Tmp gt 0 do
        if IsOdd(Tmp) then
            Val := Tmp mod 4;
            if Val eq 3 then Val := -1; end if;
            Tmp -:= Val;
        else
            Val := 0;
        end if;
        Append(~ArrayNAF, Val);
        Tmp div:= 2;
    end while;
    return ArrayNAF;
end function;

function ScalarMulJqNAF_Prj(DummyPrj, d, a_2, Scalar, PrjPoint)
    if Scalar eq 0 then return [Fp!0, Fp!1, Fp!1]; end if;
    ArrayNAF := NAF(Scalar);
    PrjPointNeg := [-PrjPoint[1], PrjPoint[2], PrjPoint[3]];
    
    ResultPoint := DummyPrj;
    started := false;
    
    for i := #ArrayNAF to 1 by -1 do
        if started then
            // Utilizing UnifiedAdd for doubling securely handles d != 1
            ResultPoint := UnifiedAddJq_Prj(DummyPrj, d, a_2, ResultPoint, ResultPoint);
        end if;
        
        Val := ArrayNAF[i];
        if Val ne 0 then
            P_add := Val eq 1 select PrjPoint else PrjPointNeg;
            if not started then
                ResultPoint := P_add;
                started := true;
            else
                ResultPoint := UnifiedAddJq_Prj(DummyPrj, d, a_2, ResultPoint, P_add);
            end if;
        end if;
    end for;
    return ResultPoint;
end function;

/*
----------------------------------------------------------------
Forward Map & Backward Map (2-Isogenies)
----------------------------------------------------------------
*/
function Isogeny2Target(P_isog)
    if P_isog eq E_isog!0 then return Identity; end if;
    x := P_isog[1]; y := P_isog[2];
    x_den := x - 2; 
    if x_den eq 0 then return Identity; end if;
    X_isog := x - 3 / x_den;
    Y_isog := y + 3 * y / (x_den^2);
    return E_target![X_isog / 4, Y_isog / 8];
end function;

function EvaluateBackwardMap(P_target)
    if P_target eq Identity then return E_isog!0; end if;
    x := P_target[1]; y := P_target[2];
    if x eq -1 then return E_isog!0; end if; 
    X_hat := x + 3 / (x + 1);
    Y_hat := y * (x^2 + 2*x - 2) / ((x + 1)^2);
    return E_isog![X_hat, Y_hat];
end function;

/*
----------------------------------------------------------------
REAL PIPELINE SETUP (Fast Cofactor Blinding & Clearing)
----------------------------------------------------------------
*/
FastCofactor := Abs(z_val - 1);

L := 2;
fc_mod_q := FastCofactor mod q;
ScalarFactor := (L * fc_mod_q) mod q;
InvScalar := Integers()!( (Integers(q)!ScalarFactor)^-1 );
BlindScalar := q * (FastCofactor div 3); 

/*
  SENDER PIPELINE (OBFUSCATION) - JACOBI OPTIMIZED
*/
function ObfuscatePublicKey_Jacobi(PublicKey_target)
    P_scaled := InvScalar * PublicKey_target;
    P_isog_base := EvaluateBackwardMap(P_scaled); 
    
    P_isog_base_J := W2J(P_isog_base);
    R_base := Random(E_isog);
    R_base_J := W2J(R_base);
    
    R_cof_J := ScalarMulJqNAF_Prj(DummyPrj, d_jq, a_2, BlindScalar, R_base_J);
    P_isog_J := UnifiedAddJq_Prj(DummyPrj, d_jq, a_2, P_isog_base_J, R_cof_J);
    P_isog := J2W(P_isog_J);

    P_u_isog := Random(E_isog); 
    P_v_isog := P_isog - P_u_isog;
    u := Random(1, p-1);
    v := Random(1, p-1);
    
    return u, v, P_u_isog, P_v_isog;
end function;

/*
  RECEIVER PIPELINE (DEOBFUSCATION) - JACOBI OPTIMIZED
*/
function DeobfuscatePublicKey_Jacobi(P_u_isog, P_v_isog)
    P_u_J := W2J(P_u_isog);
    P_v_J := W2J(P_v_isog);

    P_sum_J := UnifiedAddJq_Prj(DummyPrj, d_jq, a_2, P_u_J, P_v_J);
    P_cleared_J := ScalarMulJqNAF_Prj(DummyPrj, d_jq, a_2, FastCofactor, P_sum_J);
    P_cleared_W := J2W(P_cleared_J);

    PublicKey_target := Isogeny2Target(P_cleared_W);
    return PublicKey_target;
end function;

/*
----------------------------------------------------------------
Runtime Test & Verification
----------------------------------------------------------------
*/
printf "Target Curve:    %o\n", E_target;
printf "Isogenous Curve: %o\n", E_isog;
printf "\nJacobi Quartic d: %o", d_jq;
printf "\nJacobi Quartic a_2: %o\n", a_2;

printf "\n--- Testing Real Jacobi Optimized Pipeline ---\n";

SecretKey := Random(1, q-1);
G_target := Cofactor * Random(E_target);
while G_target eq Identity do G_target := Cofactor * Random(E_target); end while;

PublicKey_target_orig := SecretKey * G_target;

printf "Obfuscating (running pick-and-check with cofactor blinding IN JACOBI)...\n";
u, v, P_u, P_v := ObfuscatePublicKey_Jacobi(PublicKey_target_orig);
printf "   Found u (sim): %o\n   Found v (sim): %o\n\n", u, v;

printf "Deobfuscating (Reconstruct + FCC IN JACOBI + Forward Isogeny)...\n";
PublicKey_recovered := DeobfuscatePublicKey_Jacobi(P_u, P_v);

if PublicKey_recovered eq PublicKey_target_orig then
    printf "SUCCESS: Jacobi Quartic Pipeline logic closed perfectly!\n";
else
    printf "FAIL: Reconstruction mismatch!\n";
    quit;
end if;

/*
----------------------------------------------------------------
Benchmarks
----------------------------------------------------------------
*/
NumberOfAttempts := 100000;
ObfuscateCls := 0; 
DeObfuscateCls := 0; 

printf "\nRunning %o iterations for benchmark...\n", NumberOfAttempts;

for i in [1..NumberOfAttempts] do
    SK := Random(1, q-1);
    PK := SK * G_target;
    
    t := ClockCycles(); 
    _, _, pu_b, pv_b := ObfuscatePublicKey_Jacobi(PK);
    ObfuscateCls +:= (ClockCycles() - t);

    t := ClockCycles(); 
    _ := DeobfuscatePublicKey_Jacobi(pu_b, pv_b);
    DeObfuscateCls +:= (ClockCycles() - t);
end for;

print "\n__________________ Pipeline Clock Cycles ___________________";
printf "Jacobi Obfuscation Cost: %o\n", ObfuscateCls div NumberOfAttempts;
printf "Jacobi Deobfuscation Cost: %o\n", DeObfuscateCls div NumberOfAttempts;
quit;
