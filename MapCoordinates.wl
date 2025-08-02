(* ::Title:: *)
(*Verification of coordinate maps*)


(* ::Author:: *)
(*Davi C. Rodrigues (2025)*)


(* ::Section:: *)
(*SdSCosmo -> Kottler*)


(* ::Text:: *)
(*Here we use T and ell for Kottler coordinates.*)


ClearAll[c, m, f, ell, v, u, w, r, eta, t, T, A, \[CapitalLambda], a, ds2Kottler, ds2];

$Assumptions = {3m > c > 2m, m> 0, n \[Element] Integers};

ds2CosmoSdS = - A[a[t] r] Dt[t]^2 +  a[t]^2  A[a[t] r]^-1 Dt[r]^2 ; 

ds2CosmoSdSminus = - Aminus[a[t] r] Dt[t]^2 +  a[t]^2  Aminus[a[t] r]^-1 Dt[r]^2 ; 

A[x_] = 1/2 (f[x] + Sqrt[f[x]^2 + 4 x^2 \[CapitalLambda]/3]);

Aminus[x_] = 1/2 (f[x] - Sqrt[f[x]^2 + 4 x^2 \[CapitalLambda]/3]);

f[x_] = 1 - (2 m)/x - \[CapitalLambda]/3 x^2;

a[t_] = Exp[Sqrt[\[CapitalLambda]/3] t]; 

c /: Dt[c] = 0;
m /: Dt[m] = 0;
\[CapitalLambda] /: Dt[\[CapitalLambda]] = 0;

dtRule = Dt[t] -> (-1)^n Dt[T] - Sqrt[\[CapitalLambda]/3] ell/(f[ell] A[ell]) Dt[ell];

dtRuleMinus = Dt[t] -> (-1)^n Dt[T] - Sqrt[\[CapitalLambda]/3] ell/(f[ell] Aminus[ell]) Dt[ell];

rRule = r -> ell/a[t];

coefdwdw = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[T], 2] // Simplify
coefdudw = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[T]Dt[ell], 1] // Simplify
coefdudu = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[ell], 2] // Simplify

coefdwdw coefdudu // Simplify


(* ::Subsection:: *)
(*Negative branch solution verification*)


coefdwdw = Coefficient[((ds2CosmoSdSminus // Simplify) //. rRule // Simplify) //. dtRuleMinus // Simplify, Dt[T], 2] // Simplify
coefdudw = Coefficient[((ds2CosmoSdSminus // Simplify) //. rRule // Simplify) //. dtRuleMinus // Simplify, Dt[T]Dt[ell], 1] // Simplify
coefdudu = Coefficient[((ds2CosmoSdSminus // Simplify) //. rRule // Simplify) //. dtRuleMinus // Simplify, Dt[ell], 2] // Simplify

coefdwdw coefdudu // Simplify


(* ::Subsection::Closed:: *)
(*Conformal time quick verification (t here is conformal time)*)


ClearAll[c, m, f, ell, v, u, w, r, eta, t, T, A, \[CapitalLambda], a, ds2Kottler, ds2];

$Assumptions = {3m > c > 2m, m> 0, \[CapitalLambda] >0 , n \[Element] Integers};

ds2CosmoSdS = - a[t]^2 A[a[t] r] Dt[t]^2 +  a[t]^2  A[a[t] r]^-1 Dt[r]^2 ; 

A[x_] = 1/2 (f[x] + Sqrt[f[x]^2 + 4 x^2 \[CapitalLambda]/3]);

f[x_] = 1 - (2 m)/x - \[CapitalLambda]/3 x^2;

a[t_] = -Sqrt[(3/\[CapitalLambda])] 1/t; 

c /: Dt[c] = 0;
m /: Dt[m] = 0;
\[CapitalLambda] /: Dt[\[CapitalLambda]] = 0;

dtRule = Dt[t] -> Dt[T]/a[t] - Sqrt[\[CapitalLambda]/3] ell/(f[ell] A[ell]) 1/a[t] Dt[ell];

rRule = r ->  ell/a[t];

coefdwdw = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[T], 2] // Simplify
coefdudw = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[T]Dt[ell], 1] // Simplify
coefdudu = Coefficient[((ds2CosmoSdS // Simplify) //. rRule // Simplify) //. dtRule // Simplify, Dt[ell], 2] // Simplify

coefdwdw coefdudu // Simplify


(* ::Section:: *)
(*Kottler -> Lake *)


ClearAll[c, m, f, ell, v, u, w, r, eta, t, T, A, \[CapitalLambda], a, ds2Kottler, ds2];

$Assumptions = {3m > c > 2m, m> 0, n \[Element] Integers};

ds2Kottler = - f[ell] Dt[T]^2 +  f[ell]^-1 Dt[ell]^2 ; 

c /: Dt[c] = 0;
m /: Dt[m] = 0;
\[CapitalLambda] /: Dt[\[CapitalLambda]] = 0;

v[u_, w_] = u w (3 m - c)/c^2 + c; 

gg[u_,w_] = w/(3 u v[u,w] c^2) (-2 c (c-v[u,w])^2 + u w (2 c + v[u,w]));

f[x_] = 1 - (2 m)/x - \[CapitalLambda]/3 x^2;

\[CapitalLambda] = 3 (c - 2m)/c^3;

dtRule = Dt[T] -> (-1)^n (\!\(
\*SubscriptBox[\(\[PartialD]\), \(u\)]\ \(v[u, w]\)\)/f[v[u,w]] - 1/\!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\))Dt[u] + (-1)^n \!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\)/f[v[u,w]] Dt[w];
ellRule = ell -> v[u,w];

coefdwdw = Coefficient[((ds2Kottler // Simplify) //. dtRule // Simplify) //. ellRule // Simplify, Dt[w], 2] // Simplify
coefdudw = Coefficient[((ds2Kottler // Simplify) //. dtRule // Simplify) //. ellRule // Simplify, Dt[u]Dt[w], 1] // Simplify
coefdudu = Coefficient[((ds2Kottler // Simplify) //. dtRule // Simplify) //. ellRule // Simplify, Dt[u], 2] // Simplify

Simplify[coefdudu == gg[u,w]]


(* ::Section:: *)
(*SdSCosmo -> Lake*)


ClearAll[c, m, f, ell, v, u, w, r, eta, t, T, A, \[CapitalLambda], a, ds2Kottler, ds2];

$Assumptions = {3m > c > 2m, m> 0, \[CapitalLambda]>0 , n1 \[Element] Integers, n2 \[Element] Integers};

c /: Dt[c] = 0;
m /: Dt[m] = 0;
\[CapitalLambda] /: Dt[\[CapitalLambda]] = 0;

v[u_, w_] = u w (3 m - c)/c^2 + c; 

gg[u_,w_] = w/(3 u v[u,w] c^2) (-2 c (c-v[u,w])^2 + u w (2 c + v[u,w]));


A[x_] = 1/2 (f[x] + Sqrt[f[x]^2 + 4 x^2 \[CapitalLambda]/3]);

f[x_] = 1 - (2 m)/x - \[CapitalLambda]/3 x^2;

dtRuleCosmoKottler = Dt[t] -> (-1)^n1 Dt[T] - Sqrt[\[CapitalLambda]/3] ell/(f[ell] A[ell]) Dt[ell];
dtRuleKotlerLake = Dt[T] -> (-1)^n2 (\!\(
\*SubscriptBox[\(\[PartialD]\), \(u\)]\ \(v[u, w]\)\)/f[v[u,w]] - 1/\!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\))Dt[u] + (-1)^n2 \!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\)/f[v[u,w]] Dt[w];

ellRule = ell -> v[u,w];

\[CapitalLambda] = 3 (c - 2m)/c^3;

Echo["dt map from SdS-Cosmo to Lake, as explicit function of u and w: "];

((dtRuleCosmoKottler //. ellRule) //. dtRuleKotlerLake // Simplify) // Collect[#, {Dt[u], Dt[w]}] & (*FullSimplify takes a lot of time, but equal result*)


Clear[\[CapitalLambda], A, v, f];

dtRuleCosmoKottler = Dt[t] -> (-1)^n1 Dt[T] - Sqrt[\[CapitalLambda]/3] ell/(f[ell] A[ell]) Dt[ell];
dtRuleKotlerLake = Dt[T] -> (-1)^n2 (\!\(
\*SubscriptBox[\(\[PartialD]\), \(u\)]\ \(v[u, w]\)\)/f[v[u,w]] - 1/\!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\))Dt[u] + (-1)^n2 \!\(
\*SubscriptBox[\(\[PartialD]\), \(w\)]\(v[u, w]\)\)/f[v[u,w]] Dt[w];

ellRule = ell -> v[u,w];

Echo["dt map from SdS-Cosmo to Lake, implicitly defined: "];
((dtRuleCosmoKottler //. ellRule) //. dtRuleKotlerLake // Simplify) // Collect[#, {Dt[u], Dt[w]}] &



