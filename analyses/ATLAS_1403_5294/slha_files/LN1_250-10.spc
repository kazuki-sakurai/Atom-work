#  ISAJET SUSY parameters in SUSY Les Houches Accord 2 format
#  Created by ISALHA 2.0 Last revision: C. Balazs 21 Apr 2009
Block SPINFO   # Program information
     1    ISASUSY from ISAJET          # Spectrum Calculator
     2   7.80   29-OCT-2009 12:50:36   # Version number
Block MODSEL   # Model selection
     1     0   # general MSSM-24 model                             
Block SMINPUTS   # Standard Model inputs
     1     1.28000000E+02   # alpha_em^(-1)
     2     1.16570000E-05   # G_Fermi
     3     1.19999997E-01   # alpha_s(M_Z)
     4     9.11699982E+01   # m_{Z}(pole)
     5     4.51715469E+00   # m_{b}(m_{b})
     6     1.72500000E+02   # m_{top}(pole)
     7     1.77699995E+00   # m_{tau}(pole)
Block MINPAR   # SUSY breaking input parameters.
     3     1.00000000E+01   # tan(beta)
     4     2.50000000E+03   # mu
Block MSOFT   # Non-universal SUSY breaking parameters 
     1     7.50000000E+01   # M_1(Q)          
     2     2.50000000E+03   # M_2(Q)          
     3     3.00000000E+03   # M_3(Q)          
    11     0.00000000E+00   #     At          
    12     0.00000000E+00   #     Ab          
    13     0.00000000E+00   #   Atau          
    23     2.50000000E+03   #    MUE          
    25     1.00000000E+01   #    TB           
    26     2.50000000E+03   #    MA0          
    31     3.02000000E+02   # MeL(Q)          
    32     3.02000000E+02   # MmuL(Q)         
    33     2.50000000E+03   # MtauL(Q)        
    34     3.02000000E+02   # MeR(Q)          
    35     3.02000000E+02   # MmuR(Q)         
    36     2.50000000E+03   # MtauR(Q)        
    41     2.50000000E+03   # MqL1(Q)         
    42     2.50000000E+03   # MqL2(Q)         
    43     2.50000000E+03   # MqL3(Q)         
    44     2.20000000E+03   # MuR(Q)          
    45     2.20000000E+03   # McR(Q)          
    46     2.50000000E+03   # MtR(Q)          
    47     2.20000000E+03   # MdR(Q)          
    48     2.20000000E+03   # MsR(Q)          
    49     2.50000000E+03   # MbR(Q)          
Block MASS   # Scalar and gaugino mass spectrum
#  PDG code   mass                 particle
         6     1.72500000E+02   #  top
        24     8.04229965E+01   #  W^+
        25     1.13288345E+02   #  h^0            
        35     2.51651904E+03   #  H^0            
        36     2.50000000E+03   #  A^0            
        37     2.51772217E+03   #  H^+            
   1000001     2.50068872E+03   #  dnl            
   1000002     2.49943726E+03   #  upl            
   1000003     2.50068872E+03   #  stl            
   1000004     2.49943774E+03   #  chl            
   1000005     2.48715942E+03   #  b1             
   1000006     2.49648804E+03   #  t1             
   1000011     250.   #  el-            
   1000013     250.   #  mul-           
   2000011     250.   #  er-            
   2000013     250.   #  mur-           
   1000022     10.    #  z1ss              
   1000012     2.50000000E+03   #  nuel
   1000014     2.50000000E+03   #  numl
   1000015     2.49166260E+03   #  tau1           
   1000016     2.50000000E+03   #  nutl
   1000021     3.00000000E+03   #  glss           
   1000023     2.43826221E+03   #  z2ss           
   1000024     2.43827441E+03   #  w1ss           
   1000025    -2.50082446E+03   #  z3ss           
   1000035     2.56273926E+03   #  z4ss           
   1000037     2.56276270E+03   #  w2ss           
   2000001     2.20014307E+03   #  dnr            
   2000002     2.19971362E+03   #  upr            
   2000003     2.20014307E+03   #  str            
   2000004     2.19971411E+03   #  chr            
   2000005     2.51358838E+03   #  b2             
   2000006     2.51009106E+03   #  t2             
   2000015     2.50912305E+03   #  tau2           
Block ALPHA   # Effective Higgs mixing parameter
         -9.99946222E-02   # alpha
Block STOPMIX   # stop mixing matrix
  1  1     7.15108216E-01   # O_{11}
  1  2     6.99013770E-01   # O_{12}
  2  1    -6.99013770E-01   # O_{21}
  2  2     7.15108216E-01   # O_{22}
Block SBOTMIX   # sbottom mixing matrix
  1  1     6.99541688E-01   # O_{11}
  1  2     7.14591742E-01   # O_{12}
  2  1    -7.14591742E-01   # O_{21}
  2  2     6.99541688E-01   # O_{22}
Block STAUMIX   # stau mixing matrix
  1  1     7.05921233E-01   # O_{11}
  1  2     7.08290339E-01   # O_{12}
  2  1    -7.08290339E-01   # O_{21}
  2  2     7.05921233E-01   # O_{22}
Block NMIX   # neutralino mixing matrix
  1  1     9.99841452E-01   #
  1  2    -1.33783789E-04   #
  1  3     1.76630989E-02   #
  1  4    -2.28811800E-03   #
  2  1     1.01508610E-02   #
  2  2     7.06965804E-01   #
  2  3    -5.05136728E-01   #
  2  4     4.94907349E-01   #
  3  1    -1.08680278E-02   #
  3  2     1.01849828E-02   #
  3  3     7.06886351E-01   #
  3  4     7.07170606E-01   #
  4  1     9.80219152E-03   #
  4  2    -7.07174480E-01   #
  4  3    -4.94810194E-01   #
  4  4     5.04946887E-01   #
Block UMIX   # chargino U mixing matrix
  1  1    -6.99869037E-01   # U_{11}
  1  2     7.14271188E-01   # U_{12}
  2  1    -7.14271188E-01   # U_{21}
  2  2    -6.99869037E-01   # U_{22}
Block VMIX   # chargino V mixing matrix
  1  1    -7.14271009E-01   # V_{11}
  1  2     6.99869215E-01   # V_{12}
  2  1    -6.99869215E-01   # V_{21}
  2  2    -7.14271009E-01   # V_{22}
Block HMIX Q=  9.11699982E+01   # Higgs mixing parameters
     1     2.50000000E+03   # mu(Q)
     2     1.00000000E+01   # tan(beta)
Block AU   #
  3  3     0.00000000E+00   # A_t
Block AD   #
  3  3     0.00000000E+00   # A_b
Block AE   #
  3  3     0.00000000E+00   # A_tau
#  ISAJET decay tables in SUSY Les Houches accord format 
#  Created by ISALHD. Last revision: C. Balazs, 2005 May 25
#  isajet_helpers.f
Block DCINFO                           # Program information
     1   ISASUGRA from ISAJET          # Spectrum Calculator
     2   7.80   29-OCT-2009 12:50:36   # Version number
#         PDG         Width
DECAY         6  1.14680779E+00   # TP    decays
#          BR          NDA       ID1       ID2       ID3       ID4
      3.33333343E-01    3            2        -1         5             # TP     -->  UP     DB     BT          
      3.33333343E-01    3            4        -3         5             # TP     -->  CH     SB     BT          
      1.11111112E-01    3          -11        12         5             # TP     -->  E+     NUE    BT          
      1.11111112E-01    3          -13        14         5             # TP     -->  MU+    NUM    BT          
      1.11111112E-01    3          -15        16         5             # TP     -->  TAU+   NUT    BT          
#         PDG         Width
DECAY   1000021  9.86177673E+01   # GLSS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.94477399E-02    2     -1000002         2                       # GLSS   -->  UBL    UP                 
      2.94477399E-02    2      1000002        -2                       # GLSS   -->  UPL    UB                 
      2.93140132E-02    2     -1000001         1                       # GLSS   -->  DBL    DN                 
      2.93140132E-02    2      1000001        -1                       # GLSS   -->  DNL    DB                 
      6.72896355E-02    2     -2000002         2                       # GLSS   -->  UBR    UP                 
      6.72896355E-02    2      2000002        -2                       # GLSS   -->  UPR    UB                 
      6.72285408E-02    2     -2000001         1                       # GLSS   -->  DBR    DN                 
      6.72285408E-02    2      2000001        -1                       # GLSS   -->  DNR    DB                 
      2.93140095E-02    2     -1000003         3                       # GLSS   -->  SBL    ST                 
      2.93140095E-02    2      1000003        -3                       # GLSS   -->  STL    SB                 
      6.72285408E-02    2     -2000003         3                       # GLSS   -->  SBR    ST                 
      6.72285408E-02    2      2000003        -3                       # GLSS   -->  STR    SB                 
      2.94476114E-02    2     -1000004         4                       # GLSS   -->  CBL    CH                 
      2.94476114E-02    2      1000004        -4                       # GLSS   -->  CHL    CB                 
      6.72895089E-02    2     -2000004         4                       # GLSS   -->  CBR    CH                 
      6.72895089E-02    2      2000004        -4                       # GLSS   -->  CHR    CB                 
      3.04430090E-02    2     -1000005         5                       # GLSS   -->  BB1    BT                 
      3.04430090E-02    2      1000005        -5                       # GLSS   -->  BT1    BB                 
      2.82608978E-02    2     -2000005         5                       # GLSS   -->  BB2    BT                 
      2.82608978E-02    2      2000005        -5                       # GLSS   -->  BT2    BB                 
      1.77998170E-02    2     -1000006         6                       # GLSS   -->  TB1    TP                 
      1.77998170E-02    2      1000006        -6                       # GLSS   -->  TP1    TB                 
      3.69367413E-02    2     -2000006         6                       # GLSS   -->  TB2    TP                 
      3.69367413E-02    2      2000006        -6                       # GLSS   -->  TP2    TB                 
#         PDG         Width
DECAY   1000002  3.89320731E-01   # UPL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.03815806E-01    2      1000022         2                       # UPL    -->  Z1SS   UP                 
      3.17434743E-02    2      1000023         2                       # UPL    -->  Z2SS   UP                 
      6.44407347E-02    2      1000024         1                       # UPL    -->  W1SS+  DN                 
#         PDG         Width
DECAY   1000001  3.90860349E-01   # DNL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.03343320E-01    2      1000022         1                       # DNL    -->  Z1SS   DN                 
      3.25482525E-02    2      1000023         1                       # DNL    -->  Z2SS   DN                 
      6.41084686E-02    2     -1000024         2                       # DNL    -->  W1SS-  UP                 
#         PDG         Width
DECAY   1000003  3.90854329E-01   # STL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.03357208E-01    2      1000022         3                       # STL    -->  Z1SS   ST                 
      3.25487517E-02    2      1000023         3                       # STL    -->  Z2SS   ST                 
      6.40940741E-02    2     -1000024         4                       # STL    -->  W1SS-  CH                 
#         PDG         Width
DECAY   1000004  3.89318168E-01   # CHL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.03821468E-01    2      1000022         4                       # CHL    -->  Z1SS   CH                 
      3.17362323E-02    2      1000023         4                       # CHL    -->  Z2SS   CH                 
      6.44422844E-02    2      1000024         3                       # CHL    -->  W1SS+  ST                 
#         PDG         Width
DECAY   1000005  9.27733898E-01   # BT1   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.93637443E-01    2      1000022         5                       # BT1    -->  Z1SS   BT                 
      6.36259001E-03    2      1000023         5                       # BT1    -->  Z2SS   BT                 
#         PDG         Width
DECAY   1000006  2.98494864E+00   # TP1   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.82107222E-01    2      1000022         6                       # TP1    -->  Z1SS   TP                 
      1.78928263E-02    2      1000024         5                       # TP1    -->  W1SS+  BT                 
#         PDG         Width
DECAY   2000002  4.95951414E+00   # UPR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022         2                       # UPR    -->  Z1SS   UP                 
#         PDG         Width
DECAY   2000001  1.24012172E+00   # DNR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022         1                       # DNR    -->  Z1SS   DN                 
#         PDG         Width
DECAY   2000003  1.24012172E+00   # STR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022         3                       # STR    -->  Z1SS   ST                 
#         PDG         Width
DECAY   2000004  4.95951176E+00   # CHR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022         4                       # CHR    -->  Z1SS   CH                 
#         PDG         Width
DECAY   2000005  8.47375393E-01   # BT2   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.93072987E-01    2      1000022         5                       # BT2    -->  Z1SS   BT                 
      6.88519049E-03    2      1000023         5                       # BT2    -->  Z2SS   BT                 
      4.18870950E-05    2      1000025         5                       # BT2    -->  Z3SS   BT                 
#         PDG         Width
DECAY   2000006  3.01298022E+00   # TP2   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      5.74856589E-04    2      1000024         5                       # TP2    -->  W1SS+  BT                 
      9.99425113E-01    2      1000022         6                       # TP2    -->  Z1SS   TP                 
#         PDG         Width
DECAY   1000011  3.43114734E-01   # EL-   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        11                       # EL-    -->  Z1SS   E-                 
#         PDG         Width
DECAY   1000013  3.43114674E-01   # MUL-  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        13                       # MUL-   -->  Z1SS   MU-                
#         PDG         Width
DECAY   1000015  7.96046400E+00   # TAU1- decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97645140E-01    2      1000022        15                       # TAU1-  -->  Z1SS   TAU-               
      8.04700307E-04    2      1000023        15                       # TAU1-  -->  Z2SS   TAU-               
      1.55017537E-03    2     -1000024        16                       # TAU1-  -->  W1SS-  NUT                
#         PDG         Width
DECAY   1000012  3.28712136E-01   # NUEL  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        12                       # NUEL   -->  Z1SS   NUE                
#         PDG         Width
DECAY   1000014  3.28712136E-01   # NUML  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        14                       # NUML   -->  Z1SS   NUM                
#         PDG         Width
DECAY   1000016  3.20997262E+00   # NUTL  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.88395512E-01    2      1000022        16                       # NUTL   -->  Z1SS   NUT                
      3.73935141E-03    2      1000023        16                       # NUTL   -->  Z2SS   NUT                
      7.86515046E-03    2      1000024        15                       # NUTL   -->  W1SS+  TAU-               
#         PDG         Width
DECAY   2000011  1.37041867E+00   # ER-   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        11                       # ER-    -->  Z1SS   E-                 
#         PDG         Width
DECAY   2000013  1.37041831E+00   # MUR-  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022        13                       # MUR-   -->  Z1SS   MU-                
#         PDG         Width
DECAY   2000015  7.93793964E+00   # TAU2- decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97784853E-01    2      1000022        15                       # TAU2-  -->  Z1SS   TAU-               
      7.74544664E-04    2      1000023        15                       # TAU2-  -->  Z2SS   TAU-               
      1.44060096E-03    2     -1000024        16                       # TAU2-  -->  W1SS-  NUT                
#         PDG         Width
DECAY   1000023  2.06882324E+01   # Z2SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.39823733E-02    2      1000022        23                       # Z2SS   -->  Z1SS   Z0                 
      1.10382571E-05    3      1000022         2        -2             # Z2SS   -->  Z1SS   UP     UB          
      1.08865861E-05    3      1000022         1        -1             # Z2SS   -->  Z1SS   DN     DB          
      1.08865861E-05    3      1000022         3        -3             # Z2SS   -->  Z1SS   ST     SB          
      1.10382298E-05    3      1000022         4        -4             # Z2SS   -->  Z1SS   CH     CB          
      2.37160875E-05    3      1000022         5        -5             # Z2SS   -->  Z1SS   BT     BB          
      4.12496302E-05    3      1000022        15       -15             # Z2SS   -->  Z1SS   TAU-   TAU+        
      3.25276160E-05    3      1000022        16       -16             # Z2SS   -->  Z1SS   NUT    ANUT        
      2.37875786E-02    2      1000022        25                       # Z2SS   -->  Z1SS   HL0                
      1.22015417E-01    2      1000011       -11                       # Z2SS   -->  EL-    E+                 
      1.22015417E-01    2     -1000011        11                       # Z2SS   -->  EL+    E-                 
      1.22015417E-01    2      1000013       -13                       # Z2SS   -->  MUL-   MU+                
      1.22015417E-01    2     -1000013        13                       # Z2SS   -->  MUL+   MU-                
      2.99244512E-05    2      2000011       -11                       # Z2SS   -->  ER-    E+                 
      2.99244512E-05    2     -2000011        11                       # Z2SS   -->  ER+    E-                 
      2.99244512E-05    2      2000013       -13                       # Z2SS   -->  MUR-   MU+                
      2.99244512E-05    2     -2000013        13                       # Z2SS   -->  MUR+   MU-                
      1.18476830E-01    2      1000012       -12                       # Z2SS   -->  NUEL   ANUE               
      1.18476830E-01    2     -1000012        12                       # Z2SS   -->  ANUEL  NUE                
      1.18476830E-01    2      1000014       -14                       # Z2SS   -->  NUML   ANUM               
      1.18476830E-01    2     -1000014        14                       # Z2SS   -->  ANUML  NUM                
#         PDG         Width
DECAY   1000025  1.63191819E+00   # Z3SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.95263317E-04    3     -1000024         2        -1             # Z3SS   -->  W1SS-  UP     DB          
      6.50877701E-05    3     -1000024        12       -11             # Z3SS   -->  W1SS-  NUE    E+          
      6.50877701E-05    3     -1000024        14       -13             # Z3SS   -->  W1SS-  NUM    MU+         
      1.95263317E-04    3      1000024         1        -2             # Z3SS   -->  W1SS+  DN     UB          
      6.50877701E-05    3      1000024        11       -12             # Z3SS   -->  W1SS+  E-     ANUE        
      6.50877701E-05    3      1000024        13       -14             # Z3SS   -->  W1SS+  MU-    ANUM        
      1.95263317E-04    3     -1000024         4        -3             # Z3SS   -->  W1SS-  CH     SB          
      6.50877701E-05    3     -1000024        16       -15             # Z3SS   -->  W1SS-  NUT    TAU+        
      1.95263317E-04    3      1000024         3        -4             # Z3SS   -->  W1SS+  ST     CB          
      6.50877701E-05    3      1000024        15       -16             # Z3SS   -->  W1SS+  TAU-   ANUT        
      6.26748085E-01    2      1000022        23                       # Z3SS   -->  Z1SS   Z0                 
      2.65897281E-04    3      1000022         5        -5             # Z3SS   -->  Z1SS   BT     BB          
      1.17656506E-04    3      1000022        15       -15             # Z3SS   -->  Z1SS   TAU-   TAU+        
      5.01478389E-05    3      1000023         2        -2             # Z3SS   -->  Z2SS   UP     UB          
      6.46659682E-05    3      1000023         1        -1             # Z3SS   -->  Z2SS   DN     DB          
      6.46659682E-05    3      1000023         3        -3             # Z3SS   -->  Z2SS   ST     SB          
      5.01478389E-05    3      1000023         4        -4             # Z3SS   -->  Z2SS   CH     CB          
      6.06991998E-05    3      1000023         5        -5             # Z3SS   -->  Z2SS   BT     BB          
      1.46694438E-05    3      1000023        11       -11             # Z3SS   -->  Z2SS   E-     E+          
      1.46694438E-05    3      1000023        13       -13             # Z3SS   -->  Z2SS   MU-    MU+         
      1.44669639E-05    3      1000023        15       -15             # Z3SS   -->  Z2SS   TAU-   TAU+        
      2.91875785E-05    3      1000023        12       -12             # Z3SS   -->  Z2SS   NUE    ANUE        
      2.91875785E-05    3      1000023        14       -14             # Z3SS   -->  Z2SS   NUM    ANUM        
      2.91875785E-05    3      1000023        16       -16             # Z3SS   -->  Z2SS   NUT    ANUT        
      3.65796357E-01    2      1000022        25                       # Z3SS   -->  Z1SS   HL0                
      3.14360768E-05    2      2000002        -2                       # Z3SS   -->  UPR    UB                 
      3.14360768E-05    2     -2000002         2                       # Z3SS   -->  UBR    UP                 
      3.14357130E-05    2      2000004        -4                       # Z3SS   -->  CHR    CB                 
      3.14357130E-05    2     -2000004         4                       # Z3SS   -->  CBR    CH                 
      3.44251457E-05    2      1000005        -5                       # Z3SS   -->  BT1    BB                 
      3.44251457E-05    2     -1000005         5                       # Z3SS   -->  BB1    BT                 
      5.55175830E-05    2      1000011       -11                       # Z3SS   -->  EL-    E+                 
      5.55175830E-05    2     -1000011        11                       # Z3SS   -->  EL+    E-                 
      5.55175830E-05    2      1000013       -13                       # Z3SS   -->  MUL-   MU+                
      5.55175830E-05    2     -1000013        13                       # Z3SS   -->  MUL+   MU-                
      4.46733262E-04    2      2000011       -11                       # Z3SS   -->  ER-    E+                 
      4.46733262E-04    2     -2000011        11                       # Z3SS   -->  ER+    E-                 
      4.46733262E-04    2      2000013       -13                       # Z3SS   -->  MUR-   MU+                
      4.46733262E-04    2     -2000013        13                       # Z3SS   -->  MUR+   MU-                
      8.18824221E-04    2      1000012       -12                       # Z3SS   -->  NUEL   ANUE               
      8.18824221E-04    2     -1000012        12                       # Z3SS   -->  ANUEL  NUE                
      8.18824221E-04    2      1000014       -14                       # Z3SS   -->  NUML   ANUM               
      8.18824221E-04    2     -1000014        14                       # Z3SS   -->  ANUML  NUM                
#         PDG         Width
DECAY   1000035  2.30763988E+01   # Z4SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.25168653E-02    2      1000024       -24                       # Z4SS   -->  W1SS+  W-                 
      2.25168653E-02    2     -1000024        24                       # Z4SS   -->  W1SS-  W+                 
      1.38013521E-02    2      1000022        23                       # Z4SS   -->  Z1SS   Z0                 
      2.22446490E-02    2      1000022        25                       # Z4SS   -->  Z1SS   HL0                
      8.30612204E-04    2      1000002        -2                       # Z4SS   -->  UPL    UB                 
      8.30612204E-04    2     -1000002         2                       # Z4SS   -->  UBL    UP                 
      8.06641765E-04    2      1000001        -1                       # Z4SS   -->  DNL    DB                 
      8.06641765E-04    2     -1000001         1                       # Z4SS   -->  DBL    DN                 
      8.06637632E-04    2      1000003        -3                       # Z4SS   -->  STL    SB                 
      8.06637632E-04    2     -1000003         3                       # Z4SS   -->  SBL    ST                 
      8.30414996E-04    2      1000004        -4                       # Z4SS   -->  CHL    CB                 
      8.30414996E-04    2     -1000004         4                       # Z4SS   -->  CBL    CH                 
      3.55845666E-04    2      1000005        -5                       # Z4SS   -->  BT1    BB                 
      3.55845666E-04    2     -1000005         5                       # Z4SS   -->  BB1    BT                 
      3.87229084E-04    2      2000005        -5                       # Z4SS   -->  BT2    BB                 
      3.87229084E-04    2     -2000005         5                       # Z4SS   -->  BB2    BT                 
      1.11869477E-01    2      1000011       -11                       # Z4SS   -->  EL-    E+                 
      1.11869477E-01    2     -1000011        11                       # Z4SS   -->  EL+    E-                 
      1.11869477E-01    2      1000013       -13                       # Z4SS   -->  MUL-   MU+                
      1.11869477E-01    2     -1000013        13                       # Z4SS   -->  MUL+   MU-                
      2.63738439E-05    2      2000011       -11                       # Z4SS   -->  ER-    E+                 
      2.63738439E-05    2     -2000011        11                       # Z4SS   -->  ER+    E-                 
      2.63738439E-05    2      2000013       -13                       # Z4SS   -->  MUR-   MU+                
      2.63738439E-05    2     -2000013        13                       # Z4SS   -->  MUR+   MU-                
      1.24974875E-04    2      1000015       -15                       # Z4SS   -->  TAU1-  TAU+               
      1.24974875E-04    2     -1000015        15                       # Z4SS   -->  TAU1+  TAU-               
      1.33172958E-04    2      2000015       -15                       # Z4SS   -->  TAU2-  TAU+               
      1.33172958E-04    2     -2000015        15                       # Z4SS   -->  TAU2+  TAU-               
      1.15554094E-01    2      1000012       -12                       # Z4SS   -->  NUEL   ANUE               
      1.15554094E-01    2     -1000012        12                       # Z4SS   -->  ANUEL  NUE                
      1.15554094E-01    2      1000014       -14                       # Z4SS   -->  NUML   ANUM               
      1.15554094E-01    2     -1000014        14                       # Z4SS   -->  ANUML  NUM                
      2.84765731E-04    2      1000016       -16                       # Z4SS   -->  NUTL   ANUT               
      2.84765731E-04    2     -1000016        16                       # Z4SS   -->  ANUTL  NUT                
#         PDG         Width
DECAY   1000024  2.06899109E+01   # W1SS+ decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.19329631E-05    3      1000022         2        -1             # W1SS+  -->  Z1SS   UP     DB          
      2.19329377E-05    3      1000022         4        -3             # W1SS+  -->  Z1SS   CH     SB          
      7.37275768E-05    3      1000022       -15        16             # W1SS+  -->  Z1SS   TAU+   NUT         
      3.76332141E-02    2      1000022        24                       # W1SS+  -->  Z1SS   W+                 
      2.45718673E-01    2      1000012       -11                       # W1SS+  -->  NUEL   E+                 
      2.45718673E-01    2      1000014       -13                       # W1SS+  -->  NUML   MU+                
      2.35405907E-01    2     -1000011        12                       # W1SS+  -->  EL+    NUE                
      2.35405907E-01    2     -1000013        14                       # W1SS+  -->  MUL+   NUM                
#         PDG         Width
DECAY   1000037  2.29523487E+01   # W2SS+ decays
#          BR          NDA       ID1       ID2       ID3       ID4
      3.61176953E-02    2      1000022        24                       # W2SS+  -->  Z1SS   W+                 
      2.26307865E-02    2      1000023        24                       # W2SS+  -->  Z2SS   W+                 
      1.30780754E-05    3      1000025         2        -1             # W2SS+  -->  Z3SS   UP     DB          
      1.30780754E-05    3      1000025         4        -3             # W2SS+  -->  Z3SS   CH     SB          
      1.64534419E-03    2      1000002        -1                       # W2SS+  -->  UPL    DB                 
      1.64750393E-03    2     -1000001         2                       # W2SS+  -->  DBL    UP                 
      1.64530927E-03    2      1000004        -3                       # W2SS+  -->  CHL    SB                 
      1.71334704E-03    2     -1000003         4                       # W2SS+  -->  SBL    CH                 
      9.34604250E-05    2      1000006        -5                       # W2SS+  -->  TP1    BB                 
      2.80485977E-03    2      2000006        -5                       # W2SS+  -->  TP2    BB                 
      2.24147409E-01    2      1000012       -11                       # W2SS+  -->  NUEL   E+                 
      2.24147409E-01    2      1000014       -13                       # W2SS+  -->  NUML   MU+                
      5.60889544E-04    2      1000016       -15                       # W2SS+  -->  NUTL   TAU+               
      2.33016789E-01    2     -1000011        12                       # W2SS+  -->  EL+    NUE                
      2.33016789E-01    2     -1000013        14                       # W2SS+  -->  MUL+   NUM                
      2.57445528E-04    2     -1000015        16                       # W2SS+  -->  TAU1+  NUT                
      2.73301441E-04    2     -2000015        16                       # W2SS+  -->  TAU2+  NUT                
      1.62555240E-02    2      1000024        23                       # W2SS+  -->  W1SS+  Z0                 
#         PDG         Width
DECAY        25  3.96301225E-03   # HL0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.06451208E-04    2           13       -13                       # HL0    -->  MU-    MU+                
      5.90438209E-02    2           15       -15                       # HL0    -->  TAU-   TAU+               
      2.48025567E-03    2            3        -3                       # HL0    -->  ST     SB                 
      8.14806581E-01    2            5        -5                       # HL0    -->  BT     BB                 
      3.68169248E-02    2            4        -4                       # HL0    -->  CH     CB                 
      1.75749848E-03    2           22        22                       # HL0    -->  GM     GM                 
      3.86399738E-02    2           21        21                       # HL0    -->  GL     GL                 
      2.40926351E-03    3           24        11       -12             # HL0    -->  W+     E-     ANUE        
      2.40926351E-03    3           24        13       -14             # HL0    -->  W+     MU-    ANUM        
      2.40926351E-03    3           24        15       -16             # HL0    -->  W+     TAU-   ANUT        
      7.22779008E-03    3           24        -2         1             # HL0    -->  W+     UB     DN          
      7.22779008E-03    3           24        -4         3             # HL0    -->  W+     CB     ST          
      2.40926351E-03    3          -24       -11        12             # HL0    -->  W-     E+     NUE         
      2.40926351E-03    3          -24       -13        14             # HL0    -->  W-     MU+    NUM         
      2.40926351E-03    3          -24       -15        16             # HL0    -->  W-     TAU+   NUT         
      7.22779008E-03    3          -24         2        -1             # HL0    -->  W-     UP     DB          
      7.22779008E-03    3          -24         4        -3             # HL0    -->  W-     CH     SB          
      1.97093526E-04    3           23        12       -12             # HL0    -->  Z0     NUE    ANUE        
      1.97093526E-04    3           23        14       -14             # HL0    -->  Z0     NUM    ANUM        
      1.97093526E-04    3           23        16       -16             # HL0    -->  Z0     NUT    ANUT        
      9.91951965E-05    3           23        11       -11             # HL0    -->  Z0     E-     E+          
      9.91951965E-05    3           23        13       -13             # HL0    -->  Z0     MU-    MU+         
      9.91951965E-05    3           23        15       -15             # HL0    -->  Z0     TAU-   TAU+        
      3.39834311E-04    3           23         2        -2             # HL0    -->  Z0     UP     UB          
      3.39834311E-04    3           23         4        -4             # HL0    -->  Z0     CH     CB          
      4.37790324E-04    3           23         1        -1             # HL0    -->  Z0     DN     DB          
      4.37790324E-04    3           23         3        -3             # HL0    -->  Z0     ST     SB          
      4.37790324E-04    3           23         5        -5             # HL0    -->  Z0     BT     BB          
#         PDG         Width
DECAY        35  6.31993151E+00   # HH0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.85691291E-04    2           13       -13                       # HH0    -->  MU-    MU+                
      8.18261206E-02    2           15       -15                       # HH0    -->  TAU-   TAU+               
      3.32590239E-03    2            3        -3                       # HH0    -->  ST     SB                 
      7.39710867E-01    2            5        -5                       # HH0    -->  BT     BB                 
      1.73887298E-01    2            6        -6                       # HH0    -->  TP     TB                 
      3.17951708E-05    2           21        21                       # HH0    -->  GL     GL                 
      8.66344853E-05    2           24       -24                       # HH0    -->  W+     W-                 
      4.43907811E-05    2           23        23                       # HH0    -->  Z0     Z0                 
      3.03145905E-04    2      1000022   1000022                       # HH0    -->  Z1SS   Z1SS               
      3.31717834E-04    2           25        25                       # HH0    -->  HL0    HL0                
      1.59096617E-05    2      1000011  -1000011                       # HH0    -->  EL-    EL+                
      1.19236702E-05    2      2000011  -2000011                       # HH0    -->  ER-    ER+                
      1.59018273E-05    2      1000013  -1000013                       # HH0    -->  MUL-   MUL+               
      1.19168881E-05    2      2000013  -2000013                       # HH0    -->  MUR-   MUR+               
      5.54933977E-05    2      1000012  -1000012                       # HH0    -->  NUEL   ANUEL              
      5.54933977E-05    2      1000014  -1000014                       # HH0    -->  NUML   ANUML              
#         PDG         Width
DECAY        36  6.28864384E+00   # HA0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.85245769E-04    2           13       -13                       # HA0    -->  MU-    MU+                
      8.16986710E-02    2           15       -15                       # HA0    -->  TAU-   TAU+               
      3.32086813E-03    2            3        -3                       # HA0    -->  ST     SB                 
      7.39120722E-01    2            5        -5                       # HA0    -->  BT     BB                 
      1.75106540E-01    2            6        -6                       # HA0    -->  TP     TB                 
      6.20586216E-05    2           21        21                       # HA0    -->  GL     GL                 
      3.19922401E-04    2      1000022   1000022                       # HA0    -->  Z1SS   Z1SS               
      8.61386507E-05    2           25        23                       # HA0    -->  HL0    Z0                 
#         PDG         Width
DECAY        37  5.85075855E+00   # H+    decays
#          BR          NDA       ID1       ID2       ID3       ID4
      3.08765622E-04    2           14       -13                       # H+     -->  NUM    MU+                
      8.84351209E-02    2           16       -15                       # H+     -->  NUT    TAU+               
      3.33053549E-03    2            4        -3                       # H+     -->  CH     SB                 
      9.07234311E-01    2            6        -5                       # H+     -->  TP     BB                 
      4.57155169E-04    2      1000024   1000022                       # H+     -->  W1SS+  Z1SS               
      9.34312266E-05    2           25        24                       # H+     -->  HL0    W+                 
      7.03785263E-05    2     -1000011   1000012                       # H+     -->  EL+    NUEL               
      7.03664191E-05    2     -1000013   1000014                       # H+     -->  MUL+   NUML               