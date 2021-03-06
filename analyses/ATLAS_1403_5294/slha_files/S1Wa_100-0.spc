#  ISAJET SUSY parameters in SUSY Les Houches Accord 2 format
#  Created by ISALHA 2.0 Last revision: C. Balazs 21 Apr 2009
Block SPINFO   # Program information
     1   ISASUGRA from ISAJET          # Spectrum Calculator
     2   7.81   26-APR-2011 16:50:35   # Version number
Block MODSEL   # Model selection
     1     1   # Minimal supergravity (mSUGRA) model               
Block SMINPUTS   # Standard Model inputs
     1     1.27840950E+02   # alpha_em^(-1)
     2     1.16570000E-05   # G_Fermi
     3     1.17200002E-01   # alpha_s(M_Z)
     4     9.11699982E+01   # m_{Z}(pole)
     5     4.19999981E+00   # m_{b}(m_{b})
     6     1.73199997E+02   # m_{top}(pole)
     7     1.77699995E+00   # m_{tau}(pole)
Block MINPAR   # SUSY breaking input parameters
     1     3.30000000E+02   # m_0
     2     5.00000000E+02   # m_{1/2}
     3     4.00000000E+01   # tan(beta)
     4     1.00000000E+00   # sign(mu)
     5    -5.00000000E+02   # A_0
Block EXTPAR   # Non-universal SUSY breaking parameters
     0     1.76318844E+16   # Input scale
Block MASS   # Scalar and gaugino mass spectrum
#  PDG code   mass                 particle
        25     1.16870972E+02   #  h^0            
        35     6.28744934E+02   #  H^0            
        36     6.24697754E+02   #  A^0            
        37     6.34586792E+02   #  H^+            
   1000001     1.10702246E+03   #  dnl            
   1000002     1.10390466E+03   #  upl            
   1000003     1.10702246E+03   #  stl            
   1000004     1.10390552E+03   #  chl            
   1000005     9.25675659E+02   #  b1             
   1000006     7.72417542E+02   #  t1             
   1000011     4.73480316E+02   #  el-            
   1000012     4.65190002E+02   #  nuel           
   1000013     4.73480316E+02   #  mul-           
   1000014     4.65190002E+02   #  numl           
   1000015     2.02426895E+02   #  tau1           
   1000016     4.17797943E+02   #  nutl           
   1000021     1.15880176E+03   #  glss           
   1000022     1.     #  z1ss           
   1000024     100.   #  w1ss              
   1000023     3.94819366E+02   #  z2ss           
   1000025    -6.98430725E+02   #  z3ss           
   1000035     7.08134155E+02   #  z4ss           
   1000037     7.08306274E+02   #  w2ss           
   2000001     1.06494458E+03   #  dnr            
   2000002     1.06849890E+03   #  upr            
   2000003     1.06494458E+03   #  str            
   2000004     1.06849988E+03   #  chr            
   2000005     9.94880493E+02   #  b2             
   2000006     1.00032068E+03   #  t2             
   2000011     3.80521210E+02   #  er-            
   2000013     3.80521271E+02   #  mur-           
   2000015     4.47397339E+02   #  tau2           
Block ALPHA   # Effective Higgs mixing parameter
         -2.53561661E-02   # alpha
Block STOPMIX   # stop mixing matrix
  1  1     4.52761382E-01   # O_{11}
  1  2     8.91631722E-01   # O_{12}
  2  1    -8.91631722E-01   # O_{21}
  2  2     4.52761382E-01   # O_{22}
Block SBOTMIX   # sbottom mixing matrix
  1  1     7.81507432E-01   # O_{11}
  1  2     6.23895943E-01   # O_{12}
  2  1    -6.23895943E-01   # O_{21}
  2  2     7.81507432E-01   # O_{22}
Block STAUMIX   # stau mixing matrix
  1  1     3.44204724E-01   # O_{11}
  1  2     9.38894629E-01   # O_{12}
  2  1    -9.38894629E-01   # O_{21}
  2  2     3.44204724E-01   # O_{22}
Block NMIX   # neutralino mixing matrix
  1  1     9.97291148E-01   #
  1  2    -9.72984079E-03   #
  1  3     6.95981011E-02   #
  1  4    -2.17275135E-02   #
  2  1     2.28502564E-02   #
  2  2     9.82255638E-01   #
  2  3    -1.60866529E-01   #
  2  4     9.36675519E-02   #
  3  1     3.32481824E-02   #
  3  2    -4.86115925E-02   #
  3  3    -7.04129338E-01   #
  3  4    -7.07625329E-01   #
  4  1     6.15066849E-02   #
  4  2    -1.80875763E-01   #
  4  3    -6.88098967E-01   #
  4  4     7.00014949E-01   #
Block UMIX   # chargino U mixing matrix
  1  1    -9.73583758E-01   # U_{11}
  1  2     2.28330076E-01   # U_{12}
  2  1    -2.28330076E-01   # U_{21}
  2  2    -9.73583758E-01   # U_{22}
Block VMIX   # chargino V mixing matrix
  1  1    -9.90937054E-01   # V_{11}
  1  2     1.34326920E-01   # V_{12}
  2  1    -1.34326920E-01   # V_{21}
  2  2    -9.90937054E-01   # V_{22}
Block GAUGE Q=  8.46515137E+02   #
     1     3.57499272E-01   # g`
     2     6.52483046E-01   # g_2
     3     1.22072554E+00   # g_3
Block YU Q=  8.46515137E+02   #
  3  3     8.59023213E-01   # y_t
Block YD Q=  8.46515137E+02   #
  3  3     4.69071597E-01   # y_b
Block YE Q=  8.46515137E+02   #
  3  3     4.30355370E-01   # y_tau
Block HMIX Q=  8.46515137E+02   # Higgs mixing parameters
     1     6.94892578E+02   # mu(Q)
     2     4.00000000E+01   # tan(beta)(M_GUT)
     3     2.50802032E+02   # Higgs vev at Q
     4     3.90247281E+05   # m_A^2(Q)
Block MSOFT Q=  8.46515137E+02   # DRbar SUSY breaking parameters
     1     2.10815384E+02   # M_1(Q)          
     2     3.91559723E+02   # M_2(Q)          
     3     1.10285413E+03   # M_3(Q)          
    21    -9.10002344E+04   # M^2_Hd          
    22    -4.78996688E+05   # M^2_Hu          
    31     4.68364502E+02   # MeL(Q)          
    32     4.68364502E+02   # MmuL(Q)         
    33     4.22006989E+02   # MtauL(Q)        
    34     3.77536163E+02   # MeR(Q)          
    35     3.77536163E+02   # MmuR(Q)         
    36     2.47031448E+02   # MtauR(Q)        
    41     1.05035156E+03   # MqL1(Q)         
    42     1.05035156E+03   # MqL2(Q)         
    43     9.10856995E+02   # MqL3(Q)         
    44     1.01528900E+03   # MuR(Q)          
    45     1.01528900E+03   # McR(Q)          
    46     7.86718323E+02   # MtR(Q)          
    47     1.01090594E+03   # MdR(Q)          
    48     1.01090594E+03   # MsR(Q)          
    49     9.32971741E+02   # MbR(Q)          
Block AU Q=  8.46515137E+02   #
  1  1    -1.02420349E+03   # A_u
  2  2    -1.02420349E+03   # A_c
  3  3    -1.02420349E+03   # A_t
Block AD Q=  8.46515137E+02   #
  1  1    -1.47017371E+03   # A_d
  2  2    -1.47017371E+03   # A_s
  3  3    -1.47017371E+03   # A_b
Block AE Q=  8.46515137E+02   #
  1  1    -5.32889954E+02   # A_e
  2  2    -5.32889954E+02   # A_mu
  3  3    -5.32889954E+02   # A_tau
#  ISAJET decay tables in SUSY Les Houches accord format
#  Created by ISALHD. Last revision: C. Balazs, 2005 May 25
Block DCINFO                           # Program information
     1   ISASUGRA from ISAJET          # Spectrum Calculator
     2   7.81   26-APR-2011 16:50:35   # Version number
#         PDG         Width
DECAY         6  1.16585374E+00   # TP    decays
#          BR          NDA       ID1       ID2       ID3       ID4
      3.33333343E-01    3            2        -1         5             # TP     -->  UP     DB     BT          
      3.33333343E-01    3            4        -3         5             # TP     -->  CH     SB     BT          
      1.11111119E-01    3          -11        12         5             # TP     -->  E+     NUE    BT          
      1.11111119E-01    3          -13        14         5             # TP     -->  MU+    NUM    BT          
      1.11111119E-01    3          -15        16         5             # TP     -->  TAU+   NUT    BT          
#         PDG         Width
DECAY   1000021  1.30986242E+01   # GLSS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.95049255E-03    3      1000024         5        -6             # GLSS   -->  W1SS+  BT     TB          
      2.95049255E-03    3     -1000024         6        -5             # GLSS   -->  W1SS-  TP     BB          
      9.63696337E-04    3      1000037         5        -6             # GLSS   -->  W2SS+  BT     TB          
      9.63696337E-04    3     -1000037         6        -5             # GLSS   -->  W2SS-  TP     BB          
      7.89559749E-07    2      1000022        21                       # GLSS   -->  Z1SS   GL                 
      5.53237274E-04    3      1000022         6        -6             # GLSS   -->  Z1SS   TP     TB          
      7.83049109E-06    2      1000023        21                       # GLSS   -->  Z2SS   GL                 
      2.50290148E-03    3      1000023         6        -6             # GLSS   -->  Z2SS   TP     TB          
      5.75613885E-05    2      1000025        21                       # GLSS   -->  Z3SS   GL                 
      5.80269319E-04    3      1000025         6        -6             # GLSS   -->  Z3SS   TP     TB          
      4.21289660E-05    2      1000035        21                       # GLSS   -->  Z4SS   GL                 
      2.48247455E-03    3      1000035         6        -6             # GLSS   -->  Z4SS   TP     TB          
      8.58569052E-03    2     -1000002         2                       # GLSS   -->  UBL    UP                 
      8.58569052E-03    2      1000002        -2                       # GLSS   -->  UPL    UB                 
      7.65922479E-03    2     -1000001         1                       # GLSS   -->  DBL    DN                 
      7.65922479E-03    2      1000001        -1                       # GLSS   -->  DNL    DB                 
      2.25102771E-02    2     -2000002         2                       # GLSS   -->  UBR    UP                 
      2.25102771E-02    2      2000002        -2                       # GLSS   -->  UPR    UB                 
      2.42396072E-02    2     -2000001         1                       # GLSS   -->  DBR    DN                 
      2.42396072E-02    2      2000001        -1                       # GLSS   -->  DNR    DB                 
      7.65917124E-03    2     -1000003         3                       # GLSS   -->  SBL    ST                 
      7.65917124E-03    2      1000003        -3                       # GLSS   -->  STL    SB                 
      2.42395587E-02    2     -2000003         3                       # GLSS   -->  SBR    ST                 
      2.42395587E-02    2      2000003        -3                       # GLSS   -->  STR    SB                 
      8.58296175E-03    2     -1000004         4                       # GLSS   -->  CBL    CH                 
      8.58296175E-03    2      1000004        -4                       # GLSS   -->  CHL    CB                 
      2.25074943E-02    2     -2000004         4                       # GLSS   -->  CBR    CH                 
      2.25074943E-02    2      2000004        -4                       # GLSS   -->  CHR    CB                 
      1.28321588E-01    2     -1000005         5                       # GLSS   -->  BB1    BT                 
      1.28321588E-01    2      1000005        -5                       # GLSS   -->  BT1    BB                 
      7.15424567E-02    2     -2000005         5                       # GLSS   -->  BB2    BT                 
      7.15424567E-02    2      2000005        -5                       # GLSS   -->  BT2    BB                 
      1.67124167E-01    2     -1000006         6                       # GLSS   -->  TB1    TP                 
      1.67124167E-01    2      1000006        -6                       # GLSS   -->  TP1    TB                 
#         PDG         Width
DECAY   1000002  1.06098747E+01   # UPL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.21955825E-02    2      1000022         2                       # UPL    -->  Z1SS   UP                 
      3.24108571E-01    2      1000023         2                       # UPL    -->  Z2SS   UP                 
      2.84771639E-04    2      1000025         2                       # UPL    -->  Z3SS   UP                 
      4.36333008E-03    2      1000035         2                       # UPL    -->  Z4SS   UP                 
      6.53577685E-01    2      1000024         1                       # UPL    -->  W1SS+  DN                 
      5.47003094E-03    2      1000037         1                       # UPL    -->  W2SS+  DN                 
#         PDG         Width
DECAY   1000001  1.05146809E+01   # DNL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.52797932E-02    2      1000022         1                       # DNL    -->  Z1SS   DN                 
      3.22955042E-01    2      1000023         1                       # DNL    -->  Z2SS   DN                 
      4.80528222E-04    2      1000025         1                       # DNL    -->  Z3SS   DN                 
      5.71126351E-03    2      1000035         1                       # DNL    -->  Z4SS   DN                 
      6.39454246E-01    2     -1000024         2                       # DNL    -->  W1SS-  UP                 
      1.61191169E-02    2     -1000037         2                       # DNL    -->  W2SS-  UP                 
#         PDG         Width
DECAY   1000003  1.05146542E+01   # STL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.52798314E-02    2      1000022         3                       # STL    -->  Z1SS   ST                 
      3.22955847E-01    2      1000023         3                       # STL    -->  Z2SS   ST                 
      4.80529445E-04    2      1000025         3                       # STL    -->  Z3SS   ST                 
      5.71127841E-03    2      1000035         3                       # STL    -->  Z4SS   ST                 
      6.39453471E-01    2     -1000024         4                       # STL    -->  W1SS-  CH                 
      1.61190201E-02    2     -1000037         4                       # STL    -->  W2SS-  CH                 
#         PDG         Width
DECAY   1000004  1.06098747E+01   # CHL   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.21955564E-02    2      1000022         4                       # CHL    -->  Z1SS   CH                 
      3.24107736E-01    2      1000023         4                       # CHL    -->  Z2SS   CH                 
      2.84770096E-04    2      1000025         4                       # CHL    -->  Z3SS   CH                 
      4.36330587E-03    2      1000035         4                       # CHL    -->  Z4SS   CH                 
      6.53578579E-01    2      1000024         3                       # CHL    -->  W1SS+  ST                 
      5.47004770E-03    2      1000037         3                       # CHL    -->  W2SS+  ST                 
#         PDG         Width
DECAY   1000005  7.96929169E+00   # BT1   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      5.20670116E-02    2      1000022         5                       # BT1    -->  Z1SS   BT                 
      2.50477284E-01    2      1000023         5                       # BT1    -->  Z2SS   BT                 
      4.36264798E-02    2      1000025         5                       # BT1    -->  Z3SS   BT                 
      3.05872597E-02    2      1000035         5                       # BT1    -->  Z4SS   BT                 
      4.20057476E-01    2     -1000024         6                       # BT1    -->  W1SS-  TP                 
      6.35062456E-02    2     -1000037         6                       # BT1    -->  W2SS-  TP                 
      1.39678210E-01    2          -24   1000006                       # BT1    -->  W-     TP1                
#         PDG         Width
DECAY   1000006  3.14273500E+00   # TP1   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      3.61652523E-01    2      1000022         6                       # TP1    -->  Z1SS   TP                 
      1.56794369E-01    2      1000023         6                       # TP1    -->  Z2SS   TP                 
      4.17066365E-01    2      1000024         5                       # TP1    -->  W1SS+  BT                 
      6.44867718E-02    2      1000037         5                       # TP1    -->  W2SS+  BT                 
#         PDG         Width
DECAY   2000002  2.22934365E+00   # UPR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97896373E-01    2      1000022         2                       # UPR    -->  Z1SS   UP                 
      4.21777600E-04    2      1000023         2                       # UPR    -->  Z2SS   UP                 
      3.92875401E-04    2      1000025         2                       # UPR    -->  Z3SS   UP                 
      1.28896697E-03    2      1000035         2                       # UPR    -->  Z4SS   UP                 
#         PDG         Width
DECAY   2000001  5.55181503E-01   # DNR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97913539E-01    2      1000022         1                       # DNR    -->  Z1SS   DN                 
      4.21113917E-04    2      1000023         1                       # DNR    -->  Z2SS   DN                 
      3.89176304E-04    2      1000025         1                       # DNR    -->  Z3SS   DN                 
      1.27619202E-03    2      1000035         1                       # DNR    -->  Z4SS   DN                 
#         PDG         Width
DECAY   2000003  5.55181503E-01   # STR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97913539E-01    2      1000022         3                       # STR    -->  Z1SS   ST                 
      4.21113917E-04    2      1000023         3                       # STR    -->  Z2SS   ST                 
      3.89176304E-04    2      1000025         3                       # STR    -->  Z3SS   ST                 
      1.27619202E-03    2      1000035         3                       # STR    -->  Z4SS   ST                 
#         PDG         Width
DECAY   2000004  2.22933817E+00   # CHR   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.97896373E-01    2      1000022         4                       # CHR    -->  Z1SS   CH                 
      4.21777426E-04    2      1000023         4                       # CHR    -->  Z2SS   CH                 
      3.92873888E-04    2      1000025         4                       # CHR    -->  Z3SS   CH                 
      1.28896185E-03    2      1000035         4                       # CHR    -->  Z4SS   CH                 
#         PDG         Width
DECAY   2000005  8.65413189E+00   # BT2   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.76298951E-02    2      1000022         5                       # BT2    -->  Z1SS   BT                 
      8.36304203E-02    2      1000023         5                       # BT2    -->  Z2SS   BT                 
      6.94645569E-02    2      1000025         5                       # BT2    -->  Z3SS   BT                 
      7.83215240E-02    2      1000035         5                       # BT2    -->  Z4SS   BT                 
      1.40161589E-01    2     -1000024         6                       # BT2    -->  W1SS-  TP                 
      3.01571250E-01    2     -1000037         6                       # BT2    -->  W2SS-  TP                 
      2.99220800E-01    2          -24   1000006                       # BT2    -->  W-     TP1                
#         PDG         Width
DECAY   2000006  1.50752592E+01   # TP2   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.65493870E-01    2      1000024         5                       # TP2    -->  W1SS+  BT                 
      1.25092804E-01    2      1000037         5                       # TP2    -->  W2SS+  BT                 
      1.80052266E-01    2           23   1000006                       # TP2    -->  Z0     TP1                
      8.66630077E-02    2           25   1000006                       # TP2    -->  HL0    TP1                
      2.60817762E-02    2      1000022         6                       # TP2    -->  Z1SS   TP                 
      1.14099048E-01    2      1000023         6                       # TP2    -->  Z2SS   TP                 
      5.57712391E-02    2      1000025         6                       # TP2    -->  Z3SS   TP                 
      1.46746024E-01    2      1000035         6                       # TP2    -->  Z4SS   TP                 
#         PDG         Width
DECAY   1000011  9.06351805E-01   # EL-   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.16202903E-01    2      1000022        11                       # EL-    -->  Z1SS   E-                 
      2.02000082E-01    2      1000023        11                       # EL-    -->  Z2SS   E-                 
      3.81796658E-01    2     -1000024        12                       # EL-    -->  W1SS-  NUE                
      1.36210957E-07    3      1000015        11       -15             # EL-    -->  TAU1-  E-     TAU+        
      2.07975759E-07    3     -1000015        11        15             # EL-    -->  TAU1+  E-     TAU-        
#         PDG         Width
DECAY   1000013  9.06351507E-01   # MUL-  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.16202962E-01    2      1000022        13                       # MUL-   -->  Z1SS   MU-                
      2.01999933E-01    2      1000023        13                       # MUL-   -->  Z2SS   MU-                
      3.81796777E-01    2     -1000024        14                       # MUL-   -->  W1SS-  NUM                
      1.36211014E-07    3      1000015        13       -15             # MUL-   -->  TAU1-  MU-    TAU+        
      2.07975816E-07    3     -1000015        13        15             # MUL-   -->  TAU1+  MU-    TAU-        
#         PDG         Width
DECAY   1000012  8.31395209E-01   # NUEL  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.70447093E-01    2      1000022        12                       # NUEL   -->  Z1SS   NUE                
      1.73210651E-01    2      1000023        12                       # NUEL   -->  Z2SS   NUE                
      3.56338859E-01    2      1000024        11                       # NUEL   -->  W1SS+  E-                 
      3.09690165E-07    3      1000015        12       -15             # NUEL   -->  TAU1-  NUE    TAU+        
      3.14721035E-07    3     -1000015        12        15             # NUEL   -->  TAU1+  NUE    TAU-        
      2.77051549E-06    3     -1000015        11        16             # NUEL   -->  TAU1+  E-     NUT         
#         PDG         Width
DECAY   1000014  8.31394851E-01   # NUML  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.70447302E-01    2      1000022        14                       # NUML   -->  Z1SS   NUM                
      1.73210725E-01    2      1000023        14                       # NUML   -->  Z2SS   NUM                
      3.56338590E-01    2      1000024        13                       # NUML   -->  W1SS+  MU-                
      3.09690279E-07    3      1000015        14       -15             # NUML   -->  TAU1-  NUM    TAU+        
      3.14721177E-07    3     -1000015        14        15             # NUML   -->  TAU1+  NUM    TAU-        
      2.78562788E-06    3     -1000015        13        16             # NUML   -->  TAU1+  MU-    NUT         
#         PDG         Width
DECAY   1000016  2.32832766E+00   # NUTL  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.33451402E-01    2      1000022        16                       # NUTL   -->  Z1SS   NUT                
      8.12754314E-03    2      1000023        16                       # NUTL   -->  Z2SS   NUT                
      1.60994045E-02    2      1000024        15                       # NUTL   -->  W1SS+  TAU-               
      8.42320919E-01    2           24   1000015                       # NUTL   -->  W+     TAU1-              
      7.21977585E-07    3     -1000015        16        15             # NUTL   -->  TAU1+  NUT    TAU-        
      4.91183947E-08    3      1000015        16       -15             # NUTL   -->  TAU1-  NUT    TAU+        
#         PDG         Width
DECAY   2000011  9.50808167E-01   # ER-   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.99999940E-01    2      1000022        11                       # ER-    -->  Z1SS   E-                 
      3.20349400E-08    3      1000015        11       -15             # ER-    -->  TAU1-  E-     TAU+        
      3.21147517E-08    3     -1000015        11        15             # ER-    -->  TAU1+  E-     TAU-        
      3.60796332E-12    3      1000015        12       -16             # ER-    -->  TAU1-  NUE    ANUT        
#         PDG         Width
DECAY   2000013  9.50808465E-01   # MUR-  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      9.99999821E-01    2      1000022        13                       # MUR-   -->  Z1SS   MU-                
      3.20350004E-08    3      1000015        13       -15             # MUR-   -->  TAU1-  MU-    TAU+        
      3.21148157E-08    3     -1000015        13        15             # MUR-   -->  TAU1+  MU-    TAU-        
      1.53391213E-07    3      1000015        14       -16             # MUR-   -->  TAU1-  NUM    ANUT        
#         PDG         Width
DECAY   2000015  2.57927346E+00   # TAU2- decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.69659615E-01    2      1000022        15                       # TAU2-  -->  Z1SS   TAU-               
      2.82944925E-02    2      1000023        15                       # TAU2-  -->  Z2SS   TAU-               
      5.20175248E-02    2     -1000024        16                       # TAU2-  -->  W1SS-  NUT                
      4.66725141E-01    2           23   1000015                       # TAU2-  -->  Z0     TAU1-              
      2.83303201E-01    2           25   1000015                       # TAU2-  -->  HL0    TAU1-              
#         PDG         Width
DECAY   1000022  1.91913766E-03   # Z1SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      5.00000000E-01    2      1000015       -15                       # Z1SS   -->  TAU1-  TAU+               
      5.00000000E-01    2     -1000015        15                       # Z1SS   -->  TAU1+  TAU-               
#         PDG         Width
DECAY   1000023  2.19691008E-01   # Z2SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.21119777E-07    2      1000022        22                       # Z2SS   -->  Z1SS   GM                 
      3.65646160E-03    2      1000022        23                       # Z2SS   -->  Z1SS   Z0                 
      2.72281767E-07    3      1000022         2        -2             # Z2SS   -->  Z1SS   UP     UB          
      3.25441874E-07    3      1000022         1        -1             # Z2SS   -->  Z1SS   DN     DB          
      3.25441874E-07    3      1000022         3        -3             # Z2SS   -->  Z1SS   ST     SB          
      2.72280829E-07    3      1000022         4        -4             # Z2SS   -->  Z1SS   CH     CB          
      3.24125858E-06    3      1000022         5        -5             # Z2SS   -->  Z1SS   BT     BB          
      7.04530248E-05    3      1000022        11       -11             # Z2SS   -->  Z1SS   E-     E+          
      7.04530248E-05    3      1000022        13       -13             # Z2SS   -->  Z1SS   MU-    MU+         
      9.87348540E-05    3      1000022        15       -15             # Z2SS   -->  Z1SS   TAU-   TAU+        
      8.21191206E-05    3      1000022        12       -12             # Z2SS   -->  Z1SS   NUE    ANUE        
      8.21191206E-05    3      1000022        14       -14             # Z2SS   -->  Z1SS   NUM    ANUM        
      2.19217196E-04    3      1000022        16       -16             # Z2SS   -->  Z1SS   NUT    ANUT        
      2.80832015E-02    2      1000022        25                       # Z2SS   -->  Z1SS   HL0                
      1.20694722E-05    2      2000011       -11                       # Z2SS   -->  ER-    E+                 
      1.20694722E-05    2     -2000011        11                       # Z2SS   -->  ER+    E-                 
      1.20690565E-05    2      2000013       -13                       # Z2SS   -->  MUR-   MU+                
      1.20690565E-05    2     -2000013        13                       # Z2SS   -->  MUR+   MU-                
      4.83792126E-01    2      1000015       -15                       # Z2SS   -->  TAU1-  TAU+               
      4.83792126E-01    2     -1000015        15                       # Z2SS   -->  TAU1+  TAU-               
#         PDG         Width
DECAY   1000025  5.25831127E+00   # Z3SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      2.33761497E-07    2      1000022        22                       # Z3SS   -->  Z1SS   GM                 
      5.39223599E-08    2      1000023        22                       # Z3SS   -->  Z2SS   GM                 
      2.14747682E-01    2      1000024       -24                       # Z3SS   -->  W1SS+  W-                 
      2.14747682E-01    2     -1000024        24                       # Z3SS   -->  W1SS-  W+                 
      6.46357238E-02    2      1000022        23                       # Z3SS   -->  Z1SS   Z0                 
      1.90247342E-01    2      1000023        23                       # Z3SS   -->  Z2SS   Z0                 
      7.08019687E-09    3      1000022         2        -2             # Z3SS   -->  Z1SS   UP     UB          
      2.18991958E-09    3      1000022         1        -1             # Z3SS   -->  Z1SS   DN     DB          
      2.18991958E-09    3      1000022         3        -3             # Z3SS   -->  Z1SS   ST     SB          
      7.08016312E-09    3      1000022         4        -4             # Z3SS   -->  Z1SS   CH     CB          
      8.36535837E-05    3      1000022         5        -5             # Z3SS   -->  Z1SS   BT     BB          
      1.66435311E-05    3      1000022        15       -15             # Z3SS   -->  Z1SS   TAU-   TAU+        
      4.88501373E-09    3      1000023         2        -2             # Z3SS   -->  Z2SS   UP     UB          
      7.82452148E-09    3      1000023         1        -1             # Z3SS   -->  Z2SS   DN     DB          
      7.82452148E-09    3      1000023         3        -3             # Z3SS   -->  Z2SS   ST     SB          
      4.88499330E-09    3      1000023         4        -4             # Z3SS   -->  Z2SS   CH     CB          
      3.24260218E-05    3      1000023         5        -5             # Z3SS   -->  Z2SS   BT     BB          
      4.55672171E-06    3      1000023        15       -15             # Z3SS   -->  Z2SS   TAU-   TAU+        
      1.70735139E-02    2      1000022        25                       # Z3SS   -->  Z1SS   HL0                
      1.28090009E-02    2      1000023        25                       # Z3SS   -->  Z2SS   HL0                
      7.51439002E-05    2      1000011       -11                       # Z3SS   -->  EL-    E+                 
      7.51439002E-05    2     -1000011        11                       # Z3SS   -->  EL+    E-                 
      7.51438856E-05    2      1000013       -13                       # Z3SS   -->  MUL-   MU+                
      7.51438856E-05    2     -1000013        13                       # Z3SS   -->  MUL+   MU-                
      1.84628792E-04    2      2000011       -11                       # Z3SS   -->  ER-    E+                 
      1.84628792E-04    2     -2000011        11                       # Z3SS   -->  ER+    E-                 
      1.84628749E-04    2      2000013       -13                       # Z3SS   -->  MUR-   MU+                
      1.84628749E-04    2     -2000013        13                       # Z3SS   -->  MUR+   MU-                
      9.55013186E-02    2      1000015       -15                       # Z3SS   -->  TAU1-  TAU+               
      9.55013186E-02    2     -1000015        15                       # Z3SS   -->  TAU1+  TAU-               
      4.54898700E-02    2      2000015       -15                       # Z3SS   -->  TAU2-  TAU+               
      4.54898700E-02    2     -2000015        15                       # Z3SS   -->  TAU2+  TAU-               
      3.87133972E-04    2      1000012       -12                       # Z3SS   -->  NUEL   ANUE               
      3.87133972E-04    2     -1000012        12                       # Z3SS   -->  ANUEL  NUE                
      3.87133972E-04    2      1000014       -14                       # Z3SS   -->  NUML   ANUM               
      3.87133972E-04    2     -1000014        14                       # Z3SS   -->  ANUML  NUM                
      5.15718595E-04    2      1000016       -16                       # Z3SS   -->  NUTL   ANUT               
      5.15718595E-04    2     -1000016        16                       # Z3SS   -->  ANUTL  NUT                
#         PDG         Width
DECAY   1000035  5.51008034E+00   # Z4SS  decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.08327992E-07    2      1000022        22                       # Z4SS   -->  Z1SS   GM                 
      9.35497457E-09    2      1000023        22                       # Z4SS   -->  Z2SS   GM                 
      6.37252681E-11    2      1000025        22                       # Z4SS   -->  Z3SS   GM                 
      2.13380888E-01    2      1000024       -24                       # Z4SS   -->  W1SS+  W-                 
      2.13380888E-01    2     -1000024        24                       # Z4SS   -->  W1SS-  W+                 
      1.79086700E-02    2      1000022        23                       # Z4SS   -->  Z1SS   Z0                 
      1.56781562E-02    2      1000023        23                       # Z4SS   -->  Z2SS   Z0                 
      7.44877084E-08    3      1000022         2        -2             # Z4SS   -->  Z1SS   UP     UB          
      4.97084898E-08    3      1000022         1        -1             # Z4SS   -->  Z1SS   DN     DB          
      4.97084898E-08    3      1000022         3        -3             # Z4SS   -->  Z1SS   ST     SB          
      7.44873816E-08    3      1000022         4        -4             # Z4SS   -->  Z1SS   CH     CB          
      7.97867324E-05    3      1000022         5        -5             # Z4SS   -->  Z1SS   BT     BB          
      1.59519841E-05    3      1000022        15       -15             # Z4SS   -->  Z1SS   TAU-   TAU+        
      2.27225399E-07    3      1000023         2        -2             # Z4SS   -->  Z2SS   UP     UB          
      2.82327335E-07    3      1000023         1        -1             # Z4SS   -->  Z2SS   DN     DB          
      2.82327335E-07    3      1000023         3        -3             # Z4SS   -->  Z2SS   ST     SB          
      2.27224490E-07    3      1000023         4        -4             # Z4SS   -->  Z2SS   CH     CB          
      3.23678833E-05    3      1000023         5        -5             # Z4SS   -->  Z2SS   BT     BB          
      4.40232043E-06    3      1000023        15       -15             # Z4SS   -->  Z2SS   TAU-   TAU+        
      1.87236360E-09    3      1000025         2        -2             # Z4SS   -->  Z3SS   UP     UB          
      2.41535103E-09    3      1000025         1        -1             # Z4SS   -->  Z3SS   DN     DB          
      2.41535103E-09    3      1000025         3        -3             # Z4SS   -->  Z3SS   ST     SB          
      1.87236360E-09    3      1000025         4        -4             # Z4SS   -->  Z3SS   CH     CB          
      5.47803192E-10    3      1000025        11       -11             # Z4SS   -->  Z3SS   E-     E+          
      5.47803192E-10    3      1000025        13       -13             # Z4SS   -->  Z3SS   MU-    MU+         
      3.67608721E-10    3      1000025        15       -15             # Z4SS   -->  Z3SS   TAU-   TAU+        
      1.08995613E-09    3      1000025        12       -12             # Z4SS   -->  Z3SS   NUE    ANUE        
      1.08995613E-09    3      1000025        14       -14             # Z4SS   -->  Z3SS   NUM    ANUM        
      1.08995613E-09    3      1000025        16       -16             # Z4SS   -->  Z3SS   NUT    ANUT        
      6.12107255E-02    2      1000022        25                       # Z4SS   -->  Z1SS   HL0                
      1.87687099E-01    2      1000023        25                       # Z4SS   -->  Z2SS   HL0                
      1.78868964E-03    2      1000011       -11                       # Z4SS   -->  EL-    E+                 
      1.78868964E-03    2     -1000011        11                       # Z4SS   -->  EL+    E-                 
      1.78868952E-03    2      1000013       -13                       # Z4SS   -->  MUL-   MU+                
      1.78868952E-03    2     -1000013        13                       # Z4SS   -->  MUL+   MU-                
      6.25477580E-04    2      2000011       -11                       # Z4SS   -->  ER-    E+                 
      6.25477580E-04    2     -2000011        11                       # Z4SS   -->  ER+    E-                 
      6.25477405E-04    2      2000013       -13                       # Z4SS   -->  MUR-   MU+                
      6.25477405E-04    2     -2000013        13                       # Z4SS   -->  MUR+   MU-                
      7.62132257E-02    2      1000015       -15                       # Z4SS   -->  TAU1-  TAU+               
      7.62132257E-02    2     -1000015        15                       # Z4SS   -->  TAU1+  TAU-               
      5.09139150E-02    2      2000015       -15                       # Z4SS   -->  TAU2-  TAU+               
      5.09139150E-02    2     -2000015        15                       # Z4SS   -->  TAU2+  TAU-               
      4.02822904E-03    2      1000012       -12                       # Z4SS   -->  NUEL   ANUE               
      4.02822904E-03    2     -1000012        12                       # Z4SS   -->  ANUEL  NUE                
      4.02822904E-03    2      1000014       -14                       # Z4SS   -->  NUML   ANUM               
      4.02822904E-03    2     -1000014        14                       # Z4SS   -->  ANUML  NUM                
      5.29775722E-03    2      1000016       -16                       # Z4SS   -->  NUTL   ANUT               
      5.29775722E-03    2     -1000016        16                       # Z4SS   -->  ANUTL  NUT                
#         PDG         Width
DECAY   1000024  2.13363171E-01   # W1SS+ decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.00000000E+00    2      1000022         24        # W1SS+  -->  W   Z1SS          
#         PDG         Width
DECAY   1000037  5.45480204E+00   # W2SS+ decays
#          BR          NDA       ID1       ID2       ID3       ID4
      8.23240285E-08    3      1000022         2        -1             # W2SS+  -->  Z1SS   UP     DB          
      8.23239432E-08    3      1000022         4        -3             # W2SS+  -->  Z1SS   CH     SB          
      1.39341673E-05    3      1000022       -15        16             # W2SS+  -->  Z1SS   TAU+   NUT         
      7.27590919E-02    2      1000022        24                       # W2SS+  -->  Z1SS   W+                 
      2.30239564E-07    3      1000023         2        -1             # W2SS+  -->  Z2SS   UP     DB          
      2.30239493E-07    3      1000023         4        -3             # W2SS+  -->  Z2SS   CH     SB          
      5.56894202E-06    3      1000023       -15        16             # W2SS+  -->  Z2SS   TAU+   NUT         
      2.21155554E-01    2      1000023        24                       # W2SS+  -->  Z2SS   W+                 
      7.12335968E-09    3      1000025         2        -1             # W2SS+  -->  Z3SS   UP     DB          
      7.12335968E-09    3      1000025         4        -3             # W2SS+  -->  Z3SS   CH     SB          
      2.37430253E-09    3      1000025       -11        12             # W2SS+  -->  Z3SS   E+     NUE         
      2.37430253E-09    3      1000025       -13        14             # W2SS+  -->  Z3SS   MU+    NUM         
      2.37431297E-09    3      1000025       -15        16             # W2SS+  -->  Z3SS   TAU+   NUT         
      3.18924151E-03    2      1000012       -11                       # W2SS+  -->  NUEL   E+                 
      3.18924151E-03    2      1000014       -13                       # W2SS+  -->  NUML   MU+                
      1.00677975E-01    2      1000016       -15                       # W2SS+  -->  NUTL   TAU+               
      8.71899258E-03    2     -1000011        12                       # W2SS+  -->  EL+    NUE                
      8.71899258E-03    2     -1000013        14                       # W2SS+  -->  MUL+   NUM                
      1.28181159E-01    2     -1000015        16                       # W2SS+  -->  TAU1+  NUT                
      3.76218781E-02    2     -2000015        16                       # W2SS+  -->  TAU2+  NUT                
      2.11088791E-01    2      1000024        23                       # W2SS+  -->  W1SS+  Z0                 
      2.00140633E-07    3      1000024         1        -1             # W2SS+  -->  W1SS+  DN     DB          
      2.00139780E-07    3      1000024         3        -3             # W2SS+  -->  W1SS+  ST     SB          
      5.49638912E-07    3      1000024         2        -2             # W2SS+  -->  W1SS+  UP     UB          
      5.49638912E-07    3      1000024         4        -4             # W2SS+  -->  W1SS+  CH     CB          
      2.04677463E-01    2      1000024        25                       # W2SS+  -->  W1SS+  HL0                
#         PDG         Width
DECAY        25  4.28013969E-03   # HL0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      4.77537165E-09    2           11       -11                       # HL0    -->  E-     E+                 
      2.01623829E-04    2           13       -13                       # HL0    -->  MU-    MU+                
      5.76683506E-02    2           15       -15                       # HL0    -->  TAU-   TAU+               
      5.99213945E-06    2            1        -1                       # HL0    -->  DN     DB                 
      2.42112461E-03    2            3        -3                       # HL0    -->  ST     SB                 
      7.91881859E-01    2            5        -5                       # HL0    -->  BT     BB                 
      1.86300781E-06    2            2        -2                       # HL0    -->  UP     UB                 
      3.49842906E-02    2            4        -4                       # HL0    -->  CH     CB                 
      1.84560777E-03    2           22        22                       # HL0    -->  GM     GM                 
      3.96150686E-02    2           21        21                       # HL0    -->  GL     GL                 
      3.66468704E-03    3           24        11       -12             # HL0    -->  W+     E-     ANUE        
      3.66468704E-03    3           24        13       -14             # HL0    -->  W+     MU-    ANUM        
      3.66468704E-03    3           24        15       -16             # HL0    -->  W+     TAU-   ANUT        
      1.09940609E-02    3           24        -2         1             # HL0    -->  W+     UB     DN          
      1.09940609E-02    3           24        -4         3             # HL0    -->  W+     CB     ST          
      3.66468704E-03    3          -24       -11        12             # HL0    -->  W-     E+     NUE         
      3.66468704E-03    3          -24       -13        14             # HL0    -->  W-     MU+    NUM         
      3.66468704E-03    3          -24       -15        16             # HL0    -->  W-     TAU+   NUT         
      1.09940609E-02    3          -24         2        -1             # HL0    -->  W-     UP     DB          
      1.09940609E-02    3          -24         4        -3             # HL0    -->  W-     CH     SB          
      3.69976769E-04    3           23        12       -12             # HL0    -->  Z0     NUE    ANUE        
      3.69976769E-04    3           23        14       -14             # HL0    -->  Z0     NUM    ANUM        
      3.69976769E-04    3           23        16       -16             # HL0    -->  Z0     NUT    ANUT        
      1.86205580E-04    3           23        11       -11             # HL0    -->  Z0     E-     E+          
      1.86205580E-04    3           23        13       -13             # HL0    -->  Z0     MU-    MU+         
      1.86205580E-04    3           23        15       -15             # HL0    -->  Z0     TAU-   TAU+        
      6.37924590E-04    3           23         2        -2             # HL0    -->  Z0     UP     UB          
      6.37924590E-04    3           23         4        -4             # HL0    -->  Z0     CH     CB          
      8.21803987E-04    3           23         1        -1             # HL0    -->  Z0     DN     DB          
      8.21803987E-04    3           23         3        -3             # HL0    -->  Z0     ST     SB          
      8.21803987E-04    3           23         5        -5             # HL0    -->  Z0     BT     BB          
#         PDG         Width
DECAY        35  1.57665234E+01   # HH0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.08428519E-08    2           11       -11                       # HH0    -->  E-     E+                 
      4.57804708E-04    2           13       -13                       # HH0    -->  MU-    MU+                
      1.31116018E-01    2           15       -15                       # HH0    -->  TAU-   TAU+               
      1.33369213E-05    2            1        -1                       # HH0    -->  DN     DB                 
      5.38879121E-03    2            3        -3                       # HH0    -->  ST     SB                 
      8.47770452E-01    2            5        -5                       # HH0    -->  BT     BB                 
      1.71551694E-12    2            2        -2                       # HH0    -->  UP     UB                 
      2.55749324E-08    2            4        -4                       # HH0    -->  CH     CB                 
      8.74553400E-04    2            6        -6                       # HH0    -->  TP     TB                 
      4.17106527E-08    2           22        22                       # HH0    -->  GM     GM                 
      8.81763026E-06    2           21        21                       # HH0    -->  GL     GL                 
      6.07276320E-07    2           24       -24                       # HH0    -->  W+     W-                 
      3.03128559E-07    2           23        23                       # HH0    -->  Z0     Z0                 
      2.11170744E-04    2      1000022   1000022                       # HH0    -->  Z1SS   Z1SS               
      9.15425335E-05    2      1000022   1000023                       # HH0    -->  Z1SS   Z2SS               
      5.49548022E-05    2           25        25                       # HH0    -->  HL0    HL0                
      1.40115535E-02    2      1000015  -1000015                       # HH0    -->  TAU1-  TAU1+              
#         PDG         Width
DECAY        36  1.54689407E+01   # HA0   decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.09805027E-08    2           11       -11                       # HA0    -->  E-     E+                 
      4.63616598E-04    2           13       -13                       # HA0    -->  MU-    MU+                
      1.32784754E-01    2           15       -15                       # HA0    -->  TAU-   TAU+               
      1.35070286E-05    2            1        -1                       # HA0    -->  DN     DB                 
      5.45752328E-03    2            3        -3                       # HA0    -->  ST     SB                 
      8.58688116E-01    2            5        -5                       # HA0    -->  BT     BB                 
      1.68820968E-12    2            2        -2                       # HA0    -->  UP     UB                 
      2.51881627E-08    2            4        -4                       # HA0    -->  CH     CB                 
      1.19010801E-03    2            6        -6                       # HA0    -->  TP     TB                 
      4.99788477E-08    2           22        22                       # HA0    -->  GM     GM                 
      2.17849338E-05    2           21        21                       # HA0    -->  GL     GL                 
      3.88916204E-04    2      1000022   1000022                       # HA0    -->  Z1SS   Z1SS               
      9.91033274E-04    2      1000022   1000023                       # HA0    -->  Z1SS   Z2SS               
      5.67365362E-07    2           25        23                       # HA0    -->  HL0    Z0                 
#         PDG         Width
DECAY        37  1.30214748E+01   # H+    decays
#          BR          NDA       ID1       ID2       ID3       ID4
      1.32508511E-08    2           12       -11                       # H+     -->  NUE    E+                 
      5.59474807E-04    2           14       -13                       # H+     -->  NUM    MU+                
      1.60239652E-01    2           16       -15                       # H+     -->  NUT    TAU+               
      1.49208545E-05    2            2        -1                       # H+     -->  UP     DB                 
      6.02875138E-03    2            4        -3                       # H+     -->  CH     SB                 
      8.13027143E-01    2            6        -5                       # H+     -->  TP     BB                 
      1.50601589E-03    2      1000024   1000022                       # H+     -->  W1SS+  Z1SS               
      7.12423457E-07    2           25        24                       # H+     -->  HL0    W+                 
      1.86233055E-02    2     -1000015   1000016                       # H+     -->  TAU1+  NUTL               
