#
#                              ======================
#                              | THE SUSYHIT OUTPUT |
#                              ======================
#
#
#              ------------------------------------------------------
#              |     This is the output of the SUSY-HIT package     |
#              |  created by A.Djouadi, M.Muehlleitner and M.Spira. |
#              |  In case of problems with SUSY-HIT email to        |
#              |           margarete.muehlleitner@cern.ch           |
#              |           michael.spira@psi.ch                     |
#              |           abdelhak.djouadi@cern.ch                 |
#              ------------------------------------------------------
#
#              ------------------------------------------------------
#              |  SUSY Les Houches Accord - MSSM Spectrum + Decays  |
#              |              based on the decay programs           |
#              |                                                    |
#              |                     SDECAY 1.3b                    |
#              |                                                    |
#              |  Authors: M.Muhlleitner, A.Djouadi and Y.Mambrini  |
#              |  Ref.:    Comput.Phys.Commun.168(2005)46           |
#              |           [hep-ph/0311167]                         |
#              |                                                    |
#              |                     HDECAY 3.4                     |
#              |                                                    |
#              |  By: A.Djouadi,J.Kalinowski,M.Muhlleitner,M.Spira  |
#              |  Ref.:    Comput.Phys.Commun.108(1998)56           |
#              |           [hep-ph/9704448]                         |
#              |                                                    |
#              |                                                    |
#              |  If not stated otherwise all DRbar couplings and   |
#              |  soft SUSY breaking masses are given at the scale  |
#              |  Q=  0.40022015E+04
#              |                                                    |
#              ------------------------------------------------------
#
#
BLOCK DCINFO  # Decay Program information
     1   SDECAY/HDECAY # decay calculator
     2   1.3b  /3.4    # version number
#
BLOCK SPINFO  # Spectrum calculator information
     1    SOFTSUSY    # spectrum calculator                
     2    3.3.5       # version number                     
#
BLOCK MODSEL  # Model selection
     1     0   #  nonUniversal                                     
#
BLOCK SMINPUTS  # Standard Model inputs
         1     1.27934000E+02   # alpha_em^-1(M_Z)^MSbar
         2     1.16637000E-05   # G_F [GeV^-2]
         3     1.17200000E-01   # alpha_S(M_Z)^MSbar
         4     9.11876000E+01   # M_Z pole mass
         5     4.25000000E+00   # mb(mb)^MSbar
         6     1.73000000E+02   # mt pole mass
         7     1.77700000E+00   # mtau pole mass
#
BLOCK MINPAR  # Input parameters - minimal models
         3     1.00000000E+01   # tanb                
#
BLOCK EXTPAR  # Input parameters - non-minimal models
         0    -1.00000000E+00   # Set                 
         1     1.00000000E+00   # M_1(MX)             
         2     4.00000000E+03   # M_2(MX)             
         3     1.20000000E+03   # M_3(MX)             
        11     0.00000000E+00   # At(MX)              
        12     0.00000000E+00   # Ab(MX)              
        13     0.00000000E+00   # Atau(MX)            
        23     4.00000000E+03   # mu(MX)              
        26     4.00000000E+03   # mA(pole)            
        31     4.00000000E+03   # meL(MX)             
        32     4.00000000E+03   # mmuL(MX)            
        33     4.00000000E+03   # mtauL(MX)           
        34     4.00000000E+03   # meR(MX)             
        35     4.00000000E+03   # mmuR(MX)            
        36     4.00000000E+03   # mtauR(MX)           
        41     4.00000000E+03   # mqL1(MX)            
        42     4.00000000E+03   # mqL2(MX)            
        43     4.00000000E+03   # mqL3(MX)            
        44     4.00000000E+03   # muR(MX)             
        45     4.00000000E+03   # mcR(MX)             
        46     4.00000000E+03   # mtR(MX)             
        47     4.00000000E+03   # mdR(MX)             
        48     4.00000000E+03   # msR(MX)             
        49     4.00000000E+03   # mbR(MX)             
#
BLOCK MASS  # Mass Spectrum
# PDG code           mass       particle
        24     8.04095278E+01   # W+
        25     1.26000000E+02   # h
        35     4.00006519E+03   # H
        36     3.99997936E+03   # A
        37     4.00119035E+03   # H+
         5     4.87877998E+00   # b-quark pole mass calculated from mb(mb)_Msbar
   1000001     4.07266705E+03   # ~d_L
   2000001     4.04379662E+03   # ~d_R
   1000002     4.07205962E+03   # ~u_L
   2000002     4.04300039E+03   # ~u_R
   1000003     4.07266705E+03   # ~s_L
   2000003     4.04379662E+03   # ~s_R
   1000004     4.07205962E+03   # ~c_L
   2000004     4.04300039E+03   # ~c_R
   1000005     4.03811482E+03   # ~b_1
   2000005     4.07471332E+03   # ~b_2
   1000006     4.03829212E+03   # ~t_1
   2000006     4.07575571E+03   # ~t_2
   1000011     4.02986687E+03   # ~e_L
   2000011     4.00514666E+03   # ~e_R
   1000012     4.02879605E+03   # ~nu_eL
   1000013     4.02986687E+03   # ~mu_L
   2000013     4.00514666E+03   # ~mu_R
   1000014     4.02879605E+03   # ~nu_muL
   1000015     4.00169846E+03   # ~tau_1
   2000015     4.03231221E+03   # ~tau_2
   1000016     4.02843753E+03   # ~nu_tauL
   1000021     1.20000000E+03   # ~g
   1000022     1.00000000E+02   # ~chi_10
   1000023     3.95018270E+03   # ~chi_20
   1000025    -4.00007717E+03   # ~chi_30
   1000035     4.06855690E+03   # ~chi_40
   1000024     3.95143887E+03   # ~chi_1+
   1000037     4.06761725E+03   # ~chi_2+
#
BLOCK NMIX  # Neutralino Mixing Matrix
  1  1     9.99937761E-01   # N_11
  1  2    -4.28326615E-05   # N_12
  1  3     1.10948334E-02   # N_13
  1  4    -1.17323638E-03   # N_14
  2  1    -6.52575126E-03   # N_21
  2  2    -6.66316786E-01   # N_22
  2  3     5.30121526E-01   # N_23
  2  4    -5.24357248E-01   # N_24
  3  1    -7.01484196E-03   # N_31
  3  2     6.10062540E-03   # N_32
  3  3     7.07022535E-01   # N_33
  3  4     7.07129910E-01   # N_34
  4  1    -5.71666058E-03   # N_41
  4  2     7.45643830E-01   # N_42
  4  3     4.67939320E-01   # N_43
  4  4    -4.74357873E-01   # N_44
#
BLOCK UMIX  # Chargino Mixing Matrix U
  1  1    -6.42783576E-01   # U_11
  1  2     7.66047827E-01   # U_12
  2  1    -7.66047827E-01   # U_21
  2  2    -6.42783576E-01   # U_22
#
BLOCK VMIX  # Chargino Mixing Matrix V
  1  1    -6.51870003E-01   # V_11
  1  2     7.58330732E-01   # V_12
  2  1    -7.58330732E-01   # V_21
  2  2    -6.51870003E-01   # V_22
#
BLOCK STOPMIX  # Stop Mixing Matrix
  1  1     2.46543338E-01   # cos(theta_t)
  1  2     9.69131767E-01   # sin(theta_t)
  2  1    -9.69131767E-01   # -sin(theta_t)
  2  2     2.46543338E-01   # cos(theta_t)
#
BLOCK SBOTMIX  # Sbottom Mixing Matrix
  1  1     2.88543093E-01   # cos(theta_b)
  1  2     9.57466910E-01   # sin(theta_b)
  2  1    -9.57466910E-01   # -sin(theta_b)
  2  2     2.88543093E-01   # cos(theta_b)
#
BLOCK STAUMIX  # Stau Mixing Matrix
  1  1     2.92601533E-01   # cos(theta_tau)
  1  2     9.56234460E-01   # sin(theta_tau)
  2  1    -9.56234460E-01   # -sin(theta_tau)
  2  2     2.92601533E-01   # cos(theta_tau)
#
BLOCK ALPHA  # Higgs mixing
          -1.05041888E-01   # Mixing angle in the neutral Higgs boson sector
#
BLOCK HMIX Q=  4.00220151E+03  # DRbar Higgs Parameters
         1     4.00000000E+03   # mu(Q)MSSM           
         2     9.49828503E+00   # tan                 
         3     2.43271581E+02   # higgs               
         4     1.57485305E+07   # mA^2(Q)MSSM         
#
BLOCK GAUGE Q=  4.00220151E+03  # The gauge couplings
     1     3.65176479E-01   # gprime(Q) DRbar
     2     6.33291151E-01   # g(Q) DRbar
     3     1.00857801E+00   # g3(Q) DRbar
#
BLOCK AU Q=  4.00220151E+03  # The trilinear couplings
  1  1     7.64045290E-06   # A_u(Q) DRbar
  2  2     7.64050268E-06   # A_c(Q) DRbar
  3  3     1.15854252E-05   # A_t(Q) DRbar
#
BLOCK AD Q=  4.00220151E+03  # The trilinear couplings
  1  1     4.16784708E-06   # A_d(Q) DRbar
  2  2     4.16788485E-06   # A_s(Q) DRbar
  3  3     5.21980454E-06   # A_b(Q) DRbar
#
BLOCK AE Q=  4.00220151E+03  # The trilinear couplings
  1  1     0.00000000E+00   # A_e(Q) DRbar
  2  2     1.08826674E-07   # A_mu(Q) DRbar
  3  3     1.08725393E-07   # A_tau(Q) DRbar
#
BLOCK Yu Q=  4.00220151E+03  # The Yukawa couplings
  1  1     0.00000000E+00   # y_u(Q) DRbar
  2  2     0.00000000E+00   # y_c(Q) DRbar
  3  3     8.19594976E-01   # y_t(Q) DRbar
#
BLOCK Yd Q=  4.00220151E+03  # The Yukawa couplings
  1  1     0.00000000E+00   # y_d(Q) DRbar
  2  2     0.00000000E+00   # y_s(Q) DRbar
  3  3     1.25875011E-01   # y_b(Q) DRbar
#
BLOCK Ye Q=  4.00220151E+03  # The Yukawa couplings
  1  1     0.00000000E+00   # y_e(Q) DRbar
  2  2     0.00000000E+00   # y_mu(Q) DRbar
  3  3     9.94030174E-02   # y_tau(Q) DRbar
#
BLOCK MSOFT Q=  4.00220151E+03  # The soft SUSY breaking masses at the scale Q
         1     1.00000001E+00   # M_1(Q)              
         2     4.00000000E+03   # M_2(Q)              
         3     1.30000000E+03   # M_3(Q)              
        21    -4.43890407E+05   # mH1^2(Q)            
        22    -1.55418122E+07   # mH2^2(Q)            
        31     4.00000000E+03   # meL(Q)              
        32     4.00000000E+03   # mmuL(Q)             
        33     4.00000000E+03   # mtauL(Q)            
        34     4.00000000E+03   # meR(Q)              
        35     4.00000000E+03   # mmuR(Q)             
        36     4.00000000E+03   # mtauR(Q)            
        41     4.00000000E+03   # mqL1(Q)             
        42     4.00000000E+03   # mqL2(Q)             
        43     4.00000000E+03   # mqL3(Q)             
        44     4.00000000E+03   # muR(Q)              
        45     4.00000000E+03   # mcR(Q)              
        46     3.99999999E+03   # mtR(Q)              
        47     4.00000000E+03   # mdR(Q)              
        48     4.00000000E+03   # msR(Q)              
        49     4.00000000E+03   # mbR(Q)              
#
#
#
#                             =================
#                             |The decay table|
#                             =================
#
# - The QCD corrections to the decays gluino -> squark  + quark
#                                     squark -> gaugino + quark_prime
#                                     squark -> squark_prime + Higgs
#                                     squark -> gluino  + quark
#   are included.
#
# - The multi-body decays for the inos, stops and sbottoms are included.
#
# - The loop induced decays for the gluino, neutralinos and stops
#   are included.
#
# - The SUSY decays of the top quark are included.
#
#
#         PDG            Width
DECAY         6     1.40591219E+00   # top decays
#          BR         NDA      ID1       ID2
     1.00000000E+00    2           5        24   # BR(t ->  b    W+)
#
#         PDG            Width
DECAY   1000021     7.34077982E-05   # gluino decays
#           BR         NDA      ID1       ID2       ID3
     1.00000000E+00    3     1000022         6        -6   # BR(~g -> ~chi_10 t  tb)
#
#         PDG            Width
DECAY   1000006     2.30120471E+02   # stop1 decays
#          BR         NDA      ID1       ID2
     3.68301966E-02    2     1000022         6   # BR(~t_1 -> ~chi_10 t )
     3.55551081E-04    2     1000024         5   # BR(~t_1 -> ~chi_1+ b )
     9.62814252E-01    2     1000021         6   # BR(~t_1 -> ~g      t )
#
#         PDG            Width
DECAY   2000006     2.19261941E+02   # stop2 decays
#          BR         NDA      ID1       ID2
     6.01625156E-03    2     1000022         6   # BR(~t_2 -> ~chi_10 t )
     7.00411684E-05    2     1000024         5   # BR(~t_2 -> ~chi_1+ b )
     2.50012924E-06    2     1000037         5   # BR(~t_2 -> ~chi_2+ b )
     9.93911207E-01    2     1000021         6   # BR(~t_2 -> ~g      t )
#
#         PDG            Width
DECAY   1000005     2.21310072E+02   # sbottom1 decays
#          BR         NDA      ID1       ID2
     9.84458066E-03    2     1000022         5   # BR(~b_1 -> ~chi_10 b )
     1.64840167E-05    2     1000023         5   # BR(~b_1 -> ~chi_20 b )
     1.18966592E-06    2     1000025         5   # BR(~b_1 -> ~chi_30 b )
     9.90137746E-01    2     1000021         5   # BR(~b_1 -> ~g      b )
#
#         PDG            Width
DECAY   2000005     2.22289380E+02   # sbottom2 decays
#          BR         NDA      ID1       ID2
     3.20835981E-03    2     1000022         5   # BR(~b_2 -> ~chi_10 b )
     1.01936635E-04    2     1000023         5   # BR(~b_2 -> ~chi_20 b )
     4.21099824E-06    2     1000025         5   # BR(~b_2 -> ~chi_30 b )
     3.65862479E-07    2     1000035         5   # BR(~b_2 -> ~chi_40 b )
     9.96685127E-01    2     1000021         5   # BR(~b_2 -> ~g      b )
#
#         PDG            Width
DECAY   1000002     2.22310342E+02   # sup_L decays
#          BR         NDA      ID1       ID2
     2.61615666E-03    2     1000022         2   # BR(~u_L -> ~chi_10 u)
     2.46897109E-04    2     1000023         2   # BR(~u_L -> ~chi_20 u)
     6.07974246E-09    2     1000025         2   # BR(~u_L -> ~chi_30 u)
     4.92885108E-06    2     1000035         2   # BR(~u_L -> ~chi_40 u)
     4.63953863E-04    2     1000024         1   # BR(~u_L -> ~chi_1+ d)
     1.30635793E-05    2     1000037         1   # BR(~u_L -> ~chi_2+ d)
     9.96654994E-01    2     1000021         2   # BR(~u_L -> ~g      u)
#
#         PDG            Width
DECAY   2000002     2.28508603E+02   # sup_R decays
#          BR         NDA      ID1       ID2
     4.04483239E-02    2     1000022         2   # BR(~u_R -> ~chi_10 u)
     9.30169548E-09    2     1000023         2   # BR(~u_R -> ~chi_20 u)
     3.95367369E-09    2     1000025         2   # BR(~u_R -> ~chi_30 u)
     9.59551663E-01    2     1000021         2   # BR(~u_R -> ~g      u)
#
#         PDG            Width
DECAY   1000001     2.22357516E+02   # sdown_L decays
#          BR         NDA      ID1       ID2
     2.61832773E-03    2     1000022         1   # BR(~d_L -> ~chi_10 d)
     2.46717289E-04    2     1000023         1   # BR(~d_L -> ~chi_20 d)
     1.50998259E-08    2     1000025         1   # BR(~d_L -> ~chi_30 d)
     5.84388007E-06    2     1000035         1   # BR(~d_L -> ~chi_40 d)
     4.54217309E-04    2    -1000024         2   # BR(~d_L -> ~chi_1- u)
     1.52203003E-05    2    -1000037         2   # BR(~d_L -> ~chi_2- u)
     9.96659658E-01    2     1000021         1   # BR(~d_L -> ~g      d)
#
#         PDG            Width
DECAY   2000001     2.21640065E+02   # sdown_R decays
#          BR         NDA      ID1       ID2
     1.04275192E-02    2     1000022         1   # BR(~d_R -> ~chi_10 d)
     2.42554807E-09    2     1000023         1   # BR(~d_R -> ~chi_20 d)
     1.04230256E-09    2     1000025         1   # BR(~d_R -> ~chi_30 d)
     9.89572477E-01    2     1000021         1   # BR(~d_R -> ~g      d)
#
#         PDG            Width
DECAY   1000004     2.22310342E+02   # scharm_L decays
#          BR         NDA      ID1       ID2
     2.61615666E-03    2     1000022         4   # BR(~c_L -> ~chi_10 c)
     2.46897109E-04    2     1000023         4   # BR(~c_L -> ~chi_20 c)
     6.07974246E-09    2     1000025         4   # BR(~c_L -> ~chi_30 c)
     4.92885108E-06    2     1000035         4   # BR(~c_L -> ~chi_40 c)
     4.63953863E-04    2     1000024         3   # BR(~c_L -> ~chi_1+ s)
     1.30635793E-05    2     1000037         3   # BR(~c_L -> ~chi_2+ s)
     9.96654994E-01    2     1000021         4   # BR(~c_L -> ~g      c)
#
#         PDG            Width
DECAY   2000004     2.28508603E+02   # scharm_R decays
#          BR         NDA      ID1       ID2
     4.04483239E-02    2     1000022         4   # BR(~c_R -> ~chi_10 c)
     9.30169548E-09    2     1000023         4   # BR(~c_R -> ~chi_20 c)
     3.95367369E-09    2     1000025         4   # BR(~c_R -> ~chi_30 c)
     9.59551663E-01    2     1000021         4   # BR(~c_R -> ~g      c)
#
#         PDG            Width
DECAY   1000003     2.22357516E+02   # sstrange_L decays
#          BR         NDA      ID1       ID2
     2.61832773E-03    2     1000022         3   # BR(~s_L -> ~chi_10 s)
     2.46717289E-04    2     1000023         3   # BR(~s_L -> ~chi_20 s)
     1.50998259E-08    2     1000025         3   # BR(~s_L -> ~chi_30 s)
     5.84388007E-06    2     1000035         3   # BR(~s_L -> ~chi_40 s)
     4.54217309E-04    2    -1000024         4   # BR(~s_L -> ~chi_1- c)
     1.52203003E-05    2    -1000037         4   # BR(~s_L -> ~chi_2- c)
     9.96659658E-01    2     1000021         3   # BR(~s_L -> ~g      s)
#
#         PDG            Width
DECAY   2000003     2.21640065E+02   # sstrange_R decays
#          BR         NDA      ID1       ID2
     1.04275192E-02    2     1000022         3   # BR(~s_R -> ~chi_10 s)
     2.42554807E-09    2     1000023         3   # BR(~s_R -> ~chi_20 s)
     1.04230256E-09    2     1000025         3   # BR(~s_R -> ~chi_30 s)
     9.89572477E-01    2     1000021         3   # BR(~s_R -> ~g      s)
#
#         PDG            Width
DECAY   1000011     5.37494419E+00   # selectron_L decays
#          BR         NDA      ID1       ID2
     9.94268869E-01    2     1000022        11   # BR(~e_L -> ~chi_10 e-)
     2.05906095E-03    2     1000023        11   # BR(~e_L -> ~chi_20 e-)
     2.74227578E-09    2     1000025        11   # BR(~e_L -> ~chi_30 e-)
     3.67206742E-03    2    -1000024        12   # BR(~e_L -> ~chi_1- nu_e)
#
#         PDG            Width
DECAY   2000011     2.12485870E+01   # selectron_R decays
#          BR         NDA      ID1       ID2
     9.99999968E-01    2     1000022        11   # BR(~e_R -> ~chi_10 e-)
     3.16456338E-08    2     1000023        11   # BR(~e_R -> ~chi_20 e-)
     3.14986318E-10    2     1000025        11   # BR(~e_R -> ~chi_30 e-)
#
#         PDG            Width
DECAY   1000013     5.37494419E+00   # smuon_L decays
#          BR         NDA      ID1       ID2
     9.94268869E-01    2     1000022        13   # BR(~mu_L -> ~chi_10 mu-)
     2.05906095E-03    2     1000023        13   # BR(~mu_L -> ~chi_20 mu-)
     2.74227578E-09    2     1000025        13   # BR(~mu_L -> ~chi_30 mu-)
     3.67206742E-03    2    -1000024        14   # BR(~mu_L -> ~chi_1- nu_mu)
#
#         PDG            Width
DECAY   2000013     2.12485870E+01   # smuon_R decays
#          BR         NDA      ID1       ID2
     9.99999968E-01    2     1000022        13   # BR(~mu_R -> ~chi_10 mu-)
     3.16456338E-08    2     1000023        13   # BR(~mu_R -> ~chi_20 mu-)
     3.14986318E-10    2     1000025        13   # BR(~mu_R -> ~chi_30 mu-)
#
#         PDG            Width
DECAY   1000015     1.98825913E+01   # stau_1 decays
#          BR         NDA      ID1       ID2
     9.99858010E-01    2     1000022        15   # BR(~tau_1 -> ~chi_10  tau-)
     5.00936869E-05    2     1000023        15   # BR(~tau_1 -> ~chi_20  tau-)
     9.18960809E-05    2    -1000024        16   # BR(~tau_1 -> ~chi_1-  nu_tau)
#
#         PDG            Width
DECAY   2000015     6.73543812E+00   # stau_2 decays
#          BR         NDA      ID1       ID2
     9.95989145E-01    2     1000022        15   # BR(~tau_2 -> ~chi_10  tau-)
     1.46663699E-03    2     1000023        15   # BR(~tau_2 -> ~chi_20  tau-)
     1.49550155E-05    2     1000025        15   # BR(~tau_2 -> ~chi_30  tau-)
     2.52926337E-03    2    -1000024        16   # BR(~tau_2 -> ~chi_1-  nu_tau)
#
#         PDG            Width
DECAY   1000012     5.37460362E+00   # snu_eL decays
#          BR         NDA      ID1       ID2
     9.94363079E-01    2     1000022        12   # BR(~nu_eL -> ~chi_10 nu_e)
     1.96048829E-03    2     1000023        12   # BR(~nu_eL -> ~chi_20 nu_e)
     6.21205857E-08    2     1000025        12   # BR(~nu_eL -> ~chi_30 nu_e)
     3.67637073E-03    2     1000024        11   # BR(~nu_eL -> ~chi_1+ e-)
#
#         PDG            Width
DECAY   1000014     5.37460362E+00   # snu_muL decays
#          BR         NDA      ID1       ID2
     9.94363079E-01    2     1000022        14   # BR(~nu_muL -> ~chi_10 nu_mu)
     1.96048829E-03    2     1000023        14   # BR(~nu_muL -> ~chi_20 nu_mu)
     6.21205857E-08    2     1000025        14   # BR(~nu_muL -> ~chi_30 nu_mu)
     3.67637073E-03    2     1000024        13   # BR(~nu_muL -> ~chi_1+ mu-)
#
#         PDG            Width
DECAY   1000016     5.37435037E+00   # snu_tauL decays
#          BR         NDA      ID1       ID2
     9.94321445E-01    2     1000022        16   # BR(~nu_tauL -> ~chi_10 nu_tau)
     1.94308292E-03    2     1000023        16   # BR(~nu_tauL -> ~chi_20 nu_tau)
     6.05928875E-08    2     1000025        16   # BR(~nu_tauL -> ~chi_30 nu_tau)
     3.73541164E-03    2     1000024        15   # BR(~nu_tauL -> ~chi_1+ tau-)
#
#         PDG            Width
DECAY   1000024     1.40399066E+00   # chargino1+ decays
#          BR         NDA      ID1       ID2
     1.00000000E+00    2     1000022        24   # BR(~chi_1+ -> ~chi_10  W+)
#
#         PDG            Width
DECAY   1000037     1.72297055E+00   # chargino2+ decays
#          BR         NDA      ID1       ID2
     3.17665596E-03    2     1000006        -5   # BR(~chi_2+ -> ~t_1     bb)
     1.95456203E-03    2     1000012       -11   # BR(~chi_2+ -> ~nu_eL   e+  )
     1.95456203E-03    2     1000014       -13   # BR(~chi_2+ -> ~nu_muL  mu+ )
     1.99971032E-03    2     1000016       -15   # BR(~chi_2+ -> ~nu_tau1 tau+)
     1.88652919E-03    2    -1000011        12   # BR(~chi_2+ -> ~e_L+    nu_e)
     1.88652919E-03    2    -1000013        14   # BR(~chi_2+ -> ~mu_L+   nu_mu)
     1.58662615E-04    2    -1000015        16   # BR(~chi_2+ -> ~tau_1+  nu_tau)
     1.63381724E-03    2    -2000015        16   # BR(~chi_2+ -> ~tau_2+  nu_tau)
     1.29109412E-01    2     1000024        23   # BR(~chi_2+ -> ~chi_1+  Z )
     6.11157108E-01    2     1000022        24   # BR(~chi_2+ -> ~chi_10  W+)
     2.44448570E-01    2     1000023        24   # BR(~chi_2+ -> ~chi_20  W+)
     6.33881676E-04    2     1000022        37   # BR(~chi_2+ -> ~chi_10  H+)
#
#         PDG            Width
DECAY   1000022     0.00000000E+00   # neutralino1 decays
#
#         PDG            Width
DECAY   1000023     1.42565399E+00   # neutralino2 decays
#          BR         NDA      ID1       ID2
     3.83302707E-01    2     1000022        23   # BR(~chi_20 -> ~chi_10   Z )
     6.16697293E-01    2     1000022        25   # BR(~chi_20 -> ~chi_10   h )
#
#         PDG            Width
DECAY   1000025     2.58680102E+00   # neutralino3 decays
#          BR         NDA      ID1       ID2
     5.95041601E-01    2     1000022        23   # BR(~chi_30 -> ~chi_10   Z )
     4.04958399E-01    2     1000022        25   # BR(~chi_30 -> ~chi_10   h )
#
#         PDG            Width
DECAY   1000035     1.81744500E+00   # neutralino4 decays
#          BR         NDA      ID1       ID2
     2.54487515E-01    2     1000022        23   # BR(~chi_40 -> ~chi_10   Z )
     1.39061923E-06    2     1000023        23   # BR(~chi_40 -> ~chi_20   Z )
     1.68166203E-01    2     1000024       -24   # BR(~chi_40 -> ~chi_1+   W-)
     1.68166203E-01    2    -1000024        24   # BR(~chi_40 -> ~chi_1-   W+)
     3.97498873E-01    2     1000022        25   # BR(~chi_40 -> ~chi_10   h )
     2.70614745E-04    2     1000022        35   # BR(~chi_40 -> ~chi_10   H )
     4.09507484E-04    2     1000022        36   # BR(~chi_40 -> ~chi_10   A )
     4.93195305E-07    2     2000002        -2   # BR(~chi_40 -> ~u_R      ub)
     4.93195305E-07    2    -2000002         2   # BR(~chi_40 -> ~u_R*     u )
     1.19139045E-07    2     2000001        -1   # BR(~chi_40 -> ~d_R      db)
     1.19139045E-07    2    -2000001         1   # BR(~chi_40 -> ~d_R*     d )
     4.93195305E-07    2     2000004        -4   # BR(~chi_40 -> ~c_R      cb)
     4.93195305E-07    2    -2000004         4   # BR(~chi_40 -> ~c_R*     c )
     1.19139045E-07    2     2000003        -3   # BR(~chi_40 -> ~s_R      sb)
     1.19139045E-07    2    -2000003         3   # BR(~chi_40 -> ~s_R*     s )
     3.49410130E-05    2     1000005        -5   # BR(~chi_40 -> ~b_1      bb)
     3.49410130E-05    2    -1000005         5   # BR(~chi_40 -> ~b_1*     b )
     8.81674064E-04    2     1000011       -11   # BR(~chi_40 -> ~e_L-     e+)
     8.81674064E-04    2    -1000011        11   # BR(~chi_40 -> ~e_L+     e-)
     1.85653220E-07    2     2000011       -11   # BR(~chi_40 -> ~e_R-     e+)
     1.85653220E-07    2    -2000011        11   # BR(~chi_40 -> ~e_R+     e-)
     8.81674064E-04    2     1000013       -13   # BR(~chi_40 -> ~mu_L-    mu+)
     8.81674064E-04    2    -1000013        13   # BR(~chi_40 -> ~mu_L+    mu-)
     1.85653220E-07    2     2000013       -13   # BR(~chi_40 -> ~mu_R-    mu+)
     1.85653220E-07    2    -2000013        13   # BR(~chi_40 -> ~mu_R+    mu-)
     6.79402841E-05    2     1000015       -15   # BR(~chi_40 -> ~tau_1-   tau+)
     6.79402841E-05    2    -1000015        15   # BR(~chi_40 -> ~tau_1+   tau-)
     7.72958069E-04    2     2000015       -15   # BR(~chi_40 -> ~tau_2-   tau+)
     7.72958069E-04    2    -2000015        15   # BR(~chi_40 -> ~tau_2+   tau-)
     9.47515514E-04    2     1000012       -12   # BR(~chi_40 -> ~nu_eL    nu_eb)
     9.47515514E-04    2    -1000012        12   # BR(~chi_40 -> ~nu_eL*   nu_e )
     9.47515514E-04    2     1000014       -14   # BR(~chi_40 -> ~nu_muL   nu_mub)
     9.47515514E-04    2    -1000014        14   # BR(~chi_40 -> ~nu_muL*  nu_mu )
     9.64594452E-04    2     1000016       -16   # BR(~chi_40 -> ~nu_tau1  nu_taub)
     9.64594452E-04    2    -1000016        16   # BR(~chi_40 -> ~nu_tau1* nu_tau )
#
#         PDG            Width
DECAY        25     4.33749302E-03   # h decays
#          BR         NDA      ID1       ID2
     6.08319660E-01    2           5        -5   # BR(h -> b       bb     )
     6.02982239E-02    2         -15        15   # BR(h -> tau+    tau-   )
     2.13429425E-04    2         -13        13   # BR(h -> mu+     mu-    )
     4.51612143E-04    2           3        -3   # BR(h -> s       sb     )
     1.96130448E-02    2           4        -4   # BR(h -> c       cb     )
     6.48379316E-02    2          21        21   # BR(h -> g       g      )
     2.24434343E-03    2          22        22   # BR(h -> gam     gam    )
     1.57432183E-03    2          22        23   # BR(h -> Z       gam    )
     2.14757797E-01    2          24       -24   # BR(h -> W+      W-     )
     2.73068946E-02    2          23        23   # BR(h -> Z       Z      )
     3.82742023E-04    2     1000022   1000022   # BR(h -> ~chi_10 ~chi_10)
#
#         PDG            Width
DECAY        35     6.44434564E+00   # H decays
#          BR         NDA      ID1       ID2
     6.12453396E-01    2           5        -5   # BR(H -> b       bb     )
     1.16051587E-01    2         -15        15   # BR(H -> tau+    tau-   )
     4.10284123E-04    2         -13        13   # BR(H -> mu+     mu-    )
     4.84465162E-04    2           3        -3   # BR(H -> s       sb     )
     2.67200211E-06    2           4        -4   # BR(H -> c       cb     )
     2.69877731E-01    2           6        -6   # BR(H -> t       tb     )
     6.66746113E-05    2          21        21   # BR(H -> g       g      )
     3.03908953E-07    2          22        22   # BR(H -> gam     gam    )
     1.38627414E-07    2          23        22   # BR(H -> Z       gam    )
     6.94478434E-05    2          24       -24   # BR(H -> W+      W-     )
     3.46998478E-05    2          23        23   # BR(H -> Z       Z      )
     2.62194971E-04    2          25        25   # BR(H -> h       h      )
    -3.20940376E-24    2          36        36   # BR(H -> A       A      )
     1.18032802E-18    2          23        36   # BR(H -> Z       A      )
     1.79390257E-04    2     1000022   1000022   # BR(H -> ~chi_10 ~chi_10)
     1.07014106E-04    2     1000022   1000023   # BR(H -> ~chi_10 ~chi_20)
#
#         PDG            Width
DECAY        36     6.44834573E+00   # A decays
#          BR         NDA      ID1       ID2
     6.12283416E-01    2           5        -5   # BR(A -> b       bb     )
     1.15980772E-01    2         -15        15   # BR(A -> tau+    tau-   )
     4.10033443E-04    2         -13        13   # BR(A -> mu+     mu-    )
     4.84520096E-04    2           3        -3   # BR(A -> s       sb     )
     2.68745099E-06    2           4        -4   # BR(A -> c       cb     )
     2.70218852E-01    2           6        -6   # BR(A -> t       tb     )
     1.95231214E-04    2          21        21   # BR(A -> g       g      )
     5.51321920E-07    2          22        22   # BR(A -> gam     gam    )
     2.75817460E-07    2          23        22   # BR(A -> Z       gam    )
     6.92534697E-05    2          23        25   # BR(A -> Z       h      )
     1.87451267E-04    2     1000022   1000022   # BR(A -> ~chi_10 ~chi_10)
     1.66955639E-04    2     1000022   1000023   # BR(A -> ~chi_10 ~chi_20)
#
#         PDG            Width
DECAY        37     6.19219136E+00   # H+ decays
#          BR         NDA      ID1       ID2
     9.60624390E-04    2           4        -5   # BR(H+ -> c       bb     )
     1.20815152E-01    2         -15        16   # BR(H+ -> tau+    nu_tau )
     4.27124701E-04    2         -13        14   # BR(H+ -> mu+     nu_mu  )
     6.14796452E-06    2           2        -5   # BR(H+ -> u       bb     )
     2.41612986E-05    2           2        -3   # BR(H+ -> u       sb     )
     4.99624621E-04    2           4        -3   # BR(H+ -> c       sb     )
     8.76898976E-01    2           6        -5   # BR(H+ -> t       bb     )
     7.22091277E-05    2          24        25   # BR(H+ -> W+      h      )
     5.53424969E-13    2          24        36   # BR(H+ -> W+      A      )
     2.95980702E-04    2     1000024   1000022   # BR(H+ -> ~chi_1+ ~chi_10)
