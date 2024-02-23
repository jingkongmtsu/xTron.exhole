!  6-31G**  EMSL  Basis Set Exchange Library   5/6/16 7:32 PM
! Elements                             References
! --------                             ----------
! H - He: W.J. Hehre, R. Ditchfield and J.A. Pople, J. Chem. Phys. 56,
! Li - Ne: 2257 (1972).  Note: Li and B come from J.D. Dill and J.A.
! Pople, J. Chem. Phys. 62, 2921 (1975).
! Na - Ar: M.M. Francl, W.J. Petro, W.J. Hehre, J.S. Binkley, M.S. Gordon,
! D.J. DeFrees and J.A. Pople, J. Chem. Phys. 77, 3654 (1982)
! K  - Zn: V. Rassolov, J.A. Pople, M. Ratner and T.L. Windus, J. Chem. Phys.
! 109, 1223 (1998)
! Note: He and Ne are unpublished basis sets taken from the Gaussian
! program
! 


! Elements                             References
! --------                             ----------
! H,Li - Ne: P.C. Hariharan and J.A. Pople, Theoret. Chimica Acta 28, 213 (1973).
! Na - Ar  : M.M. Francl, W.J. Petro, W.J. Hehre, J.S. Binkley, M.S. Gordon, D.J.
!            DeFrees and J.A. Pople, J. Chem. Phys. 77, 3654 (1982).
! K  - Zn:   V.A. Rassolov, J.A. Pople, M.A. Ratner, and T.L. Windus
!            J. Chem. Phys. 109, 1223 (1998)
!            Note: He and Ne are unpublished basis sets taken from Gaussian.
!   

!
! we modify the Si basis set to match with GAMESS
! 

%basis
H     0 
S   3   1.00
     18.7311370              0.03349460       
      2.8253937              0.23472695       
      0.6401217              0.81375733       
S   1   1.00
      0.1612778              1.0000000        
P   1   1.00
      1.1000000              1.0000000        
****
He     0 
S   3   1.00
     38.4216340              0.0237660        
      5.7780300              0.1546790        
      1.2417740              0.4696300        
S   1   1.00
      0.2979640              1.0000000        
P   1   1.00
      1.1000000              1.0000000        
****
Li     0 
S   6   1.00
    642.4189200              0.0021426        
     96.7985150              0.0162089        
     22.0911210              0.0773156        
      6.2010703              0.2457860        
      1.9351177              0.4701890        
      0.6367358              0.3454708        
SP   3   1.00
      2.3249184             -0.0350917              0.0089415        
      0.6324306             -0.1912328              0.1410095        
      0.0790534              1.0839878              0.9453637        
SP   1   1.00
      0.0359620              1.0000000              1.0000000        
D   1   1.00    C
      0.2000000              1.0000000        
****
Be     0 
S   6   1.00
   1264.5857000              0.0019448        
    189.9368100              0.0148351        
     43.1590890              0.0720906        
     12.0986630              0.2371542        
      3.8063232              0.4691987        
      1.2728903              0.3565202        
SP   3   1.00
      3.1964631             -0.1126487              0.0559802        
      0.7478133             -0.2295064              0.2615506        
      0.2199663              1.1869167              0.7939723        
SP   1   1.00
      0.0823099              1.0000000              1.0000000        
D   1   1.00    C
      0.4000000              1.0000000        
****
B     0 
S   6   1.00
   2068.8823000              0.0018663        
    310.6495700              0.0142515        
     70.6830330              0.0695516        
     19.8610800              0.2325729        
      6.2993048              0.4670787        
      2.1270270              0.3634314        
SP   3   1.00
      4.7279710             -0.1303938              0.0745976        
      1.1903377             -0.1307889              0.3078467        
      0.3594117              1.1309444              0.7434568        
SP   1   1.00
      0.1267512              1.0000000              1.0000000        
D   1   1.00     C
      0.6000000              1.0000000        
****
C     0 
S   6   1.00
   3047.5249000              0.0018347        
    457.3695100              0.0140373        
    103.9486900              0.0688426        
     29.2101550              0.2321844        
      9.2866630              0.4679413        
      3.1639270              0.3623120        
SP   3   1.00
      7.8682724             -0.1193324              0.0689991        
      1.8812885             -0.1608542              0.3164240        
      0.5442493              1.1434564              0.7443083        
SP   1   1.00
      0.1687144              1.0000000              1.0000000        
D   1   1.00       C
      0.8000000              1.0000000        
****
N     0 
S   6   1.00
   4173.5110000              0.0018348        
    627.4579000              0.0139950        
    142.9021000              0.0685870        
     40.2343300              0.2322410        
     12.8202100              0.4690700        
      4.3904370              0.3604550        
SP   3   1.00
     11.6263580             -0.1149610              0.0675800        
      2.7162800             -0.1691180              0.3239070        
      0.7722180              1.1458520              0.7408950        
SP   1   1.00
      0.2120313              1.0000000              1.0000000        
D   1   1.00    C
      0.8000000              1.0000000        
****
O     0 
S   6   1.00
   5484.6717000              0.0018311        
    825.2349500              0.0139501        
    188.0469600              0.0684451        
     52.9645000              0.2327143        
     16.8975700              0.4701930        
      5.7996353              0.3585209        
SP   3   1.00
     15.5396160             -0.1107775              0.0708743        
      3.5999336             -0.1480263              0.3397528        
      1.0137618              1.1307670              0.7271586        
SP   1   1.00
      0.2700058              1.0000000              1.0000000        
D   1   1.00  C
      0.8000000              1.0000000        
****
F     0 
S   6   1.00
   7001.7130900              0.0018196169     
   1051.3660900              0.0139160796     
    239.2856900              0.0684053245     
     67.3974453              0.233185760      
     21.5199573              0.471267439      
      7.40310130             0.356618546      
SP   3   1.00
     20.8479528             -0.108506975            0.0716287243     
      4.80830834            -0.146451658            0.3459121030     
      1.34406986             1.128688580            0.7224699570     
SP   1   1.00
      0.358151393            1.0000000              1.0000000        
D   1   1.00   C
      0.8000000              1.0000000        
****
Ne     0 
S   6   1.00
   8425.8515300              0.0018843481     
   1268.5194000              0.0143368994     
    289.6214140              0.0701096233     
     81.8590040              0.2373732660     
     26.2515079              0.4730071260     
      9.09472051             0.3484012410     
SP   3   1.00
     26.5321310             -0.107118287            0.0719095885     
      6.10175501            -0.146163821            0.3495133720     
      1.69627153             1.127773500            0.7199405120     
SP   1   1.00
      0.44581870             1.0000000              1.0000000        
D   1   1.00    C
      0.8000000              1.0000000        
****
Na     0 
S   6   1.00
   9993.2000000              0.0019377        
   1499.8900000              0.0148070        
    341.9510000              0.0727060        
     94.6797000              0.2526290        
     29.7345000              0.4932420        
     10.0063000              0.3131690        
SP   6   1.00
    150.9630000             -0.0035421              0.0050017        
     35.5878000             -0.0439590              0.0355110        
     11.1683000             -0.1097521              0.1428250        
      3.9020100              0.1873980              0.3386200        
      1.3817700              0.6466990              0.4515790        
      0.4663820              0.3060580              0.2732710        
SP   3   1.00
      0.4979660             -0.2485030             -0.0230230        
      0.0843530             -0.1317040              0.9503590        
      0.0666350              1.2335200              0.0598580        
SP   1   1.00
      0.0259544              1.0000000              1.0000000        
D   1   1.00   C
      0.1750000              1.0000000        
****
Mg     0 
S   6   1.00
  11722.8000000              0.0019778        
   1759.9300000              0.0151140        
    400.8460000              0.0739110        
    112.8070000              0.2491910        
     35.9997000              0.4879280        
     12.1828000              0.3196620        
SP   6   1.00
    189.1800000             -0.0032372              0.0049281        
     45.2119000             -0.0410080              0.0349890        
     14.3563000             -0.1126000              0.1407250        
      5.1388600              0.1486330              0.3336420        
      1.9065200              0.6164970              0.4449400        
      0.7058870              0.3648290              0.2692540        
SP   3   1.00
      0.9293400             -0.2122900             -0.0224190        
      0.2690350             -0.1079850              0.1922700        
      0.1173790              1.1758400              0.8461810        
SP   1   1.00
      0.0421061              1.0000000              1.0000000        
D   1   1.00    C
      0.1750000              1.0000000        
****
Al     0 
S   6   1.00
  13983.1000000              0.00194267       
   2098.7500000              0.0148599        
    477.7050000              0.0728494        
    134.3600000              0.2468300        
     42.8709000              0.4872580        
     14.5189000              0.3234960        
SP   6   1.00
    239.6680000             -0.00292619             0.00460285       
     57.4419000             -0.0374080              0.0331990        
     18.2859000             -0.1144870              0.1362820        
      6.5991400              0.1156350              0.3304760        
      2.4904900              0.6125950              0.4491460        
      0.9445400              0.3937990              0.2657040        
SP   3   1.00
      1.2779000             -0.2276060             -0.0175130        
      0.3975900              0.00144583             0.2445330        
      0.1600950              1.0927900              0.8049340        
SP   1   1.00
      0.0556577              1.0000000              1.0000000        
D   1   1.00    C
      0.3250000              1.0000000        
****
Si     0 
S   6   1.00
       16192.1000000    0.001949239405
        2436.0900000    0.014855895466
         556.0010000    0.072568877851
         156.8130000    0.245654925022
          50.1692000    0.486059851646
          17.0300000    0.325719900585
SP   6   1.00
        293.3500000   -0.002829912441    0.004433341604
         70.1173000   -0.036073731114    0.032440211739
         22.4301000   -0.116808100749    0.133719048387
          8.1942500    0.093576880711    0.326780118247
          3.1476800    0.601705518979    0.451139163247
          1.2151500    0.422072364043    0.264105095568
SP   3   1.00
        1.6537000   -0.240599186106   -0.015177406538
        0.5407600    0.073795050368    0.275139118516
        0.2044060    1.040936478740    0.783008337279
SP   1   1.00
        0.0723837    1.000000000000    1.000000000000
D   1   1.00      C
      0.39500000             1.0000000        
****
P     0 
S   6   1.00
  19413.3000000              0.0018516        
   2909.4200000              0.0142062        
    661.3640000              0.0699995        
    185.7590000              0.2400790        
     59.1943000              0.4847620        
     20.0310000              0.3352000        
SP   6   1.00
    339.4780000             -0.00278217             0.00456462       
     81.0101000             -0.0360499              0.03369360       
     25.8780000             -0.1166310              0.13975500       
      9.4522100              0.0968328              0.33936200       
      3.6656600              0.6144180              0.45092100       
      1.4674600              0.4037980              0.23858600       
SP   3   1.00
      2.1562300             -0.2529230             -0.01776530       
      0.7489970              0.0328517              0.27405800       
      0.2831450              1.0812500              0.78542100       
SP   1   1.00
      0.0998317              1.0000000              1.00000000       
D   1   1.00    C
      0.5500000              1.0000000        
****
S     0 
S   6   1.00
  21917.1000000              0.0018690        
   3301.4900000              0.0142300        
    754.1460000              0.0696960        
    212.7110000              0.2384870        
     67.9896000              0.4833070        
     23.0515000              0.3380740        
SP   6   1.00
    423.7350000             -0.0023767              0.0040610        
    100.7100000             -0.0316930              0.0306810        
     32.1599000             -0.1133170              0.1304520        
     11.8079000              0.0560900              0.3272050        
      4.6311000              0.5922550              0.4528510        
      1.8702500              0.4550060              0.2560420        
SP   3   1.00
      2.6158400             -0.2503740             -0.0145110        
      0.9221670              0.0669570              0.3102630        
      0.3412870              1.0545100              0.7544830        
SP   1   1.00
      0.1171670              1.0000000              1.0000000        
D   1   1.00      C
      0.6500000              1.0000000        
****
Cl     0 
S   6   1.00
  25180.1000000              0.0018330        
   3780.3500000              0.0140340        
    860.4740000              0.0690970        
    242.1450000              0.2374520        
     77.3349000              0.4830340        
     26.2470000              0.3398560        
SP   6   1.00
    491.7650000             -0.0022974              0.0039894        
    116.9840000             -0.0307140              0.0303180        
     37.4153000             -0.1125280              0.1298800        
     13.7834000              0.0450160              0.3279510        
      5.4521500              0.5893530              0.4535270        
      2.2258800              0.4652060              0.2521540        
SP   3   1.00
      3.1864900             -0.2518300             -0.0142990        
      1.1442700              0.0615890              0.3235720        
      0.4203770              1.0601800              0.7435070        
SP   1   1.00
      0.1426570              1.0000000              1.0000000        
D   1   1.00     C
      0.7500000              1.0000000        
****
Ar     0 
S   6   1.00
  28348.3000000              0.00182526       
   4257.6200000              0.01396860       
    969.8570000              0.06870730       
    273.2630000              0.23620400       
     87.3695000              0.48221400       
     29.6867000              0.34204300       
SP   6   1.00
    575.8910000             -0.00215972             0.00380665       
    136.8160000             -0.02907750             0.02923050       
     43.8098000             -0.11082700             0.12646700       
     16.2094000              0.02769990             0.32351000       
      6.4608400              0.57761300             0.45489600       
      2.6511400              0.48868800             0.25663000       
SP   3   1.00
      3.8602800             -0.2555920             -0.01591970       
      1.4137300              0.0378066              0.32464600       
      0.5166460              1.0805600              0.74399000       
SP   1   1.00
      0.1738880              1.0000000              1.0000000        
D   1   1.00   C
      0.8500000              1.0000000        
%end



