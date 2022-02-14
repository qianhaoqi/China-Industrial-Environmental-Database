$title Cross Entropy Pollution Matrix Estimation

* Define sets
SETS
    r               'region'    / r11, r12, r13, r14, r15, r21, r22, r23, r31, r32, r33, r34, r35, r36, r37, r41, r42, r43, r44, r45, r46, r50, r51, r52, r53, r54, r61, r62, r63, r64, r65, total, total2 /
    rr(r)                       / r11, r12, r13, r14, r15, r21, r22, r23, r31, r32, r33, r34, r35, r36, r37, r41, r42, r43, r44, r45, r46, r50, r51, r52, r53, r54, r61, r62, r63, r64, r65 /
    j               'sector'    / s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38, s39, total /
    jj(j)                       / s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36, s37, s38, s39 /
    t               'year'      / 2015 /
    nonzero(rr,jj)  'matrix elements that can be nonzero';
;

* Define parameters and scalars
PARAMETERS
    mat0(r,t,j)         'raw pollution matrix'
    mat(r,j)            'intermediate pollution matrix'
    mat_balance(r,t,j)  'balanced matrix'
    coef0(r,j)          'prior coefficient matrix'
    totalr(rr)          'total value for each region from raw data'
    totalj(jj)          'total value for each sector from raw data'
    targetr(rr)         'total value for each region from yearbook'
    targetj0(jj)        'weighted adjusted total value for each sector'
    targetj(jj)         'scale adjusted total value for each sector'
    adjustr(rr)         'scale rates for each region'
    adjustj(jj)         'scale rates for each sector'
    epsilon             'tolerance to allow zero cells in the matrix'
    scale
;

* Read raw data to be banlanced, change pollutant name for own needs
*EXECUTE 'GDXXRW .\Data\so2removal.xlsx par=so2nremoval rng=Sheet1!A1:AP595 cdim=1 rdim=2';
$GDXIN so2removal.gdx
$LOAD mat0=so2removal
$GDXIN so2removal.gdx

* Initialize Parameters for the specific year that need balancing
mat(r,j) = mat0(r,"2015",j);
mat(rr,jj)$(mat(rr,jj) le 0) = 0;
mat("total",jj)$(mat("total",jj) le 0) = SUM(rr,mat(rr,jj));
totalj(jj) = SUM(rr,mat(rr,jj));
totalr(rr) = SUM(jj,mat(rr,jj));
targetr(rr) = mat(rr,"total");
adjustr(rr)$totalr(rr) = targetr(rr) / totalr(rr);
targetj0(jj) = SUM(rr,mat(rr,jj)*adjustr(rr));
adjustj(jj)$totalj(jj) = targetj0(jj) / totalj(jj);
coef0(rr,jj)$(mat(rr,jj) AND targetj0(jj)) = mat(rr,jj)*adjustr(rr) / targetj0(jj);
nonzero(rr,jj)$(coef0(rr,jj)) = YES;
epsilon = 1e-6;

* Define endogeneous variables
VARIABLES
    A(rr,jj)        'post SAM coefficient matrix'
    PMAT(rr,jj)     'post matrix of SAM transactions'
    DENTROPY        'entropy difference (objective)'
;

* Define core equations for the optimization procedure
EQUATIONS
    ENTROPY         'entropy difference definition'
    GENA(rr,jj)     'generate post coefficients'
    ROWSUM          'total target'
    COLSUM(jj)      'column target'
;

GENA(rr,jj)..   A(rr,jj) =E= PMAT(rr,jj) / targetj0(jj);
ENTROPY..       DENTROPY =E= SUM((rr,jj)$nonzero(rr,jj), A(rr,jj)*(log(A(rr,jj))- log(coef0(rr,jj))));
ROWSUM..        SUM((rr,jj), PMAT(rr,jj)) =E= SUM(rr,targetr(rr));
COLSUM(jj)..    SUM(rr, A(rr,jj)) =E= 1.0;

* Define the optimization model
MODEL MATENTROP / all /;

* Initialize endogeneous variables
A.L(rr,jj)      = coef0(rr,jj);
PMAT.L(rr,jj)   = mat(rr,jj);
DENTROPY.L      = 0;
A.LO(rr,jj)$nonzero(rr,jj)       = epsilon;
A.UP(rr,jj)$nonzero(rr,jj)       = 1;
A.FX(rr,jj)$(NOT nonzero(rr,jj)) = 0;
PMAT.FX(rr,jj)$(NOT nonzero(rr,jj)) = 0;

* Initialize option settings
OPTION iterLim = 10000, limRow = 0, limCol = 0, solPrint = on, nlp = MINOS;

* Start solving the optimization problem
SOLVE MATENTROP USING NLP MINIMIZING DENTROPY;

* Store results
mat_balance(rr,"2015",jj) = PMAT.L(rr,jj);
mat_balance("total","2015",jj) = SUM(rr,PMAT.L(rr,jj));
mat_balance(rr,"2015","total") = SUM(jj,PMAT.L(rr,jj));

* Display key results if needed
*DISPLAY adjustr, adjustj, totalj, totalr, targetr, coef0, mat, PMAT.L, mat_balance;

* Store results to dgx file and then convert to excel file
* Change the pollutant name and year for own needs
EXECUTE_UNLOAD '.\Result\so2removal_2015.gdx', mat_balance;
EXECUTE 'GDXXRW i=.\Result\so2removal_2015.gdx o=.\Result\so2removal_balance.xlsx par=mat_balance cdim=1 rdim=2 rng=2015!A1';
