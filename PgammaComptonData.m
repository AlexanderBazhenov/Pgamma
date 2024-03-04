%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PgammaCompton=[
-3.1
-0.2
-4.0
-2.1
-3.7
-1.7
-5.7
-2.8
-8.0
-2.1
-3.6
-7.2
]

OnesPgammaCompton=ones(length(PgammaCompton),1)

PgammaComptonStat=[
2.7
2.1
2.1
2.5
1.9
3.7
2.8
1.9
4.0
3.9
2.6
2.5
]


PgammaComptonMeanW = sum(PgammaCompton./PgammaComptonStat)/sum(OnesPgammaCompton./PgammaComptonStat)

PgammaComptonStd =midrad(PgammaCompton, PgammaComptonStat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

