% Pgamma1992Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phys Lett 1992
PgammaPh=[
-4.4
-3.4
-6.9
-1.2
-1.0
-10.8
-10.2
-6.3
-10.4
0.6
-1.8
-6.6
-4.9
-6.0
-4.0]

OnesPgamma=ones(length(PgammaPh),1)

PgammaPhStat=[
2.7
1.9
2.4
2.4
2.7
3.5
2.8
2
4.1
3.4
2
2.1
2.1
2.4
2.7]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PgammaBG=[
4.2
-3.2
12.1
12.4
9.4
1
-0.6
3.9
10.3
-4.8
4.6
-5.7
13
8.4
10.6
]

PgammaBGStat=[
6.7
4.8
9
7.2
5.1
12.4
6.1
4.3
10
10.6
4.2
4.6
3
4.6
5.5
]

PgammaBGMeanW = sum(PgammaBG./PgammaBGStat)/sum(OnesPgamma./PgammaBGStat)


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PgammaBG=[
4.2
-3.2
12.1
12.4
9.4
1
-0.6
3.9
10.3
-4.8
4.6
-5.7
13
8.4
10.6
]

PgammaBGStat=[
6.7
4.8
9
7.2
5.1
12.4
6.1
4.3
10
10.6
4.2
4.6
3
4.6
5.5
]

PgammaBGMeanW = sum(PgammaBG./PgammaBGStat)/sum(OnesPgamma./PgammaBGStat)


