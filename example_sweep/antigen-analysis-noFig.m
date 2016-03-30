(* ::Package:: *)

Analysis*of*Antigen*output

Basic*parameters
runDirectory = $CommandLine[[4]];
SetDirectory[runDirectory];

Off[IterationLimit::itlim];
Off[RecursionLimit::reclim];

demeCount = 3; 
demes = {"north", "tropics", "south"}; 
spatialColors = {RGBColor[0.765, 0.728, 0.274], 
   RGBColor[0.324, 0.609, 0.708], 
       RGBColor[0.857, 0.131, 0.132]}; 
spatialColors = (Lighter[#1, 0.1] & ) /@ spatialColors;
spatialColorRules = MapIndexed[#2[[1]] - 1 -> #1 & , spatialColors]; 
k = 11; 
mod = 8; 
clusterColors = 
  Table[Lighter[ColorData["Rainbow"][FractionalPart[i - 0.001]], 0.1], 
       {i, 1/k, mod, mod/k}]; 
narrowWidth = 255; 
extraNarrowWidth = 215; 
padding = {{35, 5}, {20, 15}}; 
fullPadding = {{37, 2}, {29, 5}}; 
partialPadding = {{15, 2}, {13, 2}}; 
imageSize = 500; 

Loess*fit

LoessFit[(x_)?VectorQ, data_, \[Alpha]_: 0.75, \[Lambda]_: 1] := 
 Table[LoessFit[x[[i]], data, \[Alpha], \[Lambda]], 
	   {i, Length[x]}];
LoessFit[(x_)?NumberQ, data_, \[Alpha]_: 0.75, \[Lambda]_: 1] := 
  WLSFit[data, LoessWts[x, data, \[Alpha]], \[Lambda], 
       x]; 
WLSFit[data_, wts_, ldegree_: 1, x_] := 
   Fit[Transpose[(wts*#1 & ) /@ 
     Join[{Table[1, {Length[data]}]}, Transpose[data]]], 
       Join[{u}, Table[v^i, {i, ldegree}]], {u, v}] /. {u -> 1, v -> x}
LoessWts[x_, data_, \[Alpha]_] := 
 Tricube[(x - First[Transpose[data]])/LoessDistance[x, data, \[Alpha]]]
Tricube = 
  Compile[{{x, _Real, 
     1}}, (If[Abs[#1] < 1, (1 - Abs[#1]^3)^3, 0] & ) /@ x]; 
LoessDistance[x_, data_, \[Alpha]_] := 
 Module[{A = Max[1, \[Alpha]], X = First[Transpose[data]], q}, 
     q = Min[Length[X], Ceiling[\[Alpha]*Length[X]]]; 
  A*Sort[Abs[X - x]][[q]]]

Options

SetOptions[Plot, 
  LabelStyle -> 
   Directive[Black, FontSize -> 9, FontFamily -> "Helvetica"], 
     Axes -> False, Frame -> {True, True, False, False}, 
  FrameStyle -> Black, 
     ImageSize -> imageSize, PlotRange -> {0, All}]; 
SetOptions[ListLinePlot, LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameStyle -> Black, ImageSize -> imageSize, 
  PlotRange -> {0, All}]; 
SetOptions[ListPlot, LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameStyle -> Black, ImageSize -> imageSize, 
  PlotRange -> {0, All}]; 
SetOptions[ListLogPlot, 
  LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameStyle -> Black, ImageSize -> imageSize, 
  PlotRange -> {0, All}]; 
SetOptions[ListContourPlot, 
  LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameStyle -> Black, ImageSize -> imageSize]; 
SetOptions[ContourPlot, 
  LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameStyle -> Black, ImageSize -> imageSize]; 
SetOptions[DiscretePlot, LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     ImageSize -> imageSize, FrameStyle -> Black, 
  PlotRange -> {0, All}]; 
SetOptions[Graphics, LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, Frame -> False, 
  FrameTicksStyle -> Black, 
     ImageSize -> imageSize, FrameStyle -> Black, 
  PlotRange -> {0, All}]; 
SetOptions[Histogram, LabelStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameTicksStyle -> Black, FrameStyle -> Black, 
  ImageSize -> imageSize, 
     PlotRange -> {0, All}, ChartStyle -> EdgeForm[Thin]]; 
SetOptions[SmoothDensityHistogram, 
  FrameStyle -> Directive[Black, FontSize -> 9, 
         FontFamily -> "Helvetica"], Axes -> False, 
  Frame -> {True, True, False, False}, 
     FrameTicksStyle -> Black, FrameStyle -> Black, 
  ImageSize -> imageSize, 
     PlotRange -> {0, All}]; 

Timeseries

Data*import

data = Import["out.timeseries", "Table"]; 
headerRules = MapIndexed[#1 -> #2[[1]] & , data[[1]]];
data = Drop[data, 1]; 
date = "date" /. headerRules;
startTime = data[[1, date]];
endTime = data[[-1, date]];
stepTime = data[[2, date]] - data[[1, date]];
popSize = data[[1, "totalN" /. headerRules]];
demePopSize = popSize/demeCount;

Prevalence

index = "totalI" /. headerRules;
N[Mean[data[[All, index]]]]
N[HarmonicMean[data[[All, index]]]]
prev = Transpose[{data[[All, date]], (data[[All, index]]*100000)/
     popSize}]; 

Incidence

index = "totalCases" /. headerRules;

Mean*incidence*per*year

meaninc = N[Total[data[[All, index]]]/endTime/popSize];
Export["meanAnnualIncidence.csv",meaninc,"CSV"];


Deme - specific*yearly*incidence

indices = (StringJoin[#1, "Cases"] & ) /@ demes /. headerRules;
incList = Table[data[[All, i]], {i, indices}]; 
meanIncList = (Total[#1] & ) /@ 
    incList/((endTime - startTime)*(popSize/demeCount)); 
Export["meanAnnualDemeIncidence.csv",Round[meanIncList, 0.001],"CSV"];


Deme - specific*prevalence
indices = (StringJoin[#1, "I"] & ) /@ demes /. headerRules;
prevList = 
  N[Table[Transpose[{data[[All, date]], 
      data[[All, i]]*(100000/demePopSize)}], 
         {i, indices}]]; 
Export["prevalenceSpatial.csv", Round[prevList,0.01], "CSV"];

Total*infecteds

indices = (StringJoin[#1, "I"] & ) /@ demes /. headerRules;
prevList = 
  N[Table[Transpose[{data[[All, date]], data[[All, i]]}], {i, 
     indices}]]; 

Deme - specific*incidence
indices = (StringJoin[#1, "Cases"] & ) /@ demes /. headerRules;
incList = N[Table[Transpose[{data[[All, date]], 
	data[[All, i]]*(100000/demePopSize)}], {i, indices}]];
max = Max[Max /@ incList[[All, All, 2]]];


Diversity

index = "diversity" /. headerRules;
div = data[[All, {date, index}]]; 
meandiv = Mean[div][[2]];
Export["meanDiversity.csv",Round[meandiv, 0.001],"CSV"];

Deme - specific*diversity

indices = (StringJoin[#1, "Diversity"] & ) /@ demes /. headerRules;
divList = 
  N[Table[Transpose[{data[[All, date]], data[[All, i]]}], {i, 
     indices}]]; 
meanDivList = ((Mean[#1] & ) /@ divList)[[All, 2]]; 
Export["meanDemeDiversity.csv", Round[meanDivList, 0.001], "CSV"];

TMRCA

index = "tmrca" /. headerRules;
tmrca = data[[All, {date, index}]]; 
meantmrca = Mean[tmrca][[2]];
Export["meanTMRCA.csv",meantmrca,"CSV"];

Deme - specific*TMRCA

indices = (StringJoin[#1, "Tmrca"] & ) /@ demes /. headerRules;
tmrcaList = 
  N[Table[Transpose[{data[[All, date]], data[[All, i]]}], {i, 
     indices}]]; 
meanTmrcaList = ((Mean[#1] & ) /@ tmrcaList)[[All, 2]]; 
Export["meanDemeTMRCA.csv",Round[meanTmrcaList, 0.001],"CSV"];

Ne*tau

index = "netau" /. headerRules;
netau = data[[All, {date, index}]]; 
netau = Select[netau, NumberQ[#1[[2]]] & ]; 
meannetau = Mean[netau][[2]];
Export["meanNeTau.csv",meannetau,"CSV"];


Deme - specific*Ne*tau

indices = (StringJoin[#1, "Netau"] & ) /@ demes /. headerRules;
netauList = 
  N[Table[Transpose[{data[[All, date]], data[[All, i]]}], {i, 
     indices}]]; 
netauList = (Select[#1, NumberQ[#1[[2]]] & ] & ) /@ netauList; 
meanNetauList = ((Mean[#1] & ) /@ netauList)[[All, 2]]; 
Export["meanDemeNetau.csv",Round[meanNetauList, 0.001],"CSV"];


Virus*samples

tipScale = 0.0015; 
tipData = Import["out.tips", "CSV"]; 
traitHeader = tipData[[1]]; 
TableForm[traitHeader, TableHeadings -> {Range[10]}];
tipData = Drop[tipData, 1]; 
name = 1; year = 2; trunk = 3; tip = 4; mark = 5; loc = 6; layout = \
7; ag1 = 8; ag2 = 9; 

First*tip

firstIndex = First[Ordering[tipData[[All, year]], 1]];
length = Length[tipData];

Graphics*settings

frameWindow = 1.; 
frameTicks = 
  Table[Table[{i, i, {0.01, 0}}, {i, -100, 100, 
     frameWindow}], {2}, {2}]; 

Map

noise[] := With[{theta = RandomReal[{0, 2*Pi}], 
         r = RandomVariate[ExponentialDistribution[1.9]]}, {r*
     Cos[theta], r*Sin[theta]}]; 
jitter = (#1 + noise[] & ) /@ tipData[[All, {ag1, ag2}]]; 
plotRange = {{Floor[Min[jitter[[All, 1]]]], 
    Ceiling[Max[jitter[[All, 1]]]]}, 
      {Floor[Min[jitter[[All, 2]]]], Ceiling[Max[jitter[[All, 2]]]]}};
avgDist[pointsA_, pointsB_] := Mean[MapThread[EuclideanDistance, 
       {RandomChoice[pointsA, 1000], RandomChoice[pointsB, 1000]}]]
clustered = FindClusters[jitter -> Range[length], k, 
       Method -> {"Agglomerate", "Linkage" -> "Ward"}, DistanceFunction -> EuclideanDistance]; 
cluster[i_] := jitter[[clustered[[i]]]]

Cluster*function

clusterCenters = (Mean[jitter[[#1]]] & ) /@ clustered; 
clusterFunction[p_] := First[Nearest[clusterCenters -> Range[k], p]]
points = ({clusterColors[[clusterFunction[#1]]], Point[#1]} & ) /@ 
   jitter; 
Graphics[points, AspectRatio -> Automatic, PlotRange -> plotRange, 
     GridLines -> {Range[-100, 100, frameWindow], Range[-100, 100, frameWindow]}, 
     GridLinesStyle -> LightGray, Frame -> False, ImageSize -> imageSize*0.95]; 

Distance*from*origin

	   dist[point_] := 
 EuclideanDistance[tipData[[firstIndex, {ag1, ag2}]], point]
list = ({#1[[1]], dist[#1[[{2, 3}]]]} & ) /@ 
   tipData[[All, {year, ag1, ag2}]]; 
meanFluxRate = (Max[list[[All, 2]]] - 
    Min[list[[All, 2]]])/(Max[tipData[[All, year]]] - 
        Min[tipData[[All, year]]])
Export["meanFluxRate.csv",meanFluxRate,"CSV"];

Cluster*turnover
times = tipData[[All,year]]; 
step = N[(4*7)/365]; 
step = stepTime;
clusters = (times[[#1]] & ) /@ clustered; 
datelist = Table[i, {i, step, endTime, step}]; 
totalCounts = BinCounts[Flatten[clusters], {0, endTime, step}]; 
pro = Table[DeleteCases[N[MapThread[If[#2 == 0, Null, {#3, #1/#2}] & , {BinCounts[clusters[[i]], {0, endTime, step}], 
	totalCounts, datelist}]], Null], {i, k}]; 
loess = Table[({#1, LoessFit[#1, pro[[j]], 0.05, 1]} & ) /@ Table[i, {i, incList[[1, All, 1]]}], {j, k}];


(Spatial*comparison)/antigenic

Loess*spline

data = ({#1[[1]], dist[#1[[{2, 3}]]]} & ) /@ 
   tipData[[All, {year, ag1, ag2}]]; 
fit[x_] := LoessFit[x, data, 0.2, 1]
startTime = Min[data[[All, 1]]];
endTime = Max[data[[All, 1]]];
max = Max[data[[All, 2]]];
loessList = Table[{x, fit[x]}, {x, startTime, endTime, 0.1}]; 

distance[list_] := 
 Module[{nf}, 
  nf = Nearest[loessList[[All, 1]] -> loessList[[All, 2]]]; 
      (#1[[2]] - First[nf[#1[[1]]]] & ) /@ list]
standardError[list_] := StandardDeviation[list]/Sqrt[Length[list]]
standardError[distance[data]]

spatialDists = 
  Table[With[{dist = distance[({#1[[1]], dist[#1[[{2, 3}]]]} & ) /@ 
                 
        Select[tipData, #1[[loc]] == i & ][[All, {year, ag1, ag2}]]]}, 
         {Mean[dist] - standardError[dist], Mean[dist], 
     Mean[dist] + standardError[dist]}], 
       {i, 0, demeCount - 1}]; 

{meanDistNorth, meanDistTropics, meanDistSouth} = 
     Table[
   With[{dist = distance[({#1[[1]], dist[#1[[{2, 3}]]]} & ) /@ 
                 
        Select[tipData, #1[[loc]] == i & ][[
         All, {year, ag1, ag2}]]]}, Mean[dist]], 
       {i, 0, demeCount - 1}]; 
Export["antigenicLagAG1.csv", Round[spatialDists, 0.001], "CSV"];

Tree

Data*import

branchScale = 0.001; tipScale = 0.0015; clip = 0.005; 
getThickness[count_] := 
 Thickness[Clip[branchScale*N[Sqrt[count]], {0, clip}]]
getPointSize[count_] := PointSize[N[tipScale*Sqrt[count]]]; 
name = 1; year = 2; trunk = 3; tip = 4; mark = 5; loc = 6; layout = \
7; ag1 = 8; ag2 = 9; 
branchData = 
  ToExpression[Import["out.branches", "Table"]]; 

Reduction

slimBranches = 
  Select[branchData, #1[[1, mark]] == 1 && #1[[2, mark]] == 1 & ]; 


Migration*statistics

trunkBranches = 
  Select[branchData, #1[[1, trunk]] == 1 && #1[[1, year]] < 
      endTime - 5 & ]; 
sideBranches = 
  Select[branchData, #1[[1, trunk]] == 0 && #1[[1, year]] < 
      endTime - 5 & ]; 

Proportion*of*each*deme*in*trunk

trunkLengths = 
  Table[Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ 
     Select[trunkBranches, #1[[1, loc]] == i && #1[[2, loc]] == 
         i & ]], 
       {i, 0, demeCount - 1}]; 
trunkProportions = trunkLengths/Total[trunkLengths]; 
{trunkProNorth, trunkProTropics, trunkProSouth} = 
  trunkLengths/Total[trunkLengths]; 
Export["trunkProportions.csv", Round[trunkProportions, 0.001], "CSV"];

Migration*rate

trunkOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ trunkBranches]; 
sideBranchOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ sideBranches]; 
trunkCount = 
  Total[(If[#1[[1, loc]] != #1[[2, loc]], 1, 0] & ) /@ trunkBranches]; 
sideBranchCount = 
  Total[(If[#1[[1, loc]] != #1[[2, loc]], 1, 0] & ) /@ sideBranches]; 

trunkMigRate = trunkCount/trunkOpportunity;
Export["trunkMigRate.csv", Round[trunkMigRate, 0.001], "CSV"];

sideBranchMigRate = sideBranchCount/sideBranchOpportunity;
Export["sideBranchMigRate.csv", Round[sideBranchMigRate, 0.001], "CSV"];


Migration*matrix--*trunk

trunkCounts = 
  Table[Total[(If[#1[[1, loc]] == j && #1[[2, loc]] == i, 1, 
        0] & ) /@ 
           trunkBranches], {i, 0, demeCount - 1}, {j, 0, 
    demeCount - 1}]; 
trunkOpportunies = 
  Table[Total[(If[#1[[1, loc]] == i && #1[[2, loc]] == i, 
                Abs[#1[[1, year]] - #1[[2, year]]], 0] & ) /@ 
     trunkBranches], {i, 0, demeCount - 1}, 
       {j, 0, demeCount - 1}]; 
trunkRates = (trunkCounts/trunkOpportunies)*(1 - 
     DiagonalMatrix[Table[1, {demeCount}]]); 

TableForm[Round[trunkRates, 0.001] /. 0. -> "", 
 TableHeadings -> {demes, demes}]
Export["trunkMigrationRateMatrix.csv", 
  Round[trunkRates, 0.001] /. 0. -> "", "CSV"];

Migration*matrix--*side*branches

sideBranchCounts = 
  Table[Total[(If[#1[[1, loc]] == j && #1[[2, loc]] == i, 1, 
        0] & ) /@ 
           sideBranches], {i, 0, demeCount - 1}, {j, 0, 
    demeCount - 1}]; 
sideBranchOpportunies = 
     Table[
   Total[(If[#1[[1, loc]] == i && #1[[2, loc]] == i, 
        Abs[#1[[1, year]] - #1[[2, year]]], 
                0] & ) /@ sideBranches], {i, 0, demeCount - 1}, {j, 
    0, demeCount - 1}]; 
sideBranchRates = (sideBranchCounts/sideBranchOpportunies)*
       (1 - DiagonalMatrix[Table[1, {demeCount}]]); 

TableForm[Round[sideBranchRates, 0.001] /. 0. -> "", 
 TableHeadings -> {demes, demes}]
Export["sideMigrationRateMatrix.csv", 
  Round[sideBranchRates, 0.001] /. (0. -> ""), "CSV"];

Antigenic*statistics

trunkBranches = 
  Select[branchData, #1[[1, trunk]] == 1 && #1[[1, year]] < 
      endTime - 5 & ]; 
sideBranches = 
  Select[branchData, #1[[1, trunk]] == 0 && #1[[1, year]] < 
      endTime - 5 & ]; 

Mutation*rate

trunkOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ trunkBranches]; 
sideBranchOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ sideBranches]; 
trunkCount = 
  Total[(If[
       EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[2, {ag1, ag2}]]] > 
        1.*^-7, 1, 
              0] & ) /@ trunkBranches]; 
sideBranchCount = 
  Total[(If[
       EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[2, {ag1, ag2}]]] > 
                1.*^-7, 1, 0] & ) /@ sideBranches]; 

trunkMutationRate = trunkCount/trunkOpportunity;
sideBranchMutationRate = sideBranchCount/sideBranchOpportunity;
										 
Antigenic*flux

trunkOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ trunkBranches]; 
sideBranchOpportunity = 
  Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ sideBranches]; 
trunkDistance = 
  Total[(EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
        2, {ag1, ag2}]]] & ) /@ 
         trunkBranches]; 
sideBranchDistance = 
  Total[(EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
        2, {ag1, ag2}]]] & ) /@ 
         sideBranches]; 
trunkFluxRate = trunkDistance/trunkOpportunity;
sideBranchFluxRate = sideBranchDistance/sideBranchOpportunity;

Mutation*size

trunkMutations = 
  Select[(EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
        2, {ag1, ag2}]]] & ) /@ 
         trunkBranches, #1 > 1.*^-7 & ]; 
sideBranchMutations = Select[(EuclideanDistance[#1[[1, {ag1, ag2}]], 
              #1[[2, {ag1, ag2}]]] & ) /@ 
    sideBranches, #1 > 1.*^-7 & ]; 
trunkMutationSize = Mean[trunkMutations];
sideBranchMutationSize = Mean[sideBranchMutations]

Combined

Export["mutationSizeFlux.csv",Round[{{sideBranchMutationRate, sideBranchMutationSize, 
	  sideBranchFluxRate}, {trunkMutationRate, trunkMutationSize, trunkFluxRate}}, 0.001]];

Mutation*rate - spatial

trunkRates = 
  Table[branches = 
    Select[trunkBranches, #1[[1, loc]] == i && #1[[2, loc]] == 
        i & ]; 
        opportunity = 
    Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ branches]; 
        count = 
    Total[(If[
         EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
            2, {ag1, ag2}]]] > 1.*^-7, 1, 
                   0] & ) /@ branches]; 
   rate = count/opportunity, {i, 0, demeCount - 1}]; 
sideBranchRates = Table[branches = Select[sideBranches, 
            #1[[1, loc]] == i && #1[[2, loc]] == i & ]; 
   opportunity = 
          
    Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ branches]; 
        count = 
    Total[(If[
         EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
            2, {ag1, ag2}]]] > 1.*^-7, 1, 
                   0] & ) /@ branches]; 
   rate = count/opportunity, {i, 0, demeCount - 1}]; 
Export["spatialMutationRate.csv",Round[{sideBranchRates, trunkRates}, 0.001]];						


Mutation*size - spatial

trunkSize = 
  Table[branches = 
    Select[trunkBranches, #1[[1, loc]] == i && #1[[2, loc]] == 
        i & ]; 
        Mean[
    Select[(EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
          2, {ag1, ag2}]]] & ) /@ branches, 
            #1 > 1.*^-7 & ]], {i, 0, demeCount - 1}]; 
sideBranchSize = Table[branches = Select[sideBranches, 
            #1[[1, loc]] == i && #1[[2, loc]] == i & ]; 
        Mean[
    Select[(EuclideanDistance[#1[[1, {ag1, ag2}]], #1[[
          2, {ag1, ag2}]]] & ) /@ branches, 
            #1 > 1.*^-7 & ]], {i, 0, demeCount - 1}]; 
Export["spatialMutationSize.csv",Round[{sideBranchSize, trunkSize}, 0.001]];						


Antigenic*flux - spatial

trunkRates = 
  Table[branches = 
    Select[trunkBranches, #1[[1, loc]] == i && #1[[2, loc]] == 
        i & ]; 
        opportunity = 
    Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ branches]; 
        antigenicDistance = 
    Total[(EuclideanDistance[#1[[1, {ag1, ag2}]], 
                   #1[[2, {ag1, ag2}]]] & ) /@ branches]; 
   rate = antigenicDistance/opportunity, 
       {i, 0, demeCount - 1}]; 
sideBranchRates = Table[branches = Select[sideBranches, 
            #1[[1, loc]] == i && #1[[2, loc]] == i & ]; 
   opportunity = 
          
    Total[(Abs[#1[[1, year]] - #1[[2, year]]] & ) /@ branches]; 
        antigenicDistance = 
    Total[(EuclideanDistance[#1[[1, {ag1, ag2}]], 
                   #1[[2, {ag1, ag2}]]] & ) /@ branches]; 
   rate = antigenicDistance/opportunity, 
       {i, 0, demeCount - 1}]; 
Export["antigenicFlux.csv",Round[{sideBranchRates, trunkRates}, 0.001],"CSV"];

Exit[]
