# NetworkDirectionality
Code to undertake the calculations in Keylock and Carbone (2026). Phys. Rev. E

The primary file is NetworkLaplacianDirectionality.m
This takes an adjacency matrix and two discretisation parameters as inputs and 
returns a set of summary indices and a more detailed output. Associated with this
file are three internal functions:
DGL.m
BetaTau.m
BetaTauImposeEpsTau.m

The latter can then be passed to PlottingToolsUtility.m to show some results
