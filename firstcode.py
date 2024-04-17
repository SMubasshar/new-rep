from dune.grid import structuredGrid
grid = structuredGrid([0,0],[1,1],[10,10])
print( grid.size(0) ) # should return 100
grid.plot() # requires the installation of `matplotlib`