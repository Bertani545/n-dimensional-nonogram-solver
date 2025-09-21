# n-dimensional-nonogram-solver
A c++ implementation for solving nonograms in n-dimensions


Create a class from this solver with: `Nonogram nonogram`.

You can either build it from an already existing puzzle:
`nonogram.buildFromSquares(path)`
where the file have the next format:
```
n_dims size_dim_1 size_dim_2 ...
data_of_1_and_0
```
For example,
```
2 5 3
0 1 1 1 0
1 1 0 1 1
0 0 1 0 1
```
Please notice that the information is listed starting from the first dimension.
In `3 x y z`, every `x` numbers you change 1 index in `y` and every `xy` numbers you chnage 1 index in `z`.


Or from hint lists:
`nonogram.nonogram.buildFromLists(path)`
where the file have the next format:
```
n_dims size_dim_1 size_dim_2 ...
hint_lists
```
For example, for the same nonogram as before,
```
2 5 3
3
2 2
1 1
1
2
1 1
2
2
```
Please notice that the lists are ordered by how you defined the dimensions, using the same rules as before.
