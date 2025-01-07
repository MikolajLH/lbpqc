LBPQC is build on top of `numpy` so the most common type on which operations are performed is `numpy.ndarray`.

However in order to improve clarity of the code we rely heavily on python's type aliases that allow to differentiate between e.g vectors and matrices (which both are represented as `ndarray`).


## `int`
Aliases for `int` are only used as an information about what kind of result does the functino returns.

### `ModInt`
```python3
type ModInt = int
```
Represents integers from interval $[0, q)$ for some positive integer $q$.

### `CenteredModInt`
```python3
type CenteredModInt = int
```
Represents integers from interval $[-\frac{q}{2}, \frac{q}{2})$ for some positive integer $q$.


## `Vector`

### `Vector`
```python3
type Vector = np.ndarray
```

### `VectorFloat`
```python3
type Vector = np.ndarray[float]
```

### `VectorInt`
```python3
type Vector = np.ndarray[int]
```

### `VectorModInt`
```python3
type Vector = np.ndarray[int]
```

### `VectorCenteredModInt`
```python3
type Vector = np.ndarray[int]
```


## `Matrix`

### `Matrix`
```python3
type Matrix = np.ndarray
```

### `MatrixFloat`
```python3
type Matrix = np.ndarray[float]
```

### `MatrixInt`
```python3
type Matrix = np.ndarray[int]
```

### `MatrixModInt`
```python3
type Matrix = np.ndarray[int]
```

### `MatrixCenteredModInt`
```python3
type Matrix = np.ndarray[int]
```