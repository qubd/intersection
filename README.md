intersection
===============

Compute intersection multiplicities of algebraic curves at rational points.

```
>>> from intersection import *
```

To find the intersection multiplicity of the affine curves f(x,y) = 0 and g(x,y) = 0 at P(x,y), start by specifying the polynomials f and g.

```
>>> f = AffineCurve([[0],[0,1],[0,0,0],[-1,0,0,0]])
>>> print f
- x^3 + y
>>> g = AffineCurve([[0],[0,2],[1,0,1]])
>>> print g
x^2 + y^2 + 2y
```

The argument passed to AffineCurve is a list of lists of coefficients. The entry in the i^th list, in the j^th position is the coefficient on x^(i-j)y^(j). The polynomial 1 + 2x + 3y + 4x^2 + 5xy + 6y^2, for example, would be specified with AffineCurve([[1],[2,3],[4,5,6]]).

This package implements arithmetic for such polynomials.

```
>>> print f + 2*g
- x^3 + 2x^2 + 2y^2 + 5y
>>> print f.pow(3)
- x^9 + 3x^6y - 3x^3y^2 + y^3
```

We can also check whether a rational point is on the curve f(x,y) = 0, and perform a change of coordinates by translation, so the new origin is at any point P.

```
>>> f.contains((0,0))
True
>>> f.contains((-1,0))
False
>>> f.contains((1,1))
True
>>> h = f.shift((1,0))
>>> print h
- x^3 - 3x^2 - 3x + y - 1
>>> h.contains((-1,0))
True
```

Lastly, we can compute the intersection multiplicity of two curves at the origin, or at any other rational point in the plane.

```
>>> I_0(f,g)
2
>>> I(f,g,(1,1))
0
>>> I_0(f,f)
'Infinite'
```
