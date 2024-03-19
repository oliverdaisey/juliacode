# TropicalHomotopies.jl
🚀 A Julia package for tropical homotopy continuation on a wide range of tropical spaces, utilising OSCAR for much of its functionality. Features include fast computations of intersections of tropical hypersurfaces, tropical linear spaces, and inverted tropical linear spaces.
## Usage
Start by defining the inital data for the supports of the spaces whose intersection you want to track. For a hypersurface, one can specify the tropical polynomial directly:
```julia
TT = tropical_semiring()
R, (x, y) = PolynomialRing(TT, ["x", "y"])
f = 1*x^4*y^2 + (-1)*x^2*y^4 + (-1)*x^3*y^3 + (-1)*x^2*y^2 + (1)*x^3*y + (1)*x^2*y^3 + (-1)*x*y^2 + (-1)*x*y + (1)*x + (-1)*y
fSupport = DualSupport{Hypersurface}(generate_support(f))
```
For an intersection involving a (inverted) tropical linear space, construct the vertices of the matroid polytope and specify the inital data as follows:
```julia
matroidVertices = [1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1]
linearSupport = DualSupport{Linear}(matroidVertices)
invertedLinearSupport = DualSupport{InvertedLinear}(matroidVertices)
```
Then define the path data for the homotopy. Begin by specifying the initial and target weights for the supports, and optionally add intermediate weights (this is necessary when dealing with a tropical linear space, to avoid leaving the Dressian). 
```julia
f = x*y + 1
fNodes = [TT.([0, 1]), TT.([1, 1])] # fNodes[1] is the initial dual weight vector
fPath = DualPath{Hypersurface, typeof(min)}(fNodes, fSupport)
```
Then specify the dual cells that support the initial mixed cells.
```julia
f_start_dual = DualCell{Hypersurface, typeof(min)}(fSupport, [1, 2], fNodes[1]) # [1, 2] are indices of fSupport that make up this dual cell
```
## License
This package (will be) released under the MIT License.
