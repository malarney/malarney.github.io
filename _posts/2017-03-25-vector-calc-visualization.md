---
title: "Visualizing Vector Calculus"
excerpt: "Using Mathematica to visualize various functions in vector calculus"
tags: 
  - mathematica
  - programming
  - math
---
One of the more difficult aspects of vector and multivariable calculus when transitioning from single variables is visualizing the functions and operations that you begin working with regularly. Simple 3D graphing comes native with Mathematica, and it’s an invaluable tool for many purposes. Below are some supplementary functions written in Mathematica, that can help in understanding a few ideas in vector calculus. Namely: Gradient fields, and line integrals over a function and vector field in 2-dimensions.

## Demonstration
### Gradient Field with Function
The gradient field function takes in an expression of two variables, after which it will:
* Plot a 3D graph of the function on a given range
* Color the function by height
* Place its gradient field below the minimum value
* Color the gradient field by magnitude

The fact that the gradient will point in the direction of steepest ascent becomes quite visible, with minimums and maximums becoming sources and sinks respectively in the gradient field.

This is demonstrated in Figure 1:
<figure>
    <a href="/assets/images/Sin2x+Cos2y_GradientField.svg"><img src="/assets/images/Sin2x+Cos2y_GradientField.svg"></a>
    <figcaption>Figure 1: \(\cos^{2}\left(x\right)+\sin^{2}\left(y\right)\) demonstrating vanishing gradients at extrema</figcaption>
</figure>
This usually implies that the blue regions in the slope field, which represent low gradient magnitudes, contain a critical point. This of course includes saddles, which are more properly represented in Figure 2, which also demonstrates where saddle points in slope fields get their namesake:
<figure>
    <a href="/assets/images/16y2-25x2_GradientField.svg"><img src="/assets/images/16y2-25x2_GradientField.svg"></a>
    <figcaption>Figure 2: \(16y^{2}-25x^{2}\) showing vanishing gradient at saddle</figcaption>
</figure>
Ready understanding of these diagrams leads to a fairly intuitive grasp on why vector fields which possess potential functions are conservative:

* The gradient dotted with a unit vector is the directional derivative
* A line integral is then summing directional derivatives of the potential function
* The direction of these derivatives at a point will be the unit tangent vector of the parameterization
* Directional derivatives are how "steep" the function is in that direction
* Therefore you are summing infinitesimal slopes of the hill, so the net change will be the change in height.
* The change in height is of course, independent of path

This immediately applies to energy in physics, with the small caveat that you gain energy when you go down a hill, instead of losing height. This is remedied by a single change in sign:

$$\mathbf{F}=-\nabla P$$

All the sinks become sources, and all the sources sinks, allowing us to view the path integral of the vector field as a work integral for gravitational energy when traversing along the potential function's surface.

### Line Integral of Function
Line integrals of 2-variable expressions can also be plotted easily, the process goes roughly:
* Graph the function that is to be integrated over
* Calculate and graph the parameterization on the xy-plane
* Plot upwards to the parameterization on the surface, creating a “wall”

This shows the rather simple geometric interpretation of this type of integral, as the area of a “wall” over a parameterization whose height is determined by the function.

An example is given by Figure 3:
<figure>
    <a href="/assets/images/cos2y+sin2x_LineIntegral.svg"><img src="/assets/images/cos2y+sin2x_LineIntegral.svg"></a>
    <figcaption>Figure 3: \(\sin^{2}\left(x\right)+\cos^{2}(y)\) on unit circle </figcaption>
</figure>
In the future, it may be both useful and aesthetically pleasing to apply this to 3-variable functions by plotting the 4th value as the magnitude of the normal vector along the parameterization. As well as this, graphing the dx and dy integrals as "shadows" may prove to be an interesting exercise for geometrically showing their link to line integrals over vector fields.

### Line Integral of Vector Field
Line Integrals on Vector Fields are graphed in something of a hybrid manner of the last two modules. The mechanism goes as follows:

* The vector plot of the field is put onto the xy plane
* A mesh is constructed starting at xy plane
*  The mesh is defined by the parameterization in the x and y directions
* Mesh is bounded by the tangent vector dotted with the field, in the z direction

Line integrals on vector fields don’t share as handy a geometric interpretation as their scalar counterparts, as the direction the path is taking will make a profound impact on the value a point gives out. Graphing in a similar matter as line integrals, however, does prove useful in demonstrating the integrals over closed loops of conservative vs non-conservative vector fields.

A conservative field is demonstrated in Figure 4:
<figure>
    <a href="/assets/images/x4y5Circle_VectorIntegral.svg"><img src="/assets/images/x4y5Circle_VectorIntegral.svg"></a>
    <figcaption>Figure 4: \(\left[\begin{array}{c}
4x^{3}\\
5y^{4}
\end{array}\right]\) displaying zero integral on closed loop</figcaption>
</figure>
Due both to the symmetry of the parametrization, and the field it is on, the integral has pairs of “walls” that are equal in magnitude but opposite in sign. More irregular parameterizations would of course sum to zero, but in a less straightforward manner. These pairs would clearly cancel in any full path along this circle, and so the integral along the closed loop is zero.

A non-conservative example is demonstrated in Figure 5:
<figure>
    <a href="/assets/images/ArgEpitrochoid_VectorIntegral.svg"><img src="/assets/images/ArgEpitrochoid_VectorIntegral.svg"></a>
    <figcaption>Figure 5: \(\nabla\tan^{-1}\left(\frac{y}{x}\right)\), the angle of a point, over an epitrochoid</figcaption>
</figure>
Here we have what can be called the gradient of {% raw %}\(\tan^{-1}\left(\frac{y}{x}\right)\){% endraw %}. This potential function is the standard angle when \(x>0\), and the angle with the negative-x ray everywhere else it's defined.  This is an interesting case as despite the field having a potential function, a closed positively oriented curve will yield a positive integral. Of course, what we're actually differentiating with respect to is just a branch of a multi-function, there are infinitely many functions that would satisfy:

$$f\left(\tan\left(x\right)\right)=x$$

What is really causing the non-conservative behavior is jump discontinuities. As we go in a positively oriented circular path around the origin, we will always be "walking uphill". What allows us to end up in the same spot is the jump discontinuities along the line \(x=0\). We are effectively "teleported" to a lower position, without walking downhill. This highlights the importance of having the field and potential be continuous in a domain.

This is shown in Figure 6:
<figure>
    <a href="/assets/images/arctan_GradientField.svg"><img src="/assets/images/arctan_GradientField.svg"></a>
    <figcaption>Figure 6: \(\tan^{-1}\left(\frac{y}{x}\right)\) showing how jump discontinuities affect path independence</figcaption>
</figure>

## Source Code
#### Gradient Field with Function
``` ocaml
{% raw %}GradientField[funct_, range_] :=
    Module[{functiongraphic, minimum, vectortexture, xysquare, slopegraphic, gradient},
        minimum =NMinValue[{funct, -range <= x <= range, -range <= y <=range}, {x, y}];
        gradient = D[funct, {{x, y}}];
        functiongraphic =Plot3D[funct, {x, -range, range}, {y, -range, range},ColorFunction -> "DarkRainbow",PlotStyle -> {Specularity[White, 40], Opacity[.8]},PlotRange -> Full];
        vectortexture =Texture@VectorDensityPlot[{gradient[[1]],gradient[[2]]}, {x, -1*range, range}, {y, -1*range, range}, Frame -> False, ImageSize -> Large, PlotRangePadding -> None, ColorFunction -> "DarkRainbow", VectorPoints -> 30,VectorStyle -> White];
        xysquare = Polygon[{{-range, -range, minimum}, {range, -range, minimum}, {range, range, minimum}, {-range, range, minimum}}, VertexTextureCoordinates -> {{0, 0}, {1, 0}, {1, 1}, {0, 1}}];
        slopegraphic = Graphics3D[{vectortexture, xysquare}, Lighting -> "Ambient"];
        Show[functiongraphic, slopegraphic, ImageSize -> Full,AxesLabel -> {x, y, z}]
    ]{% endraw %}
```
#### Line Integral of Function

``` ocaml
{% raw %}LineIntegralPlot[function_, param_, range_] :=
 	Module[{rx, ry, integrand, xygraphic, integralgraphic, fullgraphic},
  		rx = param[[1]];
  		ry = param[[2]];
  		integrand = (function /. {x -> rx, y -> ry});
  		integralgraphic = ParametricPlot3D[{rx, ry, z *integrand}, {t, 0, 1}, {z, 0, 1}, BoundaryStyle -> Directive[Thick, Black], BoxRatios -> {1, 1, 1/2}, PlotStyle -> Opacity[.5], PlotRange -> All];
  		xygraphic = Plot3D[0, {x, -range, range}, {y, -range, range}, PlotStyle -> Opacity[.2], MeshStyle -> Opacity[.2]];
  		fullgraphic = Plot3D[function, {x, -range, range}, {y, -range, range}, PlotStyle -> Opacity[.15], PlotStyle -> Opacity[.15]];
  		Show[integralgraphic, xygraphic, fullgraphic, ImageSize -> Full, AxesLabel -> {x, y, z}]
  ]{% endraw %}
```
#### Line Integral of Vector Field
``` ocaml
{% raw %}VectorIntegralPlot[funct_, param_, range_] :=
    Module[{integrand, fx, fy, rx, ry, xygraphic, vectorplot, slopetexture, xysquare, slopegraphic, integralgraphic},
        fx = funct[[1]];
        fy = funct[[2]];
        rx = param[[1]];
        ry = param[[2]];
        integrand = ((fx /. {x -> rx, y -> ry})*((D[rx, t])) + (fy /. {x -> rx, y -> ry})*(D[ry, t]))/Sqrt[(D[rx, t])^2 + (D[ry, t])^2];
        xygraphic = Plot3D[0, {x, -range, range}, {y, -range, range}, PlotStyle -> Opacity[.2], MeshStyle -> Opacity[0]];
        vectorplot = VectorPlot[{fx, fy}, {x, -range, range}, {y, -range, range}, Frame -> False, VectorStyle -> {Black}, Background -> None, VectorPoints -> Fine, ImageSize -> Large];
        slopetexture = Texture@ImageData@RemoveBackground@Rasterize[vectorplot];
        xysquare = Polygon[{{-range, -range, 0.001}, {range, -range, 0.001}, {range, range, 0.001}, {-range, range, 0.001}}, VertexTextureCoordinates -> {{0, 0}, {1, 0}, {1, 1}, {0, 1}}];
        slopegraphic = Graphics3D[{slopetexture, xysquare}, Lighting -> "Neutral"];
        integralgraphic = ParametricPlot3D[{rx, ry, z*integrand}, {t, 0, 1}, {z, 0, 1}, BoundaryStyle -> Directive[Thick, Black], PlotStyle -> Opacity[.5]];
        Show[xygraphic, slopegraphic, integralgraphic, PlotRange -> All, ImageSize -> Full, BoxRatios -> {1, 1, 1/2}, Axes -> True, AxesLabel -> {x, y, z}]
    ]{% endraw %}
```
