---
title: "Visualizing Vector Calculus"
excerpt: "Using Mathematica to visualize various functions in vector calculus"
tags: 
  - mathematica
  - programming
---
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

This usually implies that the blue regions, which represent low gradient magnitudes, contain a critical point. This of course includes saddles, which are more properly represented in Figure 2, which also demonstrates where saddle points in slope fields get their namesake:

<figure>
    <a href="/assets/images/16y2-25x2_GradientField.svg"><img src="/assets/images/16y2-25x2_GradientField.svg"></a>
    <figcaption>Figure 2: \(16y^{2}-25x^{2}\) showing vanishing gradient at saddle</figcaption>
</figure>

Ready understanding of these diagrams leads to a fairly intuitive grasp on why vector fields which possess gradients are conservative.

* Line integrals are work integrals.
* The gradient dotted with the direction vector
Line integrals on vector fields are work integrals, and the gradient, so by extens is effectively giving you the opposite how steep the slope is . the work it takes to go up a hill, given by the potential function, in that direction. Giving us:

$$\mathbf{F}=-\nabla P$$

Of course, were it not for friction, the amount of energy it would take to go a distance uphill would equal the amount you would gain going down a similar height. Therefore, the energy change should be proportional to the height difference, which gives a great intuition for the gradient theorem.

### Line Integral of Function
Line integrals of 2-variable expressions can also be plotted easily, the process goes roughly:
* Graph the function that is to be integrated over
* Calculate and graph the parameterization on the xy-plane
* Plot upwards to the parameterization on the surface, creating a “wall”

This shows the rather simple geometric interpretation of this type of integral, as the area of a “wall” over a parameterization whose height is determined by the function.

An example is given by Figure 3:

<figure>
    <a href="/assets/images/x2+2xy+y2Circle_LineIntegral.svg"><img src="/assets/images/x2+2xy+y2Circle_LineIntegral.svg"></a>
    <figcaption>Figure 3: \(\left(x+y\right)^{2}\) over a circle of radius 2 </figcaption>
</figure>



### Line Integral of Vector Field

<figure>
    <a href="/assets/images/x4y5Circle_VectorIntegral.svg"><img src="/assets/images/x4y5Circle_VectorIntegral.svg"></a>
    <figcaption>Figure 4: \(\left[\begin{array}{c}
x^{4}\\
y^{5}
\end{array}\right]\) displaying zero integral on closed loop</figcaption>
</figure>

<figure>
    <a href="/assets/images/ArgEpitrochoid_VectorIntegral.svg"><img src="/assets/images/ArgEpitrochoid_VectorIntegral.svg"></a>
    <figcaption>Figure 5: \(\nabla\tan^{-1}\left(\frac{y}{x}\right)\), the angle of a point, over an epitrochoid</figcaption>
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