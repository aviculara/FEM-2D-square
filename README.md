# FEM 2D Unit Square Analysis
Analyse a 2D unit square with Finite Element Method through Matlab code. 
The user can select from five different element types: D2TR3N, D2TR6N, D2QU4N, D2QU8N, D2QU9N. 
The unit square can include a void or inclusion of three different shapes: circle, square, or rhombus.
The analysis can be performed for Extension, Expansion, or Shear type deformation.

## How to use
- Run the file with .mlapp extension.
- Select the inclusion shape, and select between Void or Inclusion.
- Select the Element Type from the drop-down menu.
- Edit the dimensions of the unit square as necessary.
- Click "Generate". The generated shape can be viewed from the graph above. The generated Number of Elements and Nodes can be viewed on the left.
- Click "Reset" before generating a different shape.
- Enter the Young's Modulus and Poisson's Ratio of the structure, and the inclusion if applicable.
- Select the deformation type, and enter the deformation magnitude.
- Click "Run" and wait for the indicator to turn green.
- Move to the "Post-Process" tab to view the results.
- Select the Stress, Strain, and Displacement values to be displayed from the drop-down menu.
- You may change the Deformation Scale for clarity.
