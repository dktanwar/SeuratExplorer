# Status Update on TODO Issues:

# Overall 
- In interactive plots, added all possible selection options like lasso, pan, box, etc
-- problem: I did not wanted them to be added to right panel. I think plotly plots have these on top of plot by default. Box select is not working. It should zoom into area I select

# Dimensional Reduction plot 
- Interactive functionality works correctly with proper drag modes and controls
- Add 3D as well as an addition option on right of interactive. Should be interactive. Should be none by default and user should be able to select if they have computed 3D UMAP

# Violin plot  
- interactive plots height not adjustable

# dot plot 
- b/w theme not displaying box around plot in interactive mode

# Heatmap cell level 
- the interactuve should be plotly heatmap not interactivecomplexheatmap. it is not working

# Heatmap Group average 
- Interactive averaged heatmap error: Number of layers provided does not match number of assays
- static panel is very very small, nothing visible
- unused argument (color_palette = "viridis")

# Ridge plot
- No theme selection

# Search Features
- Feature transfer works for all plot types, except heatmap average
