library(vegan3d)
library(rgl)
open3d()
#bg3d(color = "black")
ordirgl(best, size=4, ax.col = "black", arr.col = "blue", display = "species", type="n")
orglpoints(best, display = "species", choices = 1:3, radius=0.02, col = cem_lagos$ids$FITC.H)
orglpoints(best, display = "species", choices = 1:3, radius=0.02, col = cem_lagos$ids$SSC.H)
orglpoints(best, display = "species", choices = 1:3, radius=0.02, col = cem_lagos$ids$PerCP.Cy5.5.H)

#orgltext(ord, display = "species")

#Save
#writeWebGL(width=800,height = 800)