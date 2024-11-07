import gmsh
import pdb

# Initialize Gmsh
gmsh.initialize()

gmsh.model.add("circle_in_rectangle");

size = 0.4
refined_size = 0.4
r = 0.2

def add_circle(x,y,r):
     gmsh.model.geo.addPoint(x, y, 0, size, tag=5)
     gmsh.model.geo.addPoint(x+r, y, 0, refined_size, tag=6)
     gmsh.model.geo.addPoint(x-r, y, 0, refined_size, tag=7)
     gmsh.model.geo.addCircleArc(6, 5, 7, tag=5)
     gmsh.model.geo.addCircleArc(7, 5, 6, tag=6)

gmsh.model.geo.addPoint(0, 0, 0, size, tag=1)
gmsh.model.geo.addPoint(2, 0, 0, size, tag=2)
gmsh.model.geo.addPoint(2, 1, 0, size, tag=3)
gmsh.model.geo.addPoint(0, 1, 0, size, tag=4)

bottom = gmsh.model.geo.addLine(1, 2, tag=1)
right = gmsh.model.geo.addLine(2, 3, tag=2)
top = gmsh.model.geo.addLine(3, 4, tag=3)
left = gmsh.model.geo.addLine(4, 1, tag=4)

rectangle_domain = gmsh.model.geo.addCurveLoop([1,2,3,4], tag=1)
#circle = gmsh.model.geo.addCurveLoop([5,6], tag=2)

silicone = gmsh.model.geo.addPlaneSurface([1], tag=1)
#silicone = gmsh.model.geo.addPlaneSurface([1,2], tag=1)
#blob = gmsh.model.geo.addPlaneSurface([2], tag=2)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [bottom], tag=1, name="bottom_boundary")
gmsh.model.addPhysicalGroup(1, [right], tag=2, name="right_boundary")
gmsh.model.addPhysicalGroup(1, [top], tag=3, name="top_boundary")
gmsh.model.addPhysicalGroup(1, [left], tag=4, name="left_boundary")
gmsh.model.addPhysicalGroup(2, [silicone], tag=5, name="silicone")
#gmsh.model.addPhysicalGroup(2, [2], tag=6, name="blob")

gmsh.model.mesh.generate(2)

filename = "mesh.msh"
gmsh.write(filename)

gmsh.finalize()
