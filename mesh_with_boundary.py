import gmsh
import pdb

# Initialize Gmsh
gmsh.initialize()

gmsh.model.add("circle_in_rectangle");

size = 0.1
refined_size = 0.1
r = 0.2

def add_circle(x,y,r):
     gmsh.model.geo.addPoint(x, y, 0, size, tag=5)
     gmsh.model.geo.addPoint(x+r, y, 0, refined_size, tag=6)
     gmsh.model.geo.addPoint(x-r, y, 0, refined_size, tag=7)
     gmsh.model.geo.addCircleArc(6, 5, 7, tag=5)
     gmsh.model.geo.addCircleArc(7, 5, 6, tag=6)

rectangle_bl = gmsh.model.geo.addPoint(0, 0, 0, size)
rectangle_br = gmsh.model.geo.addPoint(2, 0, 0, size)
rectangle_tr = gmsh.model.geo.addPoint(2, 1, 0, size)
rectangle_tl = gmsh.model.geo.addPoint(0, 1, 0, size)

reciever_t = gmsh.model.geo.addPoint(-0.1, 1, 0, size)
reciever_b = gmsh.model.geo.addPoint(-0.1, 0, 0, size )

boundary_bl = gmsh.model.geo.addPoint(-2, -2, 0, size)
boundary_tl = gmsh.model.geo.addPoint(-2, 3, 0, size)
boundary_tr = gmsh.model.geo.addPoint(4, 3, 0, size)
boundary_br = gmsh.model.geo.addPoint(4, -2, 0, size)

bottom = gmsh.model.geo.addLine(rectangle_bl, rectangle_br)
right = gmsh.model.geo.addLine(rectangle_br, rectangle_tr)
top = gmsh.model.geo.addLine(rectangle_tr, rectangle_tl)
left = gmsh.model.geo.addLine(rectangle_tl, rectangle_bl)

reciever_left = gmsh.model.geo.addLine(reciever_t, reciever_b)
reciever_bottom = gmsh.model.geo.addLine(reciever_b, rectangle_bl)
reciever_top = gmsh.model.geo.addLine(rectangle_tl, reciever_t)

bottom_b = gmsh.model.geo.addLine(boundary_bl, boundary_br)
right_b = gmsh.model.geo.addLine(boundary_br,boundary_tr)
top_b = gmsh.model.geo.addLine(boundary_tr,boundary_tl)
left_b = gmsh.model.geo.addLine(boundary_tl,boundary_bl)

rectangle_domain = gmsh.model.geo.addCurveLoop([bottom,right,top,left])
reciever = gmsh.model.geo.addCurveLoop([reciever_left, reciever_bottom, -left, reciever_top])
outer_rectangle = gmsh.model.geo.addCurveLoop([bottom_b, right_b, top_b, left_b])
reciever_and_rectangle = gmsh.model.geo.addCurveLoop([reciever_left, reciever_bottom, bottom, right, top, reciever_top])
#circle = gmsh.model.geo.addCurveLoop([5,6], tag=2)

outer_layer = gmsh.model.geo.addPlaneSurface([outer_rectangle, reciever_and_rectangle])
silicone = gmsh.model.geo.addPlaneSurface([rectangle_domain])
reciever_surface = gmsh.model.geo.addPlaneSurface([reciever])
#silicone = gmsh.model.geo.addPlaneSurface([1,2], tag=1)
#blob = gmsh.model.geo.addPlaneSurface([2], tag=2)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [bottom_b], name="bottom_boundary")
gmsh.model.addPhysicalGroup(1, [right_b], name="right_boundary")
gmsh.model.addPhysicalGroup(1, [top_b], name="top_boundary")
gmsh.model.addPhysicalGroup(1, [left_b], name="left_boundary")

gmsh.model.addPhysicalGroup(1, [bottom], name="bottom")
gmsh.model.addPhysicalGroup(1, [right], name="right")
gmsh.model.addPhysicalGroup(1, [top], name="top")
gmsh.model.addPhysicalGroup(1, [left], name="left")

gmsh.model.addPhysicalGroup(2, [silicone], name="silicone")
gmsh.model.addPhysicalGroup(2, [outer_layer], name="outer_layer")
gmsh.model.addPhysicalGroup(2, [reciever_surface], name="reciever")
#gmsh.model.addPhysicalGroup(2, [2], tag=6, name="blob")

gmsh.model.mesh.generate(2)

filename = "mesh.msh"
gmsh.write(filename)

gmsh.finalize()
