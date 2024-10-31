import gmsh

# Initialize Gmsh
gmsh.initialize()

gmsh.model.add("circle_in_rectangle");

size = 2
refined_size = 0.2
r = 8

gmsh.model.geo.addPoint(0, 0, 0, refined_size, tag=1);
gmsh.model.geo.addPoint(40, 0, 0, size, tag=2);
gmsh.model.geo.addPoint(40, 30, 0, size, tag=3);
gmsh.model.geo.addPoint(0, 30, 0, refined_size, tag=4);

gmsh.model.geo.addPoint(20, 15, 0, size, tag=5);
gmsh.model.geo.addPoint(20+r, 15, 0, refined_size, tag=6);
gmsh.model.geo.addPoint(20-r, 15, 0, refined_size, tag=7);

bottom = gmsh.model.geo.addLine(1, 2, tag=1)
right = gmsh.model.geo.addLine(2, 3, tag=2)
top = gmsh.model.geo.addLine(3, 4, tag=3)
left = gmsh.model.geo.addLine(4, 1, tag=4)

gmsh.model.geo.addCircleArc(6, 5, 7, tag=5)
gmsh.model.geo.addCircleArc(7, 5, 6, tag=6)

gmsh.model.geo.addCurveLoop([1,2,3,4], tag=1)
gmsh.model.geo.addCurveLoop([5,6], tag=2)

gmsh.model.geo.addPlaneSurface([1,2], tag=1)
gmsh.model.geo.addPlaneSurface([2], tag=2)

gmsh.model.addPhysicalGroup(1, [bottom], tag=1, name="bottom boundary")
gmsh.model.addPhysicalGroup(1, [right], tag=2, name="right boundary")
gmsh.model.addPhysicalGroup(1, [top], tag=3, name="top boundary")
gmsh.model.addPhysicalGroup(1, [left], tag=4, name="left boundary")
gmsh.model.addPhysicalGroup(2, [1], tag=5, name="silicone")
gmsh.model.addPhysicalGroup(2, [2], tag=6, name="blob")

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

filename = "mesh.msh"
gmsh.write(filename)

gmsh.finalize()
