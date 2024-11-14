import gmsh
import numpy as np
import pdb

# Initialize Gmsh
gmsh.initialize()

gmsh.model.add("circle_in_rectangle");

size = 0.05
refined_size = 0.005
r = 0.2

dx = 2
dy = 1
inclusion_x = 0.5
inclusion_y = 0.4
inclusion_dx = 1
inclusion_dy = 0.2
pml_thickness = 0.2
transmitter_thickness = 0.05

def add_circle(x,y,r):
     gmsh.model.geo.addPoint(x, y, 0, size)
     gmsh.model.geo.addPoint(x+r, y, 0, refined_size)
     gmsh.model.geo.addPoint(x-r, y, 0, refined_size)
     gmsh.model.geo.addCircleArc(6, 5, 7, tag=5)
     gmsh.model.geo.addCircleArc(7, 5, 6, tag=6)

# points for the computational domain
domain_bl = gmsh.model.geo.addPoint(0 + transmitter_thickness, 0, 0, size)
domain_br = gmsh.model.geo.addPoint(dx, 0, 0, size)
domain_tr = gmsh.model.geo.addPoint(dx, dy, 0, size)
domain_tl = gmsh.model.geo.addPoint(0 + transmitter_thickness, dy, 0, size)

# points for the outer PML layer
pml_bl = gmsh.model.geo.addPoint(0 - pml_thickness, 0 - pml_thickness, 0, size)
pml_br = gmsh.model.geo.addPoint(dx + pml_thickness, 0 - pml_thickness, 0, size)
pml_tr = gmsh.model.geo.addPoint(dx + pml_thickness, dy + pml_thickness, 0, size)
pml_tl = gmsh.model.geo.addPoint(0 - pml_thickness, dy + pml_thickness, 0, size)

# points for the inclusion
inclusion_bl = gmsh.model.geo.addPoint(inclusion_x, inclusion_y, 0, size)
inclusion_br = gmsh.model.geo.addPoint(inclusion_x + inclusion_dx, inclusion_y, 0, size)
inclusion_tr = gmsh.model.geo.addPoint(inclusion_x + inclusion_dx, inclusion_y + inclusion_dy, 0, size)
inclusion_tl = gmsh.model.geo.addPoint(inclusion_x, inclusion_y + inclusion_dy, 0, size)

# points for the transmitter 
transmitter_bl = gmsh.model.geo.addPoint(0, 0, 0, size)
transmitter_br = domain_bl
transmitter_tl = gmsh.model.geo.addPoint(0, dy, 0, size)
transmitter_tr = domain_tl 

# lines for the computational domain
domain_bottom = gmsh.model.geo.addLine(domain_bl, domain_br)
domain_right = gmsh.model.geo.addLine(domain_br, domain_tr)
domain_top = gmsh.model.geo.addLine(domain_tr, domain_tl)
domain_left = gmsh.model.geo.addLine(domain_tl, domain_bl)

# lines for the pml
pml_bottom = gmsh.model.geo.addLine(pml_bl, pml_br)
pml_right = gmsh.model.geo.addLine(pml_br, pml_tr)
pml_top = gmsh.model.geo.addLine(pml_tr, pml_tl)
pml_left = gmsh.model.geo.addLine(pml_tl, pml_bl)

# lines for the computational domain
inclusion_bottom = gmsh.model.geo.addLine(inclusion_bl, inclusion_br)
inclusion_right = gmsh.model.geo.addLine(inclusion_br, inclusion_tr)
inclusion_top = gmsh.model.geo.addLine(inclusion_tr, inclusion_tl)
inclusion_left = gmsh.model.geo.addLine(inclusion_tl, inclusion_bl)

# lines for the transmitter domain
transmitter_bottom = gmsh.model.geo.addLine(transmitter_bl, transmitter_br)
transmitter_right = -domain_left
transmitter_top = gmsh.model.geo.addLine(transmitter_tr, transmitter_tl)
transmitter_left = gmsh.model.geo.addLine(transmitter_tl, transmitter_bl)

# curve loops
pml_loop = gmsh.model.geo.addCurveLoop([pml_bottom, pml_right, pml_top, pml_left])
domain_loop = gmsh.model.geo.addCurveLoop([domain_bottom, domain_right, domain_top, domain_left])
inclusion_loop = gmsh.model.geo.addCurveLoop([inclusion_bottom, inclusion_right, inclusion_top, inclusion_left])
transmitter_loop = gmsh.model.geo.addCurveLoop([transmitter_bottom, transmitter_right, transmitter_top, transmitter_left])

# surfaces
pml = gmsh.model.geo.addPlaneSurface([pml_loop, domain_loop, transmitter_loop])
domain = gmsh.model.geo.addPlaneSurface([domain_loop, inclusion_loop])
inclusion = gmsh.model.geo.addPlaneSurface([inclusion_loop])
transmitter = gmsh.model.geo.addPlaneSurface([transmitter_loop])

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [pml_bottom], name="bottom_boundary")
gmsh.model.addPhysicalGroup(1, [pml_right], name="right_boundary")
gmsh.model.addPhysicalGroup(1, [pml_top], name="top_boundary")
gmsh.model.addPhysicalGroup(1, [pml_left], name="left_boundary")
gmsh.model.addPhysicalGroup(1, [pml_bottom, pml_right, pml_top, pml_left], name="boundary")
gmsh.model.addPhysicalGroup(1, [domain_left], name="left_domain_boundary")
pml_physical = gmsh.model.addPhysicalGroup(2, [pml], name="pml")
domain_physical = gmsh.model.addPhysicalGroup(2, [domain], name="domain")
inclusion_physical = gmsh.model.addPhysicalGroup(2, [inclusion],name="inclusion")
transmitter_physical = gmsh.model.addPhysicalGroup(2, [transmitter],name="transmitter")

gmsh.model.mesh.generate(2)

## Find the elements that touch the left domain edge
#domain_types, domain_element_tags, domain_node_tags = gmsh.model.mesh.getElements(dim=2, tag=domain)
#border_types, border_element_tags, border_node_tags = gmsh.model.mesh.getElements(dim=1, tag=domain_left)
#
#
## Assume you already have these arrays from your Gmsh extraction
#triangle_tags = domain_element_tags[0]  # Element tags of the triangles
#triangle_nodes = domain_node_tags[0]  # Node tags of the triangles (flattened array)
#line_nodes = border_node_tags[0]  # Node tags of the lines (flattened array)
#
## Reshape the triangle nodes to get each triangle's 3 nodes (each row is a triangle)
#triangle_nodes = triangle_nodes.reshape(-1, 3)
#
## Reshape the line nodes to get each line's 2 nodes (each row is a line edge)
#line_edges = line_nodes.reshape(-1, 2)
#
## Create all possible edges of the triangles
#triangle_edges_1 = triangle_nodes[:, [0, 1]]
#triangle_edges_2 = triangle_nodes[:, [1, 2]]
#triangle_edges_3 = triangle_nodes[:, [2, 0]]
#
## Sort the nodes in each edge (to ensure consistent ordering)
#triangle_edges_1 = np.sort(triangle_edges_1, axis=1)
#triangle_edges_2 = np.sort(triangle_edges_2, axis=1)
#triangle_edges_3 = np.sort(triangle_edges_3, axis=1)
#line_edges = np.sort(line_edges, axis=1)
#
## Stack the triangle edges vertically for easy comparison
#all_triangle_edges = np.vstack([triangle_edges_1, triangle_edges_2, triangle_edges_3])
#
## Use broadcasting to check if any line edge matches any triangle edge
#matches = np.isin(all_triangle_edges, line_edges).all(axis=1)
#
## Find indices of matching edges and map back to triangle tags
#matching_triangle_indices = np.where(matches[:len(triangle_tags)])[0]
#matching_triangle_tags = triangle_tags[matching_triangle_indices]
#
#print("Triangles containing lines:", matching_triangle_tags)

filename = "mesh.msh"
gmsh.write(filename)

gmsh.finalize()
