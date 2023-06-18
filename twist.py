import math
import numpy as np
import pymesh

ptDensity = 100

# parametric curve
# fn: float -> [float, float, float]
def simple(uu):
	return [uu + 0, uu + 1, uu + 2]

# fn: float -> [float, float, float]
def fn0(uu):
	xx = uu
	yy = math.sin(uu * math.pi) * math.cos(uu / 2 * math.pi) / (uu + 0.1)
	zz = 0
	return [xx, yy, zz]

class Affine:
	# euler: [float, float, float]
	def Rotate(euler):
		rot  = [
			[1, 0, 0],
			[0, 1, 0],
			[0, 0, 1],
		]
		if euler[0] != 0:
			ss = math.sin(math.radians(euler[0]))
			cc = math.cos(math.radians(euler[0]))
			rr = [
				[1,  0,  0],
				[0, cc, -ss],
				[0, ss,  cc],
			]
			rot = np.matmul(rot, rr)
		if euler[1] != 0:
			ss = math.sin(math.radians(euler[1]))
			cc = math.cos(math.radians(euler[1]))
			rr = [
				[cc,  0, ss],
				[ 0,  1,  0],
				[-ss, 0, cc],
			]
			rot = np.matmul(rot, rr)
		if euler[2] != 0:
			ss = math.sin(math.radians(euler[2]))
			cc = math.cos(math.radians(euler[2]))
			rr = [
				[cc, -ss, 0],
				[ss,  cc, 0],
				[ 0,   0, 1],
			]
			rot = np.matmul(rot, rr)

		return rot

	# scale: [float, float, float]
	def Scale(scale):
		return [
			[scale[0], 0, 0],
			[0, scale[1], 0],
			[0, 0, scale[2]],
		]

# A collection of points and edges describing a path through 3space,
# created from a function operating on a closed linear space [0, 1]
# 	fn: (float) -> [float, float, float]
class Wire:
	def __init__(self):
		self.ps = []		# points in the set
		self.ln = []		# pairs of connected points
		self.mesh = None	# surface around the wire

	def forFn(self, fn):
		global ptDensity
		self.ps = [fn(uu) for uu in np.linspace(0, 1, ptDensity)]
		self.ln = [ [ii, ii+1] for ii in range(len(self.ps) - 1) ]
		return self

	def forPts(self, ps):
		self.ps = ps
		self.ln = [ [ii, ii+1] for ii in range(len(self.ps) - 1) ]
		return self

	def Remesh(self, wireThickness):
		global wireDensity

		ww = pymesh.wires.WireNetwork.create_from_data(self.ps, self.ln)
		inflator = pymesh.wires.Inflator(ww)
		inflator.set_profile(wireDensity)
		inflator.inflate(wireThickness, per_vertex_thickness=True)
		self.mesh = inflator.mesh

	def Transform(self, mtx):
		self.ps = [np.matmul(mtx, pt) for pt in self.ps]

	# use a function to generate transform matricies for each point
	def TransformPerPoint(self, fn):
		self.ps = [np.matmul(fn(pt), pt) for pt in self.ps]
		
	def __str__(self):
		ss = "Wire:"
		for ii in range(len(self.ps)):
			ss = ss + f"\t[ {self.ps[ii][0]}, {self.ps[ii][1]}, {self.ps[ii][2]} ]\n"
		ss = ss + "\n"
		for ii in range(len(self.ln)):
			ss = ss + f"\t[ {self.ln[ii][0]}, {self.ln[ii][1]} ]\n"
		return ss

def Twist(theta):
	# fn: [float, float, float] -> [float, float, float]
	def fn(pt):
		return Affine.Rotate([theta + pt[0] * 360, 0, 0])

	return fn

# a gmsh unit is 1mm in meat space
# we want 2 units to be 40mm, so rescale everythng by 20
# rescale = 20
rescale = 1
scale = Affine.Scale([3 * rescale, 1 * rescale, 1 * rescale])

wCnt = 5
theta = 0
wires = []
for ii in range(wCnt):
	wire = Wire().forFn(fn0)
	wire.Transform(scale)
	wire.TransformPerPoint(Twist(theta))
	wires.append(wire)
	theta += 360 / wCnt

# add wires for the base

wireDensity = 10 	# number of points used for a wire cross section
wireThickness = 0.1 * rescale
for wire in wires:
	wire.Remesh(wireThickness)

m1 = pymesh.merge_meshes([wire.mesh for wire in wires])

# Add cyliders, one sold, one a hole, to form a base
p0 = [-0.15 * rescale, 0, 0]
p1 = [ 0.10 * rescale, 0, 0]
r0 = 1.5 * rescale
r1 = 1.5 * rescale
m2 = pymesh.generate_cylinder(p0, p1, r0, r1, 180)

p0 = [-1.00 * rescale, 0, 0]
p1 = [ 1.00 * rescale, 0, 0]
r0 = 1.0 * rescale
r1 = 1.0 * rescale
m3 = pymesh.generate_cylinder(p0, p1, r0, r1, 180)

m4 = pymesh.boolean(m1, m3, 'difference')
m5 = pymesh.boolean(m2, m3, 'difference')
m6 = pymesh.merge_meshes([m4, m5])

pymesh.save_mesh("tmp.stl", m6, ascii=True)

