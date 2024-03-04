# coding:utf-8
# Statement：This script works for Abaqus2016, other versions may have some unknown errors
import random
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#The model is in millimeter units.
# Shot radius:mm
R=0.1
# Shot velocity:mm/s
V=65000
# Target material model
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(1.5, 1.5), point2=(-1.5, -1.5))
p = mdb.models['Model-1'].Part(name='Part-bacai', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-bacai']
p.BaseSolidExtrude(sketch=s, depth=1.5)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-bacai']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[5], sketchUpEdge=e[6], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=8.48, gridSpacing=0.21, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.rectangle(point1=(0.6, 0.6), point2=(-0.6, -0.6))
p = mdb.models['Model-1'].parts['Part-bacai']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#20 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e1[6], faces=pickedFaces, sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
f1, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[6], sketchUpEdge=e[4], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=8.48, 
    gridSpacing=0.21, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.rectangle(point1=(0.3, 0.3), point2=(-0.3, -0.3))
p = mdb.models['Model-1'].parts['Part-bacai']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#40 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e1[4], faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=e[10], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=3.39, gridSpacing=0.08, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
f1, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[1], sketchUpEdge=e1[10], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=8.48, 
    gridSpacing=0.21, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.Line(point1=(-1.5, 1.5), point2=(-0.3, 0.3))
s.Line(point1=(-1.5, -1.5), point2=(-0.3, -0.3))
s.Line(point1=(1.5, 1.5), point2=(0.6, 0.6))
s.Line(point1=(1.5, -1.5), point2=(0.6, -0.6))
p = mdb.models['Model-1'].parts['Part-bacai']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#2 ]', ), )
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e[10], faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
del p.features['Partition face-3']
p = mdb.models['Model-1'].parts['Part-bacai']
f, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[1], sketchUpEdge=e1[10], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=8.48, gridSpacing=0.21, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.Line(point1=(-1.5, 1.5), point2=(-0.3, 0.3))
p = mdb.models['Model-1'].parts['Part-bacai']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#2 ]', ), )
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e[10], faces=pickedFaces, sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
del p.features['Partition face-3']
p = mdb.models['Model-1'].parts['Part-bacai']
f1, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[0], sketchUpEdge=e1[10], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=3.39, 
    gridSpacing=0.08, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-bacai']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.Line(point1=(-1.5, 1.5), point2=(-0.3, 0.3))
s.Line(point1=(1.5, 1.5), point2=(0.3, 0.3))
s.Line(point1=(-1.5, -1.5), point2=(-0.3, -0.3))
s.Line(point1=(0.3, -0.3), point2=(1.5, -1.5))
p = mdb.models['Model-1'].parts['Part-bacai']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#3 ]', ), )
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e[10], faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
pickedEdges =(e1[7], e1[15])
p.PartitionCellByExtrudeEdge(line=e1[22], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e, d1 = p.edges, p.datums
pickedEdges =(e[10], e[16])
p.PartitionCellByExtrudeEdge(line=e[24], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
pickedEdges =(e1[12], e1[18])
p.PartitionCellByExtrudeEdge(line=e1[28], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e, d1 = p.edges, p.datums
pickedEdges =(e[15], e[23])
p.PartitionCellByExtrudeEdge(line=e[32], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
pickedEdges =(e1[20], e1[22], e1[24], e1[28])
p.PartitionCellByExtrudeEdge(line=e1[11], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#10 ]', ), )
e, d1 = p.edges, p.datums
pickedEdges =(e[42], e[43], e[44], e[45])
p.PartitionCellByExtrudeEdge(line=e[35], cells=pickedCells, edges=pickedEdges, 
    sense=FORWARD)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
# Shot model
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s1.FixedConstraint(entity=g[2])
s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, R))
s1.Line(point1=(0.0, R), point2=(0.0, -R))
s1.VerticalConstraint(entity=g[4], addUndoState=False)
s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
s1.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
s1.autoTrimCurve(curve1=g[3], point1=(-0.242308974266052, -0.215989306569099))
p = mdb.models['Model-1'].Part(name='Part-ball', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-ball']
p.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-ball']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-ball']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
v1, e1, d2 = p.vertices, p.edges, p.datums
p.PartitionCellByPlaneThreePoints(point1=v1[0], point3=v1[1], 
    cells=pickedCells, point2=p.InterestingPoint(edge=e1[0], rule=MIDDLE))
p = mdb.models['Model-1'].parts['Part-ball']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#3 ]', ), )
e, v2, d1 = p.edges, p.vertices, p.datums
p.PartitionCellByPlanePointNormal(normal=d1[1], cells=pickedCells, 
    point=p.InterestingPoint(edge=e[1], rule=MIDDLE))
p = mdb.models['Model-1'].parts['Part-ball']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#f ]', ), )
e1, v1, d2 = p.edges, p.vertices, p.datums
p.PartitionCellByPlaneNormalToEdge(edge=e1[1], cells=pickedCells, 
    point=p.InterestingPoint(edge=e1[1], rule=MIDDLE))
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
# Material
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
mdb.models['Model-1'].Material(name='matBall')
mdb.models['Model-1'].materials['matBall'].Density(table=((7.85e-09, ), ))
mdb.models['Model-1'].materials['matBall'].Elastic(table=((210000.0, 0.3), ))
mdb.models['Model-1'].materials['matBall'].Plastic(table=((1550.0, 0.0), (
    1550.0, 100.0)))
mdb.models['Model-1'].Material(name='matBACAI')
mdb.models['Model-1'].materials['matBACAI'].Density(table=((7.8e-09, ), ))
mdb.models['Model-1'].materials['matBACAI'].Depvar(n=15)
mdb.models['Model-1'].materials['matBACAI'].UserMaterial(mechanicalConstants=(
    210000.0, 0.3, 792.0, 0.154, 0.078, 18.6, 32.8, 89.8, 90.3, 60.8, 
    2.48e-07, 3.06, 0.25, 82000.0, 0.25, 0.06, 10000000.0, 100.0, 1.0, 3.2, 
    0.26, 25000000.0, 50000000.0, 0.01789, 31250000.0))
mdb.models['Model-1'].HomogeneousSolidSection(name='ball', material='matBall', 
    thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='bacai', 
    material='matBACAI', thickness=None)
p = mdb.models['Model-1'].parts['Part-ball']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#ff ]', ), )
region = p.Set(cells=cells, name='Set-matball')
p = mdb.models['Model-1'].parts['Part-ball']
p.SectionAssignment(region=region, sectionName='ball', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-bacai']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1ff ]', ), )
region = p.Set(cells=cells, name='Set-matbacai')
p = mdb.models['Model-1'].parts['Part-bacai']
p.SectionAssignment(region=region, sectionName='bacai', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
session.viewports['Viewport: 1'].view.setValues(width=13.4764, height=6.53993, 
    viewOffsetX=2.21357, viewOffsetY=-0.582657)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-bacai']
# Assembly of multiple shots
a.Instance(name='Part-bacai-1', part=p, dependent=ON)
n=1         # Coverage
j=10*n      # Total shot number
i=1         # Number of shots in the assembly
while  i<j+1:
				n1 = random.uniform(-0.3,0.3)# Distribution in the X-direction 
				n2 = random.uniform(-0.3,0.3)# Distribution in the Y-direction
				n3 = -i*R                # Distribution in the Z-direction
				p = mdb.models['Model-1'].parts['Part-ball']
				a.Instance(name='Part-ball-%d' % i, part=p, dependent=ON)   
				a = mdb.models['Model-1'].rootAssembly
				a.translate(instanceList=('Part-ball-%d' % i, ), vector=(n1, n2, n3))    
				i +=1
session.viewports['Viewport: 1'].view.setValues(session.views['Bottom'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
# Numerical simulation time
ts = (j+1)*R/V
mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial',
    timePeriod=ts)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'PE', 'PEEQ', 'LE', 'U', 'V', 'A', 'E',
    'SDV'))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
mdb.models['Model-1'].ContactProperty('IntProp-contact')
mdb.models['Model-1'].interactionProperties['IntProp-contact'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
    table=((0.2, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
    fraction=0.005, elasticSlipStiffness=None)
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-bacai-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#ff000000 #20 ]', ), )
region2=a.Surface(side1Faces=side1Faces1, name='s_Surf-jiechu')
# Define contact
i=1
a = mdb.models['Model-1'].rootAssembly
while  i<j+1:
	s1 = a.instances['Part-ball-%d' % i].faces
	side1Faces1 = s1.getSequenceFromMask(mask=('[#cc330 ]', ), )
	region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-ball-%d' % i)
	mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='Int-%d' % i, 
		createStepName='Initial', master = region1, slave = region2, 
		mechanicalConstraint=PENALTY, sliding=FINITE, 
		interactionProperty='IntProp-contact', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)
	i +=1
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=101.015, 
farPlane=170.715, width=6.43051, height=3.12066, viewOffsetX=0.987603, 
viewOffsetY=0.0111254)
session.viewports['Viewport: 1'].view.setValues(nearPlane=100.898, 
farPlane=170.833, width=7.73317, height=3.75282, viewOffsetX=1.26171, 
viewOffsetY=-0.112985)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-bacai-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#0 #80034 ]', ), )
region = a.Set(edges=edges1, name='Set-dimian')
mdb.models['Model-1'].DisplacementBC(name='BC-dimian', createStepName='Initial', 
region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
bcs=OFF, predefinedFields=OFF, interactions=OFF, constraints=OFF, 
connectors=OFF, engineeringFeatures=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
meshTechnique=ON)
p = mdb.models['Model-1'].parts['Part-bacai']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
engineeringFeatures=OFF)
p = mdb.models['Model-1'].parts['Part-ball']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.70337, 
farPlane=2.74796, width=1.37325, height=0.668729, viewOffsetX=0.219957, 
viewOffsetY=-0.0882102)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.65552, 
farPlane=2.79582, width=1.93467, height=0.942125, viewOffsetX=0.314581, 
viewOffsetY=-0.118708)
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
hourglassControl=DEFAULT, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
p = mdb.models['Model-1'].parts['Part-ball']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#ff ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
elemType3))
p = mdb.models['Model-1'].parts['Part-ball']
p.seedPart(size=0.02, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-ball']
p.generateMesh()
p = mdb.models['Model-1'].parts['Part-bacai']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.78021, 
farPlane=11.2198, width=5.16868, height=2.51699, viewOffsetX=0.66688, 
viewOffsetY=-0.267957)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.31257, 
farPlane=11.534, width=4.8122, height=2.34339, cameraPosition=(1.3179, 
7.15287, -4.42236), cameraUpVector=(0.105441, -0.580596, -0.807335), 
cameraTarget=(0.0512983, -0.158748, 0.67037), viewOffsetX=0.620885, 
viewOffsetY=-0.249476)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.40102, 
farPlane=11.2401, width=4.87962, height=2.37622, cameraPosition=(
0.677186, 6.34235, -5.34605), cameraUpVector=(0.129244, -0.677647, 
-0.723941), cameraTarget=(0.0609976, -0.267437, 0.731037), 
viewOffsetX=0.629584, viewOffsetY=-0.252971)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.28133, 
farPlane=11.0789, width=4.78838, height=2.33179, cameraPosition=(
-0.541452, 5.82348, -5.66792), cameraUpVector=(0.147849, -0.708463, 
-0.690088), cameraTarget=(0.094855, -0.372132, 0.828962), 
viewOffsetX=0.617812, viewOffsetY=-0.248241)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.20608, 
farPlane=11.1542, width=5.69601, height=2.77378, viewOffsetX=0.64759, 
viewOffsetY=-0.233121)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#0 #1000 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#80000000 #4400 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
end2Edges=pickedEdges2, minSize=0.04, maxSize=0.08, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#0 #1 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#10000000 #180 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
end2Edges=pickedEdges2, minSize=0.02, maxSize=0.03, constraint=FINER)
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#100 ]', ), )
leaf = dgm.LeafFromGeometry(cellSeq=cells)
session.viewports['Viewport: 1'].partDisplay.displayGroup.replace(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#fff ]', ), )
p.seedEdgeBySize(edges=pickedEdges, size=0.02, deviationFactor=0.1, 
constraint=FINER)
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
session.viewports['Viewport: 1'].partDisplay.displayGroup.replace(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
p.seedPart(size=0.04, deviationFactor=0.1, minSizeFactor=0.1)
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
hourglassControl=DEFAULT, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1ff ]', ), )
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
elemType3))
p = mdb.models['Model-1'].parts['Part-bacai']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.28506, 
farPlane=12.1998, width=5.7685, height=2.80908, cameraPosition=(
3.85092, 8.11152, -1.45158), cameraUpVector=(0.00765493, -0.246605, 
-0.969086), cameraTarget=(-0.035026, 0.237118, 0.521517), 
viewOffsetX=0.655831, viewOffsetY=-0.236088)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.31278, 
farPlane=11.4625, width=5.79394, height=2.82147, cameraPosition=(
0.0434989, 7.97035, -3.19286), cameraUpVector=(-0.0800861, -0.424527, 
-0.901866), cameraTarget=(-0.253113, -0.158005, 0.65964), 
viewOffsetX=0.658723, viewOffsetY=-0.237129)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.2989, 
farPlane=11.5443, width=5.7812, height=2.81526, cameraPosition=(1.0657, 
5.86346, -5.90542), cameraUpVector=(0.233918, -0.717561, -0.65604), 
cameraTarget=(-0.128806, -0.364653, 0.48079), viewOffsetX=0.657274, 
viewOffsetY=-0.236608)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.69585, 
farPlane=11.1473, width=2.92478, height=1.42428, viewOffsetX=0.400464, 
viewOffsetY=-0.139545)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.00052, 
farPlane=10.8528, width=3.05786, height=1.48909, cameraPosition=(
0.369803, 4.20066, -7.1334), cameraUpVector=(0.295146, -0.825574, 
-0.480953), cameraTarget=(-0.106271, -0.449883, 0.55722), 
viewOffsetX=0.418686, viewOffsetY=-0.145894)
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#100 ]', ), )
leaf = dgm.LeafFromGeometry(cellSeq=cells)
session.viewports['Viewport: 1'].partDisplay.displayGroup.replace(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.94721, 
farPlane=9.90606, width=1.07135, height=0.521717, viewOffsetX=0.108736, 
viewOffsetY=-0.0953227)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.95687, 
farPlane=9.89641, width=1.07266, height=0.52235, viewOffsetX=0.226814, 
viewOffsetY=-0.263038)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.95687, 
farPlane=9.8964, viewOffsetX=0.354047, viewOffsetY=-0.437155)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.9911, 
farPlane=9.86217, width=0.790616, height=0.385006, 
viewOffsetX=0.287194, viewOffsetY=-0.417761)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.96966, 
farPlane=3.92382, width=1.13415, height=0.552296, 
viewOffsetX=-0.0421134, viewOffsetY=0.133089)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.59947, 
farPlane=4.57613, width=0.992772, height=0.483449, cameraPosition=(
-0.221608, 1.66541, -2.42156), cameraUpVector=(0.195673, -0.874929, 
-0.442957), cameraTarget=(-0.0312465, 0.144913, 0.665811), 
viewOffsetX=-0.0368637, viewOffsetY=0.116498)
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
session.viewports['Viewport: 1'].partDisplay.displayGroup.replace(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.73519, 
farPlane=5.44041, width=3.11255, height=1.51572, viewOffsetX=0.218931, 
viewOffsetY=-0.198109)
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.32646, 
farPlane=10.6735, width=5.74446, height=2.79738, viewOffsetX=1.11905, 
viewOffsetY=-0.413872)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.48218, 
farPlane=10.5178, width=5.51457, height=2.68543, viewOffsetX=0.964942, 
viewOffsetY=-0.359192)
p = mdb.models['Model-1'].parts['Part-ball']
# Shot boundary conditions
i=1
a = mdb.models['Model-1'].rootAssembly
while  i<j+1:
    c1 = a.instances['Part-ball-%d' % i].cells
    cells1 = c1.getSequenceFromMask(mask=('[#ff ]', ), )
    f1 = a.instances['Part-ball-%d' % i].faces
    faces1 = f1.getSequenceFromMask(mask=('[#fffff ]', ), )
    e1 = a.instances['Part-ball-%d' % i].edges
    edges1 = e1.getSequenceFromMask(mask=('[#3ffff ]', ), )
    v1 = a.instances['Part-ball-%d' % i].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#7f ]', ), )
    region = a.Set(vertices=verts1, edges=edges1, faces=faces1, cells=cells1, name='Set-%d' % i)
    mdb.models['Model-1'].DisplacementBC(name='BC-%d' % i, 
        createStepName='Initial', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    i +=1
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
session.viewports['Viewport: 1'].partDisplay.displayGroup.replace(leaf=leaf)
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=ON, 
bcs=ON, predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
meshTechnique=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=105.659, 
farPlane=164.329, width=8.09805, height=3.92989, cameraPosition=(
-70.1068, 43.9916, 74.2455), cameraUpVector=(0.177948, 0.946159, 
-0.270401), cameraTarget=(0.30472, -0.00680161, -33.3721), 
viewOffsetX=1.32124, viewOffsetY=-0.118316)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF, adaptiveMeshConstraints=ON)
mdb.models['Model-1'].steps['Step-1'].Restart(numberIntervals=20, overlay=ON, 
    timeMarks=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
# Pellet velocity setting 
# Direction: Z-axis
i=1
a = mdb.models['Model-1'].rootAssembly
while  i<j+1:
    c1 = a.instances['Part-ball-%d' % i].cells
    cells1 = c1.getSequenceFromMask(mask=('[#ff ]', ), )
    f1 = a.instances['Part-ball-%d' % i].faces
    faces1 = f1.getSequenceFromMask(mask=('[#fffff ]', ), )
    e1 = a.instances['Part-ball-%d' % i].edges
    edges1 = e1.getSequenceFromMask(mask=('[#3ffff ]', ), )
    v1 = a.instances['Part-ball-%d' % i].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#7f ]', ), )
    region = a.Set(vertices=verts1, edges=edges1, faces=faces1, cells=cells1, 
    name='v-%d' % i)
    mdb.models['Model-1'].Velocity(name='Predefined Field-%d' % i, region=region, 
        field='', distributionType=MAGNITUDE, velocity1=0.0, velocity2=0.0, 
        velocity3=V, omega=0.0)
    i +=1
# Create element set
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16326 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d0')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16326 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16377 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d20')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16377 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16429 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d40')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16429 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16480 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d60')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16480 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16531 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d80')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16531 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16582 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d100')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16582 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16634 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d120')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16634 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16685 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d140')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16685 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16736 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d160')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16736 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16787 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d180')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16787 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16839 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d200')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16839 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16890 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d220')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16890 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16941 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d240')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:16941 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16992 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d260')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:16992 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17044 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d280')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17044 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17095 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d300')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17095 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17146 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d320')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17146 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17197 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d340')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17197 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17249 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d360')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17249 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17300 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d380')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17300 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17351 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d400')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17351 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17402 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d420')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17402 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17454 #fffffffc #ffffffff:50 #3ff ]', ), )
p.Set(elements=elements, name='d440')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17454 #fffffffc #ffffffff:50 #3ff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17505 #fffffc00 #ffffffff:50 #3ffff ]', ), )
p.Set(elements=elements, name='d460')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17505 #fffffc00 #ffffffff:50 #3ffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17556 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
p.Set(elements=elements, name='d480')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=(
'[#0:17556 #fffc0000 #ffffffff:50 #3ffffff ]', ), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17607 #fc000000 #ffffffff:51 #3 ]', 
), )
p.Set(elements=elements, name='d500')
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
p = mdb.models['Model-1'].parts['Part-bacai']
e = p.elements
elements = e.getSequenceFromMask(mask=('[#0:17607 #fc000000 #ffffffff:51 #3 ]', 
), )
leaf = dgm.LeafFromMeshElementLabels(elementSeq=elements)
session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.03573, 
farPlane=7.79837, width=2.01339, height=0.980457, cameraPosition=(
0.203365, -6.44238, -2.02329), cameraUpVector=(0.144456, 0.341803, 
-0.928603), cameraTarget=(1.42824, 1.8616, 1.22383), 
viewOffsetX=0.986628, viewOffsetY=-0.132298)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.11557, 
farPlane=7.60412, width=2.04002, height=0.993426, cameraPosition=(
-0.0176458, -6.85002, -0.477915), cameraUpVector=(0.131184, 0.118366, 
-0.984266), cameraTarget=(1.48782, 1.93371, 0.779065), 
viewOffsetX=0.999679, viewOffsetY=-0.134048)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.11184, 
farPlane=7.60785, width=2.03878, height=0.992822, viewOffsetX=0.974357, 
viewOffsetY=-0.411815)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.11183, 
farPlane=7.60786, width=2.03878, height=0.992821, viewOffsetX=0.960235, 
viewOffsetY=-0.569321)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.11183, 
farPlane=7.60786, width=2.03878, height=0.992822, viewOffsetX=0.972592, 
viewOffsetY=-0.631262)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.07334, 
farPlane=7.64636, width=2.43918, height=1.1878, viewOffsetX=0.972013, 
viewOffsetY=-0.653651)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.05286, 
farPlane=7.66684, width=2.43095, height=1.1838, viewOffsetX=1.14343, 
viewOffsetY=-0.571261)