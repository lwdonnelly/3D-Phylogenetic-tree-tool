#!/usr/bin/env python

# Author: Kwasi Mensah (kmensah@andrew.cmu.edu)
# Date: 8/05/2005
#
# This demo shows how to make quasi-fractal trees in Panda.
# Its primarily meant to be a more complex example on how to use
# Panda's Geom interface.
import sys

from direct.gui.DirectLabel import DirectLabel, BillboardEffect
from direct.gui.DirectOptionMenu import DirectOptionMenu
from direct.showbase.ShowBase import ShowBase, TextFont
from panda3d.core import Filename, InternalName, ConfigVariableString, PandaNode, CollisionNode, CollisionRay, \
    CollisionHandlerEvent, CollisionTraverser, CollisionHandlerQueue, GeoMipTerrain, CardMaker
from panda3d.core import GeomVertexArrayFormat, GeomVertexFormat
from panda3d.core import Geom, GeomNode, GeomTrifans, GeomTristrips
from panda3d.core import GeomVertexReader, GeomVertexWriter
from panda3d.core import GeomVertexRewriter, GeomVertexData
from panda3d.core import PerspectiveLens, TextNode
from panda3d.core import TransformState, CullFaceAttrib
from panda3d.core import Light, AmbientLight, Spotlight
from panda3d.core import NodePath
from panda3d.core import LVector3, LMatrix4
from panda3d.core import ModifierButtons
from direct.task.Task import Task
from direct.gui.OnscreenText import OnscreenText
from direct.showbase.DirectObject import DirectObject
import math
import random
from panda3d.core import LPoint3f
from statistics import mean
from Bio import Phylo

random.seed()
base = ShowBase()
base.disableMouse()
base.camera.setPos(0, -180, 30)
numPrimitives = 0
#phyloTree = Phylo.read('tree-of-life.xml', 'phyloxml')
phyloTree = Phylo.read('reptile-tree.xml', 'phyloxml')
#phyloTree = Phylo.read('apaf.xml', 'phyloxml')
depth = max(phyloTree.depths(unit_branch_lengths=True).values())
#calculate modifier for a leaf node thickness of .025
girthModifier = pow((1.0 + (depth * 40.0)), (1.0/depth))

phyloTree.ladderize()
Phylo.draw_ascii(phyloTree)
selectedCladeName = None
title = OnscreenText(text="Click on a branch for clade name",
                     style=1, fg=(1, 1, 1, 1), parent=base.a2dBottomCenter,
                     pos=(0, 0.1), scale=.08, mayChange=True)
qEvent = OnscreenText(
    text="Q: Start Scene Over",
    parent=base.a2dTopLeft, align=TextNode.ALeft,
    style=1, fg=(1, 1, 1, 1), pos=(0.06, -0.10),
    scale=.05)
wEvent = OnscreenText(
    text="T: Add Another Tree With Selected Clade as Root",
    parent=base.a2dTopLeft, align=TextNode.ALeft,
    style=1, fg=(1, 1, 1, 1), pos=(0.06, -0.16),
    scale=.05)


# this computes the new Axis which we'll make a branch grow alowng when we
# split
def randomAxis(vecList):
    fwd = vecList[0]
    perp1 = vecList[1]
    perp2 = vecList[2]

    nfwd = fwd + perp1 * (2 * random.random() - 1) + \
        perp2 * (2 * random.random() - 1)
    nfwd.normalize()

    nperp2 = nfwd.cross(perp1)
    nperp2.normalize()

    nperp1 = nfwd.cross(nperp2)
    nperp1.normalize()

    return [nfwd, nperp1, nperp2]


# this makes smalle variations in direction when we are growing a branch
# but not splitting
def smallRandomAxis(vecList):
    fwd = vecList[0]
    perp1 = vecList[1]
    perp2 = vecList[2]

    nfwd = fwd + perp1 * (1 * random.random() - 0.5) + \
        perp2 * (1 * random.random() - 0.5)
    nfwd.normalize()

    nperp2 = nfwd.cross(perp1)
    nperp2.normalize()

    nperp1 = nfwd.cross(nperp2)
    nperp1.normalize()

    return [nfwd, nperp1, nperp2]


# this draws the body of the tree. This draws a ring of vertices and connects the rings with
# triangles to form the body.
# this keepDrawing paramter tells the function wheter or not we're at an end
# if the vertices before you were an end, dont draw branches to it
def drawBody(nodePath, vdata, cladeName, pos, vecList, radius=1, keepDrawing=True, numVertices=8):

    circleGeom = Geom(vdata)

    vertWriter = GeomVertexWriter(vdata, "vertex")
    colorWriter = GeomVertexWriter(vdata, "color")
    normalWriter = GeomVertexWriter(vdata, "normal")
    drawReWriter = GeomVertexRewriter(vdata, "drawFlag")
    texReWriter = GeomVertexRewriter(vdata, "texcoord")

    startRow = vdata.getNumRows()
    vertWriter.setRow(startRow)
    colorWriter.setRow(startRow)
    normalWriter.setRow(startRow)

    sCoord = 0

    if (startRow != 0):
        texReWriter.setRow(startRow - numVertices)
        sCoord = texReWriter.getData2f().getX() + 1

        drawReWriter.setRow(startRow - numVertices)
        if(drawReWriter.getData1f() == False):
            sCoord -= 1

    drawReWriter.setRow(startRow)
    texReWriter.setRow(startRow)

    angleSlice = 2 * math.pi / numVertices
    currAngle = 0

    #axisAdj=LMatrix4.rotateMat(45, axis)*LMatrix4.scaleMat(radius)*LMatrix4.translateMat(pos)

    perp1 = vecList[1]
    perp2 = vecList[2]

    # vertex information is written here
    for i in range(numVertices):
        adjCircle = pos + \
            (perp1 * math.cos(currAngle) + perp2 * math.sin(currAngle)) * \
            radius
        normal = perp1 * math.cos(currAngle) + perp2 * math.sin(currAngle)
        normalWriter.addData3f(normal)
        vertWriter.addData3f(adjCircle)
        texReWriter.addData2f(sCoord, (i + 0.001) / (numVertices - 1))
        colorWriter.addData4f(0.5, 0.5, 0.5, 1)
        drawReWriter.addData1f(keepDrawing)
        currAngle += angleSlice

    if startRow == 0:
        return

    drawReader = GeomVertexReader(vdata, "drawFlag")
    drawReader.setRow(startRow - numVertices)

    # we cant draw quads directly so we use Tristrips
    if drawReader.getData1i() != 0:
        lines = GeomTristrips(Geom.UHStatic)
        half = int(numVertices * 0.5)
        for i in range(numVertices):
            lines.addVertex(i + startRow)
            if i < half:
                lines.addVertex(i + startRow - half)
            else:
                lines.addVertex(i + startRow - half - numVertices)

        lines.addVertex(startRow)
        lines.addVertex(startRow - half)
        lines.closePrimitive()
        lines.decompose()
        circleGeom.addPrimitive(lines)

        circleGeomNode = GeomNode(cladeName if cladeName is not None else "")
        circleGeomNode.addGeom(circleGeom)

        # I accidentally made the front-face face inwards. Make reverse makes the tree render properly and
        # should cause any surprises to any poor programmer that tries to use
        # this code
        circleGeomNode.setAttrib(CullFaceAttrib.makeReverse(), 1)
        global numPrimitives
        numPrimitives += numVertices * 2

        nodePath.attachNewNode(circleGeomNode)


# this draws leafs when we reach an end
def drawLeaf(nodePath, vdata, speciesName, pos=LVector3(0, 0, 0), vecList=[LVector3(0, 0, 1), LVector3(1, 0, 0), LVector3(0, -1, 0)], scale=0.125):
    # use the vectors that describe the direction the branch grows to make the right
        # rotation matrix
    newCs = LMatrix4(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    newCs.setRow(0, vecList[2])  # right
    newCs.setRow(1, vecList[1])  # up
    newCs.setRow(2, vecList[0])  # forward
    newCs.setRow(3, (0, 0, 0))
    newCs.setCol(3, (0, 0, 0, 1))

    axisAdj = LMatrix4.scaleMat(scale) * newCs * LMatrix4.translateMat(pos)

    textNodeContainer = PandaNode(speciesName + "Container")
    textNodeContainer.setTransform(TransformState.makeMat(axisAdj))

    # orginlly made the leaf out of geometry but that didnt look good
    # I also think there should be a better way to handle the leaf texture other than
    # hardcoding the filename
    fontPath = ConfigVariableString("on-screen-debug-font", "Roboto-Regular.ttf").value
    font = loader.loadFont(fontPath, color=(1, 1, 1, 1), renderMode=TextFont.RMSolid)
    text = TextNode(speciesName)
    text.setText(speciesName)
    text.setFont(font)
    text.setAlign(TextNode.ACenter)
    billboard_effect = BillboardEffect.makePointEye()
    text.setEffect(billboard_effect)

    textNodeContainer.addChild(text)
    nodePath.attachNewNode(textNodeContainer)

def findCladeByName(name):
    for clade in phyloTree.find_clades():
        if clade.name == name:
            return clade

# recursive algorthim to make the tree
def makeFractalTree(bodydata, nodePath, length, cladeName=phyloTree.root.name, pos=LVector3(0, 0, 0), clade=phyloTree.root, vecList=[LVector3(0, 0, 1), LVector3(1, 0, 0), LVector3(0, -1, 0)]):
    if len(clade.clades) > 0:
        drawBody(nodePath, bodydata, cladeName, pos, vecList, length.getX())

        # move foward along the right axis
        newPos = pos + vecList[0] * length.length()

        length = LVector3(
                length.getX() / girthModifier, length.getY() / girthModifier, length.getZ() / 1.1)
        for c in clade.clades:
                makeFractalTree(bodydata, nodePath, length, clade.name, newPos, c, randomAxis(vecList))
    else:
        drawBody(nodePath, bodydata, cladeName, pos, vecList, length.getX())
        newPos = pos + vecList[0] * length.length()
        length = LVector3(
            length.getX() / girthModifier, length.getY() / girthModifier, length.getZ() / 1.1)
        drawBody(nodePath, bodydata, clade.name, newPos, randomAxis(vecList), length.getX(), False)
        drawLeaf(nodePath, bodydata, clade.name, newPos, vecList)


alight = AmbientLight('alight')
alight.setColor((0.5, 0.5, 0.5, 1))
alnp = render.attachNewNode(alight)
render.setLight(alnp)

slight = Spotlight('slight')
slight.setColor((1, 1, 1, 1))
lens = PerspectiveLens()
slight.setLens(lens)
slnp = render.attachNewNode(slight)
render.setLight(slnp)
isMenuOpen = False

slnp.setPos(0, 0, 40)

# rotating light to show that normals are calculated correctly
def updateLight(task):
    global slnp
    currPos = slnp.getPos()
    currPos.setX(100 * math.cos(task.time) / 2)
    currPos.setY(100 * math.sin(task.time) / 2)
    slnp.setPos(currPos)

    slnp.lookAt(render)
    return Task.cont

taskMgr.add(updateLight, "rotating Light")


# add some interactivity to the program
class MyTapper(DirectObject):

    def __init__(self):
        formatArray = GeomVertexArrayFormat()
        formatArray.addColumn(
            InternalName.make("drawFlag"), 1, Geom.NTUint8, Geom.COther)

        format = GeomVertexFormat(GeomVertexFormat.getV3n3cpt2())
        format.addArray(formatArray)
        self.format = GeomVertexFormat.registerFormat(format)

        bodydata = GeomVertexData("body vertices", format, Geom.UHStatic)

        self.barkTexture = loader.loadTexture("barkTexture.jpg")
        treeNodePath = NodePath("Tree Holder")
        makeFractalTree(bodydata, treeNodePath, LVector3(depth, depth, depth * 2))

        treeNodePath.setTexture(self.barkTexture, 1)
        treeNodePath.reparentTo(render)

        self.mouse_sens = 0.05
        base.mouseWatcherNode.set_modifier_buttons(ModifierButtons())
        base.buttonThrowers[0].node().set_modifier_buttons(ModifierButtons())

        # This is used to store which keys are currently pressed.
        self.keyMap = {"move-left": False, "move-right": False, "move-forward": False, "move-back": False, "move-up": False, "move-down": False}

        self.accept("q", self.regenTree)
        self.accept("t", self.addTree)
        self.accept("/", self.jumpToNode)
        self.accept("escape", self.toggleMenu)
        self.accept("mouse1", self.leftClick)
        self.accept("arrow_up", self.upIterations)
        self.accept("arrow_down", self.downIterations)
        self.accept("arrow_right", self.upCopies)
        self.accept("arrow_left", self.downCopies)
        self.accept("w", self.setKey, ["move-forward", True])
        self.accept("a", self.setKey, ["move-left", True])
        self.accept("s", self.setKey, ["move-back", True])
        self.accept("d", self.setKey, ["move-right", True])
        self.accept("space", self.setKey, ["move-up", True])
        self.accept("lshift", self.setKey, ["move-down", True])
        self.accept("w-up", self.setKey, ["move-forward", False])
        self.accept("a-up", self.setKey, ["move-left", False])
        self.accept("s-up", self.setKey, ["move-back", False])
        self.accept("d-up", self.setKey, ["move-right", False])
        self.accept("space-up", self.setKey, ["move-up", False])
        self.accept("lshift-up", self.setKey, ["move-down", False])

        taskMgr.add(self.move, "moveTask")

        handler = CollisionHandlerQueue()
        traverser = CollisionTraverser('traverser name')

        pickerNode = CollisionNode('mouseRay')
        pickerNP = camera.attachNewNode(pickerNode)
        pickerNode.setFromCollideMask(GeomNode.getDefaultCollideMask())
        pickerRay = CollisionRay()
        pickerNode.addSolid(pickerRay)
        traverser.addCollider(pickerNP, handler)

        base.cTrav = traverser
        base.cHandler = handler
        base.pickerRay = pickerRay

        cm = CardMaker("plane")
        cm.setFrame(-2560, 2560, -2560, 2560)  # set the size here
        cm.setColor(74.0, 103.0, 65.0, 1.0)
        plane = render.attachNewNode(cm.generate())
        plane.setP(-90.)
        

    # Records the state of the arrow keys
    def setKey(self, key, value):
        self.keyMap[key] = value

    def upIterations(self):
        self.numIterations += 1
        self.upDownEvent.setText(
            "Up/Down: Increase/Decrease the number of iterations (" + str(self.numIterations) + ")")

    def downIterations(self):
        self.numIterations -= 1
        self.upDownEvent.setText(
            "Up/Down: Increase/Decrease the number of iterations (" + str(self.numIterations) + ")")

    def upCopies(self):
        self.numCopies += 1
        self.leftRightEvent.setText(
            "Left/Right: Increase/Decrease branching (" + str(self.numCopies) + ")")

    def downCopies(self):
        self.numCopies -= 1
        self.leftRightEvent.setText(
            "Left/Right: Increase/Decrease branching (" + str(self.numCopies) + ")")

    def regenTree(self):
        forest = render.findAllMatches("Tree Holder")
        forest.detach()

        bodydata = GeomVertexData("body vertices", self.format, Geom.UHStatic)

        global depth
        global girthModifier
        depth = max(phyloTree.depths(unit_branch_lengths=True).values())
        # calculate modifier for a leaf node thickness of .025
        girthModifier = pow((1.0 + (depth * 40.0)), (1.0 / depth))

        treeNodePath = NodePath("Tree Holder")
        makeFractalTree(bodydata, treeNodePath, LVector3(depth, depth, depth * 2), LVector3(0, 0, 0))

        treeNodePath.setTexture(self.barkTexture, 1)
        treeNodePath.reparentTo(render)

    def addTree(self):
        bodydata = GeomVertexData("body vertices", self.format, Geom.UHStatic)

        randomPlace = LVector3(
            random.randint(-2560, 2560), random.randint(-2560, 2560), 0)
        # randomPlace.normalize()

        treeNodePath = NodePath("Tree Holder")

        rootClade = phyloTree.root
        if selectedCladeName is not None and selectedCladeName != "":
            rootClade = findCladeByName(selectedCladeName)

        global depth
        global girthModifier
        depth = max(rootClade.depths(unit_branch_lengths=True).values())
        # calculate modifier for a leaf node thickness of .025
        girthModifier = pow((1.0 + (depth * 40.0)), (1.0 / depth))

        makeFractalTree(bodydata, treeNodePath, LVector3(depth, depth, depth * 2), rootClade.name, randomPlace, rootClade)

        treeNodePath.setTexture(self.barkTexture, 1)
        treeNodePath.reparentTo(render)

    def jumpToNode(self):
        if selectedCladeName is None:
            return
        nodePath = render.find("**/" + selectedCladeName)
        tightBounds = nodePath.getTightBounds()
        centeredX = mean([tightBounds[0].getXy().getX(), tightBounds[1].getXy().getX()])
        centeredY = mean([tightBounds[0].getXy().getY(), tightBounds[1].getXy().getY()])
        centeredZ = mean([tightBounds[0].getYz().getY(), tightBounds[1].getYz().getY()])
        camera.setPos(LPoint3f(centeredX, centeredY, centeredZ))
        camera.setY(camera, -10)

    def itemSel(self, arg):
        #print("Item Selected is: " + arg)
        if arg == "exit":
            sys.exit()



    def toggleMenu(self):
        global isMenuOpen
        isMenuOpen = not isMenuOpen
        if isMenuOpen:
            menu = DirectOptionMenu(text="test", scale=0.1, command=self.itemSel,
                                items=["exit"], initialitem=-1,
                                highlightColor=(0.65, 0.65, 0.65, 1))



    def leftClick(self):
        # First we check that the mouse is not outside the screen.
        if base.mouseWatcherNode.hasMouse():
            # This gives up the screen coordinates of the mouse.
            mpos = base.mouseWatcherNode.getMouse()

        # This makes the ray's origin the camera and makes the ray point
        # to the screen coordinates of the mouse.
        base.pickerRay.setFromLens(base.camNode, mpos.x, mpos.y)

        mpos = base.mouseWatcherNode.getMouse()
        base.pickerRay.setFromLens(base.camNode, mpos.getX(), mpos.getY())

        base.cTrav.traverse(render)
        # Assume for simplicity's sake that myHandler is a CollisionHandlerQueue.
        if base.cHandler.getNumEntries() > 0:
            # This is so we get the closest object
            base.cHandler.sortEntries()
            pickedObj = base.cHandler.getEntry(0).getIntoNodePath()
            if (pickedObj.name != "plane"):
                title.text = pickedObj.name if pickedObj.name != "" else "Selected Clade Has No Name"
                global selectedCladeName
                selectedCladeName = pickedObj.name

    def move(self, task):
        if isMenuOpen:
            return
        # Get the time that elapsed since last frame.  We multiply this with
        # the desired speed in order to find out with which distance to move
        # in order to achieve that desired speed.
        dt = globalClock.getDt()

        # If the camera-left key is pressed, move camera left.
        # If the camera-right key is pressed, move camera right.

        if self.keyMap["move-left"]:
            camera.setX(camera, -150 * dt)
        if self.keyMap["move-right"]:
            camera.setX(camera, +150 * dt)
        if self.keyMap["move-forward"]:
            camera.setY(camera, +150 * dt)
        if self.keyMap["move-back"]:
            camera.setY(camera, -150 * dt)
        if self.keyMap["move-up"]:
            camera.setZ(camera, +150 * dt)
        if self.keyMap["move-down"]:
            camera.setZ(camera, -150 * dt)

        md = base.win.getPointer(0)
        x = md.getX()
        y = md.getY()
        if base.win.movePointer(0, base.win.getXSize()//2, base.win.getYSize()//2):
            camera.setH(camera.getH() - (x - base.win.getXSize()//2)*self.mouse_sens) 
            camera.setP(camera.getP() - (y - base.win.getYSize()//2)*self.mouse_sens)

        return task.cont



t = MyTapper()
print(numPrimitives)

base.run()
