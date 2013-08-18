#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import sys
import math
import openbabel as ob

def makeResidue(mol, idx, aaatoms):
    res = mol.NewResidue()
    res.SetNum(idx)
    for atom in ob.OBMolAtomIter(mol):
        if atom.GetIdx() not in aaatoms:
            res.AddAtom(atom)

mol = ob.OBMol()
conv = ob.OBConversion()
pattern = ob.OBSmartsPattern()
pattern.Init("[NX3][$([CX4H1]([*])),$([CX4H2])][CX3](=[OX1])[OX2]")
builder = ob.OBBuilder()
conv.SetInAndOutFormats("sdf", "svg")
conv.AddOption("d", ob.OBConversion.OUTOPTIONS)
conv.AddOption("b", ob.OBConversion.OUTOPTIONS, "none")
conv.ReadFile(mol, sys.argv[1])

pattern.Match(mol)
mollist = pattern.GetUMapList()[0]
oatom = mol.GetAtom(mollist[4])
catom = mol.GetAtom(mollist[2])
firstnatom = mol.GetAtom(mollist[0])
makeResidue(mol, 0, mollist)

i = 1
for aa in sys.argv[2:]:
    mol2 = ob.OBMol()
    conv.ReadFile(mol2, aa)
    pattern.Match(mol2)
    mollist = pattern.GetUMapList()[0]
    makeResidue(mol2, i, mollist)

    molnatoms = mol.NumAtoms()
    mol += mol2

    natom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0]).GetIdx())

    builder.Connect(mol, catom.GetIdx(), natom.GetIdx())

    foatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[4]).GetIdx())
    catom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[2]).GetIdx())
    mol.DeleteHydrogens(oatom)
    mol.DeleteAtom(oatom)

    natom.SetImplicitValence(3)
    mol.DeleteHydrogens(natom)
    mol.AddHydrogens(natom)
    oatom = foatom

    i += 1

nidx = firstnatom.GetIdx()
oidx =  oatom.GetIdx()
builder.Build(mol)
natom = mol.GetAtom(nidx)
oatom = mol.GetAtom(oidx)

for res in ob.OBResidueIter(mol):
    if res.GetNum() % 2 > 0:
        color = "red"
    else:
        color = "green"
    for atom in ob.OBResidueAtomIter(res):
        for bond in ob.OBAtomBondIter(atom):
            data = ob.OBPairData()
            data.SetAttribute("color")
            data.SetValue(color)
            bond.CloneData(data)
mol.DeleteHydrogens()
#mol.AddHydrogens()

gen2d = ob.OBOp.FindType("gen2d")
gen2d.Do(mol)
print "x=%f; y=%f; z=%f" % (natom.GetX(), natom.GetY(), natom.GetZ())
print "x=%f; y=%f; z=%f" % (oatom.GetX(), oatom.GetY(), oatom.GetZ())

opp = oatom.GetY() - natom.GetY()
adj = oatom.GetX() - natom.GetX()
angle = abs(math.atan(opp / adj))
if opp > 0 and adj > 0:
    pass
elif opp > 0 and adj < 0:
    angle = math.pi - angle
elif opp < 0 and adj < 0:
    angle = math.pi + angle
elif opp < 0 and adj > 0:
    angle = 2 * math.pi - angle
print math.degrees(angle)
angle = -angle
mol.Rotate(ob.double_array([math.cos(angle), -math.sin(angle), 0,
                            math.sin(angle), math.cos(angle),  0,
                            0,               0,                1]))

svg = conv.WriteString(mol)
# need to get rid of square aspect ratio
delstart = svg.find("width")
delend = svg.find("svg", delstart)
delend = svg.find("viewBox", delend)
svgend = svg.rfind("</g>")
svg = svg[0:delstart] + svg[delend:svgend]
outf = file("test.svg", "w")
outf.write(svg)
outf.close()
