#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import sys
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
conv.ReadFile(mol, sys.argv[1])

pattern.Match(mol)
mollist = pattern.GetUMapList()[0]
natom = mol.GetAtom(mollist[0])
makeResidue(mol, 0, mollist)

i = 1
for aa in sys.argv[2:]:
    mol2 = ob.OBMol()
    conv.ReadFile(mol2, aa)
    pattern.Match(mol2)
    mollist = pattern.GetUMapList()[0]
    makeResidue(mol2, i, mollist)

    catom = mol2.GetAtom(mollist[2])
    oatom = mol2.GetAtom(mollist[4])

    molnatoms = mol.NumAtoms()
    mol += mol2

    builder.Connect(mol, natom.GetIdx(), molnatoms + catom.GetIdx())
    foatom = mol.GetAtom(molnatoms + oatom.GetIdx())
    fnatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0]).GetIdx())
    mol.DeleteHydrogens(foatom)
    mol.DeleteAtom(foatom)

    natom.SetImplicitValence(3)
    mol.DeleteHydrogens(natom)
    mol.AddHydrogens(natom)
    natom = fnatom
    i += 1

builder.Build(mol)
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

conv.WriteFile(mol, "test.svg")
