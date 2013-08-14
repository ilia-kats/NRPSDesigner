#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import sys
import openbabel as ob

mol = ob.OBMol()
conv = ob.OBConversion()
pattern = ob.OBSmartsPattern()
pattern.Init("[NX3][$([CX4H1]([*])),$([CX4H2])][CX3](=[OX1])[OX2]")
builder = ob.OBBuilder()
conv.SetInAndOutFormats("sdf", "svg")
conv.ReadFile(mol, sys.argv[1])
pattern.Match(mol)
mollist = pattern.GetUMapList()
natom = mol.GetAtom(mollist[0][0])

for aa in sys.argv[2:]:
    mol2 = ob.OBMol()
    conv.ReadFile(mol2, aa)
    pattern.Match(mol2)
    mollist = pattern.GetUMapList()

    catom = mol2.GetAtom(mollist[0][2])
    oatom = mol2.GetAtom(mollist[0][4])

    molnatoms = mol.NumAtoms()
    mol += mol2

    builder.Connect(mol, natom.GetIdx(), molnatoms + catom.GetIdx())
    foatom = mol.GetAtom(molnatoms + oatom.GetIdx())
    fnatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0][0]).GetIdx())
    mol.DeleteHydrogens(foatom)
    mol.DeleteAtom(foatom)
    #mol.DeleteHydrogens()
    #mol.AddHydrogens()

    natom.SetImplicitValence(3)
    mol.DeleteHydrogens(natom)
    mol.AddHydrogens(natom)
    natom = fnatom

builder.Build(mol)
conv.WriteFile(mol, "test.svg") 
