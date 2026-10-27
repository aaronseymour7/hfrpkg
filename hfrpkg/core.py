#!/usr/bin/env python3

import sys
import argparse
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from pulp import LpProblem, LpVariable, lpSum, LpMinimize, LpInteger, LpStatus
from rdkit.Chem import Descriptors, rdMolDescriptors
from itertools import combinations
from AaronTools.geometry import Geometry
from AaronTools.atoms import Atom
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from hfrpkg.common import make_solver

def kekulize_and_make_aliphatic(mol):
    try:
        mol = Chem.Mol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in mol.GetBonds():
            bond.SetIsAromatic(False)
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol)
        return mol 
    except Chem.KekulizeException:
        return None  
def get_fragments(input_mol,level):
    wildcard_reactants = []
    wildcard_products= []
    mol= kekulize_and_make_aliphatic(input_mol)
    if level==1:
        wildcard_reactants = [
            Chem.MolFromSmiles("[H][H]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*]"),
        ]
    
    if level==2:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*]#[*]"),
        ]

    if level ==3:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*]#[*]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*][*][*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*]#[*][*]"),
            Chem.MolFromSmiles("[*][*](*)[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]"),
            Chem.MolFromSmiles("[*][*@]([*])([*])[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]"),    
        ]
    if level == 4:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*][*]=[*][*]"),
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]#[*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*]#[*][*]"),
        ]

        wildcard_products = [
            Chem.MolFromSmiles("[*][*][*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*][*]#[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]"),
            Chem.MolFromSmiles("[*][*]([*])[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]"),
            Chem.MolFromSmiles("[*][*@]([*])([*])[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]=[*]"),
            Chem.MolFromSmiles("[*]=[*][*]#[*]"),
            Chem.MolFromSmiles("[*]#[*][*]#[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]=[*]"),
            Chem.MolFromSmiles("[*]=[*][*]=[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]#[*]"),
            Chem.MolFromSmarts("[*]=[*]-[!#6]-[*]=[*]"),
            Chem.MolFromSmiles("[*]([*])[*]=[*]"),
            Chem.MolFromSmarts("[*]#[*]-[!#6]-[*]=[*]"),
            Chem.MolFromSmarts("[*]#[*]-[!#6]-[*]#[*]"),
        ]
    wildcard_reactants = sorted(wildcard_reactants, key=lambda x: x.GetNumAtoms(), reverse=True)
    wildcard_products = sorted(wildcard_products, key=lambda x: x.GetNumAtoms(), reverse=True)
    return mol, wildcard_reactants, wildcard_products
def map_atoms(matched_structure, mol):
    p = Chem.AdjustQueryParameters.NoAdjustments()
    p.makeDummiesQueries = True
    query = Chem.AdjustQueryProperties(matched_structure, p) 
    matches = mol.GetSubstructMatches(query)
    if matches:
        return matches
    else:
        print(f"No atom mapping found for substructure: {Chem.MolToSmiles(matched_structure)}")
    return matched_structure
def match_substructure(mol, sub):
    p = Chem.AdjustQueryParameters()
    p.makeDummiesQueries = True
    p.adjustDegree = False
    p.adjustRingCount = False
    p.adjustAromatic = False
    p.adjustValence = False
    substructure = Chem.AdjustQueryProperties(sub, p)


    match = mol.HasSubstructMatch(substructure)
    return match
def append_unique(mol_list, mol):
    mol_smiles = Chem.MolToSmiles(mol)
    mol_list[:] = [m for m in mol_list if Chem.MolToSmiles(m) != mol_smiles]
    mol_list.append(mol)
def product_options(mol, wildcard_products):
    # Define wildcard libraries
    # Sort the wildcard libraries by size (larger first)


    # Search for matching substructures (reactants and products)
    #mol= Chem.RemoveHs(mol)
    matched_products = []
    for product in wildcard_products:
        if match_substructure(mol, product):
            matched_products.append(product)
    # Map atom IDs from reactants/products to the molecule
    mapped_products = []
    for product in matched_products:
        matches = map_atoms(product, mol)  # Get matched atom indices    
        if not matches:
            continue  # Skip if no matches found
        for match in matches:
            # Convert to an editable molecule
            rw_product = Chem.RWMol(product)
            for wild_idx, mol_idx in enumerate(match):  
                product_atom = rw_product.GetAtomWithIdx(wild_idx)  # Get atom in product
                mol_atom = mol.GetAtomWithIdx(mol_idx)  # Get corresponding atom in mol
                if product_atom.GetAtomicNum() == 0:  # If it's a wildcard atom
                    new_atomic_num = mol_atom.GetAtomicNum()  # Get atomic number from mol
                    product_atom.SetAtomicNum(new_atomic_num)  # Assign atomic number
                    product_atom.SetNoImplicit(False)
            try:
                if Chem.SanitizeMol(rw_product):
                    Chem.SanitizeMol(rw_product)
            except Exception as e:
                continue
            mapped_products.append(rw_product.GetMol()) # Convert back to immutable Mol
    #return mapped_products
    unique_products = []
    unique_smiles = set()    
    for prod in mapped_products:
            prod = Chem.AddHs(prod)
            smiles = Chem.MolToSmiles(prod)  # Convert to SMILES for easy comparison
            if smiles not in unique_smiles:
                unique_smiles.add(smiles)
                unique_products.append(prod)
    return unique_products
def reactant_options(unique_products, wildcard_reactants):
    matched_reactants = []
    mapped_reactants = []

    for rhs_product in unique_products:  # Iterate over RHS products
        for reactant in wildcard_reactants:  # Iterate over wildcard reactants
            if match_substructure(rhs_product, reactant):  # Check if reactant matches product
                matched_reactants.append((reactant, rhs_product))  # Store pair

    for reactant, rhs_product in matched_reactants:
        matches = map_atoms(reactant, rhs_product)  # Get atom mapping
        if not matches:
            continue  # Skip unmatched reactants

        rw_reactant = Chem.RWMol(reactant)  # Editable reactant molecule
        for match in matches:
            for wild_idx, prod_idx in enumerate(match):  
                reactant_atom = rw_reactant.GetAtomWithIdx(wild_idx)
                product_atom = rhs_product.GetAtomWithIdx(prod_idx)

                if reactant_atom.GetAtomicNum() == 0:  # If it's a wildcard atom
                    new_atomic_num = product_atom.GetAtomicNum()  # Assign atomic number from product
                    reactant_atom.SetAtomicNum(new_atomic_num)
                    reactant_atom.SetNoImplicit(False)
                    new_atomic_num = product_atom.GetAtomicNum()  # Get atomic number from mol
                    reactant_atom.SetAtomicNum(new_atomic_num)  # Assign atomic number
                    reactant_atom.SetNoImplicit(False)
            Chem.SanitizeMol(rw_reactant)
            mapped_reactants.append(rw_reactant.GetMol())
              # Convert back to immutable Mol


    # Ensure uniqueness of mapped reactants
    unique_reactants = []
    unique_smiles = set()

    for reactant in mapped_reactants:
        reactant = Chem.AddHs(reactant)
        smiles = Chem.MolToSmiles(reactant)  
        if smiles not in unique_smiles:
            unique_smiles.add(smiles)
            unique_reactants.append(reactant)


    return unique_reactants

# Forbid the trivial identity reaction (input -> input). Switched off only to tell "no balanced reaction"
# apart from "only the identity reaction" (paper/feasibility_scan.py).
EXCLUDE_IDENTITY = True
# Add the minimum-species objective. Switched off (together with EXCLUDE_IDENTITY) this reproduces the
# solver formulation that generated the published benchmark reactions (paper/selection_check.py).
USE_OBJECTIVE = True


def _same_molecule(a, b):
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(a)) == Chem.MolToSmiles(Chem.RemoveHs(b))
    except Exception:
        return False


def _finish_model(prob, x_vars, y_vars, rhs_options, input_mol, rhs_required):
    """The models had no objective, so any feasible point was returned - including
    the identity reaction (input -> input), whose "dHf" is just the ATcT value
    again. Forbid the identity product (unless the user demanded it) and prefer
    the smallest reaction, which also makes the answer deterministic."""
    required = rhs_required if isinstance(rhs_required, list) else ([rhs_required] if rhs_required else [])
    for i, opt in enumerate(rhs_options):
        if EXCLUDE_IDENTITY and _same_molecule(opt, input_mol) and not any(_same_molecule(opt, r) for r in required):
            prob += x_vars[i] == 0
    if USE_OBJECTIVE:
        prob += lpSum(x_vars) + lpSum(y_vars)


def Isogyric(input_mol, lhs_required=None, rhs_required=None, Substruct=None, Replacement=None):
    original_mol= input_mol
    input_mol, wildcard_reactants, wildcard_products = get_fragments(input_mol,1)
    def count_atoms(mol):
        atom_counts = {}
        if isinstance(mol, Chem.Atom):
            key = (mol.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
            return atom_counts
        for atom in mol.GetAtoms():
        
            key = (atom.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    def Balance(input_mol):
        prob = LpProblem("Isogyric_Balancing", LpMinimize)
        
        dummy = input_mol
        if Substruct and Replacement:
            if isinstance(Substruct, list) and isinstance(Replacement, list):
                for sub, rep in zip(Substruct, Replacement):
                    dummy = Chem.ReplaceSubstructs(dummy, sub, rep, True)[0]
            elif isinstance(Substruct, Chem.Mol) and isinstance(Replacement, Chem.Mol):
                dummy = Chem.ReplaceSubstructs(dummy, Substruct, Replacement, True)[0]
        rhs_options = product_options(dummy,wildcard_products)
        for atom in rhs_options:
            Chem.AddHs(atom)
        lhs_addition_options = reactant_options(rhs_options, wildcard_reactants)
        lhs_addition_options.append(input_mol)
        lhs_addition_options.append(Chem.MolFromSmiles("[H][H]"))
        input_index = lhs_addition_options.index(input_mol)
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_addition_options.append(mol)
            elif isinstance (lhs_required, Chem.Mol):
                    lhs_addition_options.append(lhs_required)
                
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    append_unique(rhs_options, mol)
            elif isinstance (rhs_required, Chem.Mol):
                append_unique(rhs_options, rhs_required)
        
        
        x_vars = [LpVariable(f"x_{i}", 0, None, LpInteger) for i in range(len(rhs_options))]
        y_vars = [LpVariable(f"y_{j}", 0, None, LpInteger) for j in range(len(lhs_addition_options))]
        prob += y_vars[input_index] == 1       
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_r_idx = lhs_addition_options.index(mol) 
                    prob += y_vars[lhs_r_idx] == 1
            elif isinstance (lhs_required, Chem.Mol):
                lhs_r_idx = lhs_addition_options.index(lhs_required) 
                prob += y_vars[lhs_r_idx] == 1
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    rhs_r_idx = rhs_options.index(mol) 
                    prob += x_vars[rhs_r_idx] == 1
            elif isinstance (rhs_required, Chem.Mol):
                rhs_r_idx = rhs_options.index(rhs_required) 
                prob += x_vars[rhs_r_idx] == 1
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], [], "Infeasible"
    
        atom_types_rhs = [count_atoms(mol) for mol in rhs_options]
        atom_types_lhs = [count_atoms(mol) for mol in lhs_addition_options]
        
        for t in set().union(*atom_types_rhs, *atom_types_lhs):
            prob += lpSum(x_vars[i] * atom_types_rhs[i].get(t, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * atom_types_lhs[j].get(t, 0) for j in range(len(lhs_addition_options)))
    
        _finish_model(prob, x_vars, y_vars, rhs_options, input_mol, rhs_required)
        solver = make_solver()

        prob.solve(solver)

        if LpStatus[prob.status] != 'Optimal':

            return [], [], LpStatus[prob.status]

        solution_x = [int(var.varValue) for var in x_vars]
        solution_y = [int(var.varValue) for var in y_vars]
    
        balanced_rhs = [(rhs_options[i], solution_x[i]) for i in range(len(rhs_options)) if solution_x[i] > 0]
        balanced_lhs = [(lhs_addition_options[j], solution_y[j]) for j in range(len(lhs_addition_options)) if solution_y[j] > 0]
        #print("Solver status:", LpStatus[prob.status])
        return balanced_rhs, balanced_lhs, LpStatus[prob.status]
    return Balance(input_mol)
            
            
            
def Isodesmic(nput_mol, lhs_required=None, rhs_required=None, Substruct=None, Replacement=None):
    input_mol = kekulize_and_make_aliphatic(nput_mol)
    def count_bonds(mol):
        bond_counts = {}
        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetSymbol()== 'H' or bond.GetEndAtom().GetSymbol()== 'H':
                continue 
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
    
            # Create a tuple where atom order doesn't matter by sorting the atoms based on their symbol and hybridization
            atom_pair = tuple(sorted(
                [(atom1.GetSymbol()), 
                 (atom2.GetSymbol())]
            ))
    
            bond_type = (atom_pair, bond.GetBondType())
            bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1
    
        return bond_counts
    def count_atoms(mol):
        atom_counts = {}
        if isinstance(mol, Chem.Atom):
            key = (mol.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
            return atom_counts
        for atom in mol.GetAtoms():

            key = (atom.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    def product_frags(atom_pair, bond_type):
        mol = Chem.RWMol()

        a1_idx = mol.AddAtom(Chem.Atom(atom_pair[0]))
        a2_idx = mol.AddAtom(Chem.Atom(atom_pair[1]))

        mol.AddBond(a1_idx, a2_idx, bond_type)
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)       
        mol = Chem.AddHs(mol)

        return mol
    def unique_atom_hydrides(mol):
        seen_elements = set()
        hydride_mols = []

        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym == 'H' or sym in seen_elements:
                continue
            seen_elements.add(sym)
            try:
                m = Chem.MolFromSmiles(sym)
                m = Chem.AddHs(m)         
                hydride_mols.append(( m))
            except:
                print(f"Could not generate molecule for {sym}")

        return hydride_mols
    def Balance(input_mol):
        prob = LpProblem("Isodesmic_Balancing", LpMinimize)
        dummy = input_mol
        if Substruct and Replacement:
            if isinstance(Substruct, list) and isinstance(Replacement, list):
                for sub, rep in zip(Substruct, Replacement):
                    dummy = Chem.ReplaceSubstructs(dummy, sub, rep, True)[0]
            elif isinstance(Substruct, Chem.Mol) and isinstance(Replacement, Chem.Mol):
                dummy = Chem.ReplaceSubstructs(dummy, Substruct, Replacement, True)[0]
        rhs_options = []
        prod_options = count_bonds(input_mol)
        for opt in prod_options:
            rhs_options.append(product_frags(*opt))
        lhs_addition_options = []
        lhs_addition_options = unique_atom_hydrides(input_mol)
        lhs_addition_options.append(Chem.MolFromSmiles("[H][H]"))
        lhs_addition_options.append(input_mol)
        input_index = lhs_addition_options.index(input_mol)
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_addition_options.append(mol)
            elif isinstance (lhs_required, Chem.Mol):
                    lhs_addition_options.append(lhs_required)
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    append_unique(rhs_options, mol)
            elif isinstance (rhs_required, Chem.Mol):
                append_unique(rhs_options, rhs_required)
        


        x_vars = [LpVariable(f"x_{i}", 0, None, LpInteger) for i in range(len(rhs_options))]
        y_vars = [LpVariable(f"y_{j}", 0, None, LpInteger) for j in range(len(lhs_addition_options))]
        prob += y_vars[input_index] == 1       
        
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_r_idx = lhs_addition_options.index(mol) 
                    prob += y_vars[lhs_r_idx] == 1
            elif isinstance (lhs_required, Chem.Mol):
                lhs_r_idx = lhs_addition_options.index(lhs_required) 
                prob += y_vars[lhs_r_idx] == 1
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    rhs_r_idx = rhs_options.index(mol) 
                    prob += x_vars[rhs_r_idx] == 1
            elif isinstance (rhs_required, Chem.Mol):
                rhs_r_idx = rhs_options.index(rhs_required) 
                prob += x_vars[rhs_r_idx] == 1
    
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], [], "Infeasible"
    
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], [], "Infeasible"
    
        bond_types_rhs = [count_bonds(mol) for mol in rhs_options]
        bond_types_lhs = [count_bonds(mol) for mol in lhs_addition_options]
        atom_types_rhs = [count_atoms(mol) for mol in rhs_options]
        atom_types_lhs = [count_atoms(mol) for mol in lhs_addition_options]
        
        for k in set().union(*bond_types_rhs, *bond_types_lhs):
            prob += lpSum(x_vars[i] * bond_types_rhs[i].get(k, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * bond_types_lhs[j].get(k, 0) for j in range(len(lhs_addition_options)))
        for t in set().union(*atom_types_rhs, *atom_types_lhs):
            prob += lpSum(x_vars[i] * atom_types_rhs[i].get(t, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * atom_types_lhs[j].get(t, 0) for j in range(len(lhs_addition_options)))
    
        _finish_model(prob, x_vars, y_vars, rhs_options, input_mol, rhs_required)
        solver = make_solver()
        prob.solve(solver)
        if LpStatus[prob.status] != 'Optimal':
            return [], [], LpStatus[prob.status]
        
        solution_x = [int(var.varValue) for var in x_vars]
        solution_y = [int(var.varValue) for var in y_vars]
    
        balanced_rhs = [(rhs_options[i], solution_x[i]) for i in range(len(rhs_options)) if solution_x[i] > 0]
        balanced_lhs = [(lhs_addition_options[j], solution_y[j]) for j in range(len(lhs_addition_options)) if solution_y[j] > 0]
        #print("Solver status:", LpStatus[prob.status])
        return balanced_rhs, balanced_lhs, LpStatus[prob.status]
    #for m in unique_atom_hydrides(input_mol):
       # print(Chem.MolToSmiles(m))
    return Balance(input_mol)
            
            
            
            
def Hypohomodesmotic(input_mol, lhs_required=None, rhs_required=None, Substruct=None, Replacement=None):
    original_mol= input_mol
    input_mol, wildcard_reactants, wildcard_products = get_fragments(input_mol,3)    
    def count_hydrogens(mol):
        mol = Chem.AddHs(mol)
        hydrogen_counts = {}
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        for atom in mol.GetAtoms():
            key = (atom.GetSymbol(), atom.GetNumExplicitHs())
            hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
        return hydrogen_counts
    def count_atoms(mol):
        mol = Chem.RemoveHs(mol)
        atom_counts = {}
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        for atom in mol.GetAtoms():
            key = (atom.GetSymbol(), atom.GetHybridization())
            atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    def Balance(input_mol):
        prob = LpProblem("Hypohomodesmotic_Balancing", LpMinimize)
        
        dummy = kekulize_and_make_aliphatic(input_mol)
        
        if Substruct and Replacement:
            if isinstance(Substruct, list) and isinstance(Replacement, list):
                for sub, rep in zip(Substruct, Replacement):
                    dummy = Chem.ReplaceSubstructs(dummy, sub, rep, True)[0]
            elif isinstance(Substruct, Chem.Mol) and isinstance(Replacement, Chem.Mol):
                dummy = Chem.ReplaceSubstructs(dummy, Substruct, Replacement, True)[0]
        rhs_options = product_options(dummy, wildcard_products)

        lhs_addition_options = reactant_options(rhs_options, wildcard_reactants)

        for i, opt in enumerate(rhs_options):
            rhs_options[i] = Chem.AddHs(opt)

        for i, opt in enumerate(lhs_addition_options):
            lhs_addition_options[i] = Chem.AddHs(opt)
        lhs_addition_options.append(input_mol)
        input_index = lhs_addition_options.index(input_mol)
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_addition_options.append(mol)
            elif isinstance (lhs_required, Chem.Mol):
                    lhs_addition_options.append(lhs_required)
                
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    append_unique(rhs_options, mol)
            elif isinstance (rhs_required, Chem.Mol):
                append_unique(rhs_options, rhs_required)
        
        
        x_vars = [LpVariable(f"x_{i}", 0, None, LpInteger) for i in range(len(rhs_options))]
        y_vars = [LpVariable(f"y_{j}", 0, None, LpInteger) for j in range(len(lhs_addition_options))]
        prob += y_vars[input_index] == 1       
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_r_idx = lhs_addition_options.index(mol) 
                    prob += y_vars[lhs_r_idx] == 1
            elif isinstance (lhs_required, Chem.Mol):
                lhs_r_idx = lhs_addition_options.index(lhs_required) 
                prob += y_vars[lhs_r_idx] == 1
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    rhs_r_idx = rhs_options.index(mol) 
                    prob += x_vars[rhs_r_idx] == 1
            elif isinstance (rhs_required, Chem.Mol):
                rhs_r_idx = rhs_options.index(rhs_required) 
                prob += x_vars[rhs_r_idx] == 1
        
        
    

    
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], [], "Infeasible"
    
        atom_types_rhs = [count_atoms(mol) for mol in rhs_options]
        atom_types_lhs = [count_atoms(mol) for mol in lhs_addition_options]
        hydrogen_rhs = [count_hydrogens(mol) for mol in rhs_options]
        hydrogen_lhs = [count_hydrogens(mol) for mol in lhs_addition_options]
    
        for k in set().union(*atom_types_rhs, *atom_types_lhs):
            prob += lpSum(x_vars[i] * atom_types_rhs[i].get(k, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * atom_types_lhs[j].get(k, 0) for j in range(len(lhs_addition_options)))
    
        for t in set().union(*hydrogen_rhs, *hydrogen_lhs):
            prob += lpSum(x_vars[i] * hydrogen_rhs[i].get(t, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * hydrogen_lhs[j].get(t, 0) for j in range(len(lhs_addition_options)))
    
        _finish_model(prob, x_vars, y_vars, rhs_options, input_mol, rhs_required)
        solver = make_solver()
        prob.solve(solver)
        if LpStatus[prob.status] != 'Optimal':
            return [], [], LpStatus[prob.status]
        solution_x = [int(var.varValue) for var in x_vars]
        solution_y = [int(var.varValue) for var in y_vars]
    
        balanced_rhs = [(rhs_options[i], solution_x[i]) for i in range(len(rhs_options)) if solution_x[i] > 0]
        balanced_lhs = [(lhs_addition_options[j], solution_y[j]) for j in range(len(lhs_addition_options)) if solution_y[j] > 0]
        #print("Solver status:", LpStatus[prob.status])
        return balanced_rhs, balanced_lhs, LpStatus[prob.status]
    return Balance(input_mol)

            
            
            
def Homodesmotic(input_mol, lhs_required=None, rhs_required=None, Substruct=None, Replacement=None):
    original_mol= input_mol
    x = kekulize_and_make_aliphatic(input_mol)
    input_mol, wildcard_reactants, wildcard_products = get_fragments(x,4)    
    def count_bonds(mol):
        bond_counts = {}
        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetSymbol()== 'H' or bond.GetEndAtom().GetSymbol()== 'H':
                continue 
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
    
            # Create a tuple where atom order doesn't matter by sorting the atoms based on their symbol and hybridization
            atom_pair = tuple(sorted(
                [(atom1.GetSymbol(), atom1.GetHybridization()), 
                 (atom2.GetSymbol(), atom2.GetHybridization())]
            ))
    
            bond_type = (atom_pair, bond.GetBondType())
            bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1
    
        return bond_counts
    
    def count_hydrogens(mol):
        #mol = Chem.AddHs(mol)
        hydrogen_counts = {}
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        for atom in mol.GetAtoms():
            key = (atom.GetSymbol(), atom.GetHybridization(), atom.GetNumExplicitHs())
            hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
        return hydrogen_counts
    def Balance(input_mol):
        prob = LpProblem("Chemical_Balancing", LpMinimize)
        input_mol = kekulize_and_make_aliphatic(input_mol)
        dummy = input_mol
        if Substruct and Replacement:
            if isinstance(Substruct, list) and isinstance(Replacement, list):
                for sub, rep in zip(Substruct, Replacement):
                    dummy = Chem.ReplaceSubstructs(dummy, sub, rep, True)[0]
            elif isinstance(Substruct, Chem.Mol) and isinstance(Replacement, Chem.Mol):
                dummy = Chem.ReplaceSubstructs(dummy, Substruct, Replacement, True)[0]
        rhs_options = product_options(dummy, wildcard_products)
        lhs_addition_options = reactant_options(rhs_options, wildcard_reactants)
        lhs_addition_options.append(input_mol)
        input_index = lhs_addition_options.index(input_mol)
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_addition_options.append(mol)
            elif isinstance (lhs_required, Chem.Mol):
                    lhs_addition_options.append(lhs_required)
                
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    append_unique(rhs_options, mol)
            elif isinstance (rhs_required, Chem.Mol):
                append_unique(rhs_options, rhs_required)
        
        
        x_vars = [LpVariable(f"x_{i}", 0, None, LpInteger) for i in range(len(rhs_options))]
        y_vars = [LpVariable(f"y_{j}", 0, None, LpInteger) for j in range(len(lhs_addition_options))]
        prob += y_vars[input_index] == 1       
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_r_idx = lhs_addition_options.index(mol) 
                    prob += y_vars[lhs_r_idx] == 1
            elif isinstance (lhs_required, Chem.Mol):
                lhs_r_idx = lhs_addition_options.index(lhs_required) 
                prob += y_vars[lhs_r_idx] == 1
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    rhs_r_idx = rhs_options.index(mol) 
                    prob += x_vars[rhs_r_idx] == 1
            elif isinstance (rhs_required, Chem.Mol):
                rhs_r_idx = rhs_options.index(rhs_required) 
                prob += x_vars[rhs_r_idx] == 1
            
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], [], "Infeasible"
    
        bond_types_rhs = [count_bonds(mol) for mol in rhs_options]
        bond_types_lhs = [count_bonds(mol) for mol in lhs_addition_options]
        hydrogen_rhs = [count_hydrogens(mol) for mol in rhs_options]
        hydrogen_lhs = [count_hydrogens(mol) for mol in lhs_addition_options]
        for k in set().union(*bond_types_rhs, *bond_types_lhs):
            prob += lpSum(x_vars[i] * bond_types_rhs[i].get(k, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * bond_types_lhs[j].get(k, 0) for j in range(len(lhs_addition_options)))
    
        for t in set().union(*hydrogen_rhs, *hydrogen_lhs):
            prob += lpSum(x_vars[i] * hydrogen_rhs[i].get(t, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * hydrogen_lhs[j].get(t, 0) for j in range(len(lhs_addition_options)))
        
        _finish_model(prob, x_vars, y_vars, rhs_options, input_mol, rhs_required)
        solver = make_solver()
        prob.solve(solver)
        if LpStatus[prob.status] != 'Optimal':
            return [], [], LpStatus[prob.status]
        #prob.solve()
        solution_x = [int(var.varValue) for var in x_vars]
        solution_y = [int(var.varValue) for var in y_vars]
    
        balanced_rhs = [(rhs_options[i], solution_x[i]) for i in range(len(rhs_options)) if solution_x[i] > 0]
        balanced_lhs = [(lhs_addition_options[j], solution_y[j]) for j in range(len(lhs_addition_options)) if solution_y[j] > 0]
        #print("Solver status:", LpStatus[prob.status])
        return balanced_rhs, balanced_lhs, LpStatus[prob.status]
    return Balance(x)


def isogyric_count(mol):
    atom_counts = {}
    if isinstance(mol, Chem.Atom):
        key = (mol.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    for atom in mol.GetAtoms():

        key = (atom.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
    return atom_counts

def isodesmic_count(mol):
    bond_counts = {}
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        atom_pair = tuple(sorted([atom1.GetSymbol(), atom2.GetSymbol()]))

        bond_type = (atom_pair, bond.GetBondType())
        bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1


    atom_counts = {}
    if isinstance(mol, Chem.Atom):
        key = (mol.GetSymbol())
        atom_counts[key] = atom_counts.get(key, 0) + 1
        return atom_counts
    for atom in mol.GetAtoms():
        if atom.GetSymbol() !='H':
            key = (atom.GetSymbol())
            atom_counts[key] = atom_counts.get(key, 0) + 1
    return bond_counts,atom_counts

def hypohomodesmotic_count(mol):
    mol = Chem.AddHs(mol)
    hydrogen_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetNumExplicitHs())
        hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
    mol = Chem.RemoveHs(mol)
    atom_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetHybridization())
        atom_counts[key] = atom_counts.get(key, 0) + 1
    return hydrogen_counts, atom_counts

def homodesmotic_count(mol):
    bond_counts = {}
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol()== 'H' or bond.GetEndAtom().GetSymbol()== 'H':
            continue 
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        atom_pair = tuple(sorted(
            [(atom1.GetSymbol(), atom1.GetHybridization()), 
             (atom2.GetSymbol(), atom2.GetHybridization())]
        ))

        bond_type = (atom_pair, bond.GetBondType())
        bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1


    hydrogen_counts = {}
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    for atom in mol.GetAtoms():
        key = (atom.GetSymbol(), atom.GetHybridization(), atom.GetNumExplicitHs())
        hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
    return bond_counts, hydrogen_counts

def display_reaction_counts(input_mol, reaction_fn):
    def format_key(k):
        if isinstance(k, tuple):
            return "(" + ", ".join(format_key(i) for i in k) + ")"
        if hasattr(k, 'name'):  # for RDKit enums like HybridizationType
            return k.name
        return str(k)

    def print_dict(title, d):
        if not d:
            return
        print(f"\n{title}:")
        print("-" * len(title))
        for key, val in sorted(d.items()):
            print(f"{format_key(key):<60} : {val}")    

    if reaction_fn is Isogyric:
        atom_counts = isogyric_count(input_mol)
        print_dict("Isogyric Atom Counts", atom_counts)

    elif reaction_fn is Isodesmic:
        bond_counts, atom_counts = isodesmic_count(input_mol)
        print_dict("Isodesmic Bond Counts", bond_counts)
        print_dict("Isodesmic Atom Counts", atom_counts)

    elif reaction_fn is Hypohomodesmotic:
        hydrogen_counts, atom_counts = hypohomodesmotic_count(input_mol)
        print_dict("Hypohomodesmotic Hydrogen Counts", hydrogen_counts)
        print_dict("Hypohomodesmotic Atom Counts", atom_counts)

    elif reaction_fn is Homodesmotic:
        bond_counts, hydrogen_counts = homodesmotic_count(input_mol)
        print_dict("Homodesmotic Bond Counts", bond_counts)
        print_dict("Homodesmotic Hydrogen Counts", hydrogen_counts)


def geom_from_rdkit(rdkitmol):
    """
    Takes an RDKit molecule (already embedded) and returns an AaronTools Geometry along with
    the weighted adjacency matrix
    """
    result = AllChem.EmbedMolecule(rdkitmol)
    if result != 0:
        raise ValueError("Embedding failed")
        
    atom_list = []
    for i, atom in enumerate(rdkitmol.GetAtoms()):
        positions = rdkitmol.GetConformer().GetAtomPosition(i)
        atom_list.append(Atom(element=atom.GetSymbol(), coords=positions))
    return Geometry(atom_list)
