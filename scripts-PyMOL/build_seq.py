#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 

# additions:
# Copyright (c) 2005 Robert L. Campbell
from pymol import cmd,stored
#import seq_convert

QuietException = parsing.QuietException
tmp_editor = "_tmp_editor"
tmp_ed_save = "_tmp_ed_save"
tpk1 = "_tmp_tpk1"
tpk2 = "_tmp_tpk2"

# attach_amino_acid is from pymol/modules/pymol/editor.py 
def attach_amino_acid(selection,amino_acid,phi,psi):

  if not selection in cmd.get_names("selections"):
    if amino_acid in cmd.get_names("objects"):
      print " Error: an object with than name already exists"
      raise QuietException
    cmd.fragment(amino_acid)
    if cmd.get_setting_legacy("auto_remove_hydrogens"):
      cmd.remove("(hydro and %s)"%amino_acid)
    if cmd.count_atoms("((%s) and name c)"%amino_acid,quiet=1):
      cmd.edit("((%s) and name c)"%amino_acid)
  else:
    cmd.fragment(amino_acid,tmp_editor)
    if cmd.count_atoms("((%s) and elem n)"%selection,quiet=1):
      cmd.select(tmp_ed_save,"(%s)"%selection)
      cmd.iterate("(%s)"%selection,"stored.resv=resv")
      stored.resi = str(stored.resv-1)
      cmd.alter(tmp_editor,"resi=stored.resi")
      cmd.fuse("(%s and name C)"%(tmp_editor),"(pk1)",2)
      if cmd.get_setting_legacy("auto_remove_hydrogens"):
        cmd.remove("(pkmol and hydro)")
      cmd.set_dihedral("(name ca and neighbor pk2)",
                            "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
      cmd.set_geometry("pk2",3,3) # make nitrogen planer
      cmd.select(tpk1,"pk2")
      cmd.select(tpk2,"pk1")
      if amino_acid[0:3]!='pro':
        cmd.set_dihedral( # PHI
            "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
            "(name ca and neighbor "+tpk1+")", # CA 
            tpk1, # N
            tpk2, # C
            phi)
      cmd.set_dihedral( # PSI (n-1)
          tpk1, # N
          tpk2, # C
          "(name ca and neighbor "+tpk2+")", # CA
          "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
          psi)
      cmd.delete(tpk1)
      cmd.delete(tpk2)
      sele = ("(name N and (byres neighbor %s) and not (byres %s))"%
                (tmp_ed_save,tmp_ed_save))
      if cmd.count_atoms(sele,quiet=1):
        cmd.edit(sele)
      cmd.delete(tmp_ed_save)
                
    elif cmd.count_atoms("((%s) and elem c)"%selection,quiet=1):
      cmd.select(tmp_ed_save,"(%s)"%selection)
      cmd.iterate("(%s)"%selection,"stored.resv=resv")
      stored.resi = str(stored.resv+1)
      cmd.alter(tmp_editor,"resi=stored.resi")
      cmd.fuse("(%s and name N)"%(tmp_editor),"(pk1)",2)
      if cmd.get_setting_legacy("auto_remove_hydrogens"):
        cmd.remove("(pkmol and hydro)")
      cmd.set_dihedral("(name ca and neighbor pk2)",
                            "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
      cmd.set_geometry("pk1",3,3) # make nitrogen planar
      cmd.select(tpk1,"pk1")
      cmd.select(tpk2,"pk2")
      if amino_acid[0:3]!='pro':
        cmd.set_dihedral( # PHI
            tpk2, # C
            tpk1, # N
            "(name ca and neighbor "+tpk1+")", # CA 
            "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
            phi)
      cmd.set_dihedral( # PSI (n-1)
          "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
          "(name ca and neighbor "+tpk2+")", # CA
          tpk2, # C
          tpk1, # N
          psi)
      cmd.delete(tpk1)
      cmd.delete(tpk2)
      sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                (tmp_ed_save,tmp_ed_save))
      if cmd.count_atoms(sele,quiet=1):
        cmd.edit(sele)
      cmd.delete(tmp_ed_save)
    elif cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
      print " Error: please pick a nitrogen or carbonyl carbon to grow from."
      cmd.delete(tmp_editor)
      raise QuietException
  
  cmd.delete(tmp_editor)
        

def seq1_to_seq3(seq1):
  """
  Convert array of 1-letter code sequence to 3-letter code.
  seq1 can be an array or a string.
  Return seq3 as array.
  """
  res3 = ['---','ace','ala','asn','asp','arg','cys','gln','glu',
           'gly','his','ile','leu','lys','met','pro','phe','ser',
           'thr','trp','tyr','val','unk',
           'ALA','ACE','ASN','ASP','ARG','CYS','GLN','GLU',
           'GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER',
           'THR','TRP','TYR','VAL','UNK']
  res1 = '-bandrcqeghilkmpfstwyvxABNDRCQEGHILKMPFSTWYVX'

  if len(seq1) > 1:
    seq3 = []
    for a1 in seq1:
      try:
        a3 = res3[res1.index(a1)]
        seq3.append(a3)
      except ValueError, err:
        print "%s:  No match for residue: %s" % (err,a1)
  else:
    seq3 = res3[res1.index(seq1)]

  return seq3

  
def build_seq(name, seq,ss=None,phi=None,psi=None):
  """
  usage: build_seq name, seq [, ss=helix | phi=phi, psi=psi]
  example: build_seq peptide, QGAADLESLGQYFEEMKTKLIQDMTE, ss=helix 

  will build the above sequence in a helical conformation
  ss can be: helix or h or alpha  (phi=-57, psi=-47)
             antiparallel or beta (phi=-139, psi=-135)
             parallel             (phi=-119, psi=113)
             3/10 helix           (phi=-40.7,psi=30)
             polypro              (phi=-78, psi=149) (polyproline helix type II)

  Alternatively, you can specify the phi and psi angles directly:

  build_seq peptide, QGAADLESLGQ, phi=-60, psi=-40

  This will create an object with the name you provided, unless
  unless that selection already exists, then it will build onto that.

  """

  if name in cmd.get_names("selections"):
    obj='pk1'
  else:
    obj=seq[0:3]

  if phi == None or psi == None:
    if ss == None:
      phi = -139
      psi = -135
    else:
      ss = ss.lower()
      if ss[0:5] == 'helix' or ss[0] == 'h' or ss[0:5] == 'alpha':
        phi = -57
        psi = -47
      elif ss[0:5] == 'antip' or ss[0:4] == 'beta':
        phi = -139
        psi = -135
      elif ss[0:5] == 'paral':
        phi = -119
        psi = 113
      elif ss[0:4] == '3/10':
        phi = -40.7
        psi = -30
      elif ss[0:7] == 'polypro':
        phi = -78
        psi = 149

  print "Building sequence: ",seq

  ### Here seq_convert is used ###
#  seq3=seq_convert.seq1_to_seq3(seq)
  seq3=seq1_to_seq3(seq)
  attach_amino_acid(obj,seq3[0].lower(),phi,psi)
  for aa in seq3[1:]:
    aa = aa.lower()
    attach_amino_acid('pk1',aa,phi,psi)
  cmd.delete('pk1')
  cmd.delete('pkmol')
  cmd.set_name(seq3[0].lower(), name)


cmd.extend('build_seq',build_seq)
