#!/bin/bash
#Dihedral PCA GROMACS
#References: 1) http://www.gromacs.org/Documentation/How-tos/Dihedral_PCA?highlight=PCA
#            2) http://www.ncbi.nlm.nih.gov/pubmed/15521057
#Preps: 
#     1) dangle.ndx: Index of atoms making phi & psi angles
#     2) covar.ndx: Index for Covar, should contain one group
#                   which has integers ranging from 1 to int(2*N/3)
#                   where N is the number of dihedrals

echo -e "\e[31m Please ensure that you have completed the steps 1-7 in the protocol (section 3.3) \e[0m"
echo -e "\e[31m If you also want to extract the first 30 structures: \e[0m"
echo -e "\e[31m Also provide a file- peptide.ndx containing ONLY ONE index group that corresponds to all atoms in the proteins_uncapped.tpr topology \e[0m"
echo "Hit ENTER/RETURN to Continue or Ctrl-C to exit"
read var
if [[ ${var} == "" ]]; then
    gmx angle -f protein_bb.xtc -n dangle.ndx -or dangle.trr -type dihedral #Extract angles
    gmx trjconv -f dangle.trr -s protein_bb.tpr -o resized.gro -n covar.ndx -e 500 #Make Reference Structure For Covar
    gmx covar -f dangle.trr -n covar.ndx -ascii -xpm -nofit -nomwa -noref -nopbc -s resized.gro #Make Covar
    gmx anaeig -v eigenvec.trr -f dangle.trr -s resized.gro -first 1 -last 2 -2d 2dproj_1_2.xvg
    gmx sham -f 2dproj_1_2.xvg -notime -bin bindex-1_2.ndx -lp prob-1_2.xpm -ls gibbs-1_2.xpm -lsh enthalpy-1_2.xpm -g shamlog-1_2 -lss entropy-1_2.xpm
fi
echo -e "\e[32 D-PCA is finished. Do you now wish to extract 30 structures one from each of the biggest states? \e[0m"
echo "Hit ENTER/RETURN to Continue or Ctrl-C to exit"
read var

if [[ ${var} == "" ]]; then
    no_states=$(cat shamlog-$1_$2.log | grep -c ".") 
    let no_states=(states-1)/2
    i=1
    j=30
    while [ $i -le $j ]; do
        ind=$(cat shamlog-1_2.log | tail -$j | awk '{print $5}' | awk "{if (NR==$i) print \$1 }") #Seek binning index from shamlog-*-*.log
        structure=$(grep -A1 "\[ $ind \]" bindex-1_2.ndx | tail -1) #Seek the timestamp for the first structure in the bin
        let frame=structure*50
        gmx trjconv -f protein_uncapped.xtc -s protein_uncapped.tpr -dump $frame.pdb -o mintemp-$i.xtc -pbc mol -n peptide.ndx
        let i=i+1
    done
