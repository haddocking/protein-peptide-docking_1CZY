#!/bin/bash
"
README:

Run this script in the folder containing you peptides and the *.mdp files.
Change the essential parameters in the "Basic Definitions" section.
BEWARE: Changing the rest of the script would alter the protocol.

"

#Basic Definitions:
min_dist=$1
force_field="amber99sb-ildn"
water_model="tip3p"
declare -a structures
pep_cnfs=(alpha polypro extended) #names of peptide conformations
parallel="-nt 48 -ntmpi 12"

#Error Handling:
if [[ ${min_dist} == "" ]]; then
    echo -e "\e[31m ERROR: Plese provide a value for minimum distance from the box for the longest peptide (in nm) \e[0m"
    echo -e "\e[32m     HINT: Make sure the distance is at least = half the length of the longest peptide\e[0m"
    read var
    if [[ ${var} == "" ]]; then
        echo -e "\e[31m ERROR Again, will exit now... \e[0m"
        exit
    fi
    min_dist=$var
fi
i=0
while [ $i -le 2 ]; do
    if [[ $(ls | grep "${pep_cnfs[${i}]}" | grep "pdb" -c) == 0 ]]; then
        echo -e "\e[31m ERROR: Please provide ${pep_cnfs[${i}]}.pdb \e[0m"
        exit
    fi
    ((i++))
done

#Calculating box size for extended peptide (longest of the three)
gmx pdb2gmx -f ${pep_cnfs[2]}.pdb -o ${pep_cnfs[2]}.gro -ignh -ff ${force_field} -water ${water_model} -ter
box_vectors=$(gmx editconf -f ${pep_cnfs[2]}.pdb -o ${pep_cnfs[2]}_pbc.gro -bt dodecahedron -d ${min_dist} | grep "new box vectors" | awk '{print $6}')

#Calibrate the number of water molecules
declare -a sol
i=0
while [ $i -le 2 ]; do
    gmx pdb2gmx -f ${pep_cnfs[${i}]}.pdb -o ${pep_cnfs[${i}]}.gro -ignh -ff ${force_field} -water ${water_model}
    gmx editconf -f ${pep_cnfs[${i}]}.pdb -o ${pep_cnfs[${i}]}_pbc.gro -bt dodecahedron -box ${box_vectors}
    gmx grompp -f vacuum.mdp -c ${pep_cnfs[${i}]}_pbc.gro -p topol.top -o ${pep_cnfs[${i}]}_vac.tpr -maxwarn 3
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_vac
    gmx solvate -cp ${pep_cnfs[${i}]}_vac.gro -cs spc216.gro -p topol.top -o ${pep_cnfs[${i}]}_solvated.gro
    grep SOL topol.top | awk '{print $2}' >> water.info
    ((i++))
done
n_water=$(sort water.info | head -1)
rm water.info
mkdir mdouts
i=0
while [ $i -le 2 ]; do
    gmx pdb2gmx -f ${pep_cnfs[${i}]}.pdb -o ${pep_cnfs[${i}]}.gro -ignh -ff ${force_field} -water ${water_model}
    gmx editconf -f ${pep_cnfs[${i}]}.pdb -o ${pep_cnfs[${i}]}_pbc.gro -bt dodecahedron -box ${box_vectors}
    gmx grompp -f vacuum.mdp -c ${pep_cnfs[${i}]}_pbc.gro -p topol.top -o ${pep_cnfs[${i}]}_vac.tpr -maxwarn 3
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_vac $parallel
    gmx solvate -cp ${pep_cnfs[${i}]}_vac.gro -cs spc216.gro -p topol.top -o ${pep_cnfs[${i}]}_solvated.gro -maxsol $n_water
    mv "#topol.top.1#" "${pep_cnfs[${i}]}_topol_after_vac_em.top"
    gmx grompp -f ions.mdp -c ${pep_cnfs[${i}]}_solvated.gro -p topol.top -o ${pep_cnfs[${i}]}_pp_ionize.tpr
    gmx genion -s ${pep_cnfs[${i}]}_pp_ionize.tpr -p topol.top -o ${pep_cnfs[${i}]}_ionized.gro -neutral
    mv "#topol.top.1#" "${pep_cnfs[${i}]}_topol_after_solvation.top"
    gmx grompp -f ions.mdp -c ${pep_cnfs[${i}]}_ionized.gro -p topol.top -o ${pep_cnfs[${i}]}_neutral_relaxed.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_neutral_relaxed $parallel
    gmx grompp -f nvt.mdp -c ${pep_cnfs[${i}]}_neutral_relaxed.gro -p topol.top -o ${pep_cnfs[${i}]}_nvt.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_nvt $parallel
    gmx grompp -f npt.mdp -c ${pep_cnfs[${i}]}_nvt.gro -p topol.top -o ${pep_cnfs[${i}]}_npt.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_npt $parallel
    sed -e 's/1000  1000  1000/ 100   100   100/g' posre.itp > tmp.itp
    gmx grompp -f npt.mdp -c ${pep_cnfs[${i}]}_npt.gro -p topol.top -o ${pep_cnfs[${i}]}_npt_progrel100.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_npt_progrel100 $parallel
    sed -e 's/100   100   100/ 10    10    10/g' posre.itp > tmp.itp
    gmx grompp -f npt.mdp -c ${pep_cnfs[${i}]}_npt_progrel100.gro -p topol.top -o ${pep_cnfs[${i}]}_npt_progrel10.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_npt_progrel10 $parallel
    gmx grompp -f unrestrained.mdp -c ${pep_cnfs[${i}]}_npt_progrel10.gro -p topol.top -o ${pep_cnfs[${i}]}_all_set.tpr
    gmx mdrun -v -deffnm ${pep_cnfs[${i}]}_all_set $parallel
    gmx grompp -f production.mdp -c ${pep_cnfs[${i}]}_all_set.gro -p topol.top -o ${pep_cnfs[${i}]}_md.tpr
    mv topol.top ${pep_cnfs[${i}]}_topol_final.top
    mv ./"#mdout"* ./mdouts
    mv ./mdout.mdp ./mdouts
    rm \#*
    ((i++))
done

echo "Run the Production MD using:"
echo ">> gmx mdrun -v -deffnm alpha_md"
echo ">> gmx mdrun -v -deffnm polypro_md"
echo ">> gmx mdrun -v -deffnm extended_md"
