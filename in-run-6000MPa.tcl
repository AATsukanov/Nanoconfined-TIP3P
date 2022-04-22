processors	2 3 4
boundary	p p p

units           real
neigh_modify    delay 1 every 1 check yes

atom_style      full
bond_style      harmonic
angle_style     charmm
dihedral_style  charmm
improper_style  harmonic
special_bonds   charmm
pair_style      lj/charmm/coul/long 8 10
pair_modify     mix arithmetic
kspace_style    pppm 1e-4

variable prefix index tip3p-1008
variable pressureMPa equal 6000
variable temperature equal 310
variable XY equal (33.93*39.18)
variable nPistonAtoms equal (108*4)
variable zforce equal -(v_pressureMPa/v_nPistonAtoms*v_XY*6.022/4.184*1.0e-4)
variable id index ${prefix}_${pressureMPa}MPa_${temperature}K.x24
shell	mkdir ${id}
shell	cd ${id}
read_data ../tip3p-252.data
replicate 2 2 1
#       1      1.008  # HT
#       2    15.9994  # OT
#       3     55.847  # FEX
#       4    15.9994  # OHX
#       5      1.008  # HOX
#set type 3 charge 0.105
#set type 4 charge -0.095
#set type 5 charge 0.0425
#make regions:
region BEPX block INF INF INF INF 10 100
region Hu3 block INF INF INF INF -100 10
group all_BEPX	region BEPX
group all_Hu3	region Hu3
#make groups:
group   h2o type 2
group	frozen	type 3 4
group	mineral type 3:5
group   substrate intersect all_Hu3 frozen
group   piston intersect all_BEPX frozen
group	mobile	subtract all frozen
#em1:
fix	em1 mineral setforce 0.0 0.0 0.0
minimize 0.0 0.0 10000 100000
unfix em1
#em2:
fix	em2 frozen setforce 0.0 0.0 0.0
minimize 0.0 0.0 10000 100000
unfix em2
write_data tip3p-1008.data
#md:
velocity mobile create ${temperature} 25122020 mom yes rot yes dist gaussian
fix md1 substrate setforce 0.0 0.0 0.0
fix md2 mobile nvt temp ${temperature} ${temperature} 100
#pressure:
compute 1 piston com
fix 1 piston rigid single torque * off off off force * off off on
fix 2 piston addforce 0.0 0.0 ${zforce}
#MSD:
compute 2 h2o msd com yes
print '(step) (msd) #the total squared displacements, A^2' file MSD.${id}.txt
fix 3 all print 100000 "$(step) $(c_2[4])" append MSD.${id}.txt
compute 3 h2o vacf
fix     4 h2o vector 1 c_3[4]
variable  diffusion equal dt*trap(f_4)
#vacf: (A/fs)^2, diff_h2o: A^2/fs
fix     5 h2o print 10000 "$(step) $(c_3[4]) $(v_diffusion)" append VACF.${id}.txt
#RDF:                      Ow-Hw Ow-Ow Hw-OH Ow-OH
compute rdfall all rdf 100  1 2   2 2   1 4  2  5
fix rdfall all ave/time 5000 100 500000 c_rdfall[*] ave running overwrite overwrite file rdf.${id}.txt mode vector
dump 1 all dcd 300000 output.${id}.uw.dcd
dump_modify 1 unwrap yes
dump 2 all xyz 1000000 frames.${id}.wrapped.xyz
fix p1 all print 100000 "$(step) $(temp) $(c_1[1]) $(c_1[2]) $(c_1[3])" append energy.${id}.txt title '(step) (temp) (piston_com_x_y_z)'
fix p2 all print 100000 "$(step) $(c_1[3])" append zvolume.${id}.txt title '(step) (piston_com_z)'
timestep 1.0
run 5000000
write_data ${id}.*.data