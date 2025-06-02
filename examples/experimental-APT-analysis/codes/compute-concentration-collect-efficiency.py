file_name = '../example-APT-data/large-cube-GeSn-raw-APT.xyz'
infile = open(file_name,'r')
data = infile.readlines()
infile.close()
N = int(data[0].split()[0])

a = float(data[1].split()[0])
b = float(data[1].split()[4])
c = float(data[1].split()[8])

V = a*b*c
density_0 = N/V

N_Si = 0
N_Ge = 0
N_Sn = 0
for line in data[2:]:
	atom_species = line.split()[0]
	if atom_species == 'Si': N_Si += 1
	if atom_species == 'Ge': N_Ge += 1
	if atom_species == 'Sn': N_Sn += 1

conc_Si = N_Si/N
conc_Ge = N_Ge/N
conc_Sn = N_Sn/N

lc_Ge = 5.6597
lc_Sn = 6.4892 #(Small Methods)
lc_Si = 5.43 #(online)
lc = lc_Si*conc_Si + lc_Ge*conc_Ge + lc_Sn*conc_Sn

#collection efficiency
V_onecell = lc**3
N_onecell = 8
N_max_atoms = V/V_onecell*N_onecell
collect_efficiency = N/N_max_atoms*100

print("Si content (%): ", round(conc_Si*100,1))
print("Ge content (%): ", round(conc_Ge*100,1))
print("Sn content (%): ", round(conc_Sn*100,1))
print("collection efficiency (%): ", round(collect_efficiency))

