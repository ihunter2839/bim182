import sys
import scipy.misc as misc
import matplotlib.pyplot as plt

p_c_iso = .01
p_c = .99
p_o_iso = .05
p_o = .95

mass_c = 12.0107
mass_o = 15.999
mass_n = 14.0067
mass_h = 1.0078
mass_s = 32.065

base_mass = 11*mass_c + 22*mass_h + 3*mass_n + 5*mass_o + mass_s

#I'm going to write this code specifically for SAM
#but the process is valid for any peptide

res = []

#iterate through number of carbon isotopes
for c in range(0,6):
    #iterate through number of oxygen isotopes
    for o in range(0,12):
        mass = base_mass + 1*c + 2*o
        #define success to be isotope incorportation
        c_iso_prob = (misc.comb(5,c) * (p_c_iso ** c) * (p_c ** (5-c)))
        o_iso_prob = (misc.comb(11,o) * (p_o_iso ** o) * (p_o ** (11-o)))
        prob = c_iso_prob * o_iso_prob
        res.append([mass, prob])

res = sorted(res, key=lambda x: x[1], reverse=True)
masses = [x[0] for x in res]
probs = [x[1] for x in res]
masses_1 = masses [0:10]
probs_1 = probs[0:10]
masses_2 = masses[10:]
probs_2 = probs[10:]

print(sum(probs))

plt.bar(masses_1, probs_1)
plt.xlabel("Mass (Da)")
plt.ylabel("Frequency")
plt.show()

plt.bar(masses_2, probs_2)
plt.xlabel("Mass (Da)")
plt.ylabel("Frequency")
plt.show()