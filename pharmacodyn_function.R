### Pharmacodynamic function 
# As in Yu ProcRoySoc 2018

### Make pharmacodynamics function 
phi_max = 1
phi_min = -1
kappa = 1.5 # From (Yu2018): for antibiotics
trate <- seq(0,100,0.1)
vs_mic = 10

#phi = phi_max - ((phi_max - phi_min)*(trate / vs_mic)^kappa) / ((trate / vs_mic)^kappa - phi_min/phi_max)
phi = phi_max * (1 - ((1 - phi_min/phi_max)*(trate / vs_mic)^kappa) / ((trate / vs_mic)^kappa - phi_min/phi_max))

plot(trate, phi)
plot(log10(trate), phi)
