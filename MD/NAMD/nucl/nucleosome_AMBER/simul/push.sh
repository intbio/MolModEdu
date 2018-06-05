#!/bin/bash
#helix
#ssh helix_anya 'mkdir -p /data/panch/nucleosome/1kx5nt_ge_cons_simul/'
#scp  * helix_anya:/data/panch/nucleosome/1kx5nt_ge_cons_simul/
#scp -r input helix_anya:/data/panch/nucleosome/1kx5nt_ge_cons_simul/

#lomonosov
ssh lomo-ak 'mkdir -p ~/nucleosome_amber/1kx5nt_simul/'
scp  * lomo-ak:~/nucleosome_amber/1kx5nt_simul/
scp -r input lomo-ak:~/nucleosome_amber/1kx5nt_simul/
