#!/bin/bash
#helix
ssh helix 'mkdir -p /data/shaytana/1kx5_link_simul_cons_zw/'
scp  * helix:/data/shaytana/1kx5_link_simul_cons_zw/
scp -r input helix:/data/shaytana/1kx5_link_simul_cons_zw/

#lomonosov
# ssh lomo-ak 'mkdir -p ~/nucleosome/1kx5_simul_link/'
# scp  * lomo-ak:~/nucleosome/1kx5_simul_link/
# scp -r input lomo-ak:~/nucleosome/1kx5_simul_link/
