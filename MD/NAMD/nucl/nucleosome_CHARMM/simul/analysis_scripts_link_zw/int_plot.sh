#!/bin/bash
/usr/bin/R --vanilla --slave < int_plot_dna_prot.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_anya_talk.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_stat_Z.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_sym_prof.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_var.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_Z.r
/usr/bin/R --vanilla --slave < int_plot_dna_prot_Z_s2.r

/usr/bin/R --vanilla --slave < int_plot_prot_prot.r
/usr/bin/R --vanilla --slave < int_plot_prot_prot_anya_talk.r

/usr/bin/R --vanilla --slave < int_plot_wat.r
/usr/bin/R --vanilla --slave < int_plot_wat_dna.r
/usr/bin/R --vanilla --slave < int_plot_wat_dna_anya_talk.r
/usr/bin/R --vanilla --slave < int_plot_wat_sym_prof.r

/usr/bin/R --vanilla --slave < int_plot_ions.r
/usr/bin/R --vanilla --slave < int_plot_ions_sym_prof.r


