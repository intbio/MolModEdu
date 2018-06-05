set alpha_hfolds_ca "((segname CHA CHE and (resid 63 to 77 or resid 85 to 114 or resid 120 to 131)) or (segname CHB CHF and                                 (resid 30 to 41 or resid 49 to 76 or resid 82 to 93))  or                   (segname CHC CHG and (resid 27 to 37 or resid 46 to 73 or resid 79 to 89)) or                   (segname CHD CHH and (resid 34 to 46 or resid 52 to 81 or resid 87 to 99))) and name CA"
set alpha_hfolds_bb "((segname CHA CHE and (resid 63 to 77 or resid 85 to 114 or resid 120 to 131)) or (segname CHB CHF and                                 (resid 30 to 41 or resid 49 to 76 or resid 82 to 93))  or                   (segname CHC CHG and (resid 27 to 37 or resid 46 to 73 or resid 79 to 89)) or                   (segname CHD CHH and (resid 34 to 46 or resid 52 to 81 or resid 87 to 99))) and name CA C O N"

set alpha_ext_ca "((segname CHA CHE and (resid 44 to 57 or resid 63 to 77 or resid 85 to 114 or resid 120 to 131)) or (segname CHB CHF and (resid 24 to 29 or resid 30 to 41 or resid 49 to 76 or resid 82 to 93))  or (segname CHC CHG and (resid 16 to 22 or resid 27 to 37 or resid 46 to 73 or resid 79 to 89 or resid 90 to 97)) or (segname CHD CHH and (resid 34 to 46 or resid 52 to 81 or resid 87 to 99 or resid 100 to 120))) and name CA"

set core_ca "((segname CHA CHE and (resid 44 to 135)) or (segname CHB CHF and (resid 24 to 102))  or (segname CHC CHG and (resid 16 to 118)) or (segname CHD CHH and (resid 30 to 122))) and name CA"
set core_sch "((segname CHA CHE and (resid 44 to 135)) or (segname CHB CHF and (resid 24 to 102))  or (segname CHC CHG and (resid 16 to 118)) or (segname CHD CHH and (resid 30 to 122))) and not name CA C O N"

# set core_ca "((segname CHA CHE and (resid 37 to 135)) or (segname CHB CHF and (resid 16 to 102))  or (segname CHC CHG and (resid 12 to 118)) or (segname CHD CHH and (resid 21 to 122))) and name CA"
# set core_sch "((segname CHA CHE and (resid 37 to 135)) or (segname CHB CHF and (resid 16 to 102))  or (segname CHC CHG and (resid 12 to 118)) or (segname CHD CHH and (resid 12 to 122))) and not name CA C O N"


set all_ca "((segname CHA CHE ) or (segname CHB CHF )  or (segname CHC CHG ) or (segname CHD CHH )) and name CA"
set all_sch "((segname CHA CHE) or (segname CHB CHF )  or (segname CHC CHG ) or (segname CHD CHH )) and not name CA C O N and noh"
