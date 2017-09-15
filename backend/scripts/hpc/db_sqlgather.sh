#!/bin/bash
#$ -N db_sqlgather
#$ -l h_rt=72:00:00
#$ -l h_vmem=20G
#$ -M j.c.wolthuis-2@umcutrecht.nl
#$ -m beas

# -----

cd /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db
TMPDIR=/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/tmp

# OTHER OPTION - not really faster but more direct
# echo "ATTACH '/hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db/pubchem.base.db' AS tmp;" >> imports.sqlite
# echo "CREATE TABLE IF NOT EXISTS base AS SELECT * FROM tmp.base;" >> imports.sqlite
# echo "CREATE INDEX IF NOT EXISTS b_idx1 on base(baseformula, charge);" >> indexes.sqlite
# 
# echo ".mode csv" >> imports.sqlite;
# echo '.separator " "' >> imports.sqlite;
# echo "CREATE TABLE IF NOT EXISTS extended(baseformula text,
# fullformula text,
# basemz decimal(30,13),
# fullmz decimal(30,13),
# adduct text,
# basecharge int,
# totalcharge int,
# isoprevalence float,
# foundinmode text);" >> imports.sqlite
# 
# for f in /hpc/cog_bioinf/ridder/users/jwolthuis/MetaboShiny/backend/db/pubchem_csv_ext/*.csv; do t=`basename $f .csv`; \
# echo ".import $f extended" >> imports.sqlite; done

echo "create index e_idx1 on extended(baseformula, basecharge);"  >> indexes.sqlite
echo "create index e_idx2 on extended(fullmz, foundinmode);"  >> indexes.sqlite

#     do the thing
#sqlite3 pubchem.full.db < imports.sqlite
sqlite3 pubchem.full.db < indexes.sqlite