python toybands/delband.py -i all
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -gfactor 1 -cp 0.2 1
python toybands/run.py -simu -loadden output/density.csv -fnm simu_gvps1 --enrange -0.05 0.1 50 --bfrange 1 7 120  -nmax 10
python toybands/delband.py -i all
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -gfactor 10 -cp 0.2 1
python toybands/run.py -simu -loadden output/density.csv -fnm simu_gvps10 --enrange -0.05 0.1 50 --bfrange 1 7 120  -nmax 10
python toybands/delband.py -i all
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
python toybands/addband.py 1e15 -gfactor 20 -cp 0.2 1
python toybands/run.py -simu -loadden output/density.csv -fnm simu_gvps20 --enrange -0.05 0.1 50 --bfrange 1 7 120  -nmax 10

