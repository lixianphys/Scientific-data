python toybands/delband.py -i all
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 403365 5.2
python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 403365 5.2
python toybands/addband.py 1e15 -gfactor 6 -cp 0.2 1
python toybands/addband.py 1e15 -gfactor 6 -cp 0.2 -1
python toybands/run.py -simu -loadden density4band.csv -fnm test_noparity --enrange -0.05 0.1 20 --bfrange 1 7 30  -nmax 5
python toybands/run.py -simu -loadden density4band.csv -fnm test_wparity --enrange -0.05 0.1 20 --bfrange 1 7 30  -nmax 5 -parity
# python toybands/run.py -dosm -loadden density4band.csv -fnm dosm_noparity --enrange -0.05 0.1 50 --bfrange 1 7 240  -nmax 10
# python toybands/run.py -dosm -loadden density4band.csv -fnm dosm_wparity --enrange -0.05 0.1 50 --bfrange 1 7 240  -nmax 10 -parity

# python toybands/delband.py -i all
# python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
# python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
# python toybands/addband.py 1e15 -gfactor 10 -cp 0.2 1
# python toybands/run.py -simu -loadden output/density.csv -fnm simu_gvps10 --enrange -0.05 0.1 50 --bfrange 1 7 120  -nmax 10
# python toybands/delband.py -i all
# python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
# python toybands/addband.py 1e15 -is_dirac -is_cond -gfactor 1.8 -dp 0 403365
# python toybands/addband.py 1e15 -gfactor 20 -cp 0.2 1
# python toybands/run.py -simu -loadden density4band.csv -fnm simu_4band --enrange -0.05 0.1 50 --bfrange 1 7 240  -nmax 10

