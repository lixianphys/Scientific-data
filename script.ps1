python toybands/run.py -simu -loadden output/density.csv -fnm simu_ --enrange -0.2 0.1 100 --bfrange 1 7 240 -nmax 15
python toybands/run.py -dosm -loadden output/density.csv -fnm dos --enrange -0.2 0.1 100 --bfrange 1 7 240  -nmax 15
# python toybands/run.py -simu --allden "3e15 1e15 2e15 2e15 0 6e15" -nos 60 -fnm simu_3bands_gf1p8_vf4e5_6 --enrange -0.2 0.1 30 --bfrange 1 7 30 -nmax 5
# python toybands/run.py -denplot -fnm den_12 --enrange -0.2 0.1 100 --bfrange 1 7 120 -nmax 10
# python toybands/run.py -enplot -fnm en_8 --enrange -0.2 0.1 100 --bfrange 1 7 120 -nmax 10
