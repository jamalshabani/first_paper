
import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3", "rm -rf test4", "rm -rf test5",
                "python3 directFirstComputation.py -tao_monitor -tao_max_funcs 10000 -tao_max_it 5000 -tao_view -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_ls_type unit -er 1.0e0 -es 1.0e1 -lr 0.1275 -ls 0.12 -vr 0.3 -vs 0.3 -k 6.0e-4 -e 2.0e-3"]




i = 1
for program in program_list:
    print("------------------------------------------------------------------------------")
    print("")
    print("Running test #{}".format(i))
    print("")
    print(program)
    print("")
    subprocess.run(program, shell = True)
    i = i + 1
