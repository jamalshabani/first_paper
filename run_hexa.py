
import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3", "rm -rf test4", "rm -rf test5",
                "python3 hexagonal.py -tao_type bncg -tao_max_funcs 10000 -tao_monitor -tao_max_it 100 -m 'hexa.msh' -o 'test1' -er 1.0e -es 1.0e-2 -lr 0.15 -ls 0.015 -vr 0.3 -vs 0.3 -k 1.0e-4 -e 1.0e-4",
                "python3 hexagonal.py -tao_type bncg -tao_max_funcs 10000 -tao_monitor -tao_max_it 100 -m 'hexa.msh' -o 'test2' -er 1.0e1 -es 1.0e-1 -lr 0.15 -ls 0.015 -vr 0.3 -vs 0.3 -k 1.0e-4 -e 1.0e-4",
                "python3 hexagonal.py -tao_type bncg -tao_max_funcs 10000 -tao_monitor -tao_max_it 100 -m 'hexa.msh' -o 'test3' -er 1.0e2 -es 1.0e0 -lr 0.15 -ls 0.015 -vr 0.3 -vs 0.3 -k 1.0e-4 -e 1.0e-4"]



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