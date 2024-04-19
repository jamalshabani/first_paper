import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 test.py -tao_monitor -tao_max_it 500 -o 'test1' -er 1.0e0 -es 1.0e-2 -vr 0.3 -vs 0.3 -q 1.0",
                "python3 test.py -tao_monitor -tao_max_it 500 -o 'test2' -er 1.0e0 -es 1.0e-2 -vr 0.3 -vs 0.3 -q 2.0"]








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
