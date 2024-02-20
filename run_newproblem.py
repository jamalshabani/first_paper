import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 newproblem.py -n 500 -o 'test1' -er 1.0 -es 1.0e-2 -vr 0.3 -vs 0.3 -k 2.5e-3 -e 2.0e-2 -p 2.0 -q 2.0"]



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

# "python3 leastsquare.py -n 500 -o 'test2' -er 1.0 -es 1.0e-2 -vr 0.3 -vs 0.3 -k 5.0e-4 -e 5.0e-3 -p 2.0 -q 2.0"]
# "python3 leastsquare.py -n 10000 -o 'test1' -er 1.0 -es 1.0e-2 -vr 0.3 -vs 0.3 -k 5.0e-3 -e 2.5e-3 -p 2.0 -q 2.0"]

