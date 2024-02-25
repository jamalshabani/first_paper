
import subprocess

program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3", "rm -rf test4", "rm -rf test5",
                "python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 500 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test1' -er 20e-2 -es 200e-1 -lr 0.05 -ls 0.5 -vr 0.4 -vs 0.4 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0",
                "python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 500 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test2' -er 25e-2 -es 250e-1 -lr 0.05 -ls 0.5 -vr 0.4 -vs 0.4 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0",
                "python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 500 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test3' -er 13e-2 -es 130e-1 -lr 0.05 -ls 0.5 -vr 0.4 -vs 0.4 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0",
                "python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 500 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test4' -er 14e-2 -es 140e-1 -lr 0.05 -ls 0.5 -vr 0.4 -vs 0.4 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0",
                "python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 500 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test5' -er 15e-2 -es 150e-1 -lr 0.05 -ls 0.5 -vr 0.4 -vs 0.4 -k 1.0e-4 -e 4.0e-3 -p 2.0 -q 1.0"]



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




# Scripts for ratio E_r/E_s = 1.0
#"python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 100 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test2' -er 40e-1 -es 40e-1 -lr 2.0 -ls 1.0 -vr 0.4 -vs 0.4 -k 5.0e-4 -e 4.0e-3 -p 2.0 -q 1.0",
#"python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 100 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test3' -er 50e-1 -es 50e-1 -lr 2.0 -ls 1.0 -vr 0.4 -vs 0.4 -k 1.0e-3 -e 4.0e-3 -p 2.0 -q 1.0",
#"python3 3_mat_compliance_responsive.py -tao_type bncg -tao_max_funcs 10000 -tao_gatol 1.0e-7 -tao_grtol 1.0e-7 -tao_gttol 1.0e-7 -tao_converged_reason -tao_monitor -tao_max_it 100 -tao_ls_type armijo -m '1_to_3_mesh.msh' -o 'test4' -er 60e-1 -es 60e-1 -lr 2.0 -ls 1.0 -vr 0.4 -vs 0.4 -k 1.0e-3 -e 4.0e-3 -p 2.0 -q 1.0",