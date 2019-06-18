#! /usr/bin/python

import os

def is_within_tolerence(test_value, target_value, tolerence):
    return abs(test_value - target_value) < tolerence


command             = "./UnitedAtom --config-file=unittests/test1.config 2>&1"
expected_text       = "Info: Processing '1n5u' (R = 50)\nInfo: Saving map to: unittests/results/1n5u_50_0.map"
simple_average      = -26.31246
boltzman_average    = -163.42631
standard_deviation  = 2.77865
current_dir         = os.path.dirname(__file__)
root_dir            = os.path.join(current_dir, '..')

os.chdir(root_dir)
results = os.popen(command).read()
os.chdir(current_dir)

result_text = "\n".join(results.split('\n')[0:2])
assert result_text == expected_text, "Test 1 Failed.\nResults text did not match expected output\n---------- Expected ----------\n{}\n---------- Recieved ----------\n{}\n------------------------------".format(expected_text, result_text)

name, size, sa, ba, sd = filter(None, results.split('\n')[2].split(' '))

tolerence = 1.0
assert is_within_tolerence(float(sa), simple_average, tolerence), "Test 1 Failed.\nSimple average result was not within tolerence\nabs( {} - {} ) !< {}".format(sa, simple_average, tolerence)

tolerence = 5.0
assert is_within_tolerence(float(ba), boltzman_average, tolerence), "Test 1 Failed.\nBoltzmann average result was not within tolerence\nabs( {} - {} ) !< {}".format(sa, simple_average, tolerence)


print("Test 1 passed!")
