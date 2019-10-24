import subprocess, os

command = "bc-mozart4-aer general.cfg  mozart4.cfg 20140316 20140615"
print command
os.system(command)

command = "bc-mozart4-du general.cfg  mozart4.cfg 20140316 20140615"
print command
os.system(command)

command = "bc-mozart4-gas general.cfg  mozart4.cfg 20140316 20140616"
print command
os.system(command)

command = "bc-mozart4-ss general.cfg  mozart4.cfg 20140316 20140615"
print command
os.system(command)

command = "ic-mozart4-aer general.cfg  mozart4.cfg 20140316"
print command
os.system(command)

command = "ic-mozart4-du general.cfg  mozart4.cfg 20140316"
print command
os.system(command)

command = "ic-mozart4-gas general.cfg  mozart4.cfg 20140316"
print command
os.system(command)

command = "ic-mozart4-ss general.cfg  mozart4.cfg 20140316"
print command
os.system(command)
